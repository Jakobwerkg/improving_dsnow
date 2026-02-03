#-------------------------------------------------------------------------
# Parameter optimization for deltasnow model (ALT-WINTER split only)
# Optimize deltasnow model w.r.t. observed vs modeled SWE (RMSE).
# - No duplicate-timestamp handling (as requested)
# - Robust column mapping for Hobs / SWEobs variants
# - Alternating winters only (even index -> fit, odd -> val)
#
# Harald Schellander, 08.2025  |
#-------------------------------------------------------------------------

# --- Libraries -------------------------------------------------------------
suppressPackageStartupMessages({
  library(optimx)
  library(zoo)
  library(foreach)
  library(doParallel)
  library(lubridate)
  library(nixmass)
  library(tidyverse)
})

# --- Working directory + user config ---------------------------------------
setwd("/Users/jakobwerkgarner/code/Master_Delta/calibration")

# Free-form tag that gets embedded in log file names (unchanged)
calib_comment <- "Win21_Mag25_data_exeriment_tigth_rho_experimental"

# Hydrological season window (Aug 1 .. Jul 31)
season.start <- "-08-01"
season.end   <- "-07-31"

# Optional cutoff at season start (keep >= cutoff_date). Set to NA to disable.
cutoff.year <- NA_integer_       # e.g., 2020 for "2020-08-01"; NA = keep all

# --- Data loading & station filtering --------------------------------------
rda_path <- "/Users/jakobwerkgarner/code/Master_Delta/calibration_data/exports/d_obs.rda"
d_obs    <- get(load(rda_path))  # list(name -> zoo), columns should include depth+SWE

# Normalize station names a bit (trim, collapse spaces)
clean_names <- function(x) trimws(gsub("\\s+", " ", x))
names(d_obs) <- clean_names(names(d_obs))

# Remove stations entirely (typos, duplicates, or outliers)
remove_stations <- clean_names(c(
  "Adelboden",
  "Bourg_St_Pierre",
  "Davos Flueelastr", "Davos_Flueelastr",
  "Juf",
  "Muerren","Muenster","Münster","Mürren",
  "Obersaxen",
  "Pusserein",
  "St_Margrethenberg",
  "Ulrichen",
  "Weissfluh Joch","Weissfluhjoch",
  "Zuoz",
  "kuehtai"
))
d_obs[intersect(names(d_obs), remove_stations)] <- NULL

# Apply optional cutoff at start of a hydrological season
if (!is.na(cutoff.year)) {
  cutoff_date <- as.Date(paste0(cutoff.year, season.start))  # e.g., "2020-08-01"
  d_obs <- lapply(d_obs, function(z) window(z, start = cutoff_date))
  d_obs <- Filter(NROW, d_obs)
}

# --- Helpers: winters + splitting (ALT-WINTERS ONLY) ------------------------
# Build all valid winters for a station (as a list of zoo objects).
# "Valid" => >= 200 rows; if < 365 rows, reject when first/last Hobs > edge_thresh.
get_valid_winters <- function(d, season.start, season.end, edge_thresh = 0.05) {
  years <- sort(unique(lubridate::year(index(d))))
  if (length(years) < 2) return(list())
  res <- list()
  for (y in years[seq_len(length(years) - 1)]) {
    win <- subset(
      d,
      index(d) >= as.Date(paste0(y, season.start)) &
        index(d) <  as.Date(paste0(y + 1, season.end)) + 1
    )
    if (nrow(win) < 200) next()
    if ("Hobs" %in% colnames(win)) {
      if (nrow(win) < 365 && (as.numeric(win$Hobs[1]) > edge_thresh ||
                              as.numeric(win$Hobs[nrow(win)]) > edge_thresh)) next()
    }
    res[[length(res) + 1]] <- win
  }
  res
}

# Strictly alternating winters: even index -> fit, odd -> val
build_splits_alt <- function(d_obs, season.start = "-08-01", season.end = "-07-31") {
  stations    <- names(d_obs)
  d_obs_fit   <- vector("list", length(stations)); names(d_obs_fit) <- stations
  d_obs_val   <- vector("list", length(stations)); names(d_obs_val) <- stations
  
  for (n in stations) {
    d <- d_obs[[n]]
    d <- d[order(index(d)), ]  # no duplicate handling, just stable order
    
    winters <- get_valid_winters(d, season.start, season.end)
    if (length(winters) == 0) next()
    
    fit_list <- winters[which(seq_along(winters) %% 2 == 0)]
    val_list <- winters[which(seq_along(winters) %% 2 == 1)]
    
    if (length(fit_list)) d_obs_fit[[n]] <- do.call(rbind, fit_list)
    if (length(val_list)) d_obs_val[[n]] <- do.call(rbind, val_list)
  }
  
  d_obs_fit <- Filter(NROW, d_obs_fit)
  d_obs_val <- Filter(NROW, d_obs_val)
  list(fit = d_obs_fit, val = d_obs_val)
}

message("Split mode: alt_winters")
splits    <- build_splits_alt(d_obs, season.start, season.end)
d_obs_fit <- splits$fit
d_obs_val <- splits$val

# Optionally exclude Sta.Maria from the fit set (analysis-only)
d_obs_fit[["Sta.Maria"]] <- NULL

# --- Hydrological blocks ----------------------------------------------------
start_of_block <- 8  # August

# Map calendar dates to hydrological season year
set_season <- function(x, m) {
  x  <- as.Date(x)
  yr <- as.integer(format(x, "%Y"))
  mt <- as.integer(format(x, "%m"))
  ifelse(mt < m, yr - 1, yr)
}

# --- Robust column mapping for building the fit tibble ----------------------
.pick_obs_cols <- function(core_df, station_name) {
  nms <- tolower(names(core_df))
  h_candidates   <- c("hobs","hs","h","hsobs","snowdepth","snow_depth","h_cm","hs_cm")
  swe_candidates <- c("sweobs","swe","swe_obs","swe_mm","swe_kgm2","swe_kg/m2")
  
  h_idx   <- which(nms %in% h_candidates)[1]
  swe_idx <- which(nms %in% swe_candidates)[1]
  
  if (is.na(h_idx) || is.na(swe_idx)) {
    message(sprintf("Station '%s' missing Hobs/SWEobs columns. Available: %s",
                    station_name, paste(names(core_df), collapse = ", ")))
    return(NULL)
  }
  list(h_col = names(core_df)[h_idx], swe_col = names(core_df)[swe_idx])
}

# Build calibration tibble (skips stations that don't have the needed columns)
build_fit_tibble <- function(d_obs_fit, start_of_block = 8) {
  out <- lapply(names(d_obs_fit), function(nm) {
    x <- d_obs_fit[[nm]]
    if (is.null(x) || NROW(x) == 0) return(NULL)
    
    cd <- as.data.frame(coredata(x))
    colmap <- .pick_obs_cols(cd, nm)
    if (is.null(colmap)) return(NULL)
    
    hs_raw  <- cd[[colmap$h_col]]    # assume cm -> convert to m below
    swe_raw <- cd[[colmap$swe_col]]  # assume mm or kg/m^2 (treated equivalent)
    
    tibble(
      date    = index(x),
      name    = nm,
      hs      = as.numeric(hs_raw) / 100,   # cm -> m
      swe_obs = as.numeric(swe_raw)
    ) |>
      mutate(block = set_season(date, start_of_block))
  })
  bind_rows(out)
}

# Build the fit tibble
d_obs_fit_tibble <- build_fit_tibble(d_obs_fit, start_of_block)
if (is.null(d_obs_fit_tibble) || nrow(d_obs_fit_tibble) == 0) {
  stop("No usable data in fit tibble after filtering/mapping.")
}
message(sprintf("Built fit tibble with %d rows across %d stations.",
                nrow(d_obs_fit_tibble), dplyr::n_distinct(d_obs_fit_tibble$name)))

# --- Parameter handling (names, scaling, bounds) ---------------------------
# Base (UNSCALED) parameter guess with stable names (k.exp used positionally as 'k').
# --- Parameter handling (names, scaling, bounds) ---------------------------

# Base (UNSCALED) parameter guess (real-scale)
par_guess_unscaled <- c(
  rho.max   = 410,   # kg m^-3  (ρmax)
  rho.null  = 110,        # kg m^-3  (ρ0)
  c.ov      = 0.0004986025, # Pa^-1   (c_ov)
  k.ov      = 0.227786,   # –        (k_ov)
  k.exp     = 0.0290256,  # m^3 kg^-1 (k)
  tau       = 0.02489556, # m (≈2.49 cm)  (τ)
  eta.null  = 8792253     # Pa·s     (η0)
)

# Per-parameter scale for optimizer (optional; guesses remain real-scale)
par_scale <- c(
  rho.max   = 100,
  rho.null  = 100,
  c.ov      = 0.001,
  k.ov      = 1,
  k.exp     = 0.1,
  tau       = 0.1,   # in meters
  eta.null  = 1e7
)

# ---- Tight UNSCALED bounds per Winkler et al. (2021) ----
lower_unscaled <- c(
  rho.max   = 400,     # kg m^-3
  rho.null  = 100,      # kg m^-3
  c.ov      = 0.0,     # Pa^-1
  k.ov      = 0.01,    # –
  k.exp     = 0.01,    # m^3 kg^-1
  tau       = 0.01,    # m (1 cm)
  eta.null  = 1e6      # Pa·s
)
upper_unscaled <- c(
  rho.max   = 500,     # kg m^-3
  rho.null  = 200,     # kg m^-3
  c.ov      = 1e-3,    # Pa^-1
  k.ov      = 10,      # –
  k.exp     = 0.20,    # m^3 kg^-1
  tau       = 0.10,    # m (20 cm)
  eta.null  = 15e6     # Pa·s
)

# Helper to validate and scale vectors consistently
.check_and_scale <- function(vec_unscaled, ref_names, par_scale, label) {
  if (is.null(names(vec_unscaled)) || !all(ref_names %in% names(vec_unscaled))) {
    missing <- setdiff(ref_names, names(vec_unscaled))
    stop(sprintf("%s must be a NAMED numeric vector with names: %s. Missing: %s",
                 label, paste(ref_names, collapse=", "),
                 ifelse(length(missing)==0,"<none>",paste(missing, collapse=", "))))
  }
  vec_unscaled <- vec_unscaled[ref_names]
  vec_scaled   <- vec_unscaled / par_scale[ref_names]
  as.numeric(vec_scaled)
}

# Build scaled par / bounds with strict name checking
par_names <- names(par_guess_unscaled)
par_delta <- .check_and_scale(par_guess_unscaled, par_names, par_scale, "par_guess_unscaled")
lower     <- .check_and_scale(lower_unscaled,     par_names, par_scale, "lower_unscaled")
upper     <- .check_and_scale(upper_unscaled,     par_names, par_scale, "upper_unscaled")

# Sanity check: initial guess inside bounds
in_bounds <- ifelse(lower <= par_delta & upper >= par_delta, "in bounds", "ERROR")
print(cbind(lower=signif(lower,6), par=signif(par_delta,6), upper=signif(upper,6), in_bounds))

# --- Objective (RMSE) with finite penalty ----------------------------------
.big_penalty <- 1e12  # large finite value to keep optimizers alive

minimize_score <- function(data, par, scale, verbose = FALSE) {
  # par is scaled inside optim; convert back to unscaled here
  par <- par * scale
  cat(sprintf("par=%s", paste(signif(par,6), collapse=",")))
  
  ll <- foreach(
    s = unique(data$name),
    .packages = c("dplyr","zoo","nixmass")
  ) %dopar% {
    if (verbose) cat(paste0(" ", s))
    data1 <- dplyr::filter(data, name == s)
    
    yrs <- unique(data1$block)
    l <- lapply(yrs, function(y) {
      left   <- dplyr::filter(data1, block == y) |>
        dplyr::select(date, hs, swe_obs)
      if (nrow(left) == 0) return(NULL)
      
      right  <- data.frame(date = seq(min(left$date), max(left$date), by = "1 day"))
      joined <- dplyr::right_join(left, right, by = "date")
      
      swe_mod <- tryCatch({
        joined |>
          dplyr::select(date, hs) |>
          mutate(date = as.character(date)) |>
          nixmass::swe.delta.snow(
            model_opts = list(
              rho.max   = par[1],
              rho.null  = par[2],
              c.ov      = par[3],
              k.ov      = par[4],
              k         = par[5],    # k.exp positionally
              tau       = par[6],
              eta.null  = par[7]
            ),
            dyn_rho_max = FALSE
          )
      }, error = function(e) rep(NA_real_, nrow(joined)))
      
      out <- suppressWarnings(na.omit(cbind(joined, swe_mod)))
      if (!is.null(out) && nrow(out) > 0) out else NULL
    })
    if (length(l)) do.call(rbind, Filter(NROW, l)) else NULL
  }
  
  ll <- Filter(NROW, ll)
  if (!length(ll)) {
    cat(" | no data -> rmse=Penalty\n")
    return(.big_penalty)
  }
  dff  <- tryCatch(do.call(rbind, ll), error = function(e) NULL)
  if (is.null(dff) || !all(c("swe_mod","swe_obs") %in% colnames(dff))) {
    cat(" | missing columns -> rmse=Penalty\n")
    return(.big_penalty)
  }
  rmse <- with(dff, sqrt(mean((swe_mod - swe_obs)^2)))
  bias <- with(dff, abs(mean(swe_mod - swe_obs)))
  if (!is.finite(rmse)) {
    cat(" | rmse non-finite -> Penalty\n")
    return(.big_penalty)
  }
  cat(sprintf(" | |bias|=%g rmse=%g\n", bias, rmse))
  rmse
}

# --- Metrics helper (no console printing) -----------------------------------
eval_metrics_for_par <- function(data, par_unscaled) {
  ll <- foreach(
    s = unique(data$name),
    .packages = c("dplyr","zoo","nixmass")
  ) %dopar% {
    data1 <- dplyr::filter(data, name == s)
    
    yrs <- unique(data1$block)
    l <- lapply(yrs, function(y) {
      left   <- dplyr::filter(data1, block == y) |> dplyr::select(date, hs, swe_obs)
      if (nrow(left) == 0) return(NULL)
      right  <- data.frame(date = seq(min(left$date), max(left$date), by = "1 day"))
      joined <- dplyr::right_join(left, right, by = "date")
      
      swe_mod <- tryCatch({
        joined |>
          dplyr::select(date, hs) |>
          mutate(date = as.character(date)) |>
          nixmass::swe.delta.snow(
            model_opts = list(
              rho.max  = par_unscaled[1],
              rho.null = par_unscaled[2],
              c.ov     = par_unscaled[3],
              k.ov     = par_unscaled[4],
              k        = par_unscaled[5],
              tau      = par_unscaled[6],
              eta.null = par_unscaled[7]
            ),
            dyn_rho_max = FALSE
          )
      }, error = function(e) rep(NA_real_, nrow(joined)))
      out <- suppressWarnings(na.omit(cbind(joined, swe_mod)))
      if (!is.null(out) && nrow(out) > 0) out else NULL
    })
    if (length(l)) do.call(rbind, Filter(NROW, l)) else NULL
  }
  ll <- Filter(NROW, ll)
  if (!length(ll)) return(list(rmse=NA_real_, bias=NA_real_))
  dff  <- do.call(rbind, ll)
  rmse <- with(dff, sqrt(mean((swe_mod - swe_obs)^2)))
  bias <- with(dff, mean(swe_mod - swe_obs))
  list(rmse = rmse, bias = bias)
}

# --- Parallel backend -------------------------------------------------------
nc <- max(1L, parallel::detectCores(logical = TRUE))
cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl)

# Ensure workers have required packages loaded
parallel::clusterEvalQ(cl, {
  library(foreach)
  library(dplyr)
  library(zoo)
  library(nixmass)
})

on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

# --- Optimization (with graceful fallback) ----------------------------------
opt_try <- try(optimx(
  fn      = minimize_score,
  data    = d_obs_fit_tibble,
  par     = par_delta,
  scale   = as.numeric(par_scale[par_names]),
  verbose = FALSE,
  method  = c('L-BFGS-B','bobyqa'),
  control = list(trace = 4, follow.on = TRUE, badval = .big_penalty),
  lower   = lower,
  upper   = upper
), silent = TRUE)

# If optimx errored or produced unusable result, synthesize a fallback row
make_fallback_result <- function() {
  met0 <- eval_metrics_for_par(d_obs_fit_tibble, par_guess_unscaled)
  # Build a 1-row data.frame with expected columns
  out <- data.frame(matrix(nrow = 1, ncol = 0))
  out[,"value"] <- if (is.finite(met0$rmse)) met0$rmse else .big_penalty
  for (nm in par_names) out[, nm] <- par_delta[which(names(par_delta)==nm)]
  out[,"method"] <- "optimization_failed"
  out
}

use_fallback <- function(df) {
  if (inherits(opt_try, "try-error")) return(TRUE)
  if (!is.data.frame(df) || nrow(df) == 0) return(TRUE)
  # If all values are NA or equal to penalty, treat as failure
  if (!("value" %in% names(df))) return(TRUE)
  if (all(!is.finite(df$value)) || all(df$value >= .big_penalty*0.999)) return(TRUE)
  FALSE
}

opt_df <- if (use_fallback(as.data.frame(opt_try))) make_fallback_result() else as.data.frame(opt_try)
saveRDS(opt_df, file = "opt_results_orig.rds")  # keep filename the same

# --- Simple logging (same structure/labels/order as original) --------------
best_idx <- which.min(opt_df$value)

# Robust parameter extraction even if optimizer didn't return named cols
param_cols <- intersect(colnames(opt_df), par_names)
if (length(param_cols) != length(par_names)) {
  # try to guess numeric param columns (excluding known meta cols)
  numcols <- names(opt_df)[sapply(opt_df, is.numeric)]
  numcols <- setdiff(numcols, c("value","fevals","gevals","nitns","xtimes"))
  if (length(numcols) >= length(par_names)) {
    param_cols <- numcols[seq_len(length(par_names))]
  } else {
    # final fallback: start guess
    best_par_scaled <- par_delta
    names(best_par_scaled) <- par_names
  }
}
if (!exists("best_par_scaled")) {
  best_par_scaled <- as.numeric(opt_df[best_idx, param_cols, drop = TRUE])
  names(best_par_scaled) <- if (length(param_cols)==length(par_names)) par_names else param_cols
  # reorder to par_names if possible
  if (all(par_names %in% names(best_par_scaled))) {
    best_par_scaled <- best_par_scaled[par_names]
  }
}
best_par_unscaled <- best_par_scaled * as.numeric(par_scale[par_names])

# Recompute metrics at best (or fallback) parameters
metrics    <- eval_metrics_for_par(d_obs_fit_tibble, best_par_unscaled)
final_rmse <- as.numeric(metrics$rmse)
final_bias <- as.numeric(metrics$bias)

wd <- getwd()
all_stations <- paste(names(d_obs), collapse = ", ")
fit_stations <- paste(unique(d_obs_fit_tibble$name), collapse = ", ")

pkg_list <- c("optimx","zoo","foreach","doParallel","lubridate","nixmass","tidyverse")
pkg_versions <- sapply(pkg_list, function(p) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NA")
  paste0(p, "=", v)
})

os_info <- tryCatch(paste(unname(Sys.info()), collapse = " | "), error = function(e) R.version$platform)
R_str   <- R.version.string

rds_path     <- file.path(getwd(), "opt_results_orig.rds")
rda_path_str <- rda_path

base_log_dir <- "/Users/jakobwerkgarner/code/Master_Delta/calibration/calibration_logsx/A_log_files/"
if (!dir.exists(base_log_dir)) dir.create(base_log_dir, recursive = TRUE, showWarnings = FALSE)

ts_colon     <- format(Sys.time(), "%Y_%m_%d_%H%M")
safe_comment <- gsub("[^A-Za-z0-9_\\-]+", "_", calib_comment)
fname        <- paste0(ts_colon, "_calib_log_", safe_comment, ".txt")
log_file     <- "/Users/jakobwerkgarner/code/Master_Delta/calibration/calibration_logsx/A_log_files//2025_10_24_2014_calib_log_Win21_Mag25_data_exeriment_100rho.txt"
  # file.path(base_log_dir, fname)

method_str <- if ("method" %in% names(opt_df)) as.character(opt_df$method[best_idx]) else "n/a"

lines <- c(
  "===== Calibration Log =====",
  paste0("Timestamp: ", ts_colon),
  paste0("Working directory: ", wd),
  paste0("RDA path: ", rda_path_str),
  paste0("Results RDS path: ", rds_path),
  "",
  "---- Data / Stations ----",
  paste0("All stations (after initial filtering): ", all_stations),
  paste0("Fit stations actually used: ", fit_stations),
  paste0("Fit rows (n): ", nrow(d_obs_fit_tibble)),
  "",
  "---- Split configuration ----",
  paste0("split.type = ", "alt_winters"),
  paste0("split.seed = ", 12),   # kept for log continuity
  paste0("ratio.fit = ", 0.5),   # kept for log continuity
  "",
  "---- Optimization Result ----",
  paste0("Best method: ", method_str),
  paste0("Objective (RMSE from fn)       = ", signif(opt_df$value[best_idx], 6)),
  paste0("Final RMSE (recomputed)        = ", signif(final_rmse, 6)),
  paste0("Final bias (mean error)        = ", signif(final_bias, 6)),
  "",
  "Optimal parameters (UNSCALED):",
  paste(sprintf("  %s = %g", names(best_par_unscaled), signif(best_par_unscaled, 10)), collapse = "\n"),
  "",
  "---- Compute Environment ----",
  paste0("R: ", R_str),
  paste0("OS: ", os_info),
  paste0("CPU cores used: ", nc),
  paste("Packages:", paste(pkg_versions, collapse = ", ")),
  "",
  "Comment:",
  calib_comment,
  "===== End of Log ====="
)
cat(paste(lines, collapse = "\n"), file = log_file, sep = "\n")
message(sprintf("Calibration log written to: %s", log_file))
