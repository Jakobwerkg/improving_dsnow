#-------------------------------------------------------------------------
# Parameter optimization for deltasnow model
# Optimize deltasnow model w.r.t. observed vs modeled SWE (RMSE).
# - No duplicate-timestamp handling (as requested)
# - Robust column mapping for Hobs / SWEobs variants
# - Flexible winter split modes
#
# Harald Schellander, 08.2025
#-------------------------------------------------------------------------

# --- Libraries -------------------------------------------------------------
suppressPackageStartupMessages({
  library(optimx)
  library(zoo)
  library(foreach)<
  library(doParallel)
  library(lubridate)
  library(nixmass)
  library(tidyverse)
})

# --- Working directory + user config ---------------------------------------
setwd("/Users/jakobwerkgarner/code/Master_Delta/calibration")

# Free-form tag that gets embedded in log file names
calib_comment <- "dsnow_calib_Win21_Mag25_alt_winters_forced_high_rho0_drpped_stations"

# Split configuration
# 1) "alt_winters"     : alternating winters (even index -> fit, odd -> val)
# 2) "random"    : random split of winters per station (uses ratio.fit)
# 3) "random_stations" : whole stations go fit/val (uses ratio.fit)
# 4) "chronological"   : earliest X% winters -> fit, remaining -> val (uses ratio.fit)
split.type  <- "alt_winters"     # "alt_winters" | "random_50_50" | "random_stations" | "chronological"
split.seed  <- 12                # seed for randomness in modes 2–4
ratio.fit   <- 0.5              # proportion for fit set for modes 2–4

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
  "Adelboden",          # very low mean rho0 (outlier)
  "Bourg_St_Pierre",    # very low mean rho0 (outlier)
  "Davos Flueelastr",   # duplicate spelling
  "Davos_Flueelastr",   # duplicate spelling
  "Juf",                # too high altitude
  "Muerren",            # duplicate spelling
  "Muenster",           # duplicate spelling
  "Münster",            # excluded in validation
  "Mürren",             # excluded in validation
  "Obersaxen",          # excluded in validation
  "Pusserein",          # very low mean rho0 (outlier)
  "St_Margrethenberg",  # very low mean rho0 (outlier)
  "Ulrichen",           # duplicate station
  "Weissfluh Joch",     # duplicate spelling
  "Weissfluhjoch",      # too high altitude
  "Zuoz",               # very low mean rho0
  "kuehtai"             # excluded (typo/variant)
))



d_obs[intersect(names(d_obs), remove_stations)] <- NULL

# Apply optional cutoff at start of a hydrological season
if (!is.na(cutoff.year)) {
  cutoff_date <- as.Date(paste0(cutoff.year, season.start))  # e.g., "2020-08-01"
  d_obs <- lapply(d_obs, function(z) window(z, start = cutoff_date))
  d_obs <- Filter(NROW, d_obs)
}

# --- Helpers: winters + splitting ------------------------------------------
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
    if (nrow(win) < 365 && (as.numeric(win$Hobs[1]) > edge_thresh ||
                            as.numeric(win$Hobs[nrow(win)]) > edge_thresh)) next()
    res[[length(res) + 1]] <- win
  }
  res
}

# Central splitter that returns two lists (fit/val) of zoo objects per station
build_splits <- function(d_obs,
                         split.type = "alt_winters",
                         ratio.fit = 0.5,
                         split.seed = 12,
                         season.start = "-08-01",
                         season.end   = "-07-31") {
  stopifnot(split.type %in% c("alt_winters","random","random_stations","chronological"))
  set.seed(split.seed)
  
  stations    <- names(d_obs)
  d_obs_fit   <- vector("list", length(stations)); names(d_obs_fit) <- stations
  d_obs_val   <- vector("list", length(stations)); names(d_obs_val) <- stations

  # Output directory
temp_dir <- file.path(base_dir, "calibration/calibration_R/temp_data", calib_comment)
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

# Save fit + val splits
save(d_obs_fit, file = file.path(temp_dir, "d_obs_fit.rda"), compress = "xz")
save(d_obs_val, file = file.path(temp_dir, "d_obs_val.rda"), compress = "xz")

cat("Saved:\n",
    file.path(temp_dir, "d_obs_fit.rda"), "\n",
    file.path(temp_dir, "d_obs_val.rda"), "\n")
  
  # Pre-pick station assignment if mode is random_stations
  fit_station_set <- NULL
  if (split.type == "random_stations") {
    n_fit_st <- max(1, floor(length(stations) * ratio.fit))
    fit_station_set <- sample(stations, size = n_fit_st)
  }
  
  for (n in stations) {
    d <- d_obs[[n]]
    
    # No duplicate timestamp processing (requested).
    # Just sort to ensure stable order.
    d <- d[order(index(d)), ]
    
    winters <- get_valid_winters(d, season.start, season.end)
    if (length(winters) == 0) next()
    
    if (split.type == "alt_winters") {
      # Alternate winters: even index -> fit, odd -> val
      fit_list <- winters[which(seq_along(winters) %% 2 == 0)]
      val_list <- winters[which(seq_along(winters) %% 2 == 1)]
      
    } else if (split.type == "random_stations") {
      if (n %in% fit_station_set) {
        fit_list <- winters; val_list <- list()
      } else {
        fit_list <- list();  val_list <- winters
      }
      
    } else if (split.type == "random") {
      nW     <- length(winters)
      n_fit  <- max(1, floor(nW * ratio.fit))
      idxFit <- sort(sample.int(nW, n_fit))
      idxVal <- setdiff(seq_len(nW), idxFit)
      fit_list <- if (length(idxFit)) winters[idxFit] else list()
      val_list <- if (length(idxVal)) winters[idxVal] else list()
      
    } else if (split.type == "chronological") {
      nW    <- length(winters)
      n_fit <- max(1, floor(nW * ratio.fit))
      fit_list <- winters[seq_len(min(n_fit, nW))]
      val_list <- if (nW > n_fit) winters[(n_fit + 1):nW] else list()
    }
    
    if (length(fit_list)) d_obs_fit[[n]] <- do.call(rbind, fit_list)
    if (length(val_list)) d_obs_val[[n]] <- do.call(rbind, val_list)
  }
  
  # Drop stations with empty splits to keep things tidy
  d_obs_fit <- Filter(NROW, d_obs_fit)
  d_obs_val <- Filter(NROW, d_obs_val)
  list(fit = d_obs_fit, val = d_obs_val)
}

# --- Build the requested split ---------------------------------------------
message(sprintf("Split mode: %s (seed=%d, ratio.fit=%.2f)", split.type, split.seed, ratio.fit))
splits   <- build_splits(d_obs, split.type, ratio.fit, split.seed, season.start, season.end)
d_obs_fit <- splits$fit
d_obs_val <- splits$val

# Optionally exclude Sta.Maria from the fit set (analysis-only)
d_obs_fit[["Sta.Maria"]] <- NULL

# --- Hydrological blocks ----------------------------------------------------
start_of_block <- 8  # August

# Map calendar dates to hydrological season year
set_season <- function(x, m) {
  x  <- as.character(x)
  yr <- as.POSIXlt(x)$year + 1900
  mt <- as.POSIXlt(x)$mon + 1
  ifelse(mt < m, yr - 1, yr)
}

# --- Robust column mapping for building the fit tibble ----------------------
# Map plausible column name variants to the canonical fields we need
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
message(sprintf("Built fit tibble with %d rows across %d stations.",
                nrow(d_obs_fit_tibble), dplyr::n_distinct(d_obs_fit_tibble$name)))

# --- Objective (RMSE) -------------------------------------------------------
minimize_score <- function(data, par, scale, verbose = FALSE) {
  par <- par * scale
  cat(par)  # follow the optimization process
  
  ll <- foreach(
    s = unique(data$name),
    .packages = c("dplyr", "zoo", "nixmass", "foreach")
  ) %dopar% {
    if (verbose) cat(paste0(s, " ..."))
    data1 <- dplyr::filter(data, name == s)
    
    # Loop hydrological seasons serially to keep memory bounded
    l <- foreach(
      y = unique(data1$block),
      .packages = c("dplyr", "zoo", "nixmass")
    ) %do% {
      left <- data1 |>
        dplyr::filter(block == y) |>
        dplyr::select(date, hs, swe_obs)
      
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
              k         = par[5],
              tau       = par[6],
              eta.null  = par[7]
            ),
            dyn_rho_max = FALSE
          )
      }, error = function(e) rep(NA, nrow(joined)))
      
      na.omit(cbind(joined, swe_mod))
    }
    do.call(rbind, l)
  }
  
  dff  <- do.call(rbind, ll)
  rmse <- with(dff, sqrt(mean((swe_mod - swe_obs)^2)))
  bias <- with(dff, abs(mean(swe_mod - swe_obs)))
  cat(paste0(" |bias|=", bias, " rmse=", rmse, "\n"))
  rmse
}

# --- Metrics helper (no console printing) -----------------------------------
eval_metrics_for_par <- function(data, par_unscaled) {
  ll <- foreach(
    s = unique(data$name),
    .packages = c("dplyr", "zoo", "nixmass", "foreach")
  ) %dopar% {
    data1 <- dplyr::filter(data, name == s)
    l <- foreach(
      y = unique(data1$block),
      .packages = c("dplyr", "zoo", "nixmass")
    ) %do% {
      left   <- dplyr::filter(data1, block == y) |> dplyr::select(date, hs, swe_obs)
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
      }, error = function(e) rep(NA, nrow(joined)))
      na.omit(cbind(joined, swe_mod))
    }
    do.call(rbind, l)
  }
  dff  <- do.call(rbind, ll)
  rmse <- with(dff, sqrt(mean((swe_mod - swe_obs)^2)))
  bias <- with(dff, mean(swe_mod - swe_obs))
  list(rmse = rmse, bias = bias)
}

# --- Parameter setup (scaled) ----------------------------------------------
par_delta <- c(
  rho.max   = 396.9521,
  rho.null  = 100,
  c.ov      = 0.0004986025,
  k.ov      = 0.227786,
  k.exp     = 0.0290256,  # NOTE: name differs, used positionally as 'k' below
  tau       = 0.02489556,
  eta.null  = 8792253
)
par_scale <- c(1000, 1000, 0.001, 1, 0.1, 0.1, 1e7)
par_delta <- par_delta / par_scale
lower <- c(300, 90, 0, 0.01, 0.01, 0.01, 1e6) / par_scale
upper <- c(600, 200, 0.001, 10, 0.2, 0.2, 2e7) / par_scale
in_bounds <- ifelse(lower <= par_delta & upper >= par_delta, "in bounds", "ERROR")
print(cbind(lower, par_delta, upper, in_bounds))

# --- Optimization -----------------------------------------------------------
nc <- detectCores(logical = TRUE)
cl <- makeCluster(nc)
registerDoParallel(cl)
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)

opt <- optimx(
  fn      = minimize_score,
  data    = d_obs_fit_tibble,
  par     = par_delta,
  scale   = par_scale,
  verbose = FALSE,
  method = c('L-BFGS-B','bobyqa'),  # defaults used
  control = list(trace = 4, follow.on = TRUE),
  lower = lower,
  upper = upper
)
saveRDS(opt, file = "opt_results_orig.rds")

# --- Simple logging ---------------------------------------------------------
opt_df    <- as.data.frame(opt)
par_names <- names(par_delta)
best_idx  <- which.min(opt_df$value)

best_par_scaled   <- as.numeric(opt_df[best_idx, par_names, drop = TRUE])
names(best_par_scaled) <- par_names
best_par_unscaled <- best_par_scaled * par_scale

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

base_log_dir <- "calibration/calibration_R/calibration_log/"
if (!dir.exists(base_log_dir)) dir.create(base_log_dir, recursive = TRUE, showWarnings = FALSE)

ts_colon     <- format(Sys.time(), "%Y_%m_%d_%H%M")
safe_comment <- gsub("[^A-Za-z0-9_\\-]+", "_", calib_comment)
fname        <- paste0(ts_colon, "_calib_log_", safe_comment, ".txt")
log_file     <- file.path(base_log_dir, fname)

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
  paste0("split.type = ", split.type),
  paste0("split.seed = ", split.seed),
  paste0("ratio.fit = ", ratio.fit),
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
