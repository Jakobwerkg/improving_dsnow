# ============================================================================
# DeltaSnow parameter optimization (robust, reproducible, logged)
# ============================================================================
# - Station coverage filtering (keeps only stations with reliable SWE coverage)
# - Hard bounds + soft interior penalty (stable optimizer behavior)
# - Always returns finite values (no NA/Inf to optimizer)
# - Reproducible parallel foreach with doRNG
# - Per-season timeouts (bad seasons can’t stall the run)
# - Full logs + a CSV of station coverage diagnostics
# ============================================================================

suppressPackageStartupMessages({
  # Auto-install small helpers if missing
  if (!requireNamespace("doRNG", quietly = TRUE)) install.packages("doRNG")
  if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
  
  library(optimx)
  library(zoo)
  library(foreach)
  library(doParallel)
  library(doRNG)      # reproducible %dopar% streams
  library(lubridate)  # year(), month(), etc.
  library(nixmass)    # swe.delta.snow()
  library(tidyverse)  # dplyr/tibble/arrange/filter/mutate/
  library(R.utils)    # withTimeout()
})

# ------------------ Working dir + user configuration ------------------------
setwd("/Users/jakobwerkgarner/code/Master_Delta/calibration")

calib_comment <- "dsnow_calib_Win21_Mag25_cov_filter_robust"

# Split configuration
# Tip while stabilizing: start with "chronological" (no random sampling)
split.type  <- "chronological"   # "alt_winters" | "random" | "random_stations" | "chronological"
split.seed  <- 12
ratio.fit   <- 0.8               # used by modes other than "alt_winters"

# Hydrological season window (Aug 1 .. Jul 31)
season.start <- "-08-01"
season.end   <- "-07-31"

# Optional cutoff at season start (keep >= cutoff_date). Set to NA to disable.
cutoff.year <- NA_integer_       # e.g., 2020

# ---------------------- Data loading + base filtering -----------------------
rda_path <- "/Users/jakobwerkgarner/code/Master_Delta/calibration_data/exports/d_obs.rda"
d_obs    <- get(load(rda_path))      # list(station_name -> zoo)

# Normalize names (trim + collapse spaces)
clean_names <- function(x) trimws(gsub("\\s+", " ", x))
names(d_obs) <- clean_names(names(d_obs))

# Remove known dupes/outliers/typos for your dataset
remove_stations <- clean_names(c(
  "Davos Flueelastr","Davos_Flueelastr", # duplicates
  "Juf",                                  # altitude outlier
  "Weissfluh Joch","Weissfluhjoch",       # duplicates/altitude
  "Zuoz",                                 # very low mean rho0
  "kuehtai"                               # typo/variant
))
d_obs[intersect(names(d_obs), remove_stations)] <- NULL

# Optional cutoff at start of a hydrological season
if (!is.na(cutoff.year)) {
  cutoff_date <- as.Date(paste0(cutoff.year, season.start))
  d_obs <- lapply(d_obs, function(z) window(z, start = cutoff_date))
  d_obs <- Filter(NROW, d_obs)
}

# ------------------- Helpers: small utilities -------------------------------
# Coerce zoo::coredata to a data.frame regardless of it being vector/matrix
.to_df <- function(cd) {
  if (is.data.frame(cd)) return(cd)
  if (is.matrix(cd)) return(as.data.frame(cd))
  data.frame(value = cd, check.names = FALSE)
}

# Hydrological year label (e.g., Aug–Jul). m=8 means August start.
set_season <- function(x, m) {
  yr <- lubridate::year(x); mo <- lubridate::month(x)
  ifelse(mo < m, yr - 1, yr)
}

# Robustly pick HS/SWE columns by tolerant name matching
.pick_obs_cols <- function(core_df, station_name) {
  nms <- tolower(names(core_df))
  h_candidates   <- c("hobs","hs","h","hsobs","snowdepth","snow_depth","h_cm","hs_cm")
  swe_candidates <- c("sweobs","swe","swe_obs","swe_mm","swe_kgm2","swe_kg/m2","swe [mm]")
  h_idx   <- which(nms %in% h_candidates)[1]
  swe_idx <- which(nms %in% swe_candidates)[1]
  if (is.na(h_idx) || is.na(swe_idx)) {
    message(sprintf("Skipping '%s': missing HS/SWE columns. Available: %s",
                    station_name, paste(names(core_df), collapse=", ")))
    return(NULL)
  }
  list(h_col = names(core_df)[h_idx], swe_col = names(core_df)[swe_idx])
}

# ------------------- Station coverage diagnostics (subset) ------------------
coverage_by_station <- function(d_obs,
                                start_month = 8,
                                min_days_per_season = 120,
                                max_gap_days = 10) {
  out <- lapply(names(d_obs), function(nm) {
    z <- d_obs[[nm]]
    if (is.null(z) || NROW(z) == 0) return(NULL)
    df <- .to_df(zoo::coredata(z))
    colmap <- .pick_obs_cols(df, nm)
    if (is.null(colmap)) return(NULL)
    
    dat <- tibble(
      date    = as.Date(index(z)),
      hs      = suppressWarnings(as.numeric(df[[colmap$h_col]])) / 100,  # cm -> m
      swe_obs = suppressWarnings(as.numeric(df[[colmap$swe_col]]))       # mm
    ) |>
      arrange(date) |>
      mutate(hyear = ifelse(month(date) < start_month, year(date) - 1, year(date)))
    
    if (!nrow(dat)) return(NULL)
    
    per_season <- dat |>
      group_by(hyear) |>
      summarise(
        n_days         = n(),
        n_swe_obs      = sum(is.finite(swe_obs)),
        frac_swe_obs   = n_swe_obs / n_days,
        max_na_gap_swe = {
          r <- rle(!is.finite(swe_obs))
          if (length(r$lengths)) max(ifelse(r$values, r$lengths, 0)) else 0
        },
        .groups = "drop"
      )
    
    tibble(
      name = nm,
      seasons_total       = nrow(per_season),
      seasons_ok          = sum(per_season$n_swe_obs >= min_days_per_season &
                                  per_season$max_na_gap_swe <= max_gap_days),
      median_frac_swe_obs = median(per_season$frac_swe_obs, na.rm = TRUE),
      worst_max_na_gap    = max(per_season$max_na_gap_swe, na.rm = TRUE)
    )
  })
  
  ans <- dplyr::bind_rows(out)
  if (is.null(ans)) {
    ans <- tibble(
      name = character(0),
      seasons_total = integer(0),
      seasons_ok = integer(0),
      median_frac_swe_obs = numeric(0),
      worst_max_na_gap = numeric(0)
    )
  }
  ans
}

diag_tbl <- coverage_by_station(d_obs, start_month = 8,
                                min_days_per_season = 120,
                                max_gap_days = 10)

# Save diagnostics
cov_dir <- "/Users/jakobwerkgarner/code/Master_Delta/calibration/calibration_logsx/coverage"
if (!dir.exists(cov_dir)) dir.create(cov_dir, recursive = TRUE, showWarnings = FALSE)
cov_csv <- file.path(cov_dir, "station_coverage.csv")
readr::write_csv(diag_tbl, cov_csv)

# ---------------- Station selection: strict → loose → top-N fallback --------
strict_sel <- diag_tbl |>
  filter(seasons_ok >= 2,
         median_frac_swe_obs >= 0.60,
         worst_max_na_gap <= 20) |>
  arrange(desc(median_frac_swe_obs)) |>
  pull(name)

sel_names <- strict_sel
if (length(sel_names) == 0) {
  message("Strict selection found 0 stations — relaxing thresholds.")
  loose_sel <- diag_tbl |>
    filter(seasons_ok >= 1,
           median_frac_swe_obs >= 0.50,
           worst_max_na_gap <= 30) |>
    arrange(desc(median_frac_swe_obs)) |>
    pull(name)
  sel_names <- loose_sel
}
if (length(sel_names) == 0) {
  message("Loose selection also found 0 stations — using top-N by coverage.")
  N_MIN <- max(1L, min(10L, nrow(diag_tbl)))
  sel_names <- diag_tbl |>
    arrange(desc(median_frac_swe_obs), desc(seasons_ok)) |>
    slice(1:N_MIN) |>
    pull(name)
}
d_obs <- d_obs[intersect(names(d_obs), sel_names)]
if (length(d_obs) == 0) stop("After selection, no stations remain. Check HS/SWE column names.")
message(sprintf("Using %d station(s): %s", length(d_obs), paste(names(d_obs), collapse = ", ")))

# --------------- Helpers: winters + splitting for fit/val -------------------
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
    
    # Optional edge guard on HS (only if a plausible HS column exists)
    hs_vec <- tryCatch({
      cd <- zoo::coredata(win)
      df <- if (is.matrix(cd)) as.data.frame(cd) else data.frame(value = cd, check.names = FALSE)
      nms <- tolower(names(df))
      h_candidates <- c("hobs","hs","h","hsobs","snowdepth","snow_depth","h_cm","hs_cm")
      h_idx <- which(nms %in% h_candidates)[1]
      if (!is.na(h_idx)) as.numeric(df[[h_idx]]) else rep(NA_real_, nrow(df))
    }, error = function(e) rep(NA_real_, nrow(win)))
    
    if (nrow(win) < 365 && any(is.finite(hs_vec))) {
      first_hs <- suppressWarnings(as.numeric(hs_vec[1]))
      last_hs  <- suppressWarnings(as.numeric(hs_vec[length(hs_vec)]))
      if (is.finite(first_hs) && first_hs > edge_thresh) next()
      if (is.finite(last_hs)  && last_hs  > edge_thresh) next()
    }
    res[[length(res) + 1]] <- win
  }
  res
}

build_splits <- function(d_obs,
                         split.type = "alt_winters",
                         ratio.fit = 0.5,
                         split.seed = 12,
                         season.start = "-08-01",
                         season.end   = "-07-31") {
  stopifnot(split.type %in% c("alt_winters","random","random_stations","chronological"))
  set.seed(split.seed)
  
  stations <- names(d_obs)
  if (length(stations) == 0) stop("build_splits: no stations provided.")
  
  d_obs_fit <- vector("list", length(stations)); names(d_obs_fit) <- stations
  d_obs_val <- vector("list", length(stations)); names(d_obs_val) <- stations
  
  fit_station_set <- NULL
  if (split.type == "random_stations") {
    n_fit_st <- max(1L, floor(length(stations) * ratio.fit))
    n_fit_st <- min(n_fit_st, length(stations))
    fit_station_set <- if (n_fit_st > 0) sample(stations, size = n_fit_st) else character(0)
  }
  
  for (n in stations) {
    d <- d_obs[[n]]
    d <- d[order(index(d)), ]  # date sort (no duplicate handling requested)
    winters <- get_valid_winters(d, season.start, season.end)
    if (length(winters) == 0) next()
    
    if (split.type == "alt_winters") {
      fit_list <- winters[which(seq_along(winters) %% 2 == 0)]
      val_list <- winters[which(seq_along(winters) %% 2 == 1)]
    } else if (split.type == "random_stations") {
      if (n %in% fit_station_set) { fit_list <- winters; val_list <- list()
      } else { fit_list <- list();  val_list <- winters }
    } else if (split.type == "random") {
      nW <- length(winters)
      n_fit  <- max(1L, floor(nW * ratio.fit)); n_fit <- min(n_fit, nW)
      idxFit <- if (n_fit > 0) sort(sample.int(nW, n_fit)) else integer(0)
      idxVal <- setdiff(seq_len(nW), idxFit)
      fit_list <- if (length(idxFit)) winters[idxFit] else list()
      val_list <- if (length(idxVal)) winters[idxVal] else list()
    } else if (split.type == "chronological") {
      nW <- length(winters)
      n_fit <- max(1L, floor(nW * ratio.fit)); n_fit <- min(n_fit, nW)
      fit_list <- winters[seq_len(n_fit)]
      val_list <- if (nW > n_fit) winters[(n_fit + 1):nW] else list()
    }
    
    if (length(fit_list)) d_obs_fit[[n]] <- do.call(rbind, fit_list)
    if (length(val_list)) d_obs_val[[n]] <- do.call(rbind, val_list)
  }
  
  d_obs_fit <- Filter(NROW, d_obs_fit)
  d_obs_val <- Filter(NROW, d_obs_val)
  if (!length(d_obs_fit) && !length(d_obs_val)) stop("build_splits: no valid winters after filtering.")
  
  list(fit = d_obs_fit, val = d_obs_val)
}

message(sprintf("Split mode: %s (seed=%d, ratio.fit=%.2f)", split.type, split.seed, ratio.fit))
splits    <- build_splits(d_obs, split.type, ratio.fit, split.seed, season.start, season.end)
d_obs_fit <- splits$fit
d_obs_val <- splits$val

# Optional: exclude specific station from fit
d_obs_fit[["Sta.Maria"]] <- NULL

# ---------------------- Build the fit tibble (robust) -----------------------
start_of_block <- 8  # August

build_fit_tibble <- function(d_obs_fit, start_of_block = 8) {
  out <- lapply(names(d_obs_fit), function(nm) {
    x <- d_obs_fit[[nm]]
    if (is.null(x) || NROW(x) == 0) return(NULL)
    
    cd <- .to_df(zoo::coredata(x))
    colmap <- .pick_obs_cols(cd, nm)
    if (is.null(colmap)) return(NULL)
    
    hs_raw  <- suppressWarnings(as.numeric(cd[[colmap$h_col]]))    # cm
    swe_raw <- suppressWarnings(as.numeric(cd[[colmap$swe_col]]))  # mm
    
    tibble(
      date    = as.Date(index(x)),
      name    = nm,
      hs      = hs_raw / 100,                 # cm -> m
      swe_obs = swe_raw
    ) |>
      mutate(block = set_season(date, start_of_block)) |>
      filter(is.finite(hs), is.finite(swe_obs)) |>
      mutate(
        hs      = pmax(pmin(hs, 15), 0),       # hard caps (keeps solver stable)
        swe_obs = pmax(pmin(swe_obs, 5000), 0)
      )
  })
  bind_rows(out)
}

d_obs_fit_tibble <- build_fit_tibble(d_obs_fit, start_of_block)
message(sprintf("Fit tibble rows: %d across %d stations.",
                nrow(d_obs_fit_tibble), dplyr::n_distinct(d_obs_fit_tibble$name)))
if (!"name" %in% names(d_obs_fit_tibble) || nrow(d_obs_fit_tibble) == 0) {
  stop("Fit tibble missing 'name' or is empty. Loosen selection thresholds or check column names.")
}

# -------------------------- Parameter setup ---------------------------------
par_delta <- c(
  rho.max   = 396.9521,
  rho.null  = 100,
  c.ov      = 0.0004986025,
  k.ov      = 0.227786,
  k.exp     = 0.0290256,  # positional as 'k'
  tau       = 0.02489556,
  eta.null  = 8792253
)
par_scale <- c(1000, 1000, 0.001, 1, 0.1, 0.1, 1e7)
par_delta <- par_delta / par_scale

lower <- c(300, 75, 0, 0.01, 0.01, 0.01, 1e6) / par_scale
upper <- c(600, 200, 0.001, 10, 0.2, 0.2, 2e7) / par_scale

lower_unscaled <- c(300, 75, 0, 0.01, 0.01, 0.01, 1e6)
upper_unscaled <- c(600, 200, 0.001, 10, 0.2, 0.2, 2e7)

# --------------------- Parallel + reproducible RNG --------------------------
nc <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
cl <- parallel::makeCluster(nc, type = "PSOCK")
doParallel::registerDoParallel(cl)
RNGkind("L'Ecuyer-CMRG")
set.seed(split.seed)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(zoo); library(dplyr); library(nixmass) })
  NULL
})

# ---------------- Objective: robust RMSE (finite; compare on obs only) ------
minimize_score <- function(data, par, scale, verbose = FALSE) {
  if (is.null(data) || !is.data.frame(data) || !"name" %in% names(data) || nrow(data) == 0) {
    return(1e11)  # penalty
  }
  stations <- unique(data$name)
  if (length(stations) == 0) return(1e11)
  
  par_unscaled <- pmin(pmax(as.numeric(par) * scale, lower_unscaled), upper_unscaled)
  
  eps_interior <- c(10, 5, 1e-5, 0.02, 0.01, 0.01, 5e5)
  edge_dist <- pmin(par_unscaled - lower_unscaled, upper_unscaled - par_unscaled)
  soft_pen  <- sum((pmax(0, eps_interior - edge_dist) / eps_interior)^2) * 1e4
  
  if (any(!is.finite(par_unscaled))) return(1e12)
  
  `%doRNG%` <- doRNG::`%dorng%`
  
  safe_eval_one <- function(joined) {
    swe_mod <- tryCatch({
      R.utils::withTimeout({
        out <- joined |>
          dplyr::select(date, hs) |>
          dplyr::mutate(date = as.character(date)) |>
          nixmass::swe.delta.snow(
            model_opts = list(
              rho.max   = par_unscaled[1],
              rho.null  = par_unscaled[2],
              c.ov      = par_unscaled[3],
              k.ov      = par_unscaled[4],
              k         = par_unscaled[5],
              tau       = par_unscaled[6],
              eta.null  = par_unscaled[7]
            ),
            dyn_rho_max = FALSE
          )
        as.numeric(out)
      }, timeout = 15)
    }, TimeoutException = function(e) rep(NA_real_, nrow(joined)),
    error            = function(e) rep(NA_real_, nrow(joined)))
    swe_mod
  }
  
  ll <- foreach(
    s = stations,
    .packages = c("dplyr", "zoo", "nixmass")
  ) %doRNG% {
    data1 <- dplyr::filter(data, name == s)
    years <- unique(data1$block)
    l <- lapply(years, function(y) {
      left <- data1 |>
        dplyr::filter(block == y) |>
        dplyr::select(date, hs, swe_obs) |>
        dplyr::arrange(date)
      
      if (nrow(left) < 20) return(NULL)
      
      backbone <- tibble::tibble(date = seq(min(left$date), max(left$date), by = "1 day"))
      joined <- backbone |>
        dplyr::left_join(left, by = "date") |>
        dplyr::arrange(date) |>
        dplyr::mutate(
          hs = zoo::na.approx(hs, x = date, na.rm = FALSE, maxgap = 7)  # small gaps only
        )
      
      swe_mod <- safe_eval_one(joined)
      
      tibble::tibble(date = joined$date, swe_obs = joined$swe_obs, swe_mod = swe_mod) |>
        dplyr::filter(is.finite(swe_obs), is.finite(swe_mod))
    })
    dplyr::bind_rows(l)
  }
  
  dff <- dplyr::bind_rows(ll)
  if (nrow(dff) < 40) return(1e11 + soft_pen)
  
  rmse <- sqrt(mean((dff$swe_mod - dff$swe_obs)^2))
  if (!is.finite(rmse)) rmse <- 1e10
  rmse + soft_pen
}

# ---------------- Metrics helper (same evaluation policy) -------------------
eval_metrics_for_par <- function(data, par_unscaled) {
  if (is.null(data) || !is.data.frame(data) || !"name" %in% names(data) || nrow(data) == 0) {
    return(list(rmse = NA_real_, bias = NA_real_))
  }
  stations <- unique(data$name)
  if (length(stations) == 0) return(list(rmse = NA_real_, bias = NA_real_))
  
  par_unscaled <- pmin(pmax(as.numeric(par_unscaled), lower_unscaled), upper_unscaled)
  `%doRNG%` <- doRNG::`%dorng%`
  
  ll <- foreach(
    s = stations,
    .packages = c("dplyr", "zoo", "nixmass")
  ) %doRNG% {
    data1 <- dplyr::filter(data, name == s)
    l <- lapply(unique(data1$block), function(y) {
      left <- dplyr::filter(data1, block == y) |>
        dplyr::select(date, hs, swe_obs) |>
        dplyr::arrange(date)
      if (nrow(left) < 20) return(NULL)
      
      backbone <- tibble::tibble(date = seq(min(left$date), max(left$date), by = "1 day"))
      joined <- backbone |> dplyr::left_join(left, by = "date") |> dplyr::arrange(date)
      
      swe_mod <- tryCatch({
        joined |>
          dplyr::select(date, hs) |>
          dplyr::mutate(date = as.character(date)) |>
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
      
      tibble::tibble(date = joined$date,
                     swe_obs = joined$swe_obs,
                     swe_mod = as.numeric(swe_mod)) |>
        dplyr::filter(is.finite(swe_obs), is.finite(swe_mod))
    })
    dplyr::bind_rows(l)
  }
  
  dff <- dplyr::bind_rows(ll)
  rmse <- if (nrow(dff)) sqrt(mean((dff$swe_mod - dff$swe_obs)^2)) else NA_real_
  bias <- if (nrow(dff)) mean(dff$swe_mod - dff$swe_obs)          else NA_real_
  list(rmse = rmse, bias = bias)
}

# --------------------------------- Optimize --------------------------------
opt <- optimx::optimx(
  fn      = minimize_score,
  data    = d_obs_fit_tibble,
  par     = par_delta,
  scale   = par_scale,
  verbose = FALSE,
  method  = c("L-BFGS-B", "bobyqa"),
  control = list(trace = 2, follow.on = TRUE, kkt = FALSE, starttests = FALSE),
  lower   = lower,
  upper   = upper
)
saveRDS(opt, file = "opt_results_orig.rds")

# --------------------- Extract best + compute final metrics -----------------
opt_df    <- as.data.frame(opt)
par_names <- names(par_delta)
best_idx  <- which.min(opt_df$value)

best_par_scaled   <- as.numeric(opt_df[best_idx, par_names, drop = TRUE]); names(best_par_scaled) <- par_names
best_par_unscaled <- best_par_scaled * par_scale
best_val          <- opt_df$value[best_idx]
method_str        <- if ("method" %in% names(opt_df)) as.character(opt_df$method[best_idx]) else "n/a"

# Best-so-far checkpoint
saveRDS(list(
  timestamp = Sys.time(),
  par_scaled = best_par_scaled,
  par_unscaled = best_par_unscaled,
  value = best_val,
  method = method_str
), "opt_checkpoint_best.rds")

metrics    <- eval_metrics_for_par(d_obs_fit_tibble, best_par_unscaled)
final_rmse <- as.numeric(metrics$rmse)
final_bias <- as.numeric(metrics$bias)

# --------------------------------- Logging ---------------------------------
wd <- getwd()
fit_stations <- paste(unique(d_obs_fit_tibble$name), collapse = ", ")

pkg_list <- c("optimx","zoo","foreach","doParallel","doRNG","lubridate","nixmass","tidyverse","R.utils")
pkg_versions <- sapply(pkg_list, function(p) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NA")
  paste0(p, "=", v)
})

os_info <- tryCatch(paste(unname(Sys.info()), collapse = " | "), error = function(e) R.version$platform)
R_str   <- R.version.string

base_log_dir <- "/Users/jakobwerkgarner/code/Master_Delta/calibration/calibration_logsx/A_log_files/"
if (!dir.exists(base_log_dir)) dir.create(base_log_dir, recursive = TRUE, showWarnings = FALSE)

ts_colon     <- format(Sys.time(), "%Y_%m_%d_%H%M")
safe_comment <- gsub("[^A-Za-z0-9_\\-]+", "_", calib_comment)
fname        <- paste0(ts_colon, "_calib_log_", safe_comment, ".txt")
log_file     <- file.path(base_log_dir, fname)

lines <- c(
  "===== Calibration Log =====",
  paste0("Timestamp: ", ts_colon),
  paste0("Working directory: ", wd),
  paste0("RDA path: ", rda_path),
  paste0("Results RDS path: ", file.path(getwd(), "opt_results_orig.rds")),
  paste0("Coverage CSV: ", cov_csv),
  "",
  "---- Data / Stations ----",
  paste0("Stations considered (after base filtering): ", length(names(d_obs))),
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
  paste0("Objective (RMSE+penalty) = ", signif(best_val, 6)),
  paste0("Final RMSE (recomputed)  = ", signif(final_rmse, 6)),
  paste0("Final bias (mean error)  = ", signif(final_bias, 6)),
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

# Optional: go sequential again
foreach::registerDoSEQ()
# ---------------------------------- END ------------------------------------