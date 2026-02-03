#-------------------------------------------------------------------------
# parameter optimization for deltasnow model
# optimize deltasnow model with respect to a reasonable score
# of observed and modeled SWE values
# score could be e.g. minimum residual sum of squares or RMSE...
#
# Harald Schellander, 08.2025
#-------------------------------------------------------------------------

library(optimx)
library(zoo)
library(foreach)
library(doParallel)
library(lubridate)
library(nixmass)
library(tidyverse)

setwd("/Users/jakobwerkgarner/code/Master_Delta/calibration")

# ========================= USER CONFIG (added) =========================
# comment that will be embedded in the log file name
calib_comment <- "hs2swe_calib_Win21"  # <- set your comment string

# -------------------------- NEW: Split configuration --------------------------
# 1) alt_winters     : alternating winters (even index -> fit, odd -> val)  [original]
# 2) random_50_50    : random split of seasons within each station (uses ratio.fit)
# 3) random_stations : random split by stations (all winters of a station kept together; uses ratio.fit)
# 4) chronological   : earliest X% winters -> fit, remaining -> val (uses ratio.fit)
split.type  <- "alt_winters"        # <- "alt_winters" | "random_50_50" | "random_stations" | "chronological"
split.seed  <- 12                   # <- seed for the random-based modes
ratio.fit   <- 0.5                  # <- proportion for fit set (used by random_50_50, random_stations, chronological)
# ======================================================================


#---------------------------------------------------------------------------------------------------
# real obs: HD are weekly, Swiss data are bi-weekly
# split data into two halfs for fitting and verification
# use every 2nd winter for fitting, every other for validation
# create a tibble

rda_path <- "/Users/jakobwerkgarner/code/Master_Delta/Winkler21_Data/dSnow_calib_data/H_SWE_obs.Rda"
d_obs <- get(load(rda_path))




d_obs[["kuehtai"]] <- NULL # use Kühtai and Weissfluhjoch for analysis purposes
d_obs[["Weissfluhjoch"]] <- NULL
season.start <- "-08-01"
season.end <- "-07-31"
d_obs_fit <- d_obs_val <- list()
for (n in names(d_obs)) {
  print(n)
  d <- d_obs[[n]]
  years <- unique(year(index(d)))
  winters4fit <- winters4val <- zoo()
  idx_y <- 1
  for (y in years[1:(length(years) - 1)]) {
    winter <- subset(
      d,
      index(d) >= as.Date(paste0(y, season.start)) &
        index(d) < as.Date(paste0(y + 1, season.end)) + 1
    )
    if (nrow(winter) < 200) {
      print(paste0("only ", nrow(winter), " values...skipping year ", y))
      next()
    } else {
      if (
        nrow(winter) < 365 &
        (as.numeric(winter$Hobs[1]) > 0.05 |
         as.numeric(winter$Hobs[nrow(winter)]) > 0.05)
      ) {
        print(paste0(
          "only ",
          nrow(winter),
          " values and last/first value > 0.05...skipping year ",
          y
        ))
        next()
      }
      if (idx_y %% 2 == 0) {
        if (is.null(nrow(winters4fit))) {
          winters4fit <- winter
        } else {
          winters4fit <- rbind(winters4fit, winter) # even years for fitting
        }
      } else {
        if (is.null(nrow(winters4val))) {
          winters4val <- winter
        } else {
          winters4val <- rbind(winters4val, winter)
        }
      }
    }
    d_obs_fit[[n]] <- winters4fit
    d_obs_val[[n]] <- winters4val
    idx_y <- idx_y + 1
  }
}

# =================== NEW: Enhanced split modes (added) ===================
# Rebuild d_obs_fit / d_obs_val when split.type != "alt_winters"
if (!split.type %in% c("alt_winters", "random_50_50", "random_stations", "chronological")) {
  stop("Unknown split.type specified.")
}
if (split.type != "alt_winters") {
  message(sprintf("Applying enhanced split mode: %s (seed=%d, ratio.fit=%.2f)", split.type, split.seed, ratio.fit))
  set.seed(split.seed)
  
  # helper: list valid winters per station (zoo objects)
  get_valid_winters <- function(d, years, season.start, season.end) {
    res <- list()
    for (y in years[1:(length(years) - 1)]) {
      winter <- subset(
        d,
        index(d) >= as.Date(paste0(y, season.start)) &
          index(d) < as.Date(paste0(y + 1, season.end)) + 1
      )
      if (nrow(winter) < 200) next()
      if (nrow(winter) < 365 &&
          (as.numeric(winter$Hobs[1]) > 0.05 ||
           as.numeric(winter$Hobs[nrow(winter)]) > 0.05)) next()
      res[[length(res) + 1]] <- winter
    }
    res
  }
  
  # prepare station list
  stations <- names(d_obs)
  
  # random_stations: decide station membership first
  fit_station_set <- NULL
  if (split.type == "random_stations") {
    n_fit_st <- max(1, floor(length(stations) * ratio.fit))
    fit_station_set <- sample(stations, size = n_fit_st)
  }
  
  d_obs_fit <- list()
  d_obs_val <- list()
  
  for (n in stations) {
    d <- d_obs[[n]]
    years <- unique(year(index(d)))
    winters <- get_valid_winters(d, years, season.start, season.end)
    if (length(winters) == 0) next()
    
    if (split.type == "random_stations") {
      # all winters of a station to either fit or val
      if (n %in% fit_station_set) {
        d_obs_fit[[n]] <- do.call(rbind, winters)
      } else {
        d_obs_val[[n]] <- do.call(rbind, winters)
      }
      
    } else if (split.type == "random_50_50") {
      nW <- length(winters)
      n_fit <- max(1, floor(nW * ratio.fit))
      idx_fit <- sort(sample.int(nW, n_fit))
      idx_val <- setdiff(seq_len(nW), idx_fit)
      if (length(idx_fit)) d_obs_fit[[n]] <- do.call(rbind, winters[idx_fit])
      if (length(idx_val)) d_obs_val[[n]] <- do.call(rbind, winters[idx_val])
      
    } else if (split.type == "chronological") {
      nW <- length(winters)
      n_fit <- max(1, floor(nW * ratio.fit))
      fit_list <- winters[seq_len(min(n_fit, nW))]
      val_list <- if (nW > n_fit) winters[(n_fit + 1):nW] else list()
      if (length(fit_list)) d_obs_fit[[n]] <- do.call(rbind, fit_list)
      if (length(val_list)) d_obs_val[[n]] <- do.call(rbind, val_list)
    }
  }
}
# =======================================================================


#---------------------------------------------------------------------------------------------------
# prepare calibration data

# handling yearly series is more efficient
start_of_block <- 8

# x...character date
# m...integer start month of block (1-12)
# return season as year
set_season <- function(x, m) {
  x <- as.character(x)
  yr <- as.POSIXlt(x)$year + 1900
  mt <- as.POSIXlt(x)$mon + 1
  # for different block sizes adjust
  # e.g. yr - 1 to months - months_block_size
  ifelse(mt < m, yr - 1, yr)
}

d_obs_fit[["Sta.Maria"]] <- NULL # use Sta. Maria only for analysis/validation purposes
d_obs_fit_tibble <- lapply(seq_along(d_obs_fit), function(i) {
  x <- d_obs_fit[[i]]
  tibble(
    date = index(x),
    name = names(d_obs_fit)[i],
    as_tibble(coredata(x))[, c("Hobs", "SWEobs")]
  ) |>
    rename(hs = "Hobs", swe_obs = "SWEobs") |>
    mutate(hs = hs / 100, block = set_season(date, start_of_block))
})
d_obs_fit_tibble <- do.call(rbind, d_obs_fit_tibble)


#---------------------------------------------------------------------------------------------------
# function for score to be minimized
minimize_score <- function(data, par, scale, verbose = FALSE) {
  par <- par * scale
  cat(par) # this is handy to follow the optimization process
  
  # parallelized
  ll <- foreach(
    s = unique(data$name),
    .packages = c("dplyr", "zoo", "nixmass", "foreach")
  ) %dopar%
    {
      if (verbose) {
        cat(paste0(s, " ..."))
      }
      data1 <- data |>
        filter(name == s)
      # serial, no gain when parallelized
      l <- foreach(
        i = 1:length(unique(data1$block)),
        .packages = c("dplyr", "zoo", "nixmass")
      ) %do%
        {
          y <- unique(data1$block)[i]
          left <- data1 |>
            filter(block == y) |>
            dplyr::select(date, hs, swe_obs)
          right <- data.frame(
            date = seq(min(left$date), max(left$date), by = "1 day")
          )
          joined <- left |>
            right_join(right, by = "date")
          
          # catch possible problems with model results
          swe_mod <- tryCatch(
            {
              joined |>
                dplyr::select(date, hs) |>
                mutate(date = as.character(date)) |>
                nixmass::hs2(
                  model_opts = list(
                    rho.max = par[1],
                    rho.null = par[2],
                    c.ov = par[3],
                    k.ov = par[4],
                    k = par[5],
                    tau = par[6],
                    eta.null = par[7]
                  ),
                  dyn_rho_max = FALSE
                )
            },
            error = function(e) {
              return(rep(NA, nrow(joined)))
            }
          )
          na.omit(cbind(joined, swe_mod)) # remove all values without swe observation
        }
      df <- do.call(rbind, l)
      df
    }
  dff <- do.call(rbind, ll)
  rmse <- with(dff, sqrt(mean((swe_mod - swe_obs)^2)))
  bias <- with(dff, abs(mean(swe_mod - swe_obs)))
  cat(paste0(" |bias|=", bias, " rmse=", rmse, "\n"))
  rmse
}

# =================== NEW: metrics helper for logging (added) ===================
eval_metrics_for_par <- function(data, par_unscaled) {
  # same computation as in minimize_score but returns both metrics without cat()
  ll <- foreach(
    s = unique(data$name),
    .packages = c("dplyr", "zoo", "nixmass", "foreach")
  ) %dopar%
    {
      data1 <- data |> filter(name == s)
      l <- foreach(
        i = 1:length(unique(data1$block)),
        .packages = c("dplyr", "zoo", "nixmass")
      ) %do%
        {
          y <- unique(data1$block)[i]
          left <- data1 |>
            filter(block == y) |>
            dplyr::select(date, hs, swe_obs)
          right <- data.frame(date = seq(min(left$date), max(left$date), by = "1 day"))
          joined <- left |> right_join(right, by = "date")
          swe_mod <- tryCatch(
            {
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
            },
            error = function(e) rep(NA, nrow(joined))
          )
          na.omit(cbind(joined, swe_mod))
        }
      do.call(rbind, l)
    }
  dff <- do.call(rbind, ll)
  rmse <- with(dff, sqrt(mean((swe_mod - swe_obs)^2)))
  bias <- with(dff, mean(swe_mod - swe_obs))
  list(rmse = rmse, bias = bias)
}
# ==============================================================================


#-------------------------------------------------------------------------
# start values
# taken from a quick optimization
par_delta <- c(
  rho.max = 396.9521,
  rho.null = 78.88324,
  c.ov = 0.0004986025,
  k.ov = 0.227786,
  k.exp = 0.0290256,
  tau = 0.02489556,
  eta.null = 8792253
)
par_scale <- c(1000, 1000, 0.001, 1, 0.1, 0.1, 1e7)
par_delta <- par_delta / par_scale
lower <- c(300, 50, 0, 0.01, 0.01, 0.01, 1e6) / par_scale
upper <- c(600, 200, 0.001, 10, 0.2, 0.2, 2e7) / par_scale
cbind(
  lower,
  par_delta,
  upper,
  ifelse(lower <= par_delta & upper >= par_delta, "in bounds", "ERROR")
)


#-------------------------------------------------------------------------
# optimization
nc <- detectCores(logical = TRUE)
cl <- makeCluster(nc)
registerDoParallel(cl)
opt <- optimx(
  fn = minimize_score,
  data = d_obs_fit_tibble,
  par = par_delta,
  scale = par_scale,
  verbose = FALSE,
  #method = c('L-BFGS-B', 'bobyqa'),
  control = list(trace = 4, follow.on = TRUE)
  #lower = lower,
  #upper = upper
)
saveRDS(opt, file = "opt_results_orig.rds")

# ========================= STRICT: Logging block (no connections) =========================
try({
  opt_df <- as.data.frame(opt)
  par_names <- names(par_delta)
  best_idx <- which.min(opt_df$value)
  best_par_scaled <- as.numeric(opt_df[best_idx, par_names, drop = TRUE]); names(best_par_scaled) <- par_names
  best_par_unscaled <- best_par_scaled * par_scale

  # recompute metrics on fit data
  metrics <- eval_metrics_for_par(d_obs_fit_tibble, best_par_unscaled)
  final_rmse <- as.numeric(metrics$rmse)
  final_bias <- as.numeric(metrics$bias)

  wd <- getwd()
  all_stations <- paste(names(d_obs), collapse = ", ")
  fit_stations <- paste(unique(d_obs_fit_tibble$name), collapse = ", ")

  # env info
  pkg_list <- c("optimx","zoo","foreach","doParallel","lubridate","nixmass","tidyverse")
  pkg_versions <- sapply(pkg_list, function(p) {
    v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NA")
    paste0(p, "=", v)
  })
  os_info <- tryCatch(paste(unname(Sys.info()), collapse = " | "), error = function(e) R.version$platform)
  R_str <- R.version.string

  # paths
  rds_path <- tryCatch(
    normalizePath("opt_results_orig.rds", winslash = "/", mustWork = FALSE),
    error = function(e) file.path(getwd(), "opt_results_orig.rds")
  )
  rda_path_str <- tryCatch(normalizePath(rda_path, winslash = "/", mustWork = FALSE), error = function(e) rda_path)

  # EXACT directory
  base_log_dir <- "/Users/jakobwerkgarner/code/Master_Delta/calibration/calibration_logsx/A_log_files"
  if (!dir.exists(base_log_dir)) dir.create(base_log_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(base_log_dir)) stop("Target log directory does not exist: ", base_log_dir)
  if (file.access(base_log_dir, 2) != 0) stop("No write permission to: ", base_log_dir)

  # filename: YYYY_DD_HH:MM (requested); fallback to HH-MM if ':' fails
  ts_colon <- format(Sys.time(), "%Y_%d_%H:%M")
  safe_comment <- gsub("[^A-Za-z0-9_\\-]+", "_", calib_comment)
  fname <- paste0(ts_colon, "_calib_log_", safe_comment, ".txt")
  log_file <- file.path(base_log_dir, fname)

  # method may be absent depending on optimx return
  method_str <- if ("method" %in% names(opt_df)) as.character(opt_df$method[best_idx]) else "n/a"

  # compose log
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
    paste0("Objective (RMSE reported by fn) = ", signif(opt_df$value[best_idx], 6)),
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

  # write with cat() (no explicit connections)
  successful <- TRUE
  err <- try(cat(paste(lines, collapse = "\n"), file = log_file, sep = "\n"), silent = TRUE)
  if (inherits(err, "try-error")) {
    # fallback without colon in filename (some FS/locale combos dislike ':')
    fname2 <- gsub(":", "-", fname, fixed = TRUE)
    log_file <- file.path(base_log_dir, fname2)
    cat(paste(lines, collapse = "\n"), file = log_file, sep = "\n")
  }

  message(sprintf("Calibration log written to: %s", log_file))
}, silent = FALSE)
# ================================================================================================