# -----------------------------------------------------------------------------
# DELTASNOW PARAMETER OPTIMIZATION (Win21) - DIFFERENTIAL EVOLUTION
# -----------------------------------------------------------------------------
#
# Based on the first script's data (H_SWE_obs.Rda, excluding Kühtai,
# Weissfluhjoch, and Sta. Maria) but using the weighted objective function,
# DIFFERENTIAL EVOLUTION optimizer, and output formatting of the second script.
#
# Objective:
#   final_score = w1 * NRMSE_SWE + w2 * NRMSE_rho + w3 * NBIAS_SWE
#
# with NRMSE = RMSE / mean(obs) and NBIAS = |bias| / mean(obs).
# Weights are user‑adjustable below.
#
# -----------------------------------------------------------------------------

library(optimx)
library(zoo)
library(foreach)
library(doParallel)
library(lubridate)
library(nixmass)
library(tidyverse)
library(DEoptim)              # <-- NEW: Differential Evolution package

# ----------------------------------------------------------------------------
# CONFIGURATION (USER-ADJUSTABLE)
# ----------------------------------------------------------------------------

# Weights for the combined objective (should sum to 1)
WEIGHT_SWE_NRMSE  <- 0.3    # weight for normalized SWE RMSE
WEIGHT_RHO_NRMSE  <- 0.7    # weight for normalized density RMSE
WEIGHT_SWE_NBIAS  <- 0.0    # weight for normalized SWE absolute bias

# Override weights from command-line: Rscript script.R <w_swe> <w_rho> <w_bias>
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) == 3) {
  WEIGHT_SWE_NRMSE <- as.numeric(.args[1])
  WEIGHT_RHO_NRMSE <- as.numeric(.args[2])
  WEIGHT_SWE_NBIAS <- as.numeric(.args[3])
}
cat(sprintf("Weights: SWE_NRMSE=%.3f  RHO_NRMSE=%.3f  SWE_NBIAS=%.3f\n",
            WEIGHT_SWE_NRMSE, WEIGHT_RHO_NRMSE, WEIGHT_SWE_NBIAS))

# Data settings
season_start   <- "-08-01"       # start of hydrological year (month-day)
season_end     <- "-07-31"       # end of hydrological year
start_of_block <- 8              # month number for block assignment

# Model settings
EPS <- 1e-6                      # threshold for snow depth in density calculation

# ----------------------------------------------------------------------------
# LOAD OBSERVATIONAL DATA (from first script)
# ----------------------------------------------------------------------------
d_obs <- get(load("/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_Win21/data/H_SWE_obs.Rda"))

# Exclude stations as in the first script
d_obs[["kuehtai"]] <- NULL
d_obs[["Weissfluhjoch"]] <- NULL
# Sta. Maria will be removed later from the fit set (only used for validation)

# ----------------------------------------------------------------------------
# HELPER FUNCTIONS (from second script)
# ----------------------------------------------------------------------------

#' Assign a hydrological season (year) to a date
set_season <- function(date, start_month = 8) {
  date <- as.Date(date)
  yr <- year(date)
  mo <- month(date)
  ifelse(mo < start_month, yr - 1, yr)
}

#' Root Mean Square Error
rmse <- function(obs, mod) {
  ok <- is.finite(obs) & is.finite(mod)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((mod[ok] - obs[ok])^2))
}

#' Normalized RMSE
nrmse <- function(obs, mod) {
  ok <- is.finite(obs) & is.finite(mod)
  if (!any(ok)) return(NA_real_)
  rmse_val <- sqrt(mean((mod[ok] - obs[ok])^2))
  mean_obs <- mean(obs[ok], na.rm = TRUE)
  rmse_val / mean_obs
}

#' Normalized Absolute Bias
nbias <- function(obs, mod) {
  ok <- is.finite(obs) & is.finite(mod)
  if (!any(ok)) return(NA_real_)
  bias_val <- mean(mod[ok] - obs[ok], na.rm = TRUE)
  mean_obs <- mean(obs[ok], na.rm = TRUE)
  abs(bias_val) / mean_obs
}

#' Combined score: weighted sum of NRMSE_SWE, NRMSE_rho, and NBIAS_SWE
combined_score <- function(df, eps = EPS) {
  swe_obs <- df$swe_obs
  swe_mod <- df$swe_mod
  hs <- df$hs

  # SWE metrics
  nrmse_swe <- nrmse(swe_obs, swe_mod)
  nbias_swe <- nbias(swe_obs, swe_mod)

  # Bulk density metrics (only where snow depth > eps)
  rho_obs <- ifelse(is.finite(hs) & hs > eps, swe_obs / hs, NA_real_)
  rho_mod <- ifelse(is.finite(hs) & hs > eps, swe_mod / hs, NA_real_)
  nrmse_rho <- nrmse(rho_obs, rho_mod)

  # Weighted combination
  score <- WEIGHT_SWE_NRMSE * nrmse_swe +
           WEIGHT_RHO_NRMSE * nrmse_rho +
           WEIGHT_SWE_NBIAS * nbias_swe

  # Attach detailed metrics for verbose printing
  attr(score, "metrics") <- list(
    rmse_swe  = rmse(swe_obs, swe_mod),
    rmse_rho  = rmse(rho_obs, rho_mod),
    nrmse_swe = nrmse_swe,
    nrmse_rho = nrmse_rho,
    bias_swe  = mean(swe_mod - swe_obs, na.rm = TRUE),
    nbias_swe = nbias_swe,
    mean_swe  = mean(swe_obs[is.finite(swe_obs) & is.finite(swe_mod)], na.rm = TRUE),
    mean_rho  = mean(rho_obs[is.finite(rho_obs) & is.finite(rho_mod)], na.rm = TRUE)
  )

  return(score)
}

# ----------------------------------------------------------------------------
# SPLIT DATA INTO FIT AND VALIDATION (even/odd years per station)
# ----------------------------------------------------------------------------
d_obs_fit <- list()
d_obs_val <- list()

for (station in names(d_obs)) {
  cat("Processing station:", station, "\n")

  d <- d_obs[[station]]
  years <- unique(year(index(d)))

  fit_list <- list()
  val_list <- list()
  counter <- 1

  for (y in years[1:(length(years) - 1)]) {
    winter <- subset(
      d,
      index(d) >= as.Date(paste0(y, season_start)) &
        index(d) < as.Date(paste0(y + 1, season_end)) + 1
    )

    if (nrow(winter) < 200) {
      cat("  skipping year", y, "- only", nrow(winter), "values\n")
      next
    }

    if (
      nrow(winter) < 365 &&
      (as.numeric(winter$Hobs[1]) > 0.05 ||
       as.numeric(winter$Hobs[nrow(winter)]) > 0.05)
    ) {
      cat("  skipping year", y, "- incomplete winter with snow at edge\n")
      next
    }

    if (counter %% 2 == 0) {
      fit_list[[length(fit_list) + 1]] <- winter
    } else {
      val_list[[length(val_list) + 1]] <- winter
    }
    counter <- counter + 1
  }

  d_obs_fit[[station]] <- if (length(fit_list) > 0) do.call(rbind, fit_list) else NULL
  d_obs_val[[station]] <- if (length(val_list) > 0) do.call(rbind, val_list) else NULL
}

# ----------------------------------------------------------------------------
# PREPARE FIT DATA (convert zoo to tibble, remove Sta. Maria)
# ----------------------------------------------------------------------------
fit_data_list <- list()

for (station in names(d_obs_fit)) {
  # Skip Sta. Maria as in the first script
  if (station == "Sta.Maria") next

  x <- d_obs_fit[[station]]
  if (is.null(x) || length(x) == 0) next

  dat <- tryCatch(as_tibble(coredata(x)), error = function(e) NULL)
  if (is.null(dat) || ncol(dat) == 0) next

  if (ncol(dat) == 2 && all(is.na(colnames(dat)))) {
    colnames(dat) <- c("Hobs", "SWEobs")
  }

  if (!all(c("Hobs", "SWEobs") %in% colnames(dat))) {
    cat("Skipping station", station, "- missing Hobs or SWEobs\n")
    next
  }

  fit_data_list[[station]] <- tibble(
    date    = as.Date(index(x)),
    name    = station,
    hs      = dat$Hobs / 100,        # cm -> m (as in first script)
    swe_obs = dat$SWEobs,
    block   = set_season(index(x), start_of_block)
  )
}

d_obs_fit_tibble <- bind_rows(fit_data_list)

# ----------------------------------------------------------------------------
# PREPARE VAL DATA (convert zoo to tibble)
# ----------------------------------------------------------------------------
val_data_list <- list()

for (station in names(d_obs_val)) {
  if (station == "Sta.Maria") next

  x <- d_obs_val[[station]]
  if (is.null(x) || length(x) == 0) next

  dat <- tryCatch(as_tibble(coredata(x)), error = function(e) NULL)
  if (is.null(dat) || ncol(dat) == 0) next

  if (ncol(dat) == 2 && all(is.na(colnames(dat)))) {
    colnames(dat) <- c("Hobs", "SWEobs")
  }

  if (!all(c("Hobs", "SWEobs") %in% colnames(dat))) {
    cat("Skipping val station", station, "- missing Hobs or SWEobs\n")
    next
  }

  val_data_list[[station]] <- tibble(
    date    = as.Date(index(x)),
    name    = station,
    hs      = dat$Hobs / 100,
    swe_obs = dat$SWEobs,
    block   = set_season(index(x), start_of_block)
  )
}

d_obs_val_tibble <- bind_rows(val_data_list)

# ----------------------------------------------------------------------------
# OBJECTIVE FUNCTION (called by DEoptim)
# ----------------------------------------------------------------------------
minimize_score <- function(par, data, scale, verbose = FALSE) {

  result <- tryCatch({

    par_real <- par * scale
    if (verbose) {
      param_names <- names(par_real)
      if (is.null(param_names)) {
        param_names <- paste0("p", seq_along(par_real))
      }
      scaled_values <- paste(
        paste0(param_names, "=", round(par, 6)),
        collapse = ", "
      )
      unscaled_values <- paste(
        paste0(param_names, "=", round(par_real, 6)),
        collapse = ", "
      )
      cat("scaled   =", scaled_values, "\n")
      cat("unscaled =", unscaled_values, "\n")
    }

    # Run model for all stations and blocks in parallel
    station_results <- foreach(
      station = unique(data$name),
      .packages = c("dplyr", "tidyr", "nixmass"),
      .combine = bind_rows
    ) %dopar% {

      data_station <- data %>% filter(name == station)
      blocks <- unique(data_station$block)

      block_results <- lapply(blocks, function(b) {
        df_block <- data_station %>% filter(block == b)

        full_dates <- tibble(
          date = seq(min(df_block$date), max(df_block$date), by = "1 day")
        )

        joined <- df_block %>%
          select(date, hs, swe_obs) %>%
          right_join(full_dates, by = "date") %>%
          arrange(date)

        swe_mod <- tryCatch(
          {
            joined %>%
              select(date, hs) %>%
              mutate(date = as.character(date)) %>%
              nixmass::swe.delta.snow(
                model_opts = list(
                  rho.max  = par_real[1],
                  rho.null = par_real[2],
                  c.ov     = par_real[3],
                  k.ov     = par_real[4],
                  k        = par_real[5],
                  tau      = par_real[6],
                  eta.null = par_real[7]
                ),
                dyn_rho_max = FALSE
              )
          },
          error = function(e) rep(NA_real_, nrow(joined))
        )

        joined %>%
          mutate(swe_mod = swe_mod) %>%
          drop_na()
      })

      bind_rows(block_results)
    }

    dff <- station_results

    if (nrow(dff) == 0) {
      cat("No valid model output. Returning large penalty.\n")
      return(1e12)
    }

    score_with_attr <- combined_score(dff)
    metrics <- attr(score_with_attr, "metrics")
    score <- as.numeric(score_with_attr)

    if (!is.finite(score)) {
      cat("Score is not finite. Returning large penalty.\n")
      return(1e12)
    }

    if (verbose && !is.null(metrics)) {
      cat(
        "| bias_swe =", round(metrics$bias_swe, 4),
        "| RMSE_SWE =", round(metrics$rmse_swe, 4),
        "| NRMSE_SWE =", round(metrics$nrmse_swe, 4),
        "| NBIAS_SWE =", round(metrics$nbias_swe, 4),
        "| RMSE_rho =", round(metrics$rmse_rho, 4),
        "| NRMSE_rho =", round(metrics$nrmse_rho, 4),
        "| final_score =", round(score, 4),
        "\n"
      )
    } else {
      cat("Score =", round(score, 4), "\n")
    }

    return(score)

  }, error = function(e) {
    cat("Error in objective function:", e$message, "\n")
    return(1e12)
  })

  if (length(result) != 1 || !is.numeric(result)) {
    return(1e12)
  }
  return(result)
}

# ----------------------------------------------------------------------------
# START VALUES AND SCALING (from first script)
# ----------------------------------------------------------------------------
par_delta <- c(
  rho.max  = 401.2588,
  rho.null = 81.19417,
  c.ov     = 0.0005104722,
  k.ov     = 0.37856737,
  k        = 0.02993175,
  tau      = 0.02362476,
  eta.null = 8523356
)

par_scale <- c(1000, 1000, 0.001, 1, 0.1, 0.1, 1e7)
par_start <- par_delta / par_scale

cat("\nStarting parameters (scaled):\n")
print(par_start)
cat("\nStarting parameters (unscaled):\n")
print(par_start * par_scale)

# ----------------------------------------------------------------------------
# BOUNDS FOR DIFFERENTIAL EVOLUTION (physical ranges, then scaled)
# ----------------------------------------------------------------------------
# Physical (unscaled) lower and upper bounds
lower_unscaled <- c(rho.max  = 200,   rho.null = 20,   c.ov = 1e-6,
                    k.ov     = 0,     k        = 0,    tau = 0,
                    eta.null = 1e5)
upper_unscaled <- c(rho.max  = 600,   rho.null = 200,  c.ov = 0.01,
                    k.ov     = 1,     k        = 0.5,  tau = 1,
                    eta.null = 1e8)

# Scale the bounds
lower <- lower_unscaled / par_scale
upper <- upper_unscaled / par_scale

# ----------------------------------------------------------------------------
# PARALLEL SETUP (for the objective function's internal foreach)
# ----------------------------------------------------------------------------
nc <- parallel::detectCores(logical = TRUE) - 1
nc <- max(1, nc)
cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl)

# ----------------------------------------------------------------------------
# TEST OBJECTIVE FUNCTION ON START VALUES
# ----------------------------------------------------------------------------
cat("\nTesting objective function at start values...\n")
test_score <- minimize_score(
  par     = par_start,
  data    = d_obs_fit_tibble,
  scale   = par_scale,
  verbose = TRUE
)
cat("Initial score =", test_score, "\n")

# ----------------------------------------------------------------------------
# OPTIMIZATION WITH DIFFERENTIAL EVOLUTION
# ----------------------------------------------------------------------------
cat("\nStarting optimization with Differential Evolution (DEoptim)...\n")

# Wrap the objective function so that DEoptim can call it with only 'par'
obj_wrapper <- function(par) {
  minimize_score(par, data = d_obs_fit_tibble, scale = par_scale, verbose = FALSE)
}

# Set DEoptim control parameters
# NP = population size (10 * number of parameters is a common choice)
# itermax = maximum number of generations
# F = differential weighting factor
# CR = crossover probability
# trace = print progress every 'trace' generations
DE_ctrl <- DEoptim.control(
  NP       = 10 * length(par_start),   # 70
  itermax  = 200,
  F        = 0.8,
  CR       = 0.9,
  trace    = 10,                       # report every 10 generations
  parallelType = 0                     # no additional parallelisation (we use internal foreach)
)

# Run DEoptim
set.seed(123)  # for reproducibility
opt_de <- DEoptim(
  fn      = obj_wrapper,
  lower   = lower,
  upper   = upper,
  control = DE_ctrl
)

# Extract best solution
best_scaled <- opt_de$optim$bestmem
names(best_scaled) <- names(par_start)
best_unscaled <- best_scaled * par_scale
best_par   <- best_unscaled
best_value <- opt_de$optim$bestval

cat("\n--- Differential Evolution finished ---\n")
cat("Best score achieved:", best_value, "\n")
cat("Best parameters (scaled):\n")
print(best_scaled)
cat("\nBest parameters (unscaled):\n")
print(best_unscaled)

# ----------------------------------------------------------------------------
# PRINT PYTHON CALL WITH OPTIMAL PARAMETERS
# ----------------------------------------------------------------------------
fmt <- function(x) format(x, scientific = FALSE, trim = TRUE, digits = 10)

cat(
  "\nswe_results = pydeltasnow.swe_deltasnow(\n",
  "    idata,\n",
  "    rho_max   = ", fmt(best_unscaled["rho.max"]), ",\n",
  "    rho_null  = ", fmt(best_unscaled["rho.null"]), ",\n",
  "    c_ov      = ", fmt(best_unscaled["c.ov"]), ",\n",
  "    k_ov      = ", fmt(best_unscaled["k.ov"]), ",\n",
  "    k         = ", fmt(best_unscaled["k"]), ",\n",
  "    tau       = ", fmt(best_unscaled["tau"]), ",\n",
  "    eta_null  = ", fmt(best_unscaled["eta.null"]), ",\n",
  "    hs_input_unit=\"m\",\n",
  "    swe_output_unit=\"mm\",\n",
  "    output_series_name=\"SWE_mod\"\n",
  ")\n",
  sep = ""
)

# ----------------------------------------------------------------------------
# SAVE RESULTS WITH "Win21" IN FILENAME
# ----------------------------------------------------------------------------
make_weight_tag <- function(x) {
  out <- format(x, scientific = FALSE, trim = TRUE, digits = 6)
  out <- gsub("\\.", "p", out)
  out <- gsub("-", "m", out)
  out
}

weight_vals <- c(
  SWE_NRMSE = WEIGHT_SWE_NRMSE,
  RHO_NRMSE = WEIGHT_RHO_NRMSE,
  SWE_NBIAS = WEIGHT_SWE_NBIAS
)

weight_tag <- paste0(
  names(weight_vals), "_", make_weight_tag(weight_vals),
  collapse = "__"
)

# Insert "Win21" into the filename
save_dir <- "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_Win21/data/R_opt_logs"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
save_file <- file.path(save_dir, paste0("opt_results_Win21_DE__", weight_tag, ".rds"))

saveRDS(list(
  opt          = opt_de,
  fit_data     = d_obs_fit_tibble,
  val_data     = d_obs_val_tibble,
  best_par     = best_par,
  best_value   = best_value,
  weights      = weight_vals
), file = save_file)

cat("\nOptimization finished. Results saved to:\n", save_file, "\n", sep = "")

# ----------------------------------------------------------------------------
# STOP CLUSTER
# ----------------------------------------------------------------------------
parallel::stopCluster(cl)