##############################################################################
# DELTASNOW PARAMETER OPTIMIZATION WITH WEIGHT HYPERPARAMETER TUNING
##############################################################################
#
# Two-level optimization:
#   Outer: find weights (w1, w2, w3) that minimise validation score.
#   Inner: for given weights, optimise DeltaSnow parameters on fit data.
#
# Validation score uses the same functional form as fit score.
#
##############################################################################

library(optimx)
library(zoo)
library(foreach)
library(doParallel)
library(lubridate)
library(nixmass)
library(tidyverse)

# ----------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------

# Number of weight grid points per dimension (total points = n * (n+1) / 2)
GRID_SIZE <- 5   # 5 -> 15 weight combinations; increase for finer search

# Data settings
season_start   <- "-08-01"
season_end     <- "-07-31"
start_of_block <- 8

# Model settings
EPS <- 1e-6

# Optimization settings
MAXIT <- 300      # fewer iterations per inner run to save time
METHOD <- "Nelder-Mead"

# ----------------------------------------------------------------------------
# LOAD OBSERVATIONAL DATA
# ----------------------------------------------------------------------------
d_obs <- get(load(
  "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_SNOWPACK/data/d_obs_SNOWPACK.rda"
))

# ----------------------------------------------------------------------------
# HELPER FUNCTIONS
# ----------------------------------------------------------------------------
set_season <- function(date, start_month = 8) {
  date <- as.Date(date)
  yr <- year(date)
  mo <- month(date)
  ifelse(mo < start_month, yr - 1, yr)
}

rmse <- function(obs, mod) {
  ok <- is.finite(obs) & is.finite(mod)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((mod[ok] - obs[ok])^2))
}

nrmse <- function(obs, mod) {
  ok <- is.finite(obs) & is.finite(mod)
  if (!any(ok)) return(NA_real_)
  rmse_val <- sqrt(mean((mod[ok] - obs[ok])^2))
  mean_obs <- mean(obs[ok], na.rm = TRUE)
  rmse_val / mean_obs
}

nbias <- function(obs, mod) {
  ok <- is.finite(obs) & is.finite(mod)
  if (!any(ok)) return(NA_real_)
  bias_val <- mean(mod[ok] - obs[ok], na.rm = TRUE)
  mean_obs <- mean(obs[ok], na.rm = TRUE)
  abs(bias_val) / mean_obs
}

combined_score <- function(df, weights, eps = EPS) {
  swe_obs <- df$swe_obs
  swe_mod <- df$swe_mod
  hs <- df$hs
  
  nrmse_swe <- nrmse(swe_obs, swe_mod)
  nbias_swe <- nbias(swe_obs, swe_mod)
  
  rho_obs <- ifelse(is.finite(hs) & hs > eps, swe_obs / hs, NA_real_)
  rho_mod <- ifelse(is.finite(hs) & hs > eps, swe_mod / hs, NA_real_)
  nrmse_rho <- nrmse(rho_obs, rho_mod)
  
  score <- weights[1] * nrmse_swe + weights[2] * nrmse_rho + weights[3] * nbias_swe
  
  attr(score, "metrics") <- list(
    rmse_swe  = rmse(swe_obs, swe_mod),
    rmse_rho  = rmse(rho_obs, rho_mod),
    nrmse_swe = nrmse_swe,
    nrmse_rho = nrmse_rho,
    bias_swe  = mean(swe_mod - swe_obs, na.rm = TRUE),
    nbias_swe = nbias_swe
  )
  return(score)
}

# ----------------------------------------------------------------------------
# SPLIT DATA INTO FIT AND VALIDATION (same as before)
# ----------------------------------------------------------------------------
d_obs_fit <- list()
d_obs_val <- list()

for (station in names(d_obs)) {
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
    
    if (nrow(winter) < 200) next
    if (nrow(winter) < 365 && (as.numeric(winter$Hobs[1]) > 0.05 ||
                               as.numeric(winter$Hobs[nrow(winter)]) > 0.05)) next
    
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
# PREPARE DATA TIBBLES (fit and validation)
# ----------------------------------------------------------------------------
prepare_tibble <- function(data_list) {
  out_list <- list()
  for (station in names(data_list)) {
    x <- data_list[[station]]
    if (is.null(x) || length(x) == 0) next
    dat <- tryCatch(as_tibble(coredata(x)), error = function(e) NULL)
    if (is.null(dat) || ncol(dat) == 0) next
    if (ncol(dat) == 2 && all(is.na(colnames(dat)))) {
      colnames(dat) <- c("Hobs", "SWEobs")
    }
    if (!all(c("Hobs", "SWEobs") %in% colnames(dat))) next
    out_list[[station]] <- tibble(
      date    = as.Date(index(x)),
      name    = station,
      hs      = dat$Hobs,
      swe_obs = dat$SWEobs,
      block   = set_season(index(x), start_of_block)
    )
  }
  bind_rows(out_list)
}

d_obs_fit_tibble <- prepare_tibble(d_obs_fit)
d_obs_val_tibble <- prepare_tibble(d_obs_val)

# ----------------------------------------------------------------------------
# RUN DELTASNOW OVER ALL STATIONS/BLOCKS (used inside objective)
# ----------------------------------------------------------------------------
run_deltasnow_all <- function(par_real, data) {
  station_results <- foreach(
    station = unique(data$name),
    .packages = c("dplyr", "tidyr", "nixmass"),
    .combine = bind_rows
  ) %dopar% {
    data_station <- data %>% filter(name == station)
    blocks <- unique(data_station$block)
    
    block_results <- lapply(blocks, function(b) {
      df_block <- data_station %>% filter(block == b)
      full_dates <- tibble(date = seq(min(df_block$date), max(df_block$date), by = "1 day"))
      
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
              dyn_rho_max = TRUE
            )
        },
        error = function(e) rep(NA_real_, nrow(joined))
      )
      
      joined %>% mutate(swe_mod = swe_mod) %>% drop_na()
    })
    bind_rows(block_results)
  }
  return(station_results)
}

# ----------------------------------------------------------------------------
# OBJECTIVE FUNCTION FACTORY (inner optimization)
# ----------------------------------------------------------------------------
make_inner_objective <- function(weights, data_fit, scale, verbose = FALSE) {
  force(weights); force(data_fit); force(scale)
  
  function(par) {
    result <- tryCatch({
      par_real <- par * scale
      dff <- run_deltasnow_all(par_real, data_fit)
      if (nrow(dff) == 0) return(1e12)
      
      score_with_attr <- combined_score(dff, weights)
      score <- as.numeric(score_with_attr)
      if (!is.finite(score)) return(1e12)
      
      if (verbose) cat("  inner score =", round(score, 6), "\n")
      return(score)
    }, error = function(e) 1e12)
    
    if (length(result) != 1 || !is.numeric(result)) return(1e12)
    return(result)
  }
}

# ----------------------------------------------------------------------------
# OUTER OBJECTIVE: evaluate a weight vector on validation data
# ----------------------------------------------------------------------------
evaluate_weights <- function(weights, par_start, par_scale, data_fit, data_val) {
  # Inner optimization
  obj_inner <- make_inner_objective(weights, data_fit, par_scale, verbose = FALSE)
  
  opt <- tryCatch(
    optimx(
      par     = par_start,
      fn      = obj_inner,
      method  = METHOD,
      control = list(trace = 0, maxit = MAXIT)
    ),
    error = function(e) NULL
  )
  
  if (is.null(opt) || nrow(opt) == 0) return(Inf)
  
  best_par_scaled <- as.numeric(opt[1, 1:7])
  best_par_real <- best_par_scaled * par_scale
  
  # Evaluate on validation data
  dff_val <- run_deltasnow_all(best_par_real, data_val)
  if (nrow(dff_val) == 0) return(Inf)
  
  val_score <- as.numeric(combined_score(dff_val, weights))
  return(val_score)
}

# ----------------------------------------------------------------------------
# GENERATE WEIGHT GRID ON THE SIMPLEX
# ----------------------------------------------------------------------------
# Create all combinations of (w1, w2, w3) with w1 + w2 + w3 = 1, w_i >= 0
step <- 1 / (GRID_SIZE - 1)
w1_vals <- seq(0, 1, by = step)
weight_grid <- tibble()
for (w1 in w1_vals) {
  for (w2 in seq(0, 1 - w1, by = step)) {
    w3 <- 1 - w1 - w2
    if (w3 >= -1e-12) {  # tolerate tiny floating errors
      weight_grid <- bind_rows(weight_grid, tibble(w1 = w1, w2 = w2, w3 = w3))
    }
  }
}
weight_grid <- weight_grid %>% filter(w3 >= 0)

cat("Number of weight combinations to test:", nrow(weight_grid), "\n")

# ----------------------------------------------------------------------------
# PARALLEL SETUP (outer loop will be sequential, inner parallel already)
# ----------------------------------------------------------------------------
nc <- parallel::detectCores(logical = TRUE) - 1
nc <- max(1, nc)
cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl)

# ----------------------------------------------------------------------------
# START VALUES AND SCALING
# ----------------------------------------------------------------------------
par_delta <- c(
  rho.max  = 300,
  rho.null = 100,
  c.ov     = 0.0004986025,
  k.ov     = 0.227786,
  k        = 0.0290256,
  tau      = 0.02489556,
  eta.null = 8792253
)
par_scale <- c(1000, 1000, 0.001, 1, 0.1, 0.1, 1e7)
par_start <- par_delta / par_scale

# ----------------------------------------------------------------------------
# RUN WEIGHT OPTIMIZATION (GRID SEARCH)
# ----------------------------------------------------------------------------
results_weight <- list()

for (i in seq_len(nrow(weight_grid))) {
  w <- as.numeric(weight_grid[i, ])
  cat(sprintf("\n--- Weight set %d: w1=%.2f, w2=%.2f, w3=%.2f ---\n", i, w[1], w[2], w[3]))
  
  val_score <- evaluate_weights(w, par_start, par_scale, d_obs_fit_tibble, d_obs_val_tibble)
  
  cat(sprintf("Validation score: %.6f\n", val_score))
  
  results_weight[[i]] <- tibble(
    w1 = w[1], w2 = w[2], w3 = w[3],
    val_score = val_score
  )
}

results_weight_df <- bind_rows(results_weight)

# ----------------------------------------------------------------------------
# FIND BEST WEIGHTS AND RE-OPTIMIZE WITH FULL ITERATIONS
# ----------------------------------------------------------------------------
best_idx <- which.min(results_weight_df$val_score)
best_weights <- as.numeric(results_weight_df[best_idx, 1:3])
cat("\nBest weights from grid search:\n")
print(best_weights)

# Final optimization with best weights and more iterations
cat("\nRunning final optimization with best weights (maxit = 500)...\n")
obj_final <- make_inner_objective(best_weights, d_obs_fit_tibble, par_scale, verbose = TRUE)

opt_final <- optimx(
  par     = par_start,
  fn      = obj_final,
  method  = METHOD,
  control = list(trace = 1, maxit = 500)
)

best_par_scaled <- as.numeric(opt_final[1, 1:7])
best_par_real <- best_par_scaled * par_scale

final_results <- tibble(
  parameter = c("rho_max", "rho_null", "c_ov", "k_ov", "k", "tau", "eta_null"),
  value = best_par_real
)
print(final_results)

# ----------------------------------------------------------------------------
# SAVE
# ----------------------------------------------------------------------------
saveRDS(list(weight_grid_results = results_weight_df,
             best_weights = best_weights,
             final_parameters = final_results),
        file = "weight_optimization_results.rds")

# ----------------------------------------------------------------------------
# STOP CLUSTER
# ----------------------------------------------------------------------------
parallel::stopCluster(cl)

cat("\nWeight optimization complete.\n")