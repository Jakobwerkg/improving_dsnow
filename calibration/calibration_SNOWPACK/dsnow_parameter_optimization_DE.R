##############################################################################
# DELTASNOW PARAMETER OPTIMIZATION (Differential Evolution)
##############################################################################
#
# Objective:
#   final_score = w1 * NRMSE_SWE + w2 * NRMSE_rho + w3 * NBIAS_SWE
#
# where:
#   NRMSE_SWE  = RMSE_SWE  / mean_obs_SWE
#   NRMSE_rho  = RMSE_rho  / mean_obs_rho
#   NBIAS_SWE  = |bias_SWE| / mean_obs_SWE
#
# Weights (w1, w2, w3) can be set below.
#
##############################################################################

library(DEoptim)
library(zoo)
library(foreach)
library(doParallel)
library(lubridate)
library(nixmass)
library(tidyverse)

# ----------------------------------------------------------------------------
# CONFIGURATION (USER-ADJUSTABLE)
# ----------------------------------------------------------------------------

# Weights for the combined objective (should sum to 1)
WEIGHT_SWE_NRMSE  <- 0.2    # weight for normalized SWE RMSE
WEIGHT_RHO_NRMSE  <- 0.8    # weight for normalized density RMSE
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

# DE settings
DE_ITERMAX   <- 100              # maximum number of DE iterations
DE_NP        <- 70               # population size (rule of thumb: 10 * n_params)
DE_F         <- 0.8              # differential weight
DE_CR        <- 0.9              # crossover probability
DE_STRATEGY  <- 2                # 1 = DE/rand/1/bin, 2 = DE/local-to-best/1/bin

# ----------------------------------------------------------------------------
# LOAD OBSERVATIONAL DATA
# ----------------------------------------------------------------------------
# Adjust the path to your .rda file containing a list of zoo objects named 'd_obs'
d_obs <- get(load(
  "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_SNOWPACK/data/d_obs_SNOWPACK.rda"
))

# ----------------------------------------------------------------------------
# HELPER FUNCTIONS
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
# PREPARE FIT DATA (convert zoo to tibble)
# ----------------------------------------------------------------------------
fit_data_list <- list()

for (station in names(d_obs_fit)) {
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
    hs      = dat$Hobs,
    swe_obs = dat$SWEobs,
    block   = set_season(index(x), start_of_block)
  )
}

d_obs_fit_tibble <- bind_rows(fit_data_list)

# ----------------------------------------------------------------------------
# OBJECTIVE FUNCTION (called by DEoptim)
# ----------------------------------------------------------------------------
# Note: DEoptim passes parameters as a vector directly (no scaling needed internally
# if we define bounds in real space). We keep verbose output consistent.

minimize_score_de <- function(par, data, verbose = FALSE) {
  
  # Ensure a scalar numeric is returned even on error
  result <- tryCatch({
    
    if (verbose) {
      param_names <- c("rho.max", "rho.null", "c.ov", "k.ov", "k", "tau", "eta.null")
      values <- paste(
        paste0(param_names, "=", round(par, 6)),
        collapse = ", "
      )
      cat("params =", values, "\n")
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
                  rho.max  = par[1],
                  rho.null = par[2],
                  c.ov     = par[3],
                  k.ov     = par[4],
                  k        = par[5],
                  tau      = par[6],
                  eta.null = par[7]
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
    score <- as.numeric(score_with_attr)   # strip attributes
    
    if (!is.finite(score)) {
      cat("Score is not finite. Returning large penalty.\n")
      return(1e12)
    }
    
    # Verbose output
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
    } else if (verbose) {
      cat("Score =", round(score, 4), "\n")
    }
    
    return(score)
    
  }, error = function(e) {
    cat("Error in objective function:", e$message, "\n")
    return(1e12)
  })
  
  # Final safeguard
  if (length(result) != 1 || !is.numeric(result)) {
    return(1e12)
  }
  return(result)
}

# ----------------------------------------------------------------------------
# PARAMETER BOUNDS (in real/unscaled space)
# ----------------------------------------------------------------------------
# Define reasonable lower and upper bounds for each parameter

par_lower <- c(
  rho.max  = 300,
  rho.null = 60,
  c.ov     = 1e-6,
  k.ov     = 0.01,
  k        = 0.001,
  tau      = 0.001,
  eta.null = 1e5
)

par_upper <- c(
  rho.max  = 600,
  rho.null = 150,
  c.ov     = 0.01,
  k.ov     = 1.0,
  k        = 0.1,
  tau      = 0.1,
  eta.null = 1e8
)

cat("\nParameter bounds:\n")
cat("Lower:", par_lower, "\n")
cat("Upper:", par_upper, "\n")

# ----------------------------------------------------------------------------
# PARALLEL SETUP
# ----------------------------------------------------------------------------
nc <- parallel::detectCores(logical = TRUE) - 1
nc <- max(1, nc)
cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl)

# Export required objects to cluster for nested parallelism
parallel::clusterExport(cl, c("WEIGHT_SWE_NRMSE", "WEIGHT_RHO_NRMSE", 
                               "WEIGHT_SWE_NBIAS", "EPS"))

# ----------------------------------------------------------------------------
# TEST OBJECTIVE FUNCTION AT INITIAL GUESS
# ----------------------------------------------------------------------------
par_init <- c(
  rho.max  = 401.2588,
  rho.null = 81.19417,
  c.ov     = 0.0005104722,
  k.ov     = 0.37856737,
  k        = 0.02993175,
  tau      = 0.02362476,
  eta.null = 8523356
)

cat("\nTesting objective function at initial guess...\n")
test_score <- minimize_score_de(
  par     = par_init,
  data    = d_obs_fit_tibble,
  verbose = TRUE
)
cat("Initial score =", test_score, "\n")

# ----------------------------------------------------------------------------
# OPTIMIZATION (Differential Evolution)
# ----------------------------------------------------------------------------
cat("\nStarting Differential Evolution optimization...\n")
cat("  itermax =", DE_ITERMAX, "\n")
cat("  NP =", DE_NP, "\n")
cat("  F =", DE_F, "\n")
cat("  CR =", DE_CR, "\n")
cat("  strategy =", DE_STRATEGY, "\n\n")

opt <- DEoptim(
  fn      = minimize_score_de,
  lower   = par_lower,
  upper   = par_upper,
  data    = d_obs_fit_tibble,
  verbose = TRUE,
  control = DEoptim.control(
    itermax   = DE_ITERMAX,
    NP        = DE_NP,
    F         = DE_F,
    CR        = DE_CR,
    strategy  = DE_STRATEGY,
    trace     = 10,           # print progress every 10 iterations
    parallelType = 0          # we handle parallelism inside the objective
  )
)

# Extract best parameters
best_par <- opt$optim$bestmem
names(best_par) <- names(par_lower)
best_value <- opt$optim$bestval

cat("\n======================================\n")
cat("Optimization complete!\n")
cat("Best score:", best_value, "\n")
cat("Best parameters:\n")
print(best_par)

# Print final parameters as Python pydeltasnow call
fmt <- function(x) format(x, scientific = FALSE, trim = TRUE, digits = 10)

cat(
  "\nswe_results = pydeltasnow.swe_deltasnow(\n",
  "    idata,\n",
  "    rho_max   = ", fmt(best_par["rho.max"]), ",\n",
  "    rho_null  = ", fmt(best_par["rho.null"]), ",\n",
  "    c_ov      = ", fmt(best_par["c.ov"]), ",\n",
  "    k_ov      = ", fmt(best_par["k.ov"]), ",\n",
  "    k         = ", fmt(best_par["k"]), ",\n",
  "    tau       = ", fmt(best_par["tau"]), ",\n",
  "    eta_null  = ", fmt(best_par["eta.null"]), ",\n",
  "    hs_input_unit=\"m\",\n",
  "    swe_output_unit=\"mm\",\n",
  "    output_series_name=\"SWE_mod\"\n",
  ")\n",
  sep = ""
)

# ----------------------------------------------------------------------------
# SAVE RESULTS
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

save_dir <- "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_SNOWPACK/data/R_opt_logs_DE"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

save_file <- file.path(
  save_dir,
  paste0("opt_results_DE__", weight_tag, ".rds")
)

saveRDS(opt, file = save_file)

cat("\nOptimization finished. Results saved to:\n", save_file, "\n", sep = "")

# ----------------------------------------------------------------------------
# STOP CLUSTER
# ----------------------------------------------------------------------------
parallel::stopCluster(cl)
