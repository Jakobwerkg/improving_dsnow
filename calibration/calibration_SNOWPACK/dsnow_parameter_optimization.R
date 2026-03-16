##############################################################################
# DELTASNOW PARAMETER OPTIMIZATION
##############################################################################
#
# Objective:
#   final_score = 0.5 * RMSE(SWE) + 0.5 * RMSE(bulk density)
#
# bulk density = SWE / HS
#
##############################################################################

library(optimx)
library(zoo)
library(foreach)
library(doParallel)
library(lubridate)
library(nixmass)
library(tidyverse)

##############################################################################
# SETTINGS
##############################################################################

season_start   <- "-08-01"
season_end     <- "-07-31"
start_of_block <- 8

##############################################################################
# LOAD DATA
##############################################################################

d_obs <- get(load(
  "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_SNOWPACK/data/d_obs_SNOWPACK.rda"
))

##############################################################################
# HELPER FUNCTIONS
##############################################################################

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

combined_rmse <- function(df, eps = 1e-6) {
  rmse_swe <- rmse(df$swe_obs, df$swe_mod)

  rho_obs <- ifelse(is.finite(df$hs) & df$hs > eps, df$swe_obs / df$hs, NA_real_)
  rho_mod <- ifelse(is.finite(df$hs) & df$hs > eps, df$swe_mod / df$hs, NA_real_)

  rmse_rho <- rmse(rho_obs, rho_mod)

  0.5 * rmse_swe + 0.5 * rmse_rho
}

##############################################################################
# SPLIT DATA INTO FIT AND VALIDATION
##############################################################################

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

##############################################################################
# PREPARE FIT DATA
##############################################################################

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

##############################################################################
# OBJECTIVE FUNCTION
##############################################################################

minimize_score <- function(par, data, scale, verbose = FALSE) {

  par_real <- par * scale

  cat("pars =", paste(round(par_real, 6), collapse = ", "), "\n")

  station_results <- foreach(
    station = unique(data$name),
    .packages = c("dplyr", "tidyr", "nixmass")
  ) %dopar% {

    if (verbose) cat("station:", station, "\n")

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

  dff <- bind_rows(station_results)

  if (nrow(dff) == 0) {
    cat("No valid model output. Returning large penalty.\n")
    return(1e12)
  }

  score <- combined_rmse(dff)

  if (!is.finite(score)) {
    cat("Score is not finite. Returning large penalty.\n")
    return(1e12)
  }

  rmse_swe <- rmse(dff$swe_obs, dff$swe_mod)

  rho_obs <- ifelse(dff$hs > 1e-6, dff$swe_obs / dff$hs, NA_real_)
  rho_mod <- ifelse(dff$hs > 1e-6, dff$swe_mod / dff$hs, NA_real_)
  rmse_rho <- rmse(rho_obs, rho_mod)

  bias_swe <- mean(dff$swe_mod - dff$swe_obs, na.rm = TRUE)

  cat(
    "| bias_swe =", round(bias_swe, 4),
    "| rmse_swe =", round(rmse_swe, 4),
    "| rmse_density =", round(rmse_rho, 4),
    "| final_score =", round(score, 4),
    "\n"
  )

  score
}

##############################################################################
# START VALUES AND SCALING
##############################################################################

par_delta <- c(
  rho.max  = 396.9521,
  rho.null = 78.88324,
  c.ov     = 0.0004986025,
  k.ov     = 0.227786,
  k        = 0.0290256,
  tau      = 0.02489556,
  eta.null = 8792253
)

par_scale <- c(1000, 1000, 0.001, 1, 0.1, 0.1, 1e7)

par_start <- par_delta / par_scale

print(par_start)

##############################################################################
# PARALLEL SETUP
##############################################################################

nc <- parallel::detectCores(logical = TRUE) - 1
nc <- max(1, nc)

cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl)

##############################################################################
# TEST OBJECTIVE FUNCTION ON START VALUES
##############################################################################

test_score <- minimize_score(
  par   = par_start,
  data  = d_obs_fit_tibble,
  scale = par_scale,
  verbose = FALSE
)

cat("Initial score =", test_score, "\n")

##############################################################################
# OPTIMIZATION
##############################################################################

opt <- optimx(
  par     = par_start,
  fn      = minimize_score,
  data    = d_obs_fit_tibble,
  scale   = par_scale,
  verbose = FALSE,
  method  = "Nelder-Mead",
  control = list(trace = 1, maxit = 500)
)

##############################################################################
# SAVE RESULTS
##############################################################################

saveRDS(opt, file = "opt_results_combined_rmse_nelder.rds")

##############################################################################
# STOP CLUSTER
##############################################################################

parallel::stopCluster(cl)