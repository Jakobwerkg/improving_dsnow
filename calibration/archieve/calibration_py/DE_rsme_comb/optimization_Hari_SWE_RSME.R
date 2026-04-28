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


#---------------------------------------------------------------------------------------------------
# real obs: HD are weekly, Swiss data are bi-weekly
# split data into two halfs for fitting and verification
# use every 2nd winter for fitting, every other for validation
# create a tibble

d_obs <- get(load("/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_data/raw_data/dsnow/Win21_calib/H_SWE_obs.Rda"))


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
# Pre-split data by station and block for efficiency (avoids repeated filter calls)
data_split <- d_obs_fit_tibble |>
  group_by(name, block) |>
  group_split()

# pre-compute station assignment for each split element
data_split_station <- sapply(data_split, function(x) x$name[1])

# pre-compute mean observed SWE and rho_bulk for normalizing the combined score
d_obs_nona <- d_obs_fit_tibble |> filter(!is.na(swe_obs))
mean_swe_obs <- mean(d_obs_nona$swe_obs)
mean_rho_obs <- mean(d_obs_nona$swe_obs[d_obs_nona$hs > 0] /
                     d_obs_nona$hs[d_obs_nona$hs > 0])
cat(paste0("Normalization: mean_swe_obs=", round(mean_swe_obs, 4),
           " mean_rho_obs=", round(mean_rho_obs, 4), "\n"))


#---------------------------------------------------------------------------------------------------
# function for score to be minimized
# score: which metric to minimize
#   "rmse_swe"  - RMSE of SWE            (default)
#   "bias_swe"  - |bias| of SWE
#   "rmse_rho"  - RMSE of rho_bulk (SWE/HS)
#   "bias_rho"  - |bias| of rho_bulk (SWE/HS)
#   "rmse_combined" - 50/50 weighted: 0.5*(RMSE_SWE/mean_SWE) + 0.5*(RMSE_rho/mean_rho)
minimize_score <- function(data, par, scale, score = "0.5*(RMSE_rho/mean_rho)",
                           data_split, data_split_station,
                           mean_swe_obs = NULL, mean_rho_obs = NULL,
                           verbose = FALSE) {
  par <- par * scale
  cat(par)

  # build model_opts once
  mopts <- list(
    rho.max = par[1], rho.null = par[2], c.ov = par[3],
    k.ov = par[4], k = par[5], tau = par[6], eta.null = par[7]
  )

  # parallel over pre-split station x block chunks
  ll <- foreach(
    chunk = data_split,
    .packages = c("dplyr", "zoo", "nixmass"),
    .export = "mopts"
  ) %dopar% {
    left <- chunk |> dplyr::select(date, hs, swe_obs)
    all_dates <- data.frame(
      date = seq(min(left$date), max(left$date), by = "1 day")
    )
    joined <- left |> right_join(all_dates, by = "date")

    swe_mod <- tryCatch(
      {
        joined |>
          dplyr::select(date, hs) |>
          mutate(date = as.character(date)) |>
          nixmass::swe.delta.snow(model_opts = mopts, dyn_rho_max = TRUE)
      },
      error = function(e) rep(NA, nrow(joined))
    )
    na.omit(cbind(joined, swe_mod))
  }

  dff <- do.call(rbind, ll)
  swe_err <- dff$swe_mod - dff$swe_obs
  rmse_swe <- sqrt(mean(swe_err^2))
  bias_swe <- abs(mean(swe_err))

  # rho_bulk (SWE/HS) where HS > 0
  idx_pos <- dff$hs > 0
  if (any(idx_pos)) {
    rho_err <- (dff$swe_mod[idx_pos] / dff$hs[idx_pos]) -
               (dff$swe_obs[idx_pos] / dff$hs[idx_pos])
    rmse_rho <- sqrt(mean(rho_err^2))
    bias_rho <- abs(mean(rho_err))
  } else {
    rmse_rho <- NA
    bias_rho <- NA
  }

  # combined score: 50/50 normalized RMSE
  if (!is.na(rmse_rho) && !is.null(mean_swe_obs) && !is.null(mean_rho_obs)) {
    rmse_combined <- 0.5 * (rmse_swe / mean_swe_obs) + 0.5 * (rmse_rho / mean_rho_obs)
  } else {
    rmse_combined <- NA
  }

  cat(paste0(
    " |bias_swe|=", round(bias_swe, 4),
    " rmse_swe=", round(rmse_swe, 4),
    " |bias_rho|=", round(bias_rho, 4),
    " rmse_rho=", round(rmse_rho, 4),
    " rmse_comb=", round(rmse_combined, 4),
    "\n"
  ))

  switch(score,
    "rmse_swe" = rmse_swe,
    "bias_swe" = bias_swe,
    "rmse_rho" = rmse_rho,
    "bias_rho" = bias_rho,
    "rmse_combined" = rmse_combined,
    stop(paste0("Unknown score: ", score,
                ". Use rmse_swe, bias_swe, rmse_rho, bias_rho, or rmse_combined."))
  )
}


#-------------------------------------------------------------------------
# start values
# taken from a quick optimization
par_delta <- c(
  rho.max = 420,
  rho.null = 100,
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
# optimization: run all 4 score objectives
nc <- detectCores(logical = TRUE)
cl <- makeCluster(nc)
registerDoParallel(cl)

scores <- c("rmse_swe", "bias_swe", "rmse_rho", "bias_rho", "rmse_combined")
opt_results <- list()

for (sc in scores) {
  cat(paste0("\n===== Optimizing for: ", sc, " =====\n"))
  opt_results[[sc]] <- optimx(
    fn = minimize_score,
    data = d_obs_fit_tibble,
    par = par_delta,
    scale = par_scale,
    score = sc,
    data_split = data_split,
    data_split_station = data_split_station,
    mean_swe_obs = mean_swe_obs,
    mean_rho_obs = mean_rho_obs,
    verbose = FALSE,
    method = c("Nelder-Mead", "newuoa", "bobyqa"),
    control = list(trace = 4, follow.on = TRUE)
    #lower = lower,
    #upper = upper
  )
}

saveRDS(opt_results, file = "opt_results_all_scores.rds")
stopCluster(cl)
