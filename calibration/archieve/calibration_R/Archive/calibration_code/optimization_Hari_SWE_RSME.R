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
d_obs <- get(load("H_SWE_obs.Rda"))
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
                nixmass::swe.delta.snow(
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

  rho_bulk_mod <- swe_mod / data$hs
  rho_bulk_obs <- swe_obs / data$hs 


  rmse <- with(dff, sqrt(mean((rho_bulk_mod - rho_bulk_obs)^2)))
  bias <- with(dff, abs(mean(rho_bulk_mod - rho_bulk_obs)))

  # rmse <- with(dff, sqrt(mean((swe_mod - swe_obs)^2)))
  # bias <- with(dff, abs(mean(swe_mod - swe_obs)))
  cat(paste0(" |bias|=", bias, " rmse rho =", rmse, "\n"))
  rmse
}


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
stopCluster(cl)
