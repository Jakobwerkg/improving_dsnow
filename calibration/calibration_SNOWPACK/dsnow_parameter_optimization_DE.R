#-------------------------------------------------------------------------
# parameter optimization for deltasnow model using Differential Evolution
# optimize deltasnow model with respect to RMSE between observed and modeled SWE
# revised to avoid nested parallel problems
#-------------------------------------------------------------------------

library(DEoptim)
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

d_obs <- get(load(
  paste0(
    "/Users/jakobwerkgarner/code/mt_dsnow/calibration/",
    "calibration_SNOWPACK/data/d_obs_SNOWPACK.rda"
  )
))

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

    if (length(winter) == 0 || nrow(winter) < 200) {
      print(paste0("only ", ifelse(length(winter) == 0, 0, nrow(winter)), " values...skipping year ", y))
      next
    }

    dat_w <- as_tibble(coredata(winter))
    if (all(is.na(colnames(dat_w))) && ncol(dat_w) == 2) {
      colnames(dat_w) <- c("Hobs", "SWEobs")
    }

    if (!all(c("Hobs", "SWEobs") %in% colnames(dat_w))) {
      print(paste0("missing Hobs/SWEobs...skipping year ", y))
      next
    }

    if (
      nrow(winter) < 365 &&
      (as.numeric(dat_w$Hobs[1]) > 0.05 ||
         as.numeric(dat_w$Hobs[nrow(dat_w)]) > 0.05)
    ) {
      print(paste0(
        "only ", nrow(winter),
        " values and last/first value > 0.05...skipping year ", y
      ))
      next
    }

    if (idx_y %% 2 == 0) {
      if (length(winters4fit) == 0) {
        winters4fit <- winter
      } else {
        winters4fit <- rbind(winters4fit, winter)
      }
    } else {
      if (length(winters4val) == 0) {
        winters4val <- winter
      } else {
        winters4val <- rbind(winters4val, winter)
      }
    }

    idx_y <- idx_y + 1
  }

  d_obs_fit[[n]] <- winters4fit
  d_obs_val[[n]] <- winters4val
}


# prepare calibration data
start_of_block <- 8

# x...character date
# m...integer start month of block (1-12)
# return season as year
set_season <- function(x, m) {
  x <- as.character(x)
  yr <- as.POSIXlt(x)$year + 1900
  mt <- as.POSIXlt(x)$mon + 1
  ifelse(mt < m, yr - 1, yr)
}

d_obs_fit_tibble <- lapply(seq_along(d_obs_fit), function(i) {
  x <- d_obs_fit[[i]]
  station_name <- names(d_obs_fit)[i]

  if (is.null(x) || length(x) == 0) return(NULL)

  dat <- as_tibble(coredata(x))
  if (all(is.na(colnames(dat))) && ncol(dat) == 2) {
    colnames(dat) <- c("Hobs", "SWEobs")
  }

  if (!all(c("Hobs", "SWEobs") %in% colnames(dat))) {
    message("Skipping station ", station_name, ": missing Hobs/SWEobs")
    return(NULL)
  }

  tibble(
    date = index(x),
    name = station_name,
    hs = as.numeric(dat$Hobs),
    swe_obs = as.numeric(dat$SWEobs)
  ) |>
    mutate(block = set_season(date, start_of_block))
})

d_obs_fit_tibble <- bind_rows(d_obs_fit_tibble)


# function for score to be minimized
# IMPORTANT: use %do% here because DEoptim already parallelizes externally
minimize_score <- function(par, data, scale, verbose = FALSE) {
  par <- par * scale
  cat(par)

  ll <- foreach(
    s = unique(data$name),
    .packages = c("dplyr", "zoo", "nixmass", "foreach")
  ) %do% {
    if (verbose) cat(paste0(s, " ..."))

    data1 <- data |>
      dplyr::filter(name == s)

    l <- foreach(
      i = seq_along(unique(data1$block)),
      .packages = c("dplyr", "zoo", "nixmass")
    ) %do% {
      y <- unique(data1$block)[i]

      left <- data1 |>
        dplyr::filter(block == y) |>
        dplyr::select(date, hs, swe_obs)

      if (nrow(left) == 0) return(NULL)

      right <- data.frame(
        date = seq(min(left$date), max(left$date), by = "1 day")
      )

      joined <- left |>
        dplyr::right_join(right, by = "date") |>
        dplyr::arrange(date)

      swe_mod <- tryCatch(
        {
          joined |>
            dplyr::select(date, hs) |>
            dplyr::mutate(date = as.character(date)) |>
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
        error = function(e) rep(NA_real_, nrow(joined))
      )

      out <- cbind(joined, swe_mod = swe_mod)
      out <- out[is.finite(out$swe_obs) & is.finite(out$swe_mod), , drop = FALSE]
      if (nrow(out) == 0) return(NULL)
      out
    }

    l <- l[!vapply(l, is.null, logical(1))]
    if (length(l) == 0) return(NULL)
    do.call(rbind, l)
  }

  ll <- ll[!vapply(ll, is.null, logical(1))]
  if (length(ll) == 0) return(1e12)

  dff <- do.call(rbind, ll)
  if (is.null(dff) || nrow(dff) == 0) return(1e12)

  rmse <- with(dff, sqrt(mean((swe_mod - swe_obs)^2, na.rm = TRUE)))
  bias <- with(dff, abs(mean(swe_mod - swe_obs, na.rm = TRUE)))
  if (!is.finite(rmse)) rmse <- 1e12

  cat(paste0(" |bias|=", bias, " rmse=", rmse, "\n"))
  rmse
}

# DEoptim objective wrapper
de_obj <- function(par) {
  minimize_score(
    par = par,
    data = d_obs_fit_tibble,
    scale = par_scale,
    verbose = FALSE
  )
}

# start values
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


# optimization with Differential Evolution
nc <- detectCores(logical = TRUE)
cl <- makeCluster(nc)
registerDoParallel(cl)

clusterExport(
  cl,
  varlist = c("minimize_score", "de_obj", "d_obs_fit_tibble", "par_scale"),
  envir = environment()
)

ctrl <- DEoptim.control(
  trace = TRUE,
  itermax = 200,
  NP = 70,
  parallelType = 1,
  packages = c("dplyr", "zoo", "nixmass", "foreach")
)

opt <- DEoptim(
  fn = de_obj,
  lower = lower,
  upper = upper,
  control = ctrl
)

saveRDS(opt, file = "opt_results_DEoptim.rds")

cat("\nBest RMSE:", opt$optim$bestval, "\n")
cat("Best parameters (scaled):\n")
print(opt$optim$bestmem)
cat("Best parameters (original scale):\n")
print(opt$optim$bestmem * par_scale)

stopCluster(cl)