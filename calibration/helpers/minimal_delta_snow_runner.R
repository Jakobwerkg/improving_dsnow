#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  args <- commandArgs(trailingOnly = TRUE)
}))

get_arg <- function(flag, default = NULL) {
  hit <- which(args %in% flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop(paste("Missing value for", flag))
  args[hit + 1]
}

# Required args
in_path <- get_arg("--in")
if (is.null(in_path)) stop("Missing input path.")

# Model parameters
rho.max   <- as.numeric(get_arg("--rho.max", "400"))
rho.null  <- as.numeric(get_arg("--rho.null", "80"))
c.ov      <- as.numeric(get_arg("--c.ov", "0.0005"))
k.ov      <- as.numeric(get_arg("--k.ov", "0.25"))
k         <- as.numeric(get_arg("--k", "0.03"))
tau       <- as.numeric(get_arg("--tau", "0.02"))
eta.null  <- as.numeric(get_arg("--eta.null", "8500000"))

suppressPackageStartupMessages(library(nixmass))

df <- read.csv(in_path)
df$date <- as.POSIXct(df$date)
df$hs <- as.numeric(df$hs)
df <- df[order(df$date), ]
df$hs[is.na(df$hs) | df$hs < 0] <- 0

model_opts <- list(
  rho.max = rho.max,
  rho.null = rho.null,
  c.ov = c.ov,
  k.ov = k.ov,
  k = k,
  tau = tau,
  eta.null = eta.null
)



res <- swe.delta.snow(
  data = df,
  model_opts = model_opts,
  layer = FALSE
)

out <- data.frame(date = df$date, hs = df$hs, swe_mod = as.numeric(res))
write.csv(out, stdout(), row.names = FALSE)


library(nixmass)
