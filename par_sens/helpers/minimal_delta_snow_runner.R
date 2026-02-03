#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  args <- commandArgs(trailingOnly = TRUE)
}))

# ----------------------------------------
# CLI parser
# ----------------------------------------
get_arg <- function(flag, default=NULL) {
  hit <- which(args %in% flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop(paste("Missing value for", flag))
  args[hit + 1]
}

# ----------------------------------------
# Required input file
# ----------------------------------------
in_path <- get_arg("--in")
if (is.null(in_path)) stop("Missing input path.", call.=FALSE)

# ----------------------------------------
# Model parameters
# ----------------------------------------
rho.max  <- as.numeric(get_arg("--rho.max",  "400"))
rho.null <- as.numeric(get_arg("--rho.null", "80"))
c.ov     <- as.numeric(get_arg("--c.ov",     "0.0005"))
k.ov     <- as.numeric(get_arg("--k.ov",     "0.25"))
k        <- as.numeric(get_arg("--k",        "0.03"))
tau      <- as.numeric(get_arg("--tau",      "0.02"))
eta.null <- as.numeric(get_arg("--eta.null", "8500000"))

suppressPackageStartupMessages(library(nixmass))

# ----------------------------------------
# Load data
# ----------------------------------------
df <- read.csv(in_path)
df$date <- as.POSIXct(df$date)
df$hs <- as.numeric(df$hs)
df <- df[order(df$date), ]

# Replace invalid HS
df$hs[is.na(df$hs) | df$hs < 0] <- 0

# ----------------------------------------
# Prepare model
# ----------------------------------------
model_opts <- list(
  rho.max  = rho.max,
  rho.null = rho.null,
  c.ov     = c.ov,
  k.ov     = k.ov,
  k        = k,
  tau      = tau,
  eta.null = eta.null
)

# ----------------------------------------
# Run Δ-SNOW
# ----------------------------------------
res <- try(
  suppressWarnings(
    swe.delta.snow(
      data = df,
      model_opts = model_opts,
      layer = FALSE,
      dyn_rho_max = FALSE
    )
  ),
  silent = TRUE
)

# If the model fails → output NA values as CSV
if (inherits(res, "try-error")) {
  out <- data.frame(date=df$date, hs=df$hs, swe_mod=NA)
  write.csv(out, stdout(), row.names=FALSE)
  quit(save="no", status=0)
}

# ----------------------------------------
# Output clean CSV to stdout (mandatory)
# ----------------------------------------
out <- data.frame(
  date    = df$date,
  hs      = df$hs,
  swe_mod = as.numeric(res)
)

write.csv(out, stdout(), row.names=FALSE)