#!/usr/bin/env Rscript
# =============================================================================
# Script: run_delta_snow.R
# Purpose:
#   Run Δ-SNOW (swe.delta.snow) on HS data with fully configurable parameters.
#
# Input:
#   A CSV with columns:
#     - date : timestamp (parseable; any common format)
#     - hs   : total snow height in METERS (non-negative)
#
# Output:
#   A CSV with columns:
#     - date, hs, SWE
#   If --layer true is used and the function returns layerwise SWE, an
#   additional CSV is written (wide or long; see below).
#
# How the model code is loaded:
#   1) If --library nixmass (default): library(nixmass)
#   2) Else if --project <path>: devtools::load_all(<path>)
#
# CLI usage:
#   Rscript run_delta_snow.R \
#     --in hs.csv \
#     --out swe.csv \
#     [--layers-out layers.csv] \
#     [--library nixmass] [--project /path/to/Rproj] \
#     [--tz Europe/Vienna] [--pad-before-start true] [--pad-days 1] \
#     [--rho.max 401.2588] [--rho.null 81.19417] [--c.ov 0.0005104722] \
#     [--k.ov 0.37856737] [--k 0.02993175] [--tau 0.02362476] [--eta.null 8523356] \
#     [--layer false]
#
# Notes:
#   - If your HS is in centimeters, convert to meters before calling.
#   - This script pads one day with hs=0 before the first record by default
#     (Δ-SNOW requirement). Disable with --pad-before-start false.
# =============================================================================

suppressWarnings(suppressMessages({
  args <- commandArgs(trailingOnly = TRUE)
}))

# ---------- CLI helper ----------
get_arg <- function(flag, default = NULL, has_value = TRUE) {
  hit <- which(args %in% flag)
  if (length(hit) == 0) return(default)
  if (!has_value) return(TRUE)
  if (hit == length(args)) stop(paste("Missing value for", flag))
  args[hit + 1]
}
as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("true","t","1","yes","y")
}

show_help <- isTRUE(get_arg("--help", has_value = FALSE)) || isTRUE(get_arg("-h", has_value = FALSE))
if (show_help) {
  cat("Usage: see header of run_delta_snow.R\n")
  quit(status = 0)
}

# ---------- Required / optional args ----------
in_path        <- get_arg(c("--in","-i"))
out_path       <- get_arg(c("--out","-o"), "swe.csv")
layers_out     <- get_arg("--layers-out", NULL)

lib_name       <- get_arg("--library", "nixmass")
proj_path      <- get_arg("--project", NULL)

tz_arg         <- get_arg("--tz", "Europe/Vienna")
pad_before     <- as_bool(get_arg("--pad-before-start", "true"), default = TRUE)
pad_days       <- suppressWarnings(as.integer(get_arg("--pad-days", "1")))

rho.max        <- suppressWarnings(as.numeric(get_arg("--rho.max", "401.2588")))
rho.null       <- suppressWarnings(as.numeric(get_arg("--rho.null", "81.19417")))
c.ov           <- suppressWarnings(as.numeric(get_arg("--c.ov", "0.0005104722")))
k.ov           <- suppressWarnings(as.numeric(get_arg("--k.ov", "0.37856737")))
k              <- suppressWarnings(as.numeric(get_arg("--k", "0.02993175")))
tau            <- suppressWarnings(as.numeric(get_arg("--tau", "0.02362476")))
eta.null       <- suppressWarnings(as.numeric(get_arg("--eta.null", "8523356")))
layer_flag     <- as_bool(get_arg("--layer", "false"), default = FALSE)

if (is.null(in_path)) stop("Please provide --in <hs.csv>")

# ---------- Load model code ----------
suppressWarnings(suppressMessages({
  ok <- FALSE
  if (!is.null(lib_name)) {
    ok <- require(lib_name, quietly = TRUE, character.only = TRUE)
  }
  if (!ok && !is.null(proj_path)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      stop("Could not load '", lib_name, "' and 'devtools' not available to load project: ", proj_path)
    }
    devtools::load_all(proj_path)
    ok <- TRUE
  }
  if (!ok && is.null(lib_name) && is.null(proj_path)) {
    stop("No model code loaded: set --library nixmass or --project /path/to/project")
  }
  if (!exists("swe.delta.snow")) {
    stop("Function 'swe.delta.snow' not found after loading. Check --library / --project.")
  }
}))

# ---------- Read input ----------
read_csv_robust <- function(path) {
  df <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || ncol(df) == 1) {
    df <- tryCatch(read.csv2(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  }
  if (is.null(df)) stop("Could not read input CSV: ", path)
  df
}
df <- read_csv_robust(in_path)
if (!all(c("date","hs") %in% names(df))) {
  stop("Input must have columns 'date' and 'hs' (hs in meters).")
}

# Parse time (be generous with formats)
parse_time <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXt")) return(x)
  x <- trimws(as.character(x))
  fmts <- c(
    "%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S",
    "%d.%m.%Y %H:%M:%S", "%d.%m.%Y %H:%M",
    "%Y-%m-%d %H:%M",    "%Y/%m/%d %H:%M",
    "%Y-%m-%d",          "%d.%m.%Y",      "%m/%d/%Y"
  )
  for (f in fmts) {
    tt <- as.POSIXct(x, format = f, tz = tz)
    if (all(!is.na(tt))) return(tt)
  }
  tt <- as.POSIXct(x, tz = tz)
  if (any(is.na(tt))) stop("Failed to parse some 'date' values. Use an ISO-like format.")
  tt
}

df$date <- parse_time(df$date, tz = tz_arg)
df$hs   <- suppressWarnings(as.numeric(gsub(",", ".", df$hs)))
df <- df[order(df$date), ]

# Enforce non-negative numeric hs
df$hs[is.na(df$hs) | df$hs < 0] <- 0

# Optional: pad one day with hs=0 before first record
if (isTRUE(pad_before) && nrow(df) > 0 && is.finite(pad_days) && pad_days > 0) {
  first_ts <- min(df$date, na.rm = TRUE)
  pad_ts   <- first_ts - pad_days * 24 * 3600
  pad_row  <- data.frame(date = pad_ts, hs = 0.0)
  df <- rbind(pad_row, df)
  df <- df[order(df$date), ]
}

# Δ-Snow prefers evenly spaced series (often daily). We do not resample here;
# handle resampling upstream if needed.

# ---------- Build call safely ----------
# Accept only args the function actually supports
formal_nms <- names(formals(swe.delta.snow))
arglist <- list(data = df) # in many impls, swe.delta.snow(data.frame(date,hs)) works
# Some versions require named args; we add the recognized ones:
maybe_add <- function(nm, val) {
  if (nm %in% formal_nms) arglist[[nm]] <<- val
}
maybe_add("rho.max",  rho.max)
maybe_add("rho.null", rho.null)
maybe_add("c.ov",     c.ov)
maybe_add("k.ov",     k.ov)
maybe_add("k",        k)
maybe_add("tau",      tau)
maybe_add("eta.null", eta.null)
maybe_add("layer",    layer_flag)

# Some versions want hs/time vectors rather than a 'data' frame.
# If 'data' is not in formals but 'hs' and 'date' are, adapt:
if (!("data" %in% formal_nms) && all(c("hs","date") %in% names(df))) {
  if ("hs" %in% formal_nms)   arglist$hs   <- df$hs
  if ("date" %in% formal_nms) arglist$date <- df$date
}

# ---------- Run model ----------
res <- tryCatch({
  do.call(swe.delta.snow, arglist)
}, error = function(e) {
  stop("swe.delta.snow failed: ", e$message)
})

# ---------- Handle output & write ----------
# Common cases:
#  - numeric vector SWE with length nrow(df)
#  - list containing SWE and maybe layerwise outputs
write_main <- function(out_path, df, SWE) {
  out <- data.frame(date = df$date, hs = df$hs, SWE = as.numeric(SWE))
  write.csv(out, out_path, row.names = FALSE)
}

if (is.numeric(res) && length(res) == nrow(df)) {
  write_main(out_path, df, res)

} else if (is.list(res)) {
  # Try common names
  if (!is.null(res$SWE) && length(res$SWE) == nrow(df)) {
    write_main(out_path, df, res$SWE)
  } else if (!is.null(res$swe) && length(res$swe) == nrow(df)) {
    write_main(out_path, df, res$swe)
  } else {
    stop("Unrecognized result structure from swe.delta.snow.")
  }

  # Optional layerwise output if present and requested
  if (isTRUE(layer_flag) && !is.null(layers_out)) {
    # Try to find a layer matrix/time series in the result list
    lyr <- NULL
    cand_names <- c("layers","layer","SWE_layers","swe_layers","swe.layer")
    for (nm in cand_names) if (nm %in% names(res)) { lyr <- res[[nm]]; break }
    if (!is.null(lyr)) {
      # If matrix with rows = time and columns = layers:
      if (is.matrix(lyr) && nrow(lyr) == nrow(df)) {
        lyr_df <- as.data.frame(lyr)
        lyr_df <- cbind(date = df$date, lyr_df)
        write.csv(lyr_df, layers_out, row.names = FALSE)
      } else if (is.data.frame(lyr) && nrow(lyr) == nrow(df)) {
        lyr_df <- cbind(date = df$date, lyr)
        write.csv(lyr_df, layers_out, row.names = FALSE)
      } else {
        warning("Layerwise output present but not in an expected shape; layers CSV not written.")
      }
    } else {
      warning("No recognizable layerwise output in result, cannot write --layers-out.")
    }
  }

} else {
  stop("Unexpected output type from swe.delta.snow.")
}

cat(sprintf("Wrote: %s%s\n",
            normalizePath(out_path, winslash = "/"),
            if (!is.null(layers_out)) paste0(" and ", normalizePath(layers_out, winslash = "/")) else ""))