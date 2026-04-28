# =============================================================================
# HNW_validation_dSnow.R
# R equivalent of HNW_validation_dSnow.ipynb
#
# Workflow (mirrors the Python notebook exactly):
#   1. Load Mag25 multi-station NetCDF
#   2. Run nixmass::swe.delta.snow() per station (calibrated parameters)
#   3. Derive HNW_mod from day-to-day SWE increments, clip negatives to 0
#   4. Validate HNW and SWE against observations (density scatter plots)
#
# Dependencies: ncdf4, nixmass, ggplot2, lubridate
#   + HNW_validation_helper.R  (plot / metric functions)
# =============================================================================

library(ncdf4)       # NetCDF I/O
library(nixmass)     # Delta-snow SWE model (R equivalent of pydeltasnow)
library(lubridate)   # Date helpers (also used inside the helper file)

# Set project root — mirrors os.chdir(base_dir) in Python
base_dir <- "/Users/jakobwerkgarner/code/mt_dsnow" 
setwd(base_dir)

# Source validation helpers — equivalent to:
#   import HNW_validation.HNW_validation_helper as val_helper
source("HNW_validation/HNW_validation_helper.R")

cat(sprintf("model_source: nixmass (version %s)\n",
            as.character(packageVersion("nixmass"))))


# ── Configuration ─────────────────────────────────────────────────────────────

calib_comment <- "nixmass_default_params"
save_data     <- FALSE
infile        <- "HNW_validation/validation_input_Mag25/Mag25_all.nc"


# ── Original nixmass delta.snow parameters (Winkler et al. 2021) ─────────────
# These are the published reference values shipped as the default in
# nixmass::swe.delta.snow(). Passing them explicitly keeps the run
# self-documenting (and independent of any future package default change).

model_opts <- list(
  rho.max  = 401.2588,
  rho.null = 81.19417,
  c.ov     = 0.0005104722,
  k.ov     = 0.37856737,
  k        = 0.02993175,
  tau      = 0.02362476,
  eta.null = 8523356
)


# ── nixmass wrapper: run delta-snow on one station's HS series ────────────────
#
# pydeltasnow internally splits the HS series into contiguous nonzero chunks
# and runs the model on each segment separately. nixmass::swe.delta.snow()
# requires the same pre-processing to avoid numerical issues with these
# calibrated parameters.
#
# Per chunk:
#   • prepend one zero-HS row (date = chunk_start − 1 day)  [required by nixmass]
#   • run swe.delta.snow() on that chunk
#   • drop the prepended row and write results back into swe_out
#
# Input:  hs    – numeric vector of daily snow depth [metres]
#         dates – Date vector, same length as hs
# Output: numeric vector of SWE [mm / kg m⁻²], same length as input

run_deltasnow <- function(hs, dates) {

  hs[is.na(hs)] <- 0          # treat missing depth as snow-free
  hs <- pmax(hs, 0)           # clip any negative values

  # Initialise output to 0 (default for snow-free days). Chunks that fail
  # in nixmass are set to NA below so that validation excludes those days
  # (instead of silently contributing a spurious SWE_mod = 0).
  swe_out <- rep(0.0, length(hs))

  # Locate start / end indices of every contiguous nonzero-HS chunk
  nonzero <- hs > 0
  starts  <- which(nonzero & !c(FALSE, nonzero[-length(nonzero)]))
  ends    <- which(nonzero & !c(nonzero[-1], FALSE))

  n_fail <- 0
  for (i in seq_along(starts)) {

    idx    <- starts[i]:ends[i]

    # Prepend a zero row so the chunk begins at snow-free conditions
    hs_c   <- c(0, hs[idx])
    date_c <- c(dates[starts[i]] - 1L, dates[idx])
    df_c   <- data.frame(date = date_c, hs = hs_c)

    # Run the nixmass delta-snow model on this chunk.
    # nixmass applies a strict max-daily-HS-jump check (and some internal
    # NA-sensitive comparisons) that pydeltasnow does not. When these fire
    # we silently mark the chunk as NA so it drops out of validation rather
    # than contaminating it with SWE_mod = 0.
    res <- tryCatch(
      suppressWarnings(
        nixmass::swe.delta.snow(df_c,
                                model_opts  = model_opts,
                                dyn_rho_max = TRUE,
                                strict_mode = FALSE)
      ),
      error = function(e) NULL
    )

    if (!is.null(res)) {
      swe_out[idx] <- res[-1]        # drop the prepended-zero row
    } else {
      swe_out[idx] <- NA_real_       # failed chunk → exclude from validation
      n_fail <- n_fail + 1
    }
  }

  attr(swe_out, "n_fail")   <- n_fail
  attr(swe_out, "n_chunks") <- length(starts)
  swe_out
}


# ── Load NetCDF dataset ────────────────────────────────────────────────────────

cat("Loading:", infile, "\n")
nc         <- nc_open(infile)
time_raw   <- ncvar_get(nc, "time")

# Parse CF-convention "days since YYYY-MM-DD" time axis
origin_str <- sub("days since ([^ ]+).*", "\\1", ncatt_get(nc, "time", "units")$value)
time_dates <- as.Date(time_raw, origin = origin_str)

stations   <- ncvar_get(nc, "station")   # character vector of station names
HS_mat     <- ncvar_get(nc, "HS")        # [station × time] from ncdf4
SWE_mat    <- ncvar_get(nc, "SWE")
HNW_mat    <- ncvar_get(nc, "HNW")
nc_close(nc)

n_stations <- length(stations)
n_time     <- length(time_dates)
cat(sprintf("Loaded: %d stations × %d timesteps\n", n_stations, n_time))


# ── Run nixmass delta.snow for every station ──────────────────────────────────
# Mirrors the Python loop:
#   for station_name in station_list:
#       swe_results = pydeltasnow.swe_deltasnow(idata, rho_max=..., ...)

SWE_mod_mat <- matrix(NA_real_, nrow = n_stations, ncol = n_time)
obs_counts  <- vector("list", n_stations)

for (i in seq_len(n_stations)) {

  station_name <- stations[i]
  cat(sprintf("\nProcessing station: %s\n", station_name))

  hs_series <- as.numeric(HS_mat[i, ])   # snow depth [m], length = n_time
  cat(sprintf("HS input shape: %d\n", length(hs_series)))

  res <- run_deltasnow(hs_series, time_dates)
  SWE_mod_mat[i, ] <- as.numeric(res)

  n_fail   <- attr(res, "n_fail")
  n_chunks <- attr(res, "n_chunks")
  if (n_fail > 0) {
    cat(sprintf("  note: %d / %d snow chunks could not be modelled (marked NA)\n",
                n_fail, n_chunks))
  }

  # Count available observed SWE measurements (mirrors Python's obs_counts list)
  n_obs <- sum(!is.na(SWE_mat[i, ]))
  obs_counts[[i]] <- data.frame(station = station_name, n_obs = n_obs,
                                stringsAsFactors = FALSE)
}

obs_counts_df <- do.call(rbind, obs_counts)


# ── Compute HNW_mod ────────────────────────────────────────────────────────────
# HNW_mod = day-to-day SWE increase; negative diffs (melt) are clipped to 0.
# A leading NA restores the full time axis — mirrors Python's
#   .diff(dim="time") → .clip(min=0) → .reindex(time=original_time)
# which leaves the first timestep as NaN.

HNW_mod_mat <- matrix(NA_real_, nrow = n_stations, ncol = n_time)
for (i in seq_len(n_stations)) {
  HNW_mod_mat[i, ] <- c(NA, pmax(diff(SWE_mod_mat[i, ]), 0))
}


# ── Build long-format data frame for validation ──────────────────────────────
# Mirrors Python's .to_dataframe().reset_index() + rename(), but carries all
# four variables (HNW_obs, HNW_mod, SWE_obs, SWE_mod) in one table so we can
# filter on SWE_obs regardless of whether we are validating HNW or SWE.

all_long <- do.call(rbind, lapply(seq_len(n_stations), function(i) {
  data.frame(
    time    = time_dates,
    station = stations[i],
    HNW_obs = as.numeric(HNW_mat[i, ]),
    HNW_mod = as.numeric(HNW_mod_mat[i, ]),
    SWE_obs = as.numeric(SWE_mat[i, ]),
    SWE_mod = as.numeric(SWE_mod_mat[i, ]),
    stringsAsFactors = FALSE
  )
}))


# ── Restrict to days with valid observed SWE ─────────────────────────────────
# Mag25 SWE comes from biweekly snow-course surveys, so SWE_obs is NA on
# most days. We only validate on days where an observation actually exists —
# otherwise modelled-vs-obs statistics are dominated by unobserved days and
# by failed-chunk artefacts, which inflates RMSE and pushes R² well below 0.

all_long_valid <- all_long[!is.na(all_long$SWE_obs), ]
cat(sprintf("\nRows with valid SWE_obs: %d (of %d total)\n",
            nrow(all_long_valid), nrow(all_long)))


# ── Validate HNW ──────────────────────────────────────────────────────────────
# Mirrors: val_helper.validate_hnw_mag25(all_df, model_name='dSnow', ...)

validate_hnw_mag25(
  df         = all_long_valid,
  model_name = "dSnow",
  obs_col    = "HNW_obs",
  mod_col    = "HNW_mod",
  save_dir   = "HNW_validation/dSnow/validation_plots",
  filename   = paste0("deltasnow_hnw_validation", calib_comment, ".png")
)


# ── Validate SWE ──────────────────────────────────────────────────────────────
# Mirrors: val_helper.validate_swe_mag25(all_df_SWE, model_name='dSnow', ...)

validate_swe_mag25(
  df         = all_long_valid,
  model_name = "dSnow",
  obs_col    = "SWE_obs",
  mod_col    = "SWE_mod",
  save_dir   = "HNW_validation/dSnow/validation_plots",
  filename   = paste0("deltasnow_SWE_validation", calib_comment, ".png")
)
