# run_dsnow_and_save_layers_nc_REVISED_season_Nov1_Apr30.R
# Only ΔSnow (nixmass::swe.delta.snow). Writes ONE NetCDF per station containing ALL seasons.
# Saves only seasons where ΔSnow has > 3 layers.
#
# NetCDF structure (dim order):
#   season, dos (Day Of Season), layer
#
# Variables stored as:
#   HS, SWE_layers, AGE, RHO, OVB : (season, dos, layer)
#   SWE_total, HS_meas_m          : (season, dos)

suppressPackageStartupMessages({
  library(nixmass)
  library(ncdf4)
  library(lubridate)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------
# SETTINGS
# -----------------------------
indir  <- "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_data/output/HS_SWE_by_station/"
outdir <- "/Users/jakobwerkgarner/code/mt_dsnow/model_diff/R_comparisson/output/HS_SWE_by_station_dsnow_only/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

MIN_LAYERS_TO_SAVE  <- 4L   # "more than 3"
MIN_DAYS_PER_SEASON <- 10L

# -----------------------------
# Helpers
# -----------------------------
is_valid_scalar_posint <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x >= 1
}

as_num_matrix <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.matrix(x)) { storage.mode(x) <- "double"; return(x) }
  if (is.data.frame(x)) { m <- as.matrix(x); storage.mode(m) <- "double"; return(m) }
  if (is.list(x) && !is.data.frame(x)) {
    m <- tryCatch(do.call(cbind, x), error = function(e) NULL)
    if (!is.null(m)) { m <- as.matrix(m); storage.mode(m) <- "double"; return(m) }
  }
  if (is.atomic(x) && !is.list(x)) {
    m <- matrix(as.numeric(x), nrow = length(x), ncol = 1)
    storage.mode(m) <- "double"
    return(m)
  }
  NULL
}

ensure_layer_time <- function(x, Tn) {
  if (is.null(x)) return(matrix(NA_real_, nrow = 0, ncol = Tn))
  if (is.vector(x) && !is.list(x)) x <- matrix(as.numeric(x), nrow = 1)
  if (!is.matrix(x)) x <- as.matrix(x)

  Tn <- as.integer(Tn)
  if (length(Tn) != 1L || is.na(Tn) || Tn < 1L) stop("ensure_layer_time: invalid Tn")

  if (ncol(x) != Tn && nrow(x) == Tn) x <- t(x)

  if (ncol(x) < Tn) {
    x <- cbind(x, matrix(NA_real_, nrow = nrow(x), ncol = Tn - ncol(x)))
  } else if (ncol(x) > Tn) {
    x <- x[, seq_len(Tn), drop = FALSE]
  }

  storage.mode(x) <- "double"
  x
}

pad_to_layers <- function(x, target_layers, fill = 0) {
  if (!is.matrix(x)) x <- as.matrix(x)
  target_layers <- as.integer(target_layers)
  if (length(target_layers) != 1L || is.na(target_layers) || target_layers < 1L)
    stop("pad_to_layers: invalid target_layers")

  if (nrow(x) < target_layers) {
    x <- rbind(x, matrix(fill, nrow = target_layers - nrow(x), ncol = ncol(x)))
  } else if (nrow(x) > target_layers) {
    x <- x[seq_len(target_layers), , drop = FALSE]
  }
  storage.mode(x) <- "double"
  x
}

safe_density <- function(SWE_layers, H_layers) {
  if (is.null(SWE_layers) || is.null(H_layers)) return(NULL)
  SWE_layers <- as.matrix(SWE_layers)
  H_layers   <- as.matrix(H_layers)

  Tn <- min(ncol(SWE_layers), ncol(H_layers))
  SWE_layers <- SWE_layers[, seq_len(Tn), drop = FALSE]
  H_layers   <- H_layers[, seq_len(Tn), drop = FALSE]

  H_layers[H_layers <= 0] <- NA_real_
  rho <- SWE_layers / H_layers
  rho[!is.finite(rho)] <- NA_real_
  storage.mode(rho) <- "double"
  rho
}

safe_overburden <- function(SWE_layers) {
  if (is.null(SWE_layers)) return(NULL)
  SWE_layers <- as.matrix(SWE_layers)
  out <- apply(SWE_layers, 2, function(col) {
    col[is.na(col)] <- 0
    rev(cumsum(rev(col)))
  })
  out <- as.matrix(out)
  storage.mode(out) <- "double"
  out
}

# Season label = year of April (end year): Nov/Dec belong to next year's season
df_season_label <- function(d) {
  y <- year(d); m <- month(d)
  ifelse(m >= 11, y + 1L, y)
}

in_season_window <- function(d, season_label) {
  start <- as.Date(sprintf("%d-11-01", season_label - 1L))
  end   <- as.Date(sprintf("%d-04-30", season_label))
  d >= start & d <= end
}

# -----------------------------
# Helpers for DOS (Nov 1 -> Apr 30)
# -----------------------------
season_start_end <- function(season_label) {
  start <- as.Date(sprintf("%d-11-01", season_label - 1L))
  end   <- as.Date(sprintf("%d-04-30", season_label))
  list(start = start, end = end)
}

season_dos_index <- function(dates, season_label) {
  se <- season_start_end(season_label)
  as.integer(as.Date(dates) - se$start) + 1L
}

season_dos_length <- function(season_label) {
  se <- season_start_end(season_label)
  as.integer(se$end - se$start) + 1L
}

# -----------------------------
# NetCDF writer: ONE file per station
# dims: season, dos, layer
# -----------------------------
write_station_dsnow_nc <- function(ncfile, season_labels, ndos, nlayer,
                                   HS, SWE_layers, AGE, RHO, OVB,
                                   SWE_total, HS_meas_m) {

  season_labels <- as.integer(season_labels)
  Sn <- length(season_labels)
  Dn <- as.integer(ndos)
  Ln <- as.integer(nlayer)

  if (Sn < 1 || Dn < 1 || Ln < 1) stop("write_station_dsnow_nc: invalid dims")

  # dims (order matters)
  dim_season <- ncdf4::ncdim_def("season", units = "year",          vals = season_labels, create_dimvar = TRUE)
  dim_dos    <- ncdf4::ncdim_def("dos",    units = "day_of_season", vals = seq_len(Dn),  create_dimvar = TRUE)
  dim_layer  <- ncdf4::ncdim_def("layer",  units = "1",            vals = seq_len(Ln),   create_dimvar = TRUE)

  # vars: (season, dos, layer)
  vHS    <- ncdf4::ncvar_def("HS",         "m", list(dim_season, dim_dos, dim_layer), missval = NA_real_, prec = "double")
  vSWE_L <- ncdf4::ncvar_def("SWE_layers", "1", list(dim_season, dim_dos, dim_layer), missval = NA_real_, prec = "double")
  vAGE   <- ncdf4::ncvar_def("AGE",        "d", list(dim_season, dim_dos, dim_layer), missval = NA_real_, prec = "double")
  vRHO   <- ncdf4::ncvar_def("RHO",        "1", list(dim_season, dim_dos, dim_layer), missval = NA_real_, prec = "double")
  vOVB   <- ncdf4::ncvar_def("OVB",        "1", list(dim_season, dim_dos, dim_layer), missval = NA_real_, prec = "double")

  # vars: (season, dos)
  vSWE_T   <- ncdf4::ncvar_def("SWE_total", "1", list(dim_season, dim_dos), missval = NA_real_, prec = "double")
  vHS_meas <- ncdf4::ncvar_def("HS_meas_m", "m", list(dim_season, dim_dos), missval = NA_real_, prec = "double")

  nc <- ncdf4::nc_create(ncfile, vars = list(vHS, vSWE_L, vAGE, vRHO, vOVB, vSWE_T, vHS_meas))
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  # write arrays (must match dim order: season, dos, layer)
  ncdf4::ncvar_put(nc, vHS,        HS)
  ncdf4::ncvar_put(nc, vSWE_L,     SWE_layers)
  ncdf4::ncvar_put(nc, vAGE,       AGE)
  ncdf4::ncvar_put(nc, vRHO,       RHO)
  ncdf4::ncvar_put(nc, vOVB,       OVB)
  ncdf4::ncvar_put(nc, vSWE_T,     SWE_total)
  ncdf4::ncvar_put(nc, vHS_meas,   HS_meas_m)

  # metadata for reconstructing dates from (season, dos)
  ncdf4::ncatt_put(nc, 0, "season_window_start_mmdd", "11-01")
  ncdf4::ncatt_put(nc, 0, "season_window_end_mmdd",   "04-30")
  ncdf4::ncatt_put(nc, 0, "dos_definition",
                   "dos = 1 corresponds to Nov 1 of (season-1); dos increments daily until Apr 30 of (season).")

  invisible(TRUE)
}

# -----------------------------
# MAIN
# -----------------------------
csv_files <- list.files(indir, pattern = "\\.csv$", full.names = TRUE)

for (infile in csv_files) {

  df <- read.csv(infile, stringsAsFactors = FALSE)

  date_col <- intersect(c("date", "timestamp"), names(df))[1]
  if (is.na(date_col)) stop("No date/timestamp column in ", infile)

  df$date   <- as.Date(df[[date_col]])
  df$season <- df_season_label(df$date)

  station <- sub("_hs_swe_obs.*", "", basename(infile))

  seasons <- sort(unique(df$season))
  per_season <- list()

  for (this_season in seasons) {

    sdf <- df[df$season == this_season & in_season_window(df$date, this_season), ]
    if (nrow(sdf) < MIN_DAYS_PER_SEASON) next
    sdf <- sdf[order(sdf$date), ]

    timestamp <- as.character(sdf$date)
    hs_vec_m <- as.numeric(sdf$hs %||% sdf$HS_meas %||% sdf$HS %||% NA_real_)
    if (all(is.na(hs_vec_m)) || sum(hs_vec_m, na.rm = TRUE) == 0) next

    # force first entry to 0
    if (length(hs_vec_m) >= 1L) {
      if (!is.finite(hs_vec_m[1]) || hs_vec_m[1] != 0) hs_vec_m[1] <- 0
    }

    Tn <- length(timestamp)
    if (Tn < MIN_DAYS_PER_SEASON) next

    dsnow_df <- data.frame(hs = hs_vec_m, date = timestamp, stringsAsFactors = FALSE)

    # ΔSnow
    fmls <- names(formals(nixmass::swe.delta.snow))
    if (!("layers" %in% fmls)) stop("Your nixmass::swe.delta.snow does not support layers=TRUE.")

    ds_out <- tryCatch({
      if ("strict_mode" %in% fmls) {
        nixmass::swe.delta.snow(dsnow_df, layers = TRUE, strict_mode = TRUE, verbose = FALSE, dyn_rho_max = FALSE)
      } else {
        nixmass::swe.delta.snow(dsnow_df, layers = TRUE, verbose = FALSE, dyn_rho_max = FALSE)
      }
    }, error = function(e) e)

    if (inherits(ds_out, "error") || !is.list(ds_out) ||
        is.null(ds_out$h) || is.null(ds_out$swe) || is.null(ds_out$age)) {
      message(sprintf("Skipping %s season %s: delta.snow invalid output.", station, this_season))
      next
    }

    H_layers   <- ensure_layer_time(as_num_matrix(ds_out$h),   Tn)
    SWE_layers <- ensure_layer_time(as_num_matrix(ds_out$swe), Tn)
    AGE_layers <- ensure_layer_time(as_num_matrix(ds_out$age), Tn)

    ds_layers_max <- max(nrow(H_layers), nrow(SWE_layers), nrow(AGE_layers))
    if (!is_valid_scalar_posint(ds_layers_max) || ds_layers_max < MIN_LAYERS_TO_SAVE) {
      message(sprintf("Skipping %s season %s: delta.snow layers = %s (need >= %d).",
                      station, this_season,
                      ifelse(is.na(ds_layers_max), "NA", ds_layers_max),
                      MIN_LAYERS_TO_SAVE))
      next
    }

    # pad to ds_layers_max
    H_layers   <- pad_to_layers(H_layers,   ds_layers_max, fill = 0)
    SWE_layers <- pad_to_layers(SWE_layers, ds_layers_max, fill = 0)
    AGE_layers <- pad_to_layers(AGE_layers, ds_layers_max, fill = 0)

    SWE_total <- ds_out$SWE
    if (is.null(SWE_total)) SWE_total <- colSums(SWE_layers, na.rm = TRUE)
    SWE_total <- as.numeric(SWE_total)
    if (length(SWE_total) != Tn) SWE_total <- rep_len(SWE_total, Tn)

    RHO_layers <- safe_density(SWE_layers, H_layers)
    OVB_layers <- safe_overburden(SWE_layers)

    per_season[[as.character(this_season)]] <- list(
      timestamp  = timestamp,
      HS_meas_m  = hs_vec_m,
      HS         = H_layers,     # (layer, time)
      SWE_layers = SWE_layers,   # (layer, time)
      AGE        = AGE_layers,   # (layer, time)
      RHO        = RHO_layers,   # (layer, time)
      OVB        = OVB_layers,   # (layer, time)
      SWE_total  = SWE_total,    # (time)
      nlayer     = ds_layers_max
    )

    message("OK: ", station, " season ", this_season, " (layers=", ds_layers_max, ")")
  }

  if (length(per_season) == 0L) {
    message("No valid seasons for station: ", station)
    next
  }

  seasons_kept <- as.integer(names(per_season))
  S_all <- length(seasons_kept)
  if (S_all < 1L) {
    message("No valid seasons for station: ", station)
    next
  }

  # fixed DOS length for Nov 1 -> Apr 30
  D_all <- season_dos_length(seasons_kept[1])

  # max layers across kept seasons
  L_all <- max(vapply(per_season, function(x) x$nlayer, integer(1)))

  # allocate arrays: (season, dos, layer)
  HS_3d   <- array(NA_real_, dim = c(S_all, D_all, L_all))
  SWE_3d  <- array(NA_real_, dim = c(S_all, D_all, L_all))
  AGE_3d  <- array(NA_real_, dim = c(S_all, D_all, L_all))
  RHO_3d  <- array(NA_real_, dim = c(S_all, D_all, L_all))
  OVB_3d  <- array(NA_real_, dim = c(S_all, D_all, L_all))

  # allocate 2D arrays: (season, dos)
  SWE_tot <- array(NA_real_, dim = c(S_all, D_all))
  HS_meas <- array(NA_real_, dim = c(S_all, D_all))

  # fill arrays season by season
  for (sidx in seq_along(seasons_kept)) {

    ss <- as.character(seasons_kept[sidx])
    obj <- per_season[[ss]]
    this_season <- seasons_kept[sidx]

    dos_idx <- season_dos_index(obj$timestamp, this_season)

    keep <- which(is.finite(dos_idx) & dos_idx >= 1L & dos_idx <= D_all)
    if (length(keep) == 0L) next

    dos_idx <- dos_idx[keep]
    Ls <- obj$nlayer

    # obj$HS etc are (layer, time); subset time columns then transpose to (time, layer)
    HS_tmp  <- t(obj$HS[, keep, drop = FALSE])           # (time, layer)
    SWE_tmp <- t(obj$SWE_layers[, keep, drop = FALSE])   # (time, layer)
    AGE_tmp <- t(obj$AGE[, keep, drop = FALSE])          # (time, layer)

    HS_3d[sidx, dos_idx,  seq_len(Ls)] <- HS_tmp[,  seq_len(Ls), drop = FALSE]
    SWE_3d[sidx, dos_idx, seq_len(Ls)] <- SWE_tmp[, seq_len(Ls), drop = FALSE]
    AGE_3d[sidx, dos_idx, seq_len(Ls)] <- AGE_tmp[, seq_len(Ls), drop = FALSE]

    if (!is.null(obj$RHO)) {
      RHO_tmp <- t(obj$RHO[, keep, drop = FALSE])
      RHO_3d[sidx, dos_idx, seq_len(Ls)] <- RHO_tmp[, seq_len(Ls), drop = FALSE]
    }
    if (!is.null(obj$OVB)) {
      OVB_tmp <- t(obj$OVB[, keep, drop = FALSE])
      OVB_3d[sidx, dos_idx, seq_len(Ls)] <- OVB_tmp[, seq_len(Ls), drop = FALSE]
    }

    SWE_tot[sidx, dos_idx] <- obj$SWE_total[keep]
    HS_meas[sidx, dos_idx] <- obj$HS_meas_m[keep]
  }

  out_nc <- file.path(outdir, paste0(station, "_dsnow_allseasons.nc"))
  write_station_dsnow_nc(
    ncfile        = out_nc,
    season_labels = seasons_kept,
    ndos          = D_all,
    nlayer        = L_all,
    HS            = HS_3d,
    SWE_layers    = SWE_3d,
    AGE           = AGE_3d,
    RHO           = RHO_3d,
    OVB           = OVB_3d,
    SWE_total     = SWE_tot,
    HS_meas_m     = HS_meas
  )

  message("Wrote station NetCDF: ", out_nc)
}

message("Done.")