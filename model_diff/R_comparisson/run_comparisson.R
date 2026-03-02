# run_models_and_save_layers_nc.R
# Runs:
#  - hs2swe() (your adapted MR-returning version)
#  - nixmass::swe.delta.snow(..., layers=TRUE)
# Saves:
#  - CSV comparison (total SWE)
#  - NetCDF for HS2SWE layers
#  - NetCDF for delta.snow layers
#
# NOTES:
# - HS2SWE expects HS in cm -> we convert HS_meas [m] to [cm] for hs2swe only.
# - NetCDF writer expects all matrix variables to be shaped (layer x time).
# - DIA is originally (5 x time). We pad it to (layer x time) with zeros to avoid ncvar_put size errors.
# - delta.snow layer matrices may come as (time x layer); we auto-transpose to (layer x time).

suppressPackageStartupMessages({
  library(nixmass)
  library(ncdf4)
})

# ---- paths (edit if needed) ----
source("model_diff/R_comparisson/HS2SWE_adapted_R.R")

infile  <- "/Users/jakobwerkgarner/code/mt_dsnow/par_sens/SNOWPACK_data_seasons_daily/OABS2_2023_2024_daily_HS_SWE.csv"
outbase <- "model_diff/R_comparisson/output/SNOWPACK_forced/"  # no extension

# -----------------------------
# Helpers
# -----------------------------
ensure_layer_time <- function(M, Tn) {
  # target: (layer x time), so ncol == Tn
  if (!is.matrix(M)) stop("ensure_layer_time expects a matrix.")
  if (ncol(M) == Tn) return(M)
  if (nrow(M) == Tn) return(t(M))
  stop(sprintf("Cannot orient matrix to (layer x time): got %d x %d, expected time=%d",
               nrow(M), ncol(M), Tn))
}

pad_to_layers <- function(M, n_layers, fill = 0) {
  if (!is.matrix(M)) stop("pad_to_layers expects a matrix.")
  if (nrow(M) == n_layers) return(M)
  if (nrow(M) > n_layers) stop("Matrix has more rows than target layer count.")
  out <- matrix(fill, nrow = n_layers, ncol = ncol(M))
  out[1:nrow(M), ] <- M
  out
}

safe_density <- function(swe, h) {
  # swe/h where h>0 and neither is NA; returns matrix with 0 elsewhere
  if (!is.matrix(swe) || !is.matrix(h)) stop("safe_density expects matrices.")
  out <- matrix(0, nrow = nrow(h), ncol = ncol(h))
  nz <- !is.na(h) & !is.na(swe) & (h > 0)
  out[nz] <- swe[nz] / h[nz]
  out
}

safe_overburden <- function(swe) {
  # OVB_i,t = sum_{j=i..end} SWE_{j,t}, treating NA as 0
  if (!is.matrix(swe)) stop("safe_overburden expects a matrix.")
  out <- matrix(0, nrow = nrow(swe), ncol = ncol(swe))
  for (t in seq_len(ncol(swe))) {
    v <- swe[, t]
    v[is.na(v)] <- 0
    out[, t] <- rev(cumsum(rev(v)))
  }
  out
}

write_model_nc <- function(ncfile, time, vars_named_list) {
  # dims
  dim_time  <- ncdim_def("time", "", seq_along(time), unlim = TRUE)

  # layer dimension is taken from first matrix variable found
  first_mat <- NULL
  for (v in vars_named_list) {
    if (is.matrix(v)) { first_mat <- v; break }
  }
  if (is.null(first_mat)) stop("No matrix variables provided to write_model_nc().")
  dim_layer <- ncdim_def("layer", "", seq_len(nrow(first_mat)))

  # define variables
  var_defs <- list()
  for (nm in names(vars_named_list)) {
    v <- vars_named_list[[nm]]
    if (is.matrix(v)) {
      var_defs[[nm]] <- ncvar_def(nm, "", list(dim_layer, dim_time),
                                  missval = NA_real_, prec = "double")
    } else if (is.numeric(v) && length(v) == length(time)) {
      var_defs[[nm]] <- ncvar_def(nm, "", list(dim_time),
                                  missval = NA_real_, prec = "double")
    } else {
      stop(sprintf("Variable '%s' must be matrix(layer x time) or numeric(time).", nm))
    }
  }

  nc <- nc_create(ncfile, var_defs)
  on.exit(nc_close(nc), add = TRUE)

  # store timestamps as attribute
  ncatt_put(nc, 0, "timestamps", paste(time, collapse = ","))

  # write data
  for (nm in names(vars_named_list)) {
    ncvar_put(nc, nm, vars_named_list[[nm]])
  }
}

# -----------------------------
# 1) Read input data
# -----------------------------
df <- read.csv(infile, stringsAsFactors = FALSE)

if (!all(c("timestamp", "HS_meas") %in% names(df))) {
  stop("Input CSV must contain at least columns: timestamp, HS_meas")
}

timestamp <- as.character(df$timestamp)
Tn <- length(timestamp)

hs_vec_m <- as.numeric(df$HS_meas)

# HS2SWE: convert meters -> centimeters (only for hs2swe)
hs_vec_cm <- hs_vec_m * 100

# delta.snow: expects hs in meters and date as character
dsnow_df <- data.frame(
  hs = hs_vec_m,
  date = timestamp,
  stringsAsFactors = FALSE
)

if (!isTRUE(hs_vec_m[1] == 0)) {
  stop("For delta.snow with strict checks, HS_meas[1] must be 0 (meters).")
}

# -----------------------------
# 2) Run HS2SWE (MR output)
# -----------------------------
mr_hs2swe <- hs2swe(hs_vec_cm)



hs2_max_layers <- nrow(mr_hs2swe[[1]][[1]]$MR$OVB)


episodes <- mr_hs2swe[[1]]
if (length(episodes) == 0) stop("HS2SWE returned no episodes (no HS>0).")




# Allocate full-season matrices (layer x time)
HS_hs2  <- matrix(NA_real_, nrow = hs2_max_layers, ncol = Tn)
RHO_hs2 <- matrix(NA_real_, nrow = hs2_max_layers, ncol = Tn)
OVB_hs2 <- matrix(NA_real_, nrow = hs2_max_layers, ncol = Tn)
AGE_hs2 <- matrix(NA_real_, nrow = hs2_max_layers, ncol = Tn)
DIA_hs2 <- matrix(NA_real_, nrow = 5, ncol = Tn)   # will be padded later
SWE_hs2 <- rep(NA_real_, Tn)

for (ep in episodes) {
  idx <- ep$index
  MR  <- ep$MR

  # MR has leading "0-state" column -> drop it
  if (is.null(MR$HS) || !is.matrix(MR$HS) || ncol(MR$HS) < 2) next
  cols <- 2:ncol(MR$HS)

  L <- nrow(MR$HS)
  HS_hs2[1:L, idx]  <- MR$HS[, cols, drop = FALSE]
  RHO_hs2[1:L, idx] <- MR$RHO[, cols, drop = FALSE]
  OVB_hs2[1:L, idx] <- MR$OVB[, cols, drop = FALSE]
  AGE_hs2[1:L, idx] <- MR$AGE[, cols, drop = FALSE]

  DIA_hs2[, idx] <- MR$DIA[, cols, drop = FALSE]
  SWE_hs2[idx]   <- as.numeric(MR$SWE[1, cols])
}

# Set non-snow periods to 0 (optional)
nosnow <- which(hs_vec_m == 0)
HS_hs2[, nosnow]  <- 0
RHO_hs2[, nosnow] <- 0
OVB_hs2[, nosnow] <- 0
AGE_hs2[, nosnow] <- 0
DIA_hs2[, nosnow] <- 0
SWE_hs2[nosnow]   <- 0

# Ensure orientation and pad DIA to layer dimension
HS_hs2  <- ensure_layer_time(HS_hs2,  Tn)
RHO_hs2 <- ensure_layer_time(RHO_hs2, Tn)
OVB_hs2 <- ensure_layer_time(OVB_hs2, Tn)
AGE_hs2 <- ensure_layer_time(AGE_hs2, Tn)

DIA_hs2 <- ensure_layer_time(DIA_hs2, Tn)              # 5 x time
DIA_hs2 <- pad_to_layers(DIA_hs2, nrow(HS_hs2), fill=0)

# -----------------------------
# 3) Run delta.snow with layers=TRUE (nixmass)
# -----------------------------
ds_out <- NULL
fmls <- names(formals(nixmass::swe.delta.snow))

if (!("layers" %in% fmls)) {
  stop("Your installed nixmass::swe.delta.snow does not support layers=TRUE.")
}

if ("strict_mode" %in% fmls) {
  ds_out <- nixmass::swe.delta.snow(dsnow_df, layers = TRUE, strict_mode = TRUE, verbose = FALSE, dyn_rho_max = FALSE)
} else {
  ds_out <- nixmass::swe.delta.snow(dsnow_df, layers = TRUE, verbose = FALSE, dyn_rho_max = FALSE)
}

if (!(is.list(ds_out) && !is.null(ds_out$h) && !is.null(ds_out$swe) && !is.null(ds_out$age))) {
  stop("Unexpected delta.snow output structure. Inspect ds_out.")
}

H_ds_layers   <- ds_out$h
SWE_ds_layers <- ds_out$swe
AGE_ds_layers <- ds_out$age
SWE_ds_total  <- ds_out$SWE
if (is.null(SWE_ds_total)) SWE_ds_total <- colSums(SWE_ds_layers, na.rm = TRUE)

# Orient to (layer x time)
H_ds_layers   <- ensure_layer_time(H_ds_layers,   Tn)
SWE_ds_layers <- ensure_layer_time(SWE_ds_layers, Tn)
AGE_ds_layers <- ensure_layer_time(AGE_ds_layers, Tn)

# Derived variables
RHO_ds_layers <- safe_density(SWE_ds_layers, H_ds_layers)
OVB_ds_layers <- safe_overburden(SWE_ds_layers)

# DIA padded to layer x time
DIA_ds <- matrix(0, nrow = 5, ncol = Tn)
DIA_ds <- ensure_layer_time(DIA_ds, Tn)                    # 5 x time
DIA_ds <- pad_to_layers(DIA_ds, nrow(H_ds_layers), fill=0) # layer x time

# -----------------------------
# 4) Save comparison CSV (total SWE)
# -----------------------------
out_csv <- data.frame(
  timestamp = timestamp,
  HS_meas_m = hs_vec_m,
  SWE_meas = if ("SWE" %in% names(df)) as.numeric(df$SWE) else NA_real_,
  SWE_hs2swe = as.numeric(SWE_hs2),
  SWE_dsnow = as.numeric(SWE_ds_total),
  station_id = if ("station_id" %in% names(df)) df$station_id else NA,
  season = if ("season" %in% names(df)) df$season else NA,
  stringsAsFactors = FALSE
)

write.csv(out_csv, paste0(outbase, "_compare_OABS2_2023_2024.csv"), row.names = FALSE)

# -----------------------------
# 5) Save HS2SWE NetCDF
# -----------------------------
write_model_nc(
  ncfile = paste0(outbase, "_hs2swe_layers_OABS2_2023_2024.nc"),
  time = timestamp,
  vars_named_list = list(
    HS  = HS_hs2,
    RHO = RHO_hs2,
    OVB = OVB_hs2,
    AGE = AGE_hs2,
    DIA = DIA_hs2,              # padded to layer x time
    SWE = as.numeric(SWE_hs2)   # time series
  )
)

# -----------------------------
# 6) Save delta.snow NetCDF
# -----------------------------
write_model_nc(
  ncfile = paste0(outbase, "_dsnow_layers_OABS2_2023_2024.nc"),
  time = timestamp,
  vars_named_list = list(
    HS         = H_ds_layers,
    SWE_layers = SWE_ds_layers,
    AGE        = AGE_ds_layers,
    RHO        = RHO_ds_layers,
    OVB        = OVB_ds_layers,
    DIA        = DIA_ds,                 # padded to layer x time
    SWE        = as.numeric(SWE_ds_total)
  )
)

message("Done.\n",
        "Wrote: ", normalizePath(paste0(outbase, "_compare_OABS2_2023_2024.csv")), "\n",
        "Wrote: ", normalizePath(paste0(outbase, "_hs2swe_layers_OABS2_2023_2024.nc")), "\n",
        "Wrote: ", normalizePath(paste0(outbase, "_dsnow_layers_OABS2_2023_2024.nc")))