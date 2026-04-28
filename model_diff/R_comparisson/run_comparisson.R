# run_models_and_save_layers_nc.R

suppressPackageStartupMessages({
  library(nixmass)
  library(ncdf4)
})

setwd("~/code/mt_dsnow")
source("model_diff/R_comparisson/HS2SWE_adapted_R.R")

infile  <- "/Users/jakobwerkgarner/code/mt_dsnow/model_diff/moduls/Dry_settling/synthetic_snowfalls/synthetic_big_snowfalls_100days.csv"
#infile <- "/Users/jakobwerkgarner/code/mt_dsnow/model_diff/moduls/Dry_settling/synthetic_snowfalls/synthetic_many_small_snowfalls_100days.csv"

outbase <- "/Users/jakobwerkgarner/code/mt_dsnow/model_diff/R_comparisson/output/big_snowfall/"
#outbase <- "/Users/jakobwerkgarner/code/mt_dsnow/model_diff/R_comparisson/output/small_snowfall/"

# -----------------------------
# Helpers
# -----------------------------
ensure_layer_time <- function(M, Tn) {
  if (!is.matrix(M)) stop("ensure_layer_time expects a matrix.")
  if (ncol(M) == Tn) return(M)
  if (nrow(M) == Tn) return(t(M))
  stop(sprintf(
    "Cannot orient matrix to (layer x time): got %d x %d, expected time=%d",
    nrow(M), ncol(M), Tn
  ))
}

pad_to_layers <- function(M, n_layers, fill = 0) {
  if (!is.matrix(M)) stop("pad_to_layers expects a matrix.")
  if (nrow(M) == n_layers) return(M)
  if (nrow(M) > n_layers) {
    stop(sprintf("Matrix has %d rows, but target layer count is %d.", nrow(M), n_layers))
  }
  out <- matrix(fill, nrow = n_layers, ncol = ncol(M))
  out[seq_len(nrow(M)), ] <- M
  out
}

safe_density <- function(swe_mm, h_m) {
  # returns density in kg m-3 when SWE is in mm and HS in m
  if (!is.matrix(swe_mm) || !is.matrix(h_m)) stop("safe_density expects matrices.")
  if (!all(dim(swe_mm) == dim(h_m))) {
    stop(sprintf(
      "safe_density dimension mismatch: swe is %d x %d, h is %d x %d",
      nrow(swe_mm), ncol(swe_mm), nrow(h_m), ncol(h_m)
    ))
  }
  out <- matrix(0, nrow = nrow(h_m), ncol = ncol(h_m))
  nz <- !is.na(h_m) & !is.na(swe_mm) & (h_m > 0)
  out[nz] <- swe_mm[nz] / h_m[nz]
  out
}

safe_overburden <- function(swe_mm) {
  # OVB in mm
  if (!is.matrix(swe_mm)) stop("safe_overburden expects a matrix.")
  out <- matrix(0, nrow = nrow(swe_mm), ncol = ncol(swe_mm))
  for (t in seq_len(ncol(swe_mm))) {
    v <- swe_mm[, t]
    v[is.na(v)] <- 0
    out[, t] <- rev(cumsum(rev(v)))
  }
  out
}

write_model_nc <- function(ncfile, time, vars_named_list) {
  dim_time <- ncdim_def("time", "", seq_along(time), unlim = TRUE)
  
  first_mat <- NULL
  for (v in vars_named_list) {
    if (is.matrix(v)) {
      first_mat <- v
      break
    }
  }
  if (is.null(first_mat)) stop("No matrix variables provided to write_model_nc().")
  
  dim_layer <- ncdim_def("layer", "", seq_len(nrow(first_mat)))
  
  var_defs <- list()
  for (nm in names(vars_named_list)) {
    v <- vars_named_list[[nm]]
    if (is.matrix(v)) {
      if (ncol(v) != length(time)) {
        stop(sprintf(
          "Matrix variable '%s' has ncol=%d, expected %d (= length(time))",
          nm, ncol(v), length(time)
        ))
      }
      if (nrow(v) != nrow(first_mat)) {
        stop(sprintf(
          "Matrix variable '%s' has nrow=%d, expected %d (= common layer count)",
          nm, nrow(v), nrow(first_mat)
        ))
      }
      var_defs[[nm]] <- ncvar_def(
        nm, "", list(dim_layer, dim_time),
        missval = NA_real_, prec = "double"
      )
    } else if (is.numeric(v) && length(v) == length(time)) {
      var_defs[[nm]] <- ncvar_def(
        nm, "", list(dim_time),
        missval = NA_real_, prec = "double"
      )
    } else {
      stop(sprintf(
        "Variable '%s' must be matrix(layer x time) or numeric(time).", nm
      ))
    }
  }
  
  nc <- nc_create(ncfile, var_defs)
  on.exit(nc_close(nc), add = TRUE)
  
  ncatt_put(nc, 0, "timestamps", paste(time, collapse = ","))
  
  for (nm in names(vars_named_list)) {
    ncvar_put(nc, nm, vars_named_list[[nm]])
  }
}

rng <- function(x) range(x, na.rm = TRUE, finite = TRUE)

print_summary_block <- function(name, HS, SWE_layers, RHO, OVB, AGE, DIA, SWE_total) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
  cat(name, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
  
  cat(sprintf("HS         : dim = %d x %d, range = [%g, %g]\n",
              nrow(HS), ncol(HS), rng(HS)[1], rng(HS)[2]))
  cat(sprintf("SWE_layers : dim = %d x %d, range = [%g, %g]\n",
              nrow(SWE_layers), ncol(SWE_layers), rng(SWE_layers)[1], rng(SWE_layers)[2]))
  cat(sprintf("RHO        : dim = %d x %d, range = [%g, %g]\n",
              nrow(RHO), ncol(RHO), rng(RHO)[1], rng(RHO)[2]))
  cat(sprintf("OVB        : dim = %d x %d, range = [%g, %g]\n",
              nrow(OVB), ncol(OVB), rng(OVB)[1], rng(OVB)[2]))
  cat(sprintf("AGE        : dim = %d x %d, range = [%g, %g]\n",
              nrow(AGE), ncol(AGE), rng(AGE)[1], rng(AGE)[2]))
  cat(sprintf("DIA        : dim = %d x %d, range = [%g, %g]\n",
              nrow(DIA), ncol(DIA), rng(DIA)[1], rng(DIA)[2]))
  cat(sprintf("SWE total  : len = %d, range = [%g, %g]\n",
              length(SWE_total), rng(SWE_total)[1], rng(SWE_total)[2]))
  
  hs_total <- colSums(HS, na.rm = TRUE)
  swe_from_layers <- colSums(SWE_layers, na.rm = TRUE)
  
  cat(sprintf("sum(HS by time)         range = [%g, %g]\n", rng(hs_total)[1], rng(hs_total)[2]))
  cat(sprintf("sum(SWE_layers by time) range = [%g, %g]\n", rng(swe_from_layers)[1], rng(swe_from_layers)[2]))
  cat(sprintf("max |SWE_total - sum(SWE_layers)| = %g\n",
              max(abs(SWE_total - swe_from_layers), na.rm = TRUE)))
  
  # plausibility checks for current unit convention:
  # HS [m], SWE [mm], RHO [kg m-3], OVB [mm]
  if (rng(HS)[2] > 10) warning(name, ": HS max > 10 m. This looks suspicious.")
  if (rng(SWE_layers)[2] > 1000) warning(name, ": SWE_layers max > 1000 mm. This looks suspicious.")
  if (rng(RHO)[2] > 700) warning(name, ": RHO max > 700 kg m-3. This looks suspicious.")
  if (rng(SWE_total)[2] > 1000) warning(name, ": SWE total max > 1000 mm. This looks suspicious.")
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

hs_vec_m  <- as.numeric(df$HS_meas)
hs_vec_cm <- hs_vec_m * 100

dsnow_df <- data.frame(
  hs = hs_vec_m,
  date = timestamp,
  stringsAsFactors = FALSE
)

if (!isTRUE(hs_vec_m[1] == 0)) {
  stop("For delta.snow with strict checks, HS_meas[1] must be 0 (meters).")
}

nosnow <- which(hs_vec_m == 0)

# -----------------------------
# 2) Run HS2SWE (raw output)
# -----------------------------
mr_hs2swe <- hs2swe(hs_vec_cm)

episodes <- mr_hs2swe[[1]]
if (length(episodes) == 0) stop("HS2SWE returned no episodes (no HS>0).")

hs2_max_layers <- max(vapply(episodes, function(ep) {
  if (is.null(ep$MR$OVB) || !is.matrix(ep$MR$OVB)) {
    stop("HS2SWE episode missing MR$OVB matrix.")
  }
  nrow(ep$MR$OVB)
}, numeric(1)))

HS_hs2_raw_cm    <- matrix(NA_real_, nrow = hs2_max_layers, ncol = Tn)
RHO_hs2_raw_kgm3 <- matrix(NA_real_, nrow = hs2_max_layers, ncol = Tn)
AGE_hs2          <- matrix(NA_real_, nrow = hs2_max_layers, ncol = Tn)
DIA_hs2_raw      <- matrix(NA_real_, nrow = 5, ncol = Tn)

for (ep in episodes) {
  idx <- ep$index
  MR  <- ep$MR
  
  if (is.null(MR$HS) || !is.matrix(MR$HS) || ncol(MR$HS) < 2) next
  
  cols <- 2:ncol(MR$HS)
  
  if (length(idx) != length(cols)) {
    stop(sprintf(
      "Length mismatch in HS2SWE episode: length(idx)=%d, length(cols)=%d",
      length(idx), length(cols)
    ))
  }
  
  L <- nrow(MR$HS)
  
  HS_hs2_raw_cm[1:L, idx]    <- MR$HS[, cols, drop = FALSE]
  RHO_hs2_raw_kgm3[1:L, idx] <- MR$RHO[, cols, drop = FALSE]
  AGE_hs2[1:L, idx]          <- MR$AGE[, cols, drop = FALSE]
  
  if (is.null(MR$DIA) || !is.matrix(MR$DIA)) {
    stop("MR$DIA is missing or not a matrix.")
  }
  if (nrow(MR$DIA) != 5) {
    stop(sprintf("MR$DIA has %d rows, expected 5.", nrow(MR$DIA)))
  }
  DIA_hs2_raw[, idx] <- MR$DIA[, cols, drop = FALSE]
}

HS_hs2_raw_cm[, nosnow]    <- 0
RHO_hs2_raw_kgm3[, nosnow] <- 0
AGE_hs2[, nosnow]          <- 0
DIA_hs2_raw[, nosnow]      <- 0

HS_hs2_raw_cm    <- ensure_layer_time(HS_hs2_raw_cm, Tn)
RHO_hs2_raw_kgm3 <- ensure_layer_time(RHO_hs2_raw_kgm3, Tn)
AGE_hs2          <- ensure_layer_time(AGE_hs2, Tn)
DIA_hs2_raw      <- ensure_layer_time(DIA_hs2_raw, Tn)
DIA_hs2          <- pad_to_layers(DIA_hs2_raw, nrow(HS_hs2_raw_cm), fill = 0)

# -----------------------------
# 2b) Convert HS2SWE to target units
# Target units:
#   HS         [m]
#   SWE_layers [mm]
#   SWE        [mm]
#   RHO        [kg m-3]
#   OVB        [mm]
# -----------------------------
HS_hs2  <- HS_hs2_raw_cm / 100
RHO_hs2 <- RHO_hs2_raw_kgm3

# SWE(mm) = HS(m) * RHO(kg m-3)
SWE_hs2_layers <- HS_hs2 * RHO_hs2
SWE_hs2_layers[is.na(SWE_hs2_layers)] <- 0

SWE_hs2 <- colSums(SWE_hs2_layers, na.rm = TRUE)
OVB_hs2 <- safe_overburden(SWE_hs2_layers)

# -----------------------------
# 3) Run delta.snow with layers=TRUE
# Assumed returned units:
#   h   [m]
#   swe [mm]
#   SWE [mm]
# -----------------------------
fmls <- names(formals(nixmass::swe.delta.snow))

if (!("layers" %in% fmls)) {
  stop("Your installed nixmass::swe.delta.snow does not support layers=TRUE.")
}

if ("strict_mode" %in% fmls) {
  ds_out <- nixmass::swe.delta.snow(
    dsnow_df,
    layers = TRUE,
    strict_mode = TRUE,
    verbose = FALSE,
    dyn_rho_max = FALSE
  )
} else {
  ds_out <- nixmass::swe.delta.snow(
    dsnow_df,
    layers = TRUE,
    verbose = FALSE,
    dyn_rho_max = FALSE
  )
}

if (!(is.list(ds_out) && !is.null(ds_out$h) && !is.null(ds_out$swe) && !is.null(ds_out$age))) {
  stop("Unexpected delta.snow output structure. Inspect ds_out.")
}

H_ds_layers   <- ensure_layer_time(ds_out$h,   Tn)
SWE_ds_layers <- ensure_layer_time(ds_out$swe, Tn)
AGE_ds_layers <- ensure_layer_time(ds_out$age, Tn)

SWE_ds_total <- ds_out$SWE
if (is.null(SWE_ds_total)) {
  SWE_ds_total <- colSums(SWE_ds_layers, na.rm = TRUE)
}
SWE_ds_total <- as.numeric(SWE_ds_total)

# density in kg m-3 because SWE is in mm and HS in m
RHO_ds_layers <- safe_density(SWE_ds_layers, H_ds_layers)
OVB_ds_layers <- safe_overburden(SWE_ds_layers)

DIA_ds_raw <- matrix(0, nrow = 5, ncol = Tn)
DIA_ds_raw <- ensure_layer_time(DIA_ds_raw, Tn)
DIA_ds     <- pad_to_layers(DIA_ds_raw, nrow(H_ds_layers), fill = 0)

H_ds_layers[, nosnow]   <- 0
SWE_ds_layers[, nosnow] <- 0
AGE_ds_layers[, nosnow] <- 0
RHO_ds_layers[, nosnow] <- 0
OVB_ds_layers[, nosnow] <- 0
DIA_ds[, nosnow]        <- 0
SWE_ds_total[nosnow]    <- 0

# -----------------------------
# 4) Print summary before saving
# -----------------------------
cat("\nInput HS_meas summary:\n")
cat(sprintf("HS_meas [m] range = [%g, %g]\n", rng(hs_vec_m)[1], rng(hs_vec_m)[2]))

cat("\nRaw HS2SWE summary before conversion:\n")
cat(sprintf("HS raw [cm]          range = [%g, %g]\n", rng(HS_hs2_raw_cm)[1], rng(HS_hs2_raw_cm)[2]))
cat(sprintf("RHO raw [kg m-3]     range = [%g, %g]\n", rng(RHO_hs2_raw_kgm3)[1], rng(RHO_hs2_raw_kgm3)[2]))

print_summary_block(
  "HS2SWE exported values",
  HS = HS_hs2,
  SWE_layers = SWE_hs2_layers,
  RHO = RHO_hs2,
  OVB = OVB_hs2,
  AGE = AGE_hs2,
  DIA = DIA_hs2,
  SWE_total = SWE_hs2
)

print_summary_block(
  "delta.snow exported values",
  HS = H_ds_layers,
  SWE_layers = SWE_ds_layers,
  RHO = RHO_ds_layers,
  OVB = OVB_ds_layers,
  AGE = AGE_ds_layers,
  DIA = DIA_ds,
  SWE_total = SWE_ds_total
)

cat("\nCross-model total comparison:\n")
cat(sprintf("max abs diff total HS  = %g\n",
            max(abs(colSums(HS_hs2, na.rm = TRUE) - colSums(H_ds_layers, na.rm = TRUE)), na.rm = TRUE)))
cat(sprintf("max abs diff total SWE = %g\n",
            max(abs(SWE_hs2 - SWE_ds_total), na.rm = TRUE)))

# -----------------------------
# 5) Save comparison CSV
# -----------------------------
out_csv <- data.frame(
  timestamp  = timestamp,
  HS_meas_m  = hs_vec_m,
  SWE_meas   = if ("SWE" %in% names(df)) as.numeric(df$SWE) else NA_real_,
  SWE_hs2swe = as.numeric(SWE_hs2),
  SWE_dsnow  = as.numeric(SWE_ds_total),
  station_id = if ("station_id" %in% names(df)) df$station_id else NA,
  season     = if ("season" %in% names(df)) df$season else NA,
  stringsAsFactors = FALSE
)

write.csv(
  out_csv,
  paste0(outbase, "_compare_synthetic_big_snowfall_100days.csv"),
  row.names = FALSE
)

# -----------------------------
# 6) Save HS2SWE NetCDF
# -----------------------------
write_model_nc(
  ncfile = paste0(outbase, "_hs2swe_layers_synthetic_big_snowfall_100days.nc"),
  time   = timestamp,
  vars_named_list = list(
    HS         = HS_hs2,
    SWE_layers = SWE_hs2_layers,
    RHO        = RHO_hs2,
    OVB        = OVB_hs2,
    AGE        = AGE_hs2,
    DIA        = DIA_hs2,
    SWE        = as.numeric(SWE_hs2)
  )
)

# -----------------------------
# 7) Save dSnow NetCDF
# -----------------------------
write_model_nc(
  ncfile = paste0(outbase, "_dsnow_layers_synthetic_big_snowfall_100days.nc"),
  time   = timestamp,
  vars_named_list = list(
    HS         = H_ds_layers,
    SWE_layers = SWE_ds_layers,
    RHO        = RHO_ds_layers,
    OVB        = OVB_ds_layers,
    AGE        = AGE_ds_layers,
    DIA        = DIA_ds,
    SWE        = as.numeric(SWE_ds_total)
  )
)

message("Done")