library(ncdf4)
library(nixmass)

infile   <- "/Users/jakobwerkgarner/code/mt_dsnow/HNW_validation/validation_input_Mag25/Mag25_all.nc"
rds_out  <- "/Users/jakobwerkgarner/code/mt_dsnow/HNW_validation/dSnow/validated_data/Mag25_SWE_mod_from_R_dyn_rho_max.rds"

nc <- nc_open(infile)

HS      <- ncvar_get(nc, "HS")
SWE     <- ncvar_get(nc, "SWE")
time    <- ncvar_get(nc, "time")
station <- ncvar_get(nc, "station")

nc_close(nc)

ntime <- dim(HS)[2]
nstat <- dim(HS)[1]

station <- as.character(station)[seq_len(nstat)]
dates <- as.Date(time, origin = "2016-09-01")

SWE_mod <- matrix(NA_real_, nrow = nstat, ncol = ntime)
obs_counts <- data.frame(station = station, n_obs = NA_integer_)

for (i in seq_len(nstat)) {
  st <- station[i]
  cat("Processing:", st, "\n")
  
  hs_i <- HS[i, ]
  ok <- !is.na(hs_i)
  
  dat_i <- data.frame(
    date = dates[ok],
    hs   = hs_i[ok]
  )
  
  swe_i <- nixmass(
    data = dat_i,
    model = "delta.snow.dyn_rho_max",
    verbose = FALSE
  )
  
  SWE_mod[i, ok] <- swe_i$swe[[1]]$SWE
  obs_counts$n_obs[i] <- sum(!is.na(SWE[i, ]))
}

print(obs_counts)



csv_out <- "/Users/jakobwerkgarner/code/mt_dsnow/HNW_validation/dSnow/validated_data/Mag25_SWE_mod_from_R_dyn_rho_max.csv"
write.csv(SWE_mod, csv_out, row.names = FALSE)


cat("Saved:", csv_out, "\n")


