# Load required packages
library(nixmass)
library(zoo)
library(ncdf4)

# Load the snow depth data
hs_data <- read.csv("calibration/calibration_data/output/HS_SWE_by_station/Obernberg_hs_swe_obs.csv")

# Convert date to Date type and filter for the specific season
hs_data$date <- as.Date(hs_data$date)
data_season <- hs_data[hs_data$date >= as.Date("1998-11-14") & hs_data$date <= as.Date("1999-04-07"), ]

# Keep only date and hs columns
data_season <- data_season[, c("date", "hs")]

# Run delta.snow model with layers
result <- swe.delta.snow(data_season, dyn_rho_max = FALSE, layers = TRUE)




# Save result as NetCDF file
# Define dimensions
n_days <- ncol(result$h)
n_layers <- nrow(result$h)

# Handle processes - replace NA with empty string
processes_clean <- ifelse(is.na(result$processes), "", as.character(result$processes))
max_strlen <- max(c(nchar(processes_clean), nchar(as.character(data_season$date))))

dim_time <- ncdim_def("time", "days", 1:n_days)
dim_layer <- ncdim_def("layer", "layer_number", 1:n_layers)
dim_str <- ncdim_def("string_length", "", 1:max_strlen, create_dimvar = FALSE)

# Define variables
var_swe <- ncvar_def("SWE", "mm", list(dim_time), -999, "Snow Water Equivalent")
var_h <- ncvar_def("h", "m", list(dim_layer, dim_time), -999, "Snow depth by layer")
var_swe_layer <- ncvar_def("swe_layer", "mm", list(dim_layer, dim_time), -999, "SWE by layer")
var_age <- ncvar_def("age", "days", list(dim_layer, dim_time), -999, "Age by layer")
var_rho <- ncvar_def("rho", "kg/m3", list(dim_layer, dim_time), -999, "Density by layer")
var_processes <- ncvar_def("processes", "", list(dim_str, dim_time), prec = "char", 
                           longname = "Active model process for each timestep")
var_dates <- ncvar_def("dates", "", list(dim_str, dim_time), prec = "char",
                       longname = "Date for each timestep")

# Create NetCDF file
nc_file <- nc_create("model_diff/dsnow/layer_model_run/delta_snow_result.nc", 
                     list(var_swe, var_h, var_swe_layer, var_age, var_rho, var_processes, var_dates))

# Write data
ncvar_put(nc_file, var_swe, result$SWE)
ncvar_put(nc_file, var_h, result$h)
ncvar_put(nc_file, var_swe_layer, result$swe)
ncvar_put(nc_file, var_age, result$age)
ncvar_put(nc_file, var_rho, rho_layers)

# Write processes as variable (pad strings to same length)
processes_padded <- sprintf(paste0("%-", max_strlen, "s"), processes_clean)
ncvar_put(nc_file, var_processes, processes_padded)

# Write dates as variable (pad strings to same length)
dates_str <- as.character(data_season$date)
dates_padded <- sprintf(paste0("%-", max_strlen, "s"), dates_str)
ncvar_put(nc_file, var_dates, dates_padded)

# Add time dimension attributes
ncatt_put(nc_file, "time", "units", "timestep")
ncatt_put(nc_file, "time", "description", "Sequential timestep index")

# Add global attributes
ncatt_put(nc_file, 0, "title", "Delta Snow Model Results")
ncatt_put(nc_file, 0, "model", "swe.delta.snow")
ncatt_put(nc_file, 0, "model_version", "nixmass package")
ncatt_put(nc_file, 0, "dyn_rho_max", "False")
ncatt_put(nc_file, 0, "start_date", as.character(min(data_season$date)))
ncatt_put(nc_file, 0, "end_date", as.character(max(data_season$date)))
ncatt_put(nc_file, 0, "station", "Obernberg")
ncatt_put(nc_file, 0, "creation_date", as.character(Sys.time()))
ncatt_put(nc_file, 0, "description", "Layer-wise delta.snow model output with dynamic rho_max")
ncatt_put(nc_file, 0, "n_days", n_days)
ncatt_put(nc_file, 0, "n_layers", n_layers)

# Close the file
nc_close(nc_file)

print("Results saved to model_diff/dsnow/layer_model_run/delta_snow_result.nc")
print(paste("Saved", length(result$processes), "process records"))
