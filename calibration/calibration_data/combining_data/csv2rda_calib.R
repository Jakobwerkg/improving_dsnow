library(zoo)

base_dir <- "/Users/jakobwerkgarner/code/mt_dsnow"

out_file_name <- "d_obs_WIN_MAG.rda"


# ---- Stations to ignore ----
ignore_stations <- c(
  "Davos Flueelastr.",
  "Sta. Maria"
)
# set input and output dir
indir   <- file.path(base_dir, "calibration/calibration_data/output/HS_SWE_by_station")
out_dir <- file.path(base_dir, "calibration/calibration_data/output/calibration_rda_files")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


#---code 
files <- list.files(indir, pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No CSV files found in: ", indir)

read_one <- function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(df)) return(NULL)
  
  names(df) <- trimws(names(df))
  if (!all(c("date", "hs", "swe_obs") %in% names(df))) {
    stop("Unexpected columns in: ", basename(f))
  }
  
  dates <- as.Date(df$date)
  if (all(is.na(dates))) stop("Could not parse dates in: ", basename(f))
  
  df <- data.frame(
    Date   = dates,
    Hobs   = as.numeric(df$hs),
    SWEobs = as.numeric(df$swe_obs),
    check.names = FALSE
  )
  
  o <- order(df$Date)
  df <- df[o, , drop = FALSE]
  
  if (any(duplicated(df$Date))) {
    df <- df[!duplicated(df$Date), , drop = FALSE]
  }
  
  z <- zoo(as.matrix(df[, c("Hobs", "SWEobs"), drop = FALSE]),
           order.by = df$Date)
  storage.mode(z) <- "numeric"
  z
}

zlist <- lapply(files, read_one)
ok <- !vapply(zlist, is.null, logical(1))
zlist <- zlist[ok]
files <- files[ok]
if (!length(zlist)) stop("No valid series built from CSVs in: ", indir)

base_names <- tools::file_path_sans_ext(basename(files))
base_names <- sub("_hs_swe_obs$", "", base_names)
station_names <- gsub("_+", " ", base_names)
station_names <- trimws(station_names)

# Clean ignore names same way as station_names
clean_ignore <- trimws(gsub("_+", " ", ignore_stations))

# Filter out ignored stations
keep <- !(station_names %in% clean_ignore)
station_names <- station_names[keep]
zlist <- zlist[keep]


names(zlist) <- station_names

drop_na_Hobs <- function(z) {
  if (!("Hobs" %in% colnames(z))) return(z)
  z[!is.na(z[, "Hobs"]), , drop = FALSE]
}
zlist <- lapply(zlist, drop_na_Hobs)

d_obs <- zlist

out_file <- file.path(out_dir, out_file_name)
save(d_obs, file = out_file, compress = "xz")

str(d_obs, max.level = 1)
cat("Saved:", out_file, "\n")