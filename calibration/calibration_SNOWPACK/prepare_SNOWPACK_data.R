library(zoo)

base_dir <- "/Users/jakobwerkgarner/code/mt_dsnow"

# Input from SNOWPACK_data
indir <- file.path(base_dir, "par_sens/SNOWPACK_data")

# Output (change if you want)
out_dir <- file.path(
    "/Users/jakobwerkgarner/code/mt_dsnow/calibration/",
    "calibration_SNOWPACK/data"
)
out_file_name <- "d_obs_SNOWPACK.rda"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Optional: stations to ignore
ignore_stations <- character(0)
# e.g. ignore_stations <- c("Davos Flueelastr.", "Sta. Maria")


# ----------------------------
# Read one .smet file
# ----------------------------
read_one_smet <- function(f) {
    lines <- readLines(f, warn = FALSE)
    if (!length(lines)) return(NULL)

    tl <- trimws(lines)
    data_idx <- which(tolower(tl) == "[data]")
    if (!length(data_idx)) stop("No [DATA] section in: ", basename(f))
    data_idx <- data_idx[1]

    # Parse header key=value
    header_lines <- lines[seq_len(data_idx - 1)]
    header_lines <- header_lines[grepl("=", header_lines)]
    kv <- strsplit(header_lines, "=", fixed = TRUE)

    keys <- trimws(tolower(vapply(kv, `[`, "", 1)))
    vals <- trimws(
        vapply(kv, function(x) paste(x[-1], collapse = "="), "")
    )
    hdr <- setNames(vals, keys)

    # Station name: prefer station_name, then station_id, then filename
    station <- hdr[["station_name"]]
    if (is.null(station) || station == "") {
        station <- hdr[["station_id"]]
    }
    if (is.null(station) || station == "") {
        station <- tools::file_path_sans_ext(basename(f))
    }
    station <- trimws(gsub("_+", " ", station))

    # Fields
    fields <- hdr[["fields"]]
    if (is.null(fields) || fields == "") {
        stop("No 'fields' in header: ", basename(f))
    }
    field_names <- strsplit(fields, "[[:space:]]+")[[1]]
    field_names <- field_names[field_names != ""]

    # Data lines (skip comments/empty)
    dlines <- lines[(data_idx + 1):length(lines)]
    dlines <- dlines[!grepl("^\\s*(#|$)", dlines)]
    if (!length(dlines)) return(NULL)

    df <- read.table(
        text = paste(dlines, collapse = "\n"),
        stringsAsFactors = FALSE,
        fill = TRUE
    )

    if (ncol(df) < length(field_names)) {
        stop("Data columns < fields in: ", basename(f))
    }
    if (ncol(df) > length(field_names)) {
        df <- df[, seq_along(field_names), drop = FALSE]
    }
    names(df) <- field_names

    # Time column
    time_col <- names(df)[tolower(names(df)) %in%
        c("timestamp", "date", "datetime", "time")]
    if (!length(time_col)) {
        stop("No timestamp/date column in: ", basename(f))
    }
    time_col <- time_col[1]

    ts_raw <- as.character(df[[time_col]])
    ts <- as.POSIXct(
        ts_raw,
        tz = "UTC",
        tryFormats = c(
            "%Y-%m-%dT%H:%M:%S",
            "%Y-%m-%dT%H:%M",
            "%Y-%m-%d %H:%M:%S",
            "%Y-%m-%d %H:%M",
            "%Y-%m-%d"
        )
    )
    if (all(is.na(ts))) stop("Could not parse timestamps in: ", basename(f))
    dates <- as.Date(ts)

    # HS column: prefer hs_mod, then hs, then anything starting with hs
    nms_low <- tolower(names(df))
    hs_col <- names(df)[nms_low == "hs_mod"]
    if (!length(hs_col)) hs_col <- names(df)[nms_low == "hs"]
    if (!length(hs_col)) {
        hs_col <- grep("^hs", names(df), ignore.case = TRUE, value = TRUE)
    }
    if (!length(hs_col)) {
        stop("No HS/hs_mod column found in: ", basename(f))
    }
    hs_col <- hs_col[1]

    # SWE column: prefer swe_obs, then swe
    swe_col <- names(df)[nms_low == "swe_obs"]
    if (!length(swe_col)) swe_col <- names(df)[nms_low == "swe"]
    if (!length(swe_col)) {
        swe_col <- grep("^swe", names(df), ignore.case = TRUE, value = TRUE)
    }

    hs <- suppressWarnings(as.numeric(df[[hs_col]]))
    swe <- if (length(swe_col)) {
        suppressWarnings(as.numeric(df[[swe_col[1]]]))
    } else {
        rep(NA_real_, nrow(df))
    }

    # Handle nodata value from header if present
    nodata <- suppressWarnings(as.numeric(hdr[["nodata"]]))
    if (length(nodata) == 1 && !is.na(nodata)) {
        hs[hs == nodata] <- NA_real_
        swe[swe == nodata] <- NA_real_
    }

    out <- data.frame(
        Date = dates,
        Hobs = as.numeric(hs),
        SWEobs = as.numeric(swe),
        check.names = FALSE
    )

    out <- out[!is.na(out$Date), , drop = FALSE]
    out <- out[order(out$Date), , drop = FALSE]
    out <- out[!duplicated(out$Date), , drop = FALSE]

    z <- zoo(as.matrix(out[, c("Hobs", "SWEobs"), drop = FALSE]), order.by = out$Date)
    storage.mode(z) <- "numeric"

    attr(z, "station_name") <- station
    z
}


# ----------------------------
# Build d_obs from all .smet files
# ----------------------------
files <- list.files(
    indir,
    pattern = "\\.smet$",
    full.names = TRUE,
    recursive = TRUE
)
if (!length(files)) stop("No .smet files found in: ", indir)

zlist <- lapply(files, read_one_smet)
ok <- !vapply(zlist, is.null, logical(1))
zlist <- zlist[ok]
files <- files[ok]
if (!length(zlist)) {
    stop("No valid series built from .smet in: ", indir)
}

station_names <- vapply(
    zlist,
    function(z) {
        s <- attr(z, "station_name")
        if (is.null(s) || !nzchar(s)) return(NA_character_)
        trimws(gsub("_+", " ", s))
    },
    character(1)
)

# Fallback to filename if no station attr
fallback <- tools::file_path_sans_ext(basename(files))
fallback <- trimws(gsub("_+", " ", fallback))
station_names[is.na(station_names) | station_names == ""] <-
    fallback[is.na(station_names) | station_names == ""]

# Ignore stations
clean_ignore <- trimws(gsub("_+", " ", ignore_stations))
keep <- !(station_names %in% clean_ignore)
station_names <- station_names[keep]
zlist <- zlist[keep]

# Ensure unique names
station_names <- make.unique(station_names, sep = "_")
names(zlist) <- station_names

# Drop NA rows in either Hobs or SWEobs
drop_na_rows <- function(z) {
    z[!(is.na(z[, "Hobs"]) | is.na(z[, "SWEobs"])), , drop = FALSE]
}
zlist <- lapply(zlist, drop_na_rows)

d_obs <- zlist

out_file <- file.path(out_dir, out_file_name)
save(d_obs, file = out_file, compress = "xz")

str(d_obs, max.level = 1)
cat("Saved:", out_file, "\n")