#!/usr/bin/env Rscript
# ============================================================
# Δ-SNOW MORRIS (MEAS CSV DIRECTORY) — ROBUST + DEBUGGABLE
#
# Input CSV per season:
#   time, hs, swe_obs
#   filename: <station>_<YYYY>_<YYYY>.csv
#
# Key points:
# - NO HS unit conversion (assume hs already in meters)
# - RMSE computed ONLY on dates where swe_obs exists
# - Morris per file; then aggregate (mean ± sd) across files
# - Writes diagnostics with first Δ-SNOW error per file
# ============================================================

options(error = function() traceback(2))

suppressPackageStartupMessages({
  library(sensitivity)
  library(ggplot2)
  library(dplyr)
  library(nixmass)   # for swe.delta.snow()
})

# -------------------------------
# USER SETTINGS
# -------------------------------
CSV_DIR <- "/Users/jakobwerkgarner/code/mt_dsnow/par_sens/seasonal_measurments/data"
OUTPUT_DIR <- "/Users/jakobwerkgarner/code/mt_dsnow/par_sens/dsnow_sensitivity_R"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

MORRIS_N      <- 20
MORRIS_LEVELS <- 4

PARAM_NAMES <- c("rho.max", "rho.null", "c.ov", "k.ov", "k", "tau", "eta.null")
BOUNDS <- rbind(
  c(300, 600),        # rho.max
  c(50,  200),        # rho.null
  c(1e-4, 1e-3),      # c.ov
  c(0.01, 1),         # k.ov
  c(0.01, 0.2),       # k
  c(0.001, 0.20),     # tau
  c(1e6,  2e7)        # eta.null
)

# scoring / robustness
MIN_MATCHED_SWE <- 5         # minimum SWE_obs days per file (after winter window)
FAIL_PENALTY    <- 1e9       # used if a run fails

# ============================================================
# Helpers
# ============================================================

safe_numeric <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "nan", "NULL", "-")] <- NA
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- gsub("[^0-9.+-eE]", "", x)
  suppressWarnings(as.numeric(x))
}

parse_meas_filename <- function(path) {
  base <- sub("\\.csv$", "", basename(path))
  m <- regexec("^(.*)_([0-9]{4})_([0-9]{4})$", base)
  reg <- regmatches(base, m)[[1]]
  if (length(reg) != 4) stop("Bad filename (need <station>_YYYY_YYYY.csv): ", basename(path))
  list(
    station = reg[2],
    season  = paste0(reg[3], "_", reg[4]),
    y0      = as.integer(reg[3]),
    y1      = as.integer(reg[4])
  )
}

get_swe_vector <- function(res) {
  # handle common return types:
  # - list with $SWE
  # - numeric vector directly
  if (is.list(res) && "SWE" %in% names(res)) {
    return(suppressWarnings(as.numeric(res$SWE)))
  }
  if (is.numeric(res)) {
    return(as.numeric(res))
  }
  NULL
}

# ============================================================
# READ + PREPARE CSV (WINTER WINDOW) for a single season file
# ============================================================

read_prepare_one_csv <- function(csv_file) {
  
  meta <- parse_meas_filename(csv_file)
  
  df <- read.csv(csv_file, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(NULL)
  
  if (!all(c("time", "hs", "swe_obs") %in% names(df))) {
    stop("Missing required columns in ", basename(csv_file), " (need time, hs, swe_obs)")
  }
  
  day <- as.Date(df$time)
  hs  <- safe_numeric(df$hs)       # <-- FIXED: hs is actually created
  swe <- safe_numeric(df$swe_obs)
  
  # Winter window for that season
  t_start <- as.Date(sprintf("%d-11-01", meta$y0))
  t_end   <- as.Date(sprintf("%d-05-01", meta$y1))
  
  keep <- day >= t_start & day < t_end
  if (!any(keep)) return(NULL)
  
  day <- day[keep]; hs <- hs[keep]; swe <- swe[keep]
  
  # unique daily (keep first occurrence)
  ord <- order(day)
  day <- day[ord]; hs <- hs[ord]; swe <- swe[ord]
  dup <- duplicated(day)
  if (any(dup)) {
    day <- day[!dup]; hs <- hs[!dup]; swe <- swe[!dup]
  }
  
  # Δ-SNOW input (daily)
  run_data <- data.frame(
    date = as.POSIXct(day, tz = "UTC"),
    hs   = as.numeric(hs),
    stringsAsFactors = FALSE
  )
  
  # match your pattern
  run_data$hs[!is.finite(run_data$hs)] <- 0
  run_data$hs[run_data$hs < 0] <- 0
  run_data$hs[1] <- 0
  
  # SWE observations (daily) for scoring
  obs_df <- data.frame(
    day     = as.Date(day),
    swe_obs = as.numeric(swe),
    stringsAsFactors = FALSE
  )
  obs_df <- obs_df[is.finite(obs_df$swe_obs), , drop = FALSE]
  
  if (nrow(obs_df) < MIN_MATCHED_SWE) return(NULL)
  
  list(
    file    = basename(csv_file),
    station = meta$station,
    season  = meta$season,
    run_data = run_data,
    obs_df   = obs_df
  )
}

# ============================================================
# MODEL → RMSE only on SWE_obs dates
# Returns list(rmse=..., err=...)
# ============================================================

model_rmse <- function(theta, prep) {
  
  params <- as.list(setNames(as.numeric(theta), PARAM_NAMES))
  
  res <- tryCatch(
    swe.delta.snow(
      data = prep$run_data,
      model_opts = params,
      layer = TRUE,
      dyn_rho_max = FALSE
    ),
    error = function(e) e
  )
  
  if (inherits(res, "error")) {
    return(list(rmse = FAIL_PENALTY, err = conditionMessage(res)))
  }
  
  SWE_mod <- get_swe_vector(res)
  if (is.null(SWE_mod)) {
    return(list(rmse = FAIL_PENALTY, err = "Model returned SWE in an unexpected format"))
  }
  
  # align model SWE with run_data length
  n <- min(length(SWE_mod), nrow(prep$run_data))
  SWE_mod <- SWE_mod[seq_len(n)]
  mod_df <- data.frame(
    day     = as.Date(prep$run_data$date[seq_len(n)]),
    swe_mod = as.numeric(SWE_mod),
    stringsAsFactors = FALSE
  )
  mod_df <- mod_df[is.finite(mod_df$swe_mod), , drop = FALSE]
  
  score_df <- merge(prep$obs_df, mod_df, by = "day", all = FALSE)
  if (nrow(score_df) < MIN_MATCHED_SWE) {
    return(list(rmse = FAIL_PENALTY, err = "Too few matched SWE days after merge"))
  }
  
  rmse <- sqrt(mean((score_df$swe_mod - score_df$swe_obs)^2))
  if (!is.finite(rmse)) {
    return(list(rmse = FAIL_PENALTY, err = "Non-finite RMSE"))
  }
  
  list(rmse = rmse, err = NA_character_)
}

# ============================================================
# MORRIS DESIGN
# ============================================================

cat("Generating Morris design...\n")
m_template <- morris(
  model   = NULL,
  factors = PARAM_NAMES,
  r       = MORRIS_N,
  design  = list(type = "oat", levels = MORRIS_LEVELS, grid.jump = 1),
  binf    = BOUNDS[, 1],
  bsup    = BOUNDS[, 2]
)

X <- m_template$X
n_runs <- nrow(X)
cat("Morris runs per file:", n_runs, "\n")

# ============================================================
# RUN OVER FILES
# ============================================================

CSV_FILES <- list.files(CSV_DIR, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(CSV_FILES) > 0)
cat("Found", length(CSV_FILES), "CSV files.\n")

all_results <- list()
diag_list <- list()

for (f in CSV_FILES) {
  
  prep <- tryCatch(read_prepare_one_csv(f), error = function(e) NULL)
  if (is.null(prep)) next
  
  cat("\n============================================================\n")
  cat("File:", prep$file, "|", prep$station, prep$season, "\n")
  cat("SWE obs days:", nrow(prep$obs_df), "\n")
  
  Y <- numeric(n_runs)
  first_err <- NA_character_
  
  for (i in seq_len(n_runs)) {
    if (i %% 20 == 0) cat(sprintf("  [%3d / %3d]\n", i, n_runs))
    tmp <- model_rmse(X[i, ], prep)
    Y[i] <- tmp$rmse
    if (is.na(first_err) && is.character(tmp$err) && !is.na(tmp$err)) {
      first_err <- tmp$err
    }
  }
  
  ok_mask <- is.finite(Y) & Y < FAIL_PENALTY
  n_ok <- sum(ok_mask)
  n_fail <- n_runs - n_ok
  
  diag_list[[length(diag_list) + 1L]] <- data.frame(
    file = prep$file, station = prep$station, season = prep$season,
    n_runs = n_runs, n_ok = n_ok, n_fail = n_fail,
    first_error = ifelse(is.na(first_err), "", first_err),
    stringsAsFactors = FALSE
  )
  
  if (n_ok == 0) {
    cat("  -> all runs failed. First error:\n")
    cat("     ", first_err, "\n")
    next
  }
  
  # Fill fails with mean of successful runs (keeps Morris stable)
  mean_ok <- mean(Y[ok_mask], na.rm = TRUE)
  Y[!ok_mask] <- mean_ok
  
  m <- m_template
  m <- tell(m, Y)
  EE <- m$ee
  stopifnot(is.matrix(EE))
  
  mu_star <- apply(EE, 2, function(x) mean(abs(x), na.rm = TRUE))
  sigma   <- apply(EE, 2, function(x) sd(x, na.rm = TRUE))
  
  res_one <- data.frame(
    file = prep$file, station = prep$station, season = prep$season,
    parameter = colnames(EE),
    mu_star = as.numeric(mu_star),
    sigma   = as.numeric(sigma),
    stringsAsFactors = FALSE
  )
  
  all_results[[length(all_results) + 1L]] <- res_one
  cat("  -> OK. Successful runs:", n_ok, "/", n_runs, "\n")
}

diag_df <- if (length(diag_list) > 0) do.call(rbind, diag_list) else data.frame()
write.csv(diag_df, file.path(OUTPUT_DIR, "series_diagnostics.csv"), row.names = FALSE)

if (length(all_results) == 0) {
  stop("No file produced valid Morris indices. Open series_diagnostics.csv and look at first_error.")
}

results_df <- do.call(rbind, all_results)
write.csv(results_df, file.path(OUTPUT_DIR, "morris_indices_by_file.csv"), row.names = FALSE)

# ============================================================
# AGGREGATE ACROSS FILES
# ============================================================

overall_stats <- results_df %>%
  group_by(parameter) %>%
  summarize(
    mu_star_mean = mean(mu_star, na.rm = TRUE),
    mu_star_sd   = sd(mu_star, na.rm = TRUE),
    sigma_mean   = mean(sigma, na.rm = TRUE),
    sigma_sd     = sd(sigma, na.rm = TRUE),
    n_files      = n_distinct(file),
    .groups      = "drop"
  ) %>%
  arrange(desc(mu_star_mean))

write.csv(overall_stats, file.path(OUTPUT_DIR, "morris_indices_overall.csv"), row.names = FALSE)

# ============================================================
# PLOTS
# ============================================================

theme_white_snow <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", linewidth = 0.5),
    axis.text          = element_text(color = "black"),
    axis.title         = element_text(color = "black"),
    plot.title         = element_text(face = "bold")
  )

p_bar <- ggplot(overall_stats, aes(x = reorder(parameter, mu_star_mean), y = mu_star_mean)) +
  geom_col(fill = "grey80") +
  geom_errorbar(
    aes(ymin = pmax(0, mu_star_mean - mu_star_sd),
        ymax = mu_star_mean + mu_star_sd),
    width = 0.15,
    color = "grey40"
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = expression(mu^"*"),
    title = "Δ-SNOW Morris Sensitivity (Overall)",
    subtitle = "Mean μ* across files ± SD"
  ) +
  theme_white_snow

ggsave(file.path(OUTPUT_DIR, "morris_mu_star_overall_barplot.png"),
       p_bar, width = 9, height = 5, dpi = 300)

p_scatter <- ggplot(overall_stats, aes(mu_star_mean, sigma_mean, label = parameter)) +
  geom_point(size = 2.7, color = "grey30") +
  geom_text(nudge_y = 0.02 * max(overall_stats$sigma_mean, na.rm = TRUE),
            size = 3, color = "grey20") +
  labs(
    x = expression(mu^"*"),
    y = expression(sigma),
    title = "Morris Scatter (Overall)"
  ) +
  theme_white_snow

ggsave(file.path(OUTPUT_DIR, "morris_overall_scatter.png"),
       p_scatter, width = 8, height = 6, dpi = 300)

cat("\nDONE. Results written to:", OUTPUT_DIR, "\n")
cat("Files used:", n_distinct(results_df$file), "\n")