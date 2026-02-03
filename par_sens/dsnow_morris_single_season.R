# ============================================================
# READ + PREPARE SMET DATA (WINTER WINDOW).  for the signle season
# ============================================================



cat("Reading SMET file...\n")

lines <- readLines(SMET_FILE)

# --- header ---
fields <- strsplit(sub("fields = ", "", grep("^fields ", lines, value = TRUE)), " ")[[1]]
data_i <- grep("^\\[DATA\\]", lines) + 1
stopifnot(length(fields) > 0, data_i > 1)

# --- data ---
df_smet <- as.data.frame(
  do.call(rbind, strsplit(lines[data_i:length(lines)], "[[:space:]]+")),
  stringsAsFactors = FALSE
)
colnames(df_smet) <- fields
df_smet[names(df_smet) != "timestamp"] <- lapply(df_smet[names(df_smet) != "timestamp"], as.numeric)

stopifnot(all(c("HS_mod", "SWE", "timestamp") %in% names(df_smet)))

# ============================================================
# DAILY RESAMPLING (HOURLY SMET → DAILY Δ-SNOW INPUT)
# ============================================================

# Parse timestamp (SMET uses ISO strings)
date_hourly <- as.POSIXct(df_smet$timestamp, tz = "UTC")
stopifnot(!anyNA(date_hourly))

# Winter window
t_start <- as.POSIXct("2023-11-01", tz = "UTC")
t_end   <- as.POSIXct("2024-05-01", tz = "UTC")

keep <- date_hourly >= t_start & date_hourly < t_end
stopifnot(sum(keep) > 24)   # at least one day

# Convert to Date for daily aggregation
day <- as.Date(date_hourly[keep])

# --- DAILY MEANS ---
HS_day  <- tapply(df_smet$HS_mod[keep], day, mean, na.rm = TRUE)
SWE_day <- tapply(df_smet$SWE[keep],    day, mean, na.rm = TRUE)

# Δ-SNOW input (DAILY)
run_data <- data.frame(
  date = as.POSIXct(names(HS_day), tz = "UTC"),
  hs   = as.numeric(HS_day)
)

run_data$hs[is.na(run_data$hs)] <- 0
run_data$hs[1] <- 0

# Reference SWE (DAILY)
SWE_ref   <- as.numeric(SWE_day)
valid_swe <- is.finite(SWE_ref)
stopifnot(sum(valid_swe) > 10)

cat("Daily period:", format(range(run_data$date)), "\n")
cat("Number of days:", nrow(run_data), "\n")
cat("Valid SWE days:", sum(valid_swe), "\n")

# ============================================================
# MODEL → RMSE (against SMET SWE)
# ============================================================

model_rmse <- function(theta) {

  params <- as.list(setNames(as.numeric(theta), PARAM_NAMES))

  res <- tryCatch(
    swe.delta.snow(
      data = run_data,
      model_opts = params,
      layer = TRUE,
      dyn_rho_max = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(res)) return(1e9)

  SWE_mod <- as.numeric(res$SWE)

  sqrt(mean((SWE_mod[valid_swe] - SWE_ref[valid_swe])^2))
}

# ============================================================
# MORRIS DESIGN
# ============================================================

cat("Generating Morris design...\n")

m <- morris(
  model   = NULL,
  factors = PARAM_NAMES,
  r       = MORRIS_N,
  design  = list(type = "oat", levels = MORRIS_LEVELS, grid.jump = 1),
  binf    = BOUNDS[, 1],
  bsup    = BOUNDS[, 2]
)

# ============================================================
# MODEL EVALUATION
# ============================================================

cat("Evaluating Δ-SNOW...\n")

n_runs <- nrow(m$X)
Y <- numeric(n_runs)

for (i in seq_len(n_runs)) {
  cat(sprintf("  [%3d / %3d]\n", i, n_runs))
  Y[i] <- tryCatch(model_rmse(m$X[i, ]), error = function(e) 1e9)
}
if (all(Y == 1e+09)) {
  warning("All Morris sensitivities shit")
}
cat("RMSE summary:\n")
print(summary(Y))
cat("Unique RMSE values:", length(unique(Y)), "\n")


# ============================================================
# MORRIS ANALYSIS
# ============================================================

# tell Morris the outputs
m <- tell(m, Y)

# ============================================================
# MORRIS INDICES (ROBUST, VERSION-INDEPENDENT)
# ============================================================

EE <- m$ee               # matrix: (r × p)
stopifnot(is.matrix(EE))

mu_star <- apply(EE, 2, function(x) mean(abs(x), na.rm = TRUE))
sigma   <- apply(EE, 2, function(x) sd(x, na.rm = TRUE))

res <- data.frame(
  parameter = colnames(EE),
  mu_star   = mu_star,
  sigma     = sigma
)

print(res[order(-res$mu_star), ], row.names = FALSE)

# ============================================================
# PLOTS
# ============================================================

p1 <- ggplot(res, aes(x = reorder(parameter, mu_star), y = mu_star)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(
    aes(ymin = pmax(0, mu_star - sigma), ymax = mu_star + sigma),
    width = 0.2
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Morris μ* (importance)",
    title = "Δ-SNOW Morris Sensitivity (RMSE vs SMET SWE)"
  ) +
  theme_minimal()

ggsave(
  file.path(OUTPUT_DIR, "morris_mu_star_barplot.png"),
  p1, width = 9, height = 5, dpi = 300
)

p2 <- ggplot(res, aes(mu_star, sigma, label = parameter)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.6, size = 3.5) +
  labs(
    x = "μ*",
    y = "σ",
    title = "Δ-SNOW Morris Scatter (SMET reference)"
  ) +
  theme_minimal()

ggsave(
  file.path(OUTPUT_DIR, "morris_mu_star_sigma_scatter.png"),
  p2, width = 8, height = 6, dpi = 300
)

cat("\nDONE. Results written to:", OUTPUT_DIR, "\n")