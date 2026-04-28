# =============================================================================
# HNW_validation_helper.R
# R equivalent of HNW_validation_helper.py
#
# Provides three internal helper functions and two public validation functions:
#   validate_hnw_mag25()  – density scatter for HNW  (lim 0–100 mm)
#   validate_swe_mag25()  – density scatter for SWE  (lim 0–1000 mm)
#
# Dependencies: ggplot2, lubridate
# =============================================================================

library(ggplot2)
library(lubridate)


# ── Internal helpers (prefixed with . to signal private use) ──────────────────

#' Filter a data frame to the snow season (Nov 1 – Apr 30).
#' Passes through unchanged when full_season = TRUE.
#' Equivalent to _filter_season() in Python.
.filter_season <- function(df, full_season = FALSE) {
  df$time <- as.Date(df$time)
  if (full_season) return(df)
  m <- lubridate::month(df$time)
  df[m >= 11 | m <= 4, ]
}


#' Calculate standard validation metrics.
#' x = observed values, y = modelled values.
#' Returns a named list: RMSE, Bias, Rel_BIAS (relative / percent bias), R2.
#' Equivalent to _calculate_metrics() in Python.
.calc_metrics <- function(x, y) {
  res    <- y - x
  rmse   <- sqrt(mean(res^2))
  bias   <- mean(res)
  pbias  <- if (sum(x) != 0) sum(res) / sum(x) else NA_real_   # Rel_BIAS
  ss_res <- sum((x - y)^2)
  ss_tot <- sum((x - mean(x))^2)
  r2     <- if (ss_tot != 0) 1 - ss_res / ss_tot else NA_real_
  list(RMSE = rmse, Bias = bias, Rel_BIAS = pbias, R2 = r2)
}


#' Build a 2-D density (histogram) validation plot with log-scaled colour.
#' Reproduces matplotlib's hist2d + LogNorm + viridis used in Python.
#'
#' @param x          Observed values (numeric vector).
#' @param y          Modelled values (numeric vector, same length as x).
#' @param stats      Named list from .calc_metrics().
#' @param title      Plot title / model name.
#' @param lim        Length-2 numeric: shared axis limits c(lo, hi).
#' @param xlabel     x-axis label string.
#' @param ylabel     y-axis label string.
#' @return A ggplot2 object.
.plot_validation <- function(x, y, stats, title, lim, xlabel, ylabel) {

  # vmax: upper colour-scale limit = max(1, n / 10)  — mirrors Python's vmax
  vmax <- max(1, length(x) / 10)

  # Colour-bar ticks: powers of 10 that fit ≤ vmax, then always append vmax
  brk_cands <- c(1, 10, 100, 1000, 10000)
  brk       <- brk_cands[brk_cands <= vmax]
  if (!length(brk)) brk <- 1
  if (tail(brk, 1) != vmax) brk <- c(brk, vmax)
  brk_labels <- c(as.character(as.integer(head(brk, -1))),
                  as.character(as.integer(vmax)))

  # 5 evenly spaced axis ticks (mirrors np.linspace(lim[0], lim[1], 5))
  axis_brk <- seq(lim[1], lim[2], length.out = 5)

  # Stats annotation text placed at axes coordinates (0.03, 0.97) = top-left
  label_text <- sprintf(
    "R\u00b2: %.2f\nBias: %.2f\nRMSE: %.1f\nRel_BIAS: %.1f%%",
    stats$R2, stats$Bias, stats$RMSE, stats$Rel_BIAS * 100
  )

  ggplot(data.frame(x = x, y = y), aes(x, y)) +

    # 2-D histogram, log-scaled fill — matches hist2d(..., norm=LogNorm, cmap="viridis")
    stat_bin2d(bins = 50, aes(fill = after_stat(count))) +
    scale_fill_viridis_c(
      name   = "Number of observations",
      trans  = "log10",
      limits = c(1, vmax),
      breaks = brk,
      labels = brk_labels
    ) +

    # 1:1 reference line — matches plt.plot(lim, lim, "--", color="gray")
    geom_abline(linetype = "dashed", colour = "grey50", linewidth = 0.8) +

    # Axis limits and ticks — matches plt.xticks / plt.yticks + plt.xlim / plt.ylim
    scale_x_continuous(limits = lim, breaks = axis_brk) +
    scale_y_continuous(limits = lim, breaks = axis_brk) +

    labs(title = title, x = xlabel, y = ylabel) +

    # Stats annotation box — matches ax.text(0.03, 0.97, ..., va="top")
    annotate("label",
             x     = lim[1] + diff(lim) * 0.03,
             y     = lim[2] - diff(lim) * 0.03,
             label = label_text,
             hjust = 0, vjust = 1,
             size  = 4,
             fill  = "white", alpha = 0.8, label.size = 0.3) +

    # Minimal theme, no grid — matches plt.grid(False)
    theme_minimal(base_size = 13) +
    theme(
      panel.grid  = element_blank(),
      plot.title  = element_text(size = 15)
    )
}


# ── Public validation functions ───────────────────────────────────────────────

#' Validate modelled HNW against observations from the Mag25 dataset.
#' Exact R equivalent of validate_hnw_mag25() in HNW_validation_helper.py.
#'
#' @param df         Data frame with columns: time, <obs_col>, <mod_col>.
#' @param model_name String used as the plot title.
#' @param obs_col    Name of the observed HNW column  (default "HNW_obs").
#' @param mod_col    Name of the modelled HNW column  (default "HNW_mod").
#' @param save_dir   Directory for saving the PNG; NULL = do not save.
#' @param filename   PNG filename (used when save_dir is not NULL).
#' @param full_season If FALSE (default) restrict to the Nov–Apr snow season.
#' @return Named list of metrics (invisibly).
validate_hnw_mag25 <- function(df,
                               model_name,
                               obs_col     = "HNW_obs",
                               mod_col     = "HNW_mod",
                               save_dir    = NULL,
                               filename    = "hnw_validation.png",
                               full_season = FALSE) {

  df <- .filter_season(df, full_season)

  # Drop NA rows, keep only non-negative observations, remove non-finite values
  # Caution: >= 0 (not > 0) is intentional — matches the Python comment exactly
  df_valid <- df[!is.na(df[[obs_col]]) & !is.na(df[[mod_col]]), ]
  df_valid <- df_valid[df_valid[[obs_col]] >= 0, ]
  df_valid <- df_valid[is.finite(df_valid[[obs_col]]) & is.finite(df_valid[[mod_col]]), ]

  x <- df_valid[[obs_col]]
  y <- df_valid[[mod_col]]

  stats   <- .calc_metrics(x, y)
  stats$N <- nrow(df_valid)
  print(stats)

  p <- .plot_validation(x, y, stats,
                        title  = model_name,
                        lim    = c(0, 100),
                        xlabel = "Observed HNW (mm)",
                        ylabel = "Modeled HNW (mm)")

  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    save_path <- file.path(save_dir, filename)
    ggsave(save_path, plot = p, dpi = 300, width = 8, height = 7,
           units = "in", bg = "white")
    cat(sprintf("Plot saved to: %s\n", save_path))
  }

  print(p)
  invisible(stats)
}


#' Validate modelled SWE against observations from the Mag25 dataset.
#' Exact R equivalent of validate_swe_mag25() in HNW_validation_helper.py.
#'
#' @inheritParams validate_hnw_mag25
#' @param obs_col Name of the observed SWE column (default "SWE_obs").
#' @param mod_col Name of the modelled SWE column (default "SWE_mod").
validate_swe_mag25 <- function(df,
                               model_name,
                               obs_col     = "SWE_obs",
                               mod_col     = "SWE_mod",
                               save_dir    = NULL,
                               filename    = "swe_validation.png",
                               full_season = FALSE) {

  df <- .filter_season(df, full_season)

  df_valid <- df[!is.na(df[[obs_col]]) & !is.na(df[[mod_col]]), ]
  df_valid <- df_valid[df_valid[[obs_col]] >= 0, ]

  # Mirrors Python's explicit print before the finite filter
  cat(sprintf("Number of valid observations after filtering: %d\n", nrow(df_valid)))

  df_valid <- df_valid[is.finite(df_valid[[obs_col]]) & is.finite(df_valid[[mod_col]]), ]

  x <- df_valid[[obs_col]]
  y <- df_valid[[mod_col]]

  stats   <- .calc_metrics(x, y)
  stats$N <- nrow(df_valid)
  print(stats)

  p <- .plot_validation(x, y, stats,
                        title  = model_name,
                        lim    = c(0, 1000),
                        xlabel = "Observed SWE (mm)",
                        ylabel = "Modeled SWE (mm)")

  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    save_path <- file.path(save_dir, filename)
    ggsave(save_path, plot = p, dpi = 300, width = 8, height = 7,
           units = "in", bg = "white")
    cat(sprintf("Plot saved to: %s\n", save_path))
  }

  print(p)
  invisible(stats)
}
