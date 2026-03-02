# ============================================================
# ΔSNOW diagnostics plot (FULL WORKING SCRIPT)
# - Full season: HS (m) + ΔHS exceedances (m/day) on sec axis
# - Drenching shading (red)
# - New layer creation vlines (green)
# - Window around k-th largest |ΔHS| (default k=3 below)
# - Window: HS_obs vs HS_model + stacked layer bars (1 bar/day)
# - Top-layer thickness time series in window
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(nixmass)

# ----------------------------
# 0) User settings
# ----------------------------
csv_path <- "/Users/jakobwerkgarner/code/mt_dsnow/par_sens/SNOWPACK_data_seasons_daily/AXLIZ1_2024_2025_daily_HS_SWE.csv"
tau <- 0.024   # [m] = 2.4 cm
k_event <- 5   # k-th largest |ΔHS| event (2 = second, 3 = third, ...)

# ----------------------------
# 1) Load data (HS in cm in CSV; model expects cm)
# ----------------------------
df <- read.csv(csv_path)

run_data <- data.frame(
  date = as.POSIXct(df$timestamp , tz = "UTC"),
  hs   = as.numeric(df$HS_mod)   # m
)

# ----------------------------
# 2) Run ΔSNOW
# ----------------------------
res_dsnow <- swe.delta.snow(
  data = run_data,
  layer = TRUE,
  dyn_rho_max = FALSE
)

# ----------------------------
# 3) Build plot dataframe (convert HS and ΔHS to meters)
# ----------------------------
plot_df <- data.frame(
  date    = as.POSIXct(run_data$date),
  hs_cm   = as.numeric(run_data$hs),
  SWE     = as.numeric(res_dsnow$SWE),
  process = as.character(res_dsnow$processes)
) %>%
  mutate(
    process = ifelse(is.na(process), "", process),
    proc_l  = tolower(process),
    is_drench   = grepl("drench", proc_l),
    is_newlayer = grepl("^create new layer", proc_l),
    
    hs  = hs_cm / 100.0,                    # [m]
    dHS = (hs_cm - lag(hs_cm)) / 100.0      # [m/day]
  ) %>%
  mutate(
    dHS = ifelse(is.na(dHS), 0, dHS),
    outside_tau = abs(dHS) > tau
  )

# ----------------------------
# 4) Helper: contiguous TRUE-intervals (for drenching shading)
# ----------------------------
get_intervals <- function(dates, mask) {
  if (length(mask) == 0 || !any(mask)) {
    return(data.frame(
      start = as.POSIXct(character()),
      end   = as.POSIXct(character())
    ))
  }
  r <- rle(mask)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  data.frame(
    start = dates[starts[r$values]],
    end   = dates[ends[r$values]]
  )
}

drench_intervals <- get_intervals(plot_df$date, plot_df$is_drench)

# ----------------------------
# 5) Secondary axis transform (ΔHS mapped onto HS axis)
# ----------------------------
dHS_lim <- max(tau, max(abs(plot_df$dHS), na.rm = TRUE))
hs_min <- min(plot_df$hs, na.rm = TRUE)
hs_max <- max(plot_df$hs, na.rm = TRUE)

scale_fac <- (hs_max - hs_min) / (2 * dHS_lim)
shift_fac <- (hs_max + hs_min) / 2

plot_df <- plot_df %>%
  mutate(dHS_on_HS = shift_fac + dHS * scale_fac)

tau_plus_on_HS  <- shift_fac + tau * scale_fac
tau_minus_on_HS <- shift_fac - tau * scale_fac

# ----------------------------
# 6) Full-season plot
# ----------------------------
p_full <- ggplot(plot_df, aes(x = date)) +
  geom_rect(
    data = drench_intervals,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "red",
    alpha = 0.08
  ) +
  geom_vline(
    data = subset(plot_df, is_newlayer),
    aes(xintercept = date),
    color = "green3",
    linewidth = 0.6,
    alpha = 0.6
  ) +
  geom_line(aes(y = hs), color = "black", linewidth = 1) +
  geom_point(
    data = subset(plot_df, outside_tau),
    aes(y = dHS_on_HS),
    color = "steelblue",
    size = 1.3,
    alpha = 0.85
  ) +
  geom_hline(yintercept = tau_plus_on_HS,  linetype = "dashed", color = "steelblue", linewidth = 0.6) +
  geom_hline(yintercept = tau_minus_on_HS, linetype = "dashed", color = "steelblue", linewidth = 0.6) +
  labs(
    title = "ΔSNOW: HS (links) und ΔHS (rechts) — alles in Metern",
    subtitle = "Grün: create new layer; Rot: Drenching; Punkte: |ΔHS| > τ; τ = 0.024 m",
    x = "Datum",
    y = "Schneehöhe HS [m]"
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(
      trans = ~ (. - shift_fac) / scale_fac,
      name  = "ΔHS [m pro Tag]"
    )
  ) +
  theme_minimal()

print(p_full)

# ============================================================
# PART 2: Window around the k-th largest |ΔHS| + layer plot
# ============================================================

# ----------------------------
# 7) Find index of k-th largest |ΔHS|
# ----------------------------
abs_dHS <- abs(plot_df$dHS)
abs_dHS[1] <- NA  # ignore first timestep (lag undefined)

ord <- order(abs_dHS, decreasing = TRUE, na.last = NA)
if (length(ord) < k_event) stop("Not enough valid ΔHS values for chosen k_event.")

idx_evt <- ord[k_event]

# ----------------------------
# 8) Window: -4 to +9 days around event index
# ----------------------------
i1 <- max(1, idx_evt - 4)
i2 <- min(nrow(plot_df), idx_evt + 9)
win_idx <- i1:i2

# ----------------------------
# 9) Window dataframe with process flags (all in meters)
# ----------------------------
hs_tot_m <- colSums(res_dsnow$h, na.rm = TRUE) / 100.0  # [m] modeled total HS (sum layers)

plot_win <- data.frame(
  idx     = win_idx,
  date    = as.POSIXct(run_data$date[win_idx]),
  hs_obs  = as.numeric(run_data$hs[win_idx]) / 100.0,    # [m]
  hs_mod  = as.numeric(hs_tot_m[win_idx]),               # [m]
  dHS     = as.numeric(plot_df$dHS[win_idx]),            # [m/day]
  process = as.character(res_dsnow$processes[win_idx])
) %>%
  mutate(
    process = ifelse(is.na(process), "", process),
    proc_l  = tolower(process),
    is_drench   = grepl("drench", proc_l),
    is_newlayer = grepl("^create new layer", proc_l)
  )

drench_int_win <- get_intervals(plot_win$date, plot_win$is_drench)
newlayer_dates_win <- plot_win$date[plot_win$is_newlayer]

# ----------------------------
# 10) Window plot: HS obs vs HS model + markers
# ----------------------------
p_win <- ggplot(plot_win, aes(x = date)) +
  geom_rect(
    data = drench_int_win,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "red",
    alpha = 0.10
  ) +
  geom_vline(
    xintercept = as.numeric(newlayer_dates_win),
    color = "green3",
    linewidth = 0.7,
    alpha = 0.6
  ) +
  geom_vline(
    xintercept = as.numeric(run_data$date[idx_evt]),
    color = "orange",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  geom_line(aes(y = hs_obs, color = "HS beobachtet"), linewidth = 1) +
  geom_line(aes(y = hs_mod, color = "HS Modell (Summe Layer)"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("HS beobachtet" = "black",
                                "HS Modell (Summe Layer)" = "steelblue")) +
  labs(
    title = sprintf("Fenster um das %d.-größte |ΔHS|-Ereignis (−4 bis +9 Tage)", k_event),
    subtitle = sprintf("Ereignisindex = %d (orange). |ΔHS| = %.3f m/Tag. Grün: create new layer. Rot: Drenching.",
                       idx_evt, abs(plot_df$dHS[idx_evt])),
    x = "Datum",
    y = "Schneehöhe [m]",
    color = ""
  ) +
  theme_minimal() +
  theme(legend.position = "top")

print(p_win)

# ----------------------------
# 11) Layer plot in window: stacked bars (1 bar/day), colored by layer
#     Highest layer index on top
# ----------------------------
H_win_m <- res_dsnow$h[, win_idx, drop = FALSE] / 100.0  # [m]
layer_df <- as.data.frame(H_win_m)
colnames(layer_df) <- as.character(win_idx)

layer_long <- layer_df %>%
  mutate(layer = row_number()) %>%
  pivot_longer(
    cols = -layer,
    names_to = "idx",
    values_to = "hs_layer"
  ) %>%
  mutate(
    idx = as.integer(idx),
    date = as.POSIXct(run_data$date[idx]),
    date_day = as.Date(date),
    hs_layer = ifelse(is.na(hs_layer), 0, hs_layer)
  )

max_layer <- max(layer_long$layer, na.rm = TRUE)
layer_long <- layer_long %>%
  mutate(layer_rev = factor(layer, levels = rev(seq_len(max_layer))))

p_layers <- ggplot(layer_long, aes(x = date_day, y = hs_layer, fill = layer_rev)) +
  geom_col(
    position = "stack",
    width = 1.0,
    color = "black",
    linewidth = 0.15
  ) +
  scale_x_date(date_breaks = "1 day", date_labels = "%b %d") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_viridis_d(option = "turbo", direction = 1) +
  labs(
    title = "Schichtaufbau im Fenster (gestapelte Balken)",
    subtitle = sprintf("Fenster um Index %d; höchster Layerindex oben", idx_evt),
    x = "Datum",
    y = "Schichtdicke hs_i [cm]",
    fill = "Layer"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p_layers)

# ----------------------------
# 12) Choose which layer to plot in the last figure
#     - If NULL: use "highest existing layer per day" (as before)
#     - If integer: plot thickness of that fixed layer index
# ----------------------------
top_layer_id <- 13  # e.g. 12  (set NULL for dynamic top layer)

if (is.null(top_layer_id)) {
  # dynamic: highest existing layer per day
  top_layer_ts <- layer_long %>%
    group_by(date_day) %>%
    summarise(
      layer_id = if (any(hs_layer > 0)) max(layer[hs_layer > 0]) else NA_integer_,
      hs_top_m = if (any(hs_layer > 0)) hs_layer[layer == max(layer[hs_layer > 0])] else 0,
      .groups = "drop"
    )
  
  plot_title <- "Dicke der obersten Schicht (höchster Layerindex pro Tag)"
} else {
  # fixed: a user-chosen layer index
  top_layer_ts <- layer_long %>%
    filter(layer == top_layer_id) %>%
    group_by(date_day) %>%
    summarise(
      layer_id = top_layer_id,
      hs_top_m = sum(hs_layer, na.rm = TRUE),  # should be single value; sum is robust
      .groups = "drop"
    ) %>%
    complete(
      date_day = seq(min(layer_long$date_day), max(layer_long$date_day), by = "day"),
      fill = list(hs_top_m = 0, layer_id = top_layer_id)
    )
  
  plot_title <- sprintf("Dicke von Layer %d (fester Layerindex)", top_layer_id)
}

p_toplayer <- ggplot(top_layer_ts, aes(x = date_day, y = hs_top_m)) +
  geom_line(color = "purple4", linewidth = 1) +
  geom_point(color = "purple4", size = 2) +
  scale_x_date(date_breaks = "1 day", date_labels = "%b %d") +
  scale_y_continuous(labels = function(x) x * 100) +
  labs(
    title = plot_title,
    subtitle = sprintf("Fenster um Index %d", idx_evt),
    x = "Datum",
    y = "Schichtdicke [cm]"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p_toplayer)
print(p_toplayer)