library(ggplot2)
library(nixmass)

  
library(ggplot2)
library(tidyr)

# Load your data (replace with your CSV path)
df <- read.csv("HNW_validation/validation_input/1617/Gadmen.csv")



run_data <- data.frame(
  date = as.POSIXct(df$date),
  hs = as.numeric(df$HS)
)



# Set model parameters (adjust as needed)
model_opts <- list(
    rho.max   = 451.6977806582531,
    rho.null  = 90.0,
    c.ov      = 2.746976858091589e-0,
    k.ov      = 0.38,
    k         = 0.020385468323087456,
    tau       = 0.01,
    eta.null  = 8.5e6
)



# Run the model
res_dsnow <- swe.delta.snow(
  data = run_data,
  model_opts = model_opts,
  layer = TRUE,
  dyn_rho_max = FALSE
)

res_hs2swe <- nixmass(
  data = run_data,
  model = "hs2swe"
)

# Output results
out <- data.frame(
  date       = df$date,
  hs         = df$hs,
  swe_dsnow  = as.numeric(res_dsnow$SWE),
  swe_hs2swe = as.numeric(res_hs2swe$SWE))




# Create plotting data
plot_df <- out %>%
  mutate(date = as.POSIXct(date)) %>%
  select(date, hs, swe_dsnow, swe_hs2swe) %>%
  pivot_longer(
    cols = c(hs, swe_dsnow, swe_hs2swe),
    names_to = "variable",
    values_to = "value"
  )

# Plot HS and both SWE models
ggplot(plot_df, aes(x = date, y = value, color = variable)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c(
      hs = "black",
      swe_dsnow = "steelblue",
      swe_hs2swe = "darkred"
    ),
    labels = c(
      hs = "HS",
      swe_dsnow = "SWE (ΔSnow)",
      swe_hs2swe = "SWE (HS → SWE)"
    )
  ) +
  labs(
    title = "HS and SWE Model Comparison",
    x = "Date",
    y = "HS / SWE",
    color = ""
  ) +
  theme_minimal() +
  theme(legend.position = "top")