"""
Plot LWD station locations on a map of Austria, coloured by elevation.
Aspect ratio is corrected so the map is not distorted.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH   = os.path.join(SCRIPT_DIR, "LWD_all", "station_summary.csv")
OUT_PATH   = os.path.join(SCRIPT_DIR, "LWD_all", "station_map.png")

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv(CSV_PATH, sep=";", decimal=",", encoding="latin-1")
df["latitude"]   = pd.to_numeric(df["latitude"],   errors="coerce")
df["longitude"]  = pd.to_numeric(df["longitude"],  errors="coerce")
df["altitude_m"] = pd.to_numeric(df["altitude_m"], errors="coerce")
df = df.dropna(subset=["latitude", "longitude", "altitude_m"])

# ── Map extent: full Tirol region with small margin ───────────────────────────
margin = 0.3
lon_min = df["longitude"].min() - margin
lon_max = df["longitude"].max() + margin
lat_min = df["latitude"].min()  - margin * 0.6
lat_max = df["latitude"].max()  + margin * 0.6

lon_range   = lon_max - lon_min
lat_range   = lat_max - lat_min
central_lat = (lat_min + lat_max) / 2

# In PlateCarree a longitude degree is shorter than a latitude degree by cos(lat).
# Correct the figure width so 1° lon == 1° lat on screen.
fig_height = 7
fig_width  = fig_height * (lon_range * np.cos(np.radians(central_lat))) / lat_range

# ── Colour scale ──────────────────────────────────────────────────────────────
cmap = plt.get_cmap("plasma")
norm = Normalize(vmin=df["altitude_m"].min(), vmax=df["altitude_m"].max())

# ── Figure ────────────────────────────────────────────────────────────────────
proj = ccrs.PlateCarree()
fig, ax = plt.subplots(figsize=(fig_width, fig_height),
                       subplot_kw={"projection": proj})
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=proj)

# Background (50 m NaturalEarth — already cached after first run)
RES = "50m"
ax.add_feature(cfeature.NaturalEarthFeature(
    "physical", "land", RES, facecolor="#f0ede8"), zorder=0)
ax.add_feature(cfeature.NaturalEarthFeature(
    "physical", "ocean", RES, facecolor="#d0e8f5"), zorder=0)
ax.add_feature(cfeature.NaturalEarthFeature(
    "physical", "lakes", RES,
    facecolor="#d0e8f5", edgecolor="#90b8d0", linewidth=0.4), zorder=1)
ax.add_feature(cfeature.NaturalEarthFeature(
    "physical", "rivers_lake_centerlines", RES,
    facecolor="none", edgecolor="#a0c8e0", linewidth=0.5), zorder=1)
ax.add_feature(cfeature.NaturalEarthFeature(
    "cultural", "admin_0_countries", RES,
    facecolor="none", edgecolor="#777777", linewidth=0.8), zorder=2)
ax.add_feature(cfeature.NaturalEarthFeature(
    "cultural", "admin_1_states_provinces", RES,
    facecolor="none", edgecolor="#aaaaaa", linewidth=0.4), zorder=2)

# Gridlines
gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="gray",
                  alpha=0.5, linestyle="--", zorder=3)
gl.top_labels   = False
gl.right_labels = False
gl.xlabel_style = {"size": 8}
gl.ylabel_style = {"size": 8}

# ── Stations ──────────────────────────────────────────────────────────────────
ax.scatter(
    df["longitude"], df["latitude"],
    c=df["altitude_m"], cmap=cmap, norm=norm,
    s=60, zorder=5, edgecolors="white", linewidths=0.6,
    transform=proj,
)

outline = [pe.withStroke(linewidth=2.5, foreground="white")]
for _, row in df.iterrows():
    ax.text(
        row["longitude"] + 0.015, row["latitude"] + 0.01,
        row["station_id"],
        fontsize=5.5, color="#222222",
        path_effects=outline,
        transform=proj, zorder=6,
        ha="left", va="bottom",
    )

# ── Colourbar ─────────────────────────────────────────────────────────────────
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, orientation="vertical",
                    fraction=0.025, pad=0.02, aspect=30)
cbar.set_label("Elevation (m a.s.l.)", fontsize=10)
cbar.ax.tick_params(labelsize=9)

# ── Title & stats box ─────────────────────────────────────────────────────────
ax.set_title(
    f"LWD Tirol – SNOWPACK forcing stations  (n = {len(df)})\n"
    "Coloured by elevation",
    fontsize=12, fontweight="bold", pad=10,
)
stats = (
    f"Elevation range:  {df['altitude_m'].min():.0f} – {df['altitude_m'].max():.0f} m\n"
    f"Median elevation: {df['altitude_m'].median():.0f} m\n"
    f"Merged pairs:     {df['was_merged'].eq('yes').sum()} / {len(df)} stations"
)
ax.text(
    0.01, 0.02, stats,
    transform=ax.transAxes, fontsize=8,
    va="bottom", ha="left", zorder=7,
    bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
              alpha=0.85, edgecolor="#cccccc"),
)

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=180, bbox_inches="tight")
print(f"Saved  {OUT_PATH}")
plt.show()
