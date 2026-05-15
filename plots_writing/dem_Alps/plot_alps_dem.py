"""
Hypsometric elevation map of the Alps with conical axes (curved frame).
Projection: Lambert Conformal Conic – no edge distortion.
The map boundary follows the projection's natural curved shape.
"""

import math
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LightSource
from matplotlib.path import Path as MplPath
import rasterio
from rasterio.enums import Resampling
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.geoaxes import GeoAxes
from shapely.geometry import Polygon

HERE = Path(__file__).parent
DEM_PATH = HERE / "alps_dem_30m.tif"
if not DEM_PATH.exists():
    DEM_PATH = HERE.parent / "alps_dem_30m.tif"
OUT_PATH = HERE / "alps_elevation_map_conical.png"

RESAMPLE_FACTOR = 10

# ── 1. Load DEM ───────────────────────────────────────────────────────────────
print("Loading DEM …")
with rasterio.open(DEM_PATH) as src:
    nodata = src.nodata
    out_h = src.height // RESAMPLE_FACTOR
    out_w = src.width  // RESAMPLE_FACTOR
    dem = src.read(1, out_shape=(out_h, out_w),
                   resampling=Resampling.average).astype(np.float32)
    bounds = src.bounds   # (left, bottom, right, top) in WGS84

if nodata is not None:
    dem[dem == nodata] = np.nan
dem = np.clip(dem, 0, None)

lon_min, lat_min, lon_max, lat_max = bounds
central_lon = (lon_min + lon_max) / 2
central_lat = (lat_min + lat_max) / 2

# ── 2. Projection (no distortion) ─────────────────────────────────────────────
proj = ccrs.LambertConformal(
    central_longitude=central_lon,
    central_latitude=central_lat,
    standard_parallels=(44.0, 48.0),
)

# ── 3. Hypsometric colormap ───────────────────────────────────────────────────
hyp_colors = [
    (0.00, "#4a7c59"), (0.08, "#7da87b"), (0.18, "#c8b866"),
    (0.32, "#a07850"), (0.52, "#7a5c40"), (0.68, "#6e5c50"),
    (0.80, "#a0968c"), (0.92, "#d4cfc9"), (1.00, "#f5f5f5"),
]
cmap = mcolors.LinearSegmentedColormap.from_list("hypsometric", hyp_colors)
norm = mcolors.Normalize(vmin=0, vmax=4500)

# ── 4. Hillshade ──────────────────────────────────────────────────────────────
ls = LightSource(azdeg=315, altdeg=35)
intensity = ls.hillshade(dem, vert_exag=3, dx=300, dy=300)
rgb = cmap(norm(dem))[:, :, :3]
blended = np.clip(rgb * intensity[:, :, np.newaxis], 0, 1)

# ── 5. Create figure with a conical (curved) boundary ─────────────────────────
fig = plt.figure(figsize=(16, 12), dpi=200)
ax = fig.add_subplot(1, 1, 1, projection=proj)

# Set the extent in geographic coordinates (this defines the curved region)
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

# --- This is the key: make the axes follow the projection's natural curve ---
# We compute the boundary polygon of the projection for the given extent
# and set it as the map's clipping path.
def get_projection_boundary(ax, lon_min, lon_max, lat_min, lat_max, n_points=100):
    """Return a Shapely polygon of the projected boundary (curved)."""
    # Create a grid of points on the geographic rectangle's border
    lons = []
    lats = []
    # bottom edge
    for i in range(n_points):
        lons.append(lon_min + (lon_max - lon_min) * i / (n_points - 1))
        lats.append(lat_min)
    # right edge
    for i in range(n_points):
        lons.append(lon_max)
        lats.append(lat_min + (lat_max - lat_min) * i / (n_points - 1))
    # top edge
    for i in range(n_points):
        lons.append(lon_max - (lon_max - lon_min) * i / (n_points - 1))
        lats.append(lat_max)
    # left edge
    for i in range(n_points):
        lons.append(lon_min)
        lats.append(lat_max - (lat_max - lat_min) * i / (n_points - 1))
    # Transform to projection coordinates
    x, y = ax.projection.transform_points(ccrs.PlateCarree(),
                                          np.array(lons), np.array(lats))[:, :2].T
    return Polygon(list(zip(x, y)))

# Get the curved boundary polygon
boundary_poly = get_projection_boundary(ax, lon_min, lon_max, lat_min, lat_max)
# Convert to matplotlib path for clipping
boundary_path = MplPath(np.array(boundary_poly.exterior.coords))

# Apply clipping to the axes
ax.set_boundary(boundary_path, transform=ax.transData)

# ── 6. Draw DEM (reprojected on the fly) ──────────────────────────────────────
print("Rendering …")
ax.imshow(blended, origin="upper",
          extent=[lon_min, lon_max, lat_min, lat_max],
          transform=ccrs.PlateCarree(), interpolation="lanczos", zorder=1)

# ── 7. Borders, coastline, gridlines ──────────────────────────────────────────
ax.add_feature(cfeature.BORDERS.with_scale("10m"),
               linewidth=0.7, edgecolor="#111111", zorder=2)
ax.add_feature(cfeature.COASTLINE.with_scale("10m"),
               linewidth=0.6, edgecolor="#111111", zorder=2)

gl = ax.gridlines(draw_labels=True, linewidth=0.35,
                  color="white", alpha=0.55, linestyle="--",
                  x_inline=False, y_inline=False,
                  crs=ccrs.PlateCarree())
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {"size": 8, "color": "#333333"}
gl.ylabel_style = {"size": 8, "color": "#333333"}

# ── 8. Colorbar ───────────────────────────────────────────────────────────────
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, orientation="vertical",
                    fraction=0.022, pad=0.02, shrink=0.9)
cbar.set_label("Elevation (m a.s.l.)", fontsize=10)

ax.set_title("Alps — Conical axes (curved frame) with Lambert Conformal Conic\n(no edge distortion)",
             fontsize=12, fontweight="bold")

plt.savefig(OUT_PATH, dpi=200, bbox_inches="tight", facecolor="white")
print(f"Saved: {OUT_PATH}")