"""
Download Copernicus DEM GLO-30 (30 m) tiles for the whole Alps and merge into one GeoTIFF.
Tiles come from the public AWS S3 bucket (no auth needed).
Coverage: lon 4–17 E, lat 43–49 N  (whole Alpine arc)
"""

import os
import time
import requests
import numpy as np
from pathlib import Path
import rasterio
from rasterio.merge import merge
from rasterio.crs import CRS

TILE_DIR = Path("dem_tiles")
OUTPUT = Path("alps_dem_30m.tif")

# Alps bounding box (1°×1° tiles)
LAT_MIN, LAT_MAX = 43, 48   # inclusive; tile N43 covers 43–44°N, etc.
LON_MIN, LON_MAX = 4, 16    # inclusive; tile E004 covers 4–5°E, etc.

BASE_URL = (
    "https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com/"
    "Copernicus_DSM_COG_10_{lat_hem}{lat:02d}_00_{lon_hem}{lon:03d}_00_DEM/"
    "Copernicus_DSM_COG_10_{lat_hem}{lat:02d}_00_{lon_hem}{lon:03d}_00_DEM.tif"
)


def tile_url(lat: int, lon: int) -> str:
    lat_hem = "N" if lat >= 0 else "S"
    lon_hem = "E" if lon >= 0 else "W"
    return BASE_URL.format(lat=abs(lat), lon=abs(lon), lat_hem=lat_hem, lon_hem=lon_hem)


def download_tile(lat: int, lon: int) -> Path | None:
    url = tile_url(lat, lon)
    name = f"COP30_N{lat:02d}_E{lon:03d}.tif"
    out = TILE_DIR / name
    if out.exists():
        return out

    for attempt in range(1, 4):
        try:
            resp = requests.get(url, stream=True, timeout=300)
            if resp.status_code == 404:
                return None  # ocean/no-data tile
            resp.raise_for_status()
            tmp = out.with_suffix(".tmp")
            with open(tmp, "wb") as f:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
            tmp.rename(out)
            return out
        except Exception as e:
            if attempt == 3:
                raise
            print(f"    retry {attempt}/3 after error: {e}")
            time.sleep(5)


def main():
    TILE_DIR.mkdir(exist_ok=True)
    tiles_needed = [
        (lat, lon)
        for lat in range(LAT_MIN, LAT_MAX + 1)
        for lon in range(LON_MIN, LON_MAX + 1)
    ]
    print(f"Downloading {len(tiles_needed)} tiles ...")

    tile_paths = []
    for i, (lat, lon) in enumerate(tiles_needed, 1):
        t0 = time.time()
        path = download_tile(lat, lon)
        elapsed = time.time() - t0
        if path:
            size_mb = path.stat().st_size / 1e6
            status = f"{size_mb:.1f} MB" if elapsed > 0.1 else "cached"
            print(f"  [{i:2d}/{len(tiles_needed)}] N{lat:02d} E{lon:03d}  {status}")
            tile_paths.append(path)
        else:
            print(f"  [{i:2d}/{len(tiles_needed)}] N{lat:02d} E{lon:03d}  (no data – ocean/missing)")

    print(f"\nMerging {len(tile_paths)} tiles into {OUTPUT} ...")
    datasets = [rasterio.open(p) for p in sorted(tile_paths)]
    mosaic, transform = merge(datasets)
    meta = datasets[0].meta.copy()
    meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": transform,
        "compress": "lzw",
        "tiled": True,
        "blockxsize": 512,
        "blockysize": 512,
    })
    with rasterio.open(OUTPUT, "w", **meta) as dst:
        dst.write(mosaic)
    for ds in datasets:
        ds.close()

    size_gb = OUTPUT.stat().st_size / 1e9
    print(f"Done: {OUTPUT}  ({size_gb:.2f} GB)")
    print(f"CRS: {meta['crs']}  |  shape: {mosaic.shape[1]} × {mosaic.shape[2]} px")


if __name__ == "__main__":
    main()
