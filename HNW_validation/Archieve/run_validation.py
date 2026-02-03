#!/usr/bin/env python3
"""
===============================================================================
ΔSNOW FULL PIPELINE — CLEAN, MODULAR, CONSISTENT NAMING

Pipeline:
---------
1. Load Mag25_all.nc
2. Extract seasonal station data (Nov–Apr)
3. For each station × season:
       → Run ΔSNOW R model
       → Merge SWE_mod, HNW_mod, HS_mod
4. Build ONE large DataFrame (no unnecessary CSV I/O)
5. Validate SWE + HNW using clean, well-documented routines
6. Optionally write outputs to disk

Author: Jakob Werkgarner
===============================================================================
"""

# =============================================================================
# Imports
# =============================================================================

import os
from pathlib import Path
import pandas as pd
import numpy as np
import xarray as xr
import subprocess
from io import StringIO
import tempfile
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


# =============================================================================
# Configuration
# =============================================================================

BASE_DIR = Path("/Users/jakobwerkgarner/code/mt_dsnow")
os.chdir(BASE_DIR)

NC_FILE = BASE_DIR / "calibration/calibration_data/raw_data/Mag25/SLF_dataset/Mag25_all.nc"

R_BIN = "/usr/local/bin/Rscript"
R_RUNNER = BASE_DIR / "calibration/helpers/minimal_delta_snow_runner.R"

OUT_DIR = BASE_DIR / "HNW_validation/validation_results"
OUT_DIR.mkdir(exist_ok=True, parents=True)

PARAMS = {
    "rho.max": 401,
    "rho.null": 81,
    "c.ov": 0.0005,
    "k.ov": 0.25,
    "k": 0.03,
    "tau": 0.024,
    "eta.null": 8500000.0,
}

SEASONS = [f"{y:02d}{(y+1)%100:02d}" for y in range(16, 22)]


# =============================================================================
# Helper: Assign season (Nov–Apr)
# =============================================================================

def assign_season(ts):
    """Return season string YYZZ for dates in Nov–Apr."""
    m, y = ts.month, ts.year

    if m in (11, 12):
        start, end = y, y + 1
    elif m in (1, 2, 3, 4):
        start, end = y - 1, y
    else:
        return None

    return f"{start % 100:02d}{end % 100:02d}"


# =============================================================================
# Step 1 — Extract seasonal data for all stations
# =============================================================================

def load_seasonal_data():
    """Load Mag25 file and return a clean pandas DataFrame."""

    ds = xr.open_dataset(NC_FILE).rename({"time": "date"})
    dates = pd.to_datetime(ds.date.values)

    ds = ds.assign_coords(season=("date", [assign_season(t) for t in dates]))
    ds = ds.where(ds.season.isin(SEASONS), drop=True)

    df = ds.to_dataframe().reset_index()
    df["date"] = pd.to_datetime(df["date"]).dt.strftime("%Y-%m-%d")

    # Renaming here makes naming consistent across entire pipeline
    df = df.rename(columns={
        "HS": "HS_obs",
        "HNW": "HNW_obs",
        "SWE": "SWE_obs",
    })

    return df[["date", "station", "season", "HS_obs", "HNW_obs", "SWE_obs"]]


# =============================================================================
# Step 2 — Run R ΔSNOW model
# =============================================================================

def run_r_model(df):
    """
    Run the ΔSNOW R model for one season × station.
    R EXPECTS a column named 'hs'.
    We convert HS_obs → hs ONLY for the R input file.
    """

    # Build R input
    df_r = df[["date", "HS_obs"]].copy()
    df_r = df_r.rename(columns={"HS_obs": "hs"})
    df_r.at[df_r.index[0], "hs"] = 0

    # Temporary CSV
    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp:
        tmp_path = Path(tmp.name)
        df_r.to_csv(tmp_path, index=False)

    # Build command
    cmd = [R_BIN, str(R_RUNNER), "--in", str(tmp_path)]
    for k, v in PARAMS.items():
        cmd += [f"--{k}", str(v)]

    # Run R
    proc = subprocess.run(cmd, capture_output=True, text=True)

    os.remove(tmp_path)

    if proc.returncode != 0:
        raise RuntimeError(f"R model failed:\nSTDERR:\n{proc.stderr}")

    # Parse output
    df_out = pd.read_csv(StringIO(proc.stdout))
    df_out["date"] = pd.to_datetime(df_out["date"]).dt.strftime("%Y-%m-%d")

    # Rename outputs consistently
    df_out = df_out.rename(columns={
        "hs": "HS_mod",
        "swe_mod": "SWE_mod",
    })

    return df_out[["date", "HS_mod", "SWE_mod"]]


# =============================================================================
# Step 3 — Process all stations × seasons
# =============================================================================

def run_all_delta_snow(df):
    """Core processing loop — everything stays in memory."""
    results = []

    for station in sorted(df.station.unique()):
        df_station = df[df.station == station]

        for season in SEASONS:
            sub = df_station[df_station.season == season]

            if len(sub) == 0:
                continue

            df_r = run_r_model(sub)
            merged = sub.merge(df_r, on="date", how="left")

            # The correct definition: ΔSWE
            merged["HNW_mod"] = merged["SWE_mod"].diff()

            results.append(merged)

    return pd.concat(results, ignore_index=True)


# =============================================================================
# Step 4 — Validation Plots
# =============================================================================

def validate_swe(df):
    df = df.dropna(subset=["SWE_obs", "SWE_mod"])
    df = df[(df["SWE_obs"] > 0) & (df["SWE_mod"] > 0)]

    x = df["SWE_obs"].values
    y = df["SWE_mod"].values

    rmse = np.sqrt(np.mean((y - x) ** 2))
    pbias = 100 * np.sum((y - x)) / np.sum(x)
    r2 = np.corrcoef(x, y)[0, 1] ** 2
    slope, intercept = np.polyfit(x, y, 1)

    plt.figure(figsize=(7, 7))
    plt.hist2d(x, y, bins=60, range=[[0, 1000], [0, 1000]], norm=LogNorm(), cmap="jet")
    plt.colorbar(label="Count")
    plt.plot([0, 1000], [0, 1000], "--", color="gray")
    plt.plot([0, 1000], slope * np.array([0, 1000]) + intercept, "--", color="red")
    plt.xlabel("Observed SWE [mm]")
    plt.ylabel("Modeled SWE [mm]")
    plt.title("ΔSNOW SWE Validation")
    plt.tight_layout()
    plt.show()


def validate_hnw(df):
    df = df.dropna(subset=["HNW_obs", "HNW_mod"])
    df = df[df["HNW_obs"] > 0]

    x = df["HNW_mod"].values
    y = df["HNW_obs"].values

    rmse = np.sqrt(np.mean((x - y) ** 2))
    pbias = 100 * np.sum((x - y)) / np.sum(y)
    r2 = np.corrcoef(x, y)[0, 1] ** 2
    slope, intercept = np.polyfit(x, y, 1)

    plt.figure(figsize=(7, 7))
    plt.hist2d(x, y, bins=60, range=[[0, 100], [0, 100]], norm=LogNorm(), cmap="jet")
    plt.colorbar(label="Count")
    plt.plot([0, 100], [0, 100], "--", color="gray")
    plt.plot([0, 100], slope * np.array([0, 100]) + intercept, "--", color="red")
    plt.xlabel("Modeled HNW [mm]")
    plt.ylabel("Observed HNW [mm]")
    plt.title("ΔSNOW HNW Validation")
    plt.tight_layout()
    plt.show()


# =============================================================================
# Main
# =============================================================================

def main():
    print("Loading seasonal data ...")
    df = load_seasonal_data()

    print("Running ΔSNOW for all stations ...")
    df_out = run_all_delta_snow(df)

    print("Saving merged dataset ...")
    df_out.to_csv(OUT_DIR / "ALL_STATIONS_ALL_YEARS.csv", index=False)

    print("Validating SWE ...")
    validate_swe(df_out)

    print("Validating HNW ...")
    validate_hnw(df_out)


if __name__ == "__main__":
    main()