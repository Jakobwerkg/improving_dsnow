#!/usr/bin/env python3
"""
Minimal SWE model runner + SWE/HNW validation plots.

Pipeline:
---------
1. Reads seasonal_data/<year>/<station>.csv
2. Calls the ΔSNOW R model (minimal_delta_snow_runner.R)
3. Merges SWE_mod + HS_mod back into the input file
4. Writes per-season and per-station outputs
5. Compiles a global ALL_STATIONS_ALL_YEARS.csv
6. Validates SWE and HNW with scatter plots + metrics

Author: Jakob Werkgarner
"""

# =============================================================================
# Imports
# =============================================================================

from pathlib import Path
import subprocess
import pandas as pd
import numpy as np
from io import StringIO
import tempfile
import os
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# =============================================================================
# CONFIGURATION
# =============================================================================

PARAMS = {
    "rho.max": 401.0,       # kg m^-3
    "rho.null": 81.0,       # kg m^-3
    "c.ov": 5.1e-4,         # Pa^-1
    "k.ov": 0.38,           # -
    "k": 0.030,             # m^3 kg^-1
    "tau": 0.024,           # m  (2.4 cm)
    "eta.null": 8.5e6       # Pa s
}

BASE_DIR = Path("/Users/jakobwerkgarner/code/mt_dsnow")
os.chdir(BASE_DIR)

DATA_DIR = Path("HNW_validation/validation_input")
OUT_DIR = Path("HNW_validation/validation_output")
OUT_DIR.mkdir(parents=True, exist_ok=True)

R_BIN = "/usr/local/bin/Rscript"
R_RUNNER = BASE_DIR / "calibration/helpers/minimal_delta_snow_runner.R"

YEARS = [1617, 1718, 1819, 1920, 2021, 2122]

STATIONS = [
    "Adelboden", "Gadmen", "Grindelwald_Bort", "Gsteig", "Gantrisch",
    "Leysin", "Muerren", "Saanenmoeser", "Wengen", "Sorenberg", "Stoos",
    "Braunwald", "Malbun", "St_Margrethenberg", "Binn", "Bourg_St_Pierre",
    "Fionnay", "Grimentz", "Lauchernalp", "Montana", "Muenster",
    "Saas_Fee", "Simplon_Dorf", "Ulrichen", "Wiler", "Bivio",
    "Davos_Flueelastr", "Juf", "Obersaxen", "Pusserein", "St_Antoenien",
    "Sedrun", "Spluegen", "Vals", "Weisfluh_Joch", "Bosco_Gurin",
    "San_Bernadino", "Maloja", "Sankt_Moritz", "Samnaun", "Zuoz"
]

# =============================================================================
# Helper Functions
# =============================================================================

def load_input(csv_path: Path) -> pd.DataFrame:
    """Load one input CSV and standardize columns."""
    df = pd.read_csv(csv_path)

    if "date" not in df.columns:
        raise ValueError(f"{csv_path} missing required column 'date'")

    df["date"] = pd.to_datetime(df["date"], errors="coerce")

    if df["date"].isna().any():
        raise ValueError(f"{csv_path} contains invalid date values")

    if "HS" in df.columns:
        df = df.rename(columns={"HS": "HS_obs"})
    if "HNW" in df.columns:
        df = df.rename(columns={"HNW": "HNW_obs"})

    df = df.sort_values("date").reset_index(drop=True)
    df["date"] = df["date"].dt.strftime("%Y-%m-%d")

    return df


def run_r_model(csv_path: Path) -> pd.DataFrame:
    """
    Runs ΔSNOW R model.
    R cannot read stdin → write temporary CSV.
    """
    df = pd.read_csv(csv_path)

    if "date" not in df.columns:
        raise ValueError(f"{csv_path} missing required column 'date'")

    if "HS" in df.columns:
        hs_col = "HS"
    elif "HS_obs" in df.columns:
        hs_col = "HS_obs"
    else:
        raise ValueError(f"{csv_path} missing HS or HS_obs")

    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    if df["date"].isna().any():
        raise ValueError(f"{csv_path} contains invalid date values")

    df = df.sort_values("date").reset_index(drop=True)

    df_r = pd.DataFrame({
        "date": df["date"].dt.strftime("%Y-%m-%d"),
        "hs": df[hs_col]
    })

    if df_r.empty:
        raise ValueError(f"{csv_path} is empty")

    # Force first value to zero as required by model a bit shit but keep it (until better option found)
    df_r.at[df_r.index[0], "hs"] = 0

    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp:
        tmp_path = Path(tmp.name)
        df_r.to_csv(tmp_path, index=False)

    cmd = [R_BIN, str(R_RUNNER), "--in", str(tmp_path)]
    for k, v in PARAMS.items():
        cmd.extend([f"--{k}", str(v)])

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
    finally:
        if tmp_path.exists():
            os.remove(tmp_path)

    if proc.returncode != 0:
        print("\n[R STDOUT]")
        print(proc.stdout)
        print("\n[R STDERR]")
        print(proc.stderr)
        raise RuntimeError(f"R model failed for {csv_path}")

    if not proc.stdout.strip():
        raise RuntimeError(f"R model returned empty output for {csv_path}")

    df_out = pd.read_csv(StringIO(proc.stdout))
    df_out["date"] = pd.to_datetime(df_out["date"], errors="coerce")

    if df_out["date"].isna().any():
        raise ValueError(f"R output contains invalid dates for {csv_path}")

    df_out = df_out.sort_values("date").reset_index(drop=True)
    df_out["date"] = df_out["date"].dt.strftime("%Y-%m-%d")

    rename_map = {}
    if "hs" in df_out.columns:
        rename_map["hs"] = "HS_mod"
    if "swe_mod" in df_out.columns:
        rename_map["swe_mod"] = "SWE_mod"
    elif "swe" in df_out.columns:
        rename_map["swe"] = "SWE_mod"

    df_out = df_out.rename(columns=rename_map)

    required = ["date", "SWE_mod"]
    missing = [c for c in required if c not in df_out.columns]
    if missing:
        raise ValueError(f"R output missing required columns {missing} for {csv_path}")

    return df_out


# =============================================================================
# Data Processing
# =============================================================================

def process_all_stations() -> pd.DataFrame | None:
    """Run model for all station-season files and compile merged dataset."""
    all_rows = []

    for station in STATIONS:
        print(f"\n=== {station} ===")
        station_rows = []

        for yr in YEARS:
            in_csv = DATA_DIR / str(yr) / f"{station}.csv"
            print(in_csv)

            if not in_csv.exists():
                print(f"[SKIP] Missing {in_csv}")
                continue

            print(f"[PROCESS] {station} {yr}")

            try:
                df_in = load_input(in_csv)
                df_r = run_r_model(in_csv)

                merged = df_in.merge(df_r, on="date", how="left")
                merged["date"] = pd.to_datetime(merged["date"])
                merged = merged.sort_values("date").reset_index(drop=True)

                # HNW from day-to-day SWE difference within this file only
                merged["HNW_mod"] = merged["SWE_mod"].diff().clip(lower=0)

                merged["station"] = station
                merged["season"] = yr

                merged["date"] = merged["date"].dt.strftime("%Y-%m-%d")

                out_path = OUT_DIR / f"{yr}_{station}_with_model.csv"
                merged.to_csv(out_path, index=False)
                print(f"  → saved: {out_path}")

                station_rows.append(merged)
                all_rows.append(merged)

            except Exception as e:
                print(f"[ERROR] {station} {yr}: {e}")
                continue

        # Optional per-station combined file
        if station_rows:
            st_all = pd.concat(station_rows, ignore_index=True)
            st_out = OUT_DIR / f"{station}_all_seasons.csv"
            st_all.to_csv(st_out, index=False)
            print(f"[WRITE] all seasons for {station}: {st_out}")

    # Write full merged dataset AFTER all stations are processed
    if all_rows:
        all_df = pd.concat(all_rows, ignore_index=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        all_out = OUT_DIR / f"{timestamp}_ALL_STATIONS_ALL_YEARS.csv"
        all_df.to_csv(all_out, index=False)
        print(f"[WRITE] full merged dataset: {all_out}")
        return all_df

    return None


# =============================================================================
# Validation
# =============================================================================

def validate_swe(df: pd.DataFrame) -> None:
    print("\n=== SWE VALIDATION ===\n")

    obs_col = "SWE"
    if obs_col not in df.columns:
        print(f"Column '{obs_col}' not found — skipping SWE validation.")
        return

    if "SWE_mod" not in df.columns:
        print("Column 'SWE_mod' not found — skipping SWE validation.")
        return

    df_valid = df.dropna(subset=[obs_col, "SWE_mod"]).copy()

    if len(df_valid) == 0:
        print("No valid SWE observations found — skipping SWE validation.")
        return

    residuals = df_valid["SWE_mod"] - df_valid[obs_col]

    rmse = np.sqrt(np.mean(residuals**2))
    bias = np.mean(residuals)

    # avoid division by zero in relative bias
    nonzero = df_valid[obs_col] != 0
    rel_bias = np.mean(residuals[nonzero] / df_valid.loc[nonzero, obs_col]) * 100 if nonzero.any() else np.nan

    corr = np.corrcoef(df_valid["SWE_mod"], df_valid[obs_col])[0, 1]
    r2 = corr**2

    print(f"N = {len(df_valid)}")
    print(f"RMSE       = {rmse:.2f} mm")
    print(f"Bias       = {bias:.2f} mm")
    print(f"Rel. Bias  = {rel_bias:.1f}%")
    print(f"R²         = {r2:.3f}")

    plt.figure(figsize=(7, 7))
    x = df_valid["SWE_mod"]
    y = df_valid[obs_col]

    plt.scatter(x, y, color="black", alpha=0.5, s=22)
    lim_max = max(x.max(), y.max()) * 1.05
    plt.plot([0, lim_max], [0, lim_max], "--", color="red")

    plt.xlabel("Modeled SWE [mm]")
    plt.ylabel("Observed SWE [mm]")
    plt.title("Modeled vs Observed SWE")
    plt.grid(True, linestyle=":", alpha=0.5)

    textstr = (
        f"N = {len(df_valid)}\n"
        f"RMSE: {rmse:.1f}\n"
        f"Bias: {bias:.1f}\n"
        f"Rel. Bias: {rel_bias:.1f}%\n"
        f"R²: {r2:.3f}"
    )

    plt.text(
        0.05, 0.95, textstr,
        transform=plt.gca().transAxes,
        fontsize=11,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.9)
    )

    plt.tight_layout()
    plt.show()


def validate_hnw(df: pd.DataFrame) -> None:
    print("\n=== HNW VALIDATION ===\n")

    required = ["HNW_obs", "HNW_mod"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"Missing columns {missing} — skipping HNW validation.")
        return

    df_valid = df.dropna(subset=["HNW_obs", "HNW_mod"]).copy()

    # Keep only positive observed HNW
    df_valid = df_valid[df_valid["HNW_obs"] > 0]

    if len(df_valid) == 0:
        print("No valid HNW observations found — skipping HNW validation.")
        return

    y = df_valid["HNW_obs"].values
    x = df_valid["HNW_mod"].values

    residuals = x - y

    rmse = np.sqrt(np.mean(residuals**2))
    bias = np.mean(residuals)
    pbias = 100 * np.sum(residuals) / np.sum(y)

    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept

    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot != 0 else np.nan

    print(f"N = {len(df_valid)}")
    print(f"RMSE  = {rmse:.3f}")
    print(f"Bias  = {bias:.3f}")
    print(f"PBIAS = {pbias:.2f}%")
    print(f"R²    = {r2:.3f}")

    plt.figure(figsize=(8, 7))

    plt.hist2d(
        x, y,
        bins=50,
        range=[[0, 100], [0, 100]],
        norm=LogNorm(vmin=1, vmax=1000),
        cmap="jet"
    )

    plt.colorbar(label="Number of observations")

    lim = [0, 100]
    plt.plot(lim, lim, "--", color="gray")
    plt.plot(lim, slope * np.array(lim) + intercept, "--", color="red")

    plt.xlabel("Modeled HNW [mm]")
    plt.ylabel("Observed HNW [mm]")
    plt.title("ΔSNOW HNW Validation")

    textstr = (
        f"$R^2$: {r2:.2f}\n"
        f"Rel. bias: {pbias:.2f}%\n"
        f"RMSE: {rmse:.1f}"
    )

    plt.text(
        0.03, 0.97, textstr,
        transform=plt.gca().transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
    )

    plt.xlim(lim)
    plt.ylim(lim)
    plt.tight_layout()
    plt.show()


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    all_df = process_all_stations()

    if all_df is None:
        print("No data processed.")
        return

    validate_swe(all_df)
    validate_hnw(all_df)


if __name__ == "__main__":
    main()