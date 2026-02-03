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
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# =============================================================================
# CONFIGURATION
# =============================================================================



PARAMS = {
    "rho.max": 451.6977806582531,
    "rho.null": 90.0,
    "c.ov": 2.746976858091589e-0,
    "k.ov": 0.38,
    "k": 0.020385468323087456,
    "tau":  0.000012,
    "eta.null": 8.5e6
}





BASE_DIR = Path("/Users/jakobwerkgarner/code/mt_dsnow")
os.chdir(BASE_DIR)

DATA_DIR = Path("HNW_validation/validation_input")

OUT_DIR = Path('HNW_validation/validation_output')
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

def load_input(csv_path: Path):
    df = pd.read_csv(csv_path)
    df["date"] = pd.to_datetime(df["date"]).dt.strftime("%Y-%m-%d")

    if "HS" in df.columns:
        df = df.rename(columns={"HS": "HS_obs"})
    if "HNW" in df.columns:
        df = df.rename(columns={"HNW": "HNW_obs"})

    return df


def run_r_model(csv_path: Path):
    """
    Runs ΔSNOW R model. R cannot read stdin → we write a temporary file.
    """
    df = pd.read_csv(csv_path)

    if "HS" in df.columns:
        hs_col = "HS"
    elif "HS_obs" in df.columns:
        hs_col = "HS_obs"
    else:
        raise ValueError(f"{csv_path} missing HS or HS_obs")

    df_r = pd.DataFrame({
        "date": df["date"],
        "hs": df[hs_col]
    })

    df_r.at[df_r.index[0], "hs"] = 0

    # temporary CSV
    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp:
        tmp_path = Path(tmp.name)
        df_r.to_csv(tmp_path, index=False)

    # R command
    cmd = [R_BIN, str(R_RUNNER), "--in", str(tmp_path)]
    for k, v in PARAMS.items():
        cmd.extend([f"--{k}", str(v)])

    # run process
    proc = subprocess.run(cmd, capture_output=True, text=True)

    os.remove(tmp_path)

    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr)
        raise RuntimeError("R model failed.")

    df_out = pd.read_csv(StringIO(proc.stdout))
    df_out["date"] = pd.to_datetime(df_out["date"]).dt.strftime("%Y-%m-%d")

    return df_out.rename(columns={"hs": "hs_mod", "swe_mod": "SWE_mod"})



# =============================================================================
# Data Processing
# =============================================================================

def process_all_stations():
    all_rows = []

    for station in STATIONS:
        print(f"\n=== {station} ===")
        station_rows = []

        for yr in YEARS:
            in_csv = DATA_DIR / str(yr) / f"{station}.csv"
            print(str(in_csv))

            if not in_csv.exists():
                print(f"[SKIP] Missing {in_csv}")
                continue

            print(f"[PROCESS] {station} {yr}")

            df_in = load_input(in_csv)
            df_r = run_r_model(in_csv)

            merged = df_in.merge(df_r, on="date", how="left")
            merged["HNW_mod"] = merged["SWE_mod"].diff()
            merged["station"] = station
            merged["season"] = yr

            out_path = OUT_DIR / f"{yr}_{station}_with_model.csv"
            merged.to_csv(out_path, index=False)
            print(f"  → saved: {out_path}")

            station_rows.append(merged)
            all_rows.append(merged)

        # # Write per-station all seasons
        # if station_rows:
        #     st_all = pd.concat(station_rows, ignore_index=True)
        #     st_out = OUT_DIR / f"{station}_all_seasons.csv"
        #     st_all.to_csv(st_out, index=False)
        #     print(f"[WRITE] all seasons for {station}")

        from datetime import datetime
        from pathlib import Path
        import pandas as pd

        # global file
        if all_rows:
            all_df = pd.concat(all_rows, ignore_index=True)

            # Create timestamp string
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

            # Build full output path safely
            all_out = OUT_DIR / f"{timestamp}_ALL_STATIONS_ALL_YEARS.csv"

            all_df.to_csv(all_out, index=False)
            print(f"[WRITE] full merged dataset: {all_out}")

            return all_df
    


    return None


# =============================================================================
# Validation
# =============================================================================

def validate_swe(df):
    print("\n=== SWE VALIDATION ===\n")

    obs_col = "SWE"
    df_valid = df[df[obs_col].notna()]

    residuals = df_valid["SWE_mod"] - df_valid[obs_col]

    rmse = np.sqrt(np.mean(residuals**2))
    bias = np.mean(residuals)
    rel_bias = np.mean(residuals / df_valid[obs_col]) * 100

    corr = np.corrcoef(df_valid["SWE_mod"], df_valid[obs_col])[0, 1]
    r2 = corr**2

    print(f"N = {len(df_valid)}")
    print(f"RMSE       = {rmse:.2f} mm")
    print(f"Bias       = {bias:.2f} mm")
    print(f"Rel. Bias  = {rel_bias:.1f}%")
    print(f"R²         = {r2:.3f}")

    # scatter plot
    plt.figure(figsize=(7,7))
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


def validate_hnw(df):
    print("\n=== HNW VALIDATION ===\n")

    df_valid = df.dropna(subset=["HNW_obs", "HNW_mod"])

    # Keep only positive observed HNW
    df_valid = df_valid[df_valid["HNW_obs"] > 0]

    # Safety check
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
    r2 = 1 - ss_res / ss_tot

    print(f"N = {len(df_valid)}")
    print(f"RMSE  = {rmse:.3f}")
    print(f"Bias  = {bias:.3f}")
    print(f"PBIAS = {pbias:.2f}%")
    print(f"R²    = {r2:.3f}")

    # ---------------------------------------------------------
    # 2D Histogram Plot
    # ---------------------------------------------------------
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

def main():
    all_df = process_all_stations()
    if all_df is None:
        print("No data processed.")
        return

    validate_swe(all_df)
    validate_hnw(all_df)


if __name__ == "__main__":
    main()