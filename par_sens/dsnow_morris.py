#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Δ-SNOW Morris Sensitivity Analysis
----------------------------------

This script:
  - Loads seasonal HS + SWE observations
  - Runs Δ-SNOW (R implementation) for many sampled parameter vectors
  - Computes RMSE between SWE_mod and SWE_obs
  - Runs Morris global sensitivity analysis (SALib)
  - Plots μ* (mu_star) as a bar plot with uncertainty (σ)
"""

# ============================================================
# USER SETTINGS — EDIT THESE ONLY
# ============================================================


# These paramters maily influence the sampling density and therfore the runnig time of the programm

MORRIS_N = 20          # number of Morris trajectories
MORRIS_LEVELS = 4      # number of grid levels per parameter





R_BIN = "/usr/local/bin/Rscript"
R_RUNNER = (
    "/Users/jakobwerkgarner/code/mt_dsnow/calibration/helpers/"
    "minimal_delta_snow_runner.R"
)
DATA_DIR = (
    "/Users/jakobwerkgarner/code/mt_dsnow/HNW_validation/validation_input"
)

SEASONS = [
    1617,
    1718,
    1819,
    1920,
    2021,
    2122,
]

STATIONS = [
    "Adelboden", "Gadmen", "Grindelwald_Bort", "Gsteig", "Gantrisch",
    "Leysin", "Muerren", "Saanenmoeser", "Wengen", "Stoos",
    "Braunwald", "Malbun", "St_Margrethenberg", "Binn", "Bourg_St_Pierre",
    "Fionnay", "Grimentz", "Lauchernalp", "Montana", "Muenster",
    "Saas_Fee", "Simplon_Dorf", "Ulrichen", "Wiler", "Bivio",
    "Davos_Flueelastr", "Juf", "Obersaxen", "Pusserein", "St_Antoenien",
    "Sedrun", "Spluegen", "Vals", "Weisfluh_Joch", "Bosco_Gurin",
    "San_Bernadino", "Maloja", "Sankt_Moritz", "Samnaun", "Zuoz",
]

PARAM_NAMES = ["rho.max", "rho.null", "c.ov", "k.ov", "k", "tau", "eta.null"]
BOUNDS = [
    [300, 600],      # rho.max
    [70, 130],       # rho.null
    [1e-9, 1e-3],    # c.ov
    [0.01, 1.0],     # k.ov
    [0.01, 0.2],     # k
    [0.01, 0.2],     # tau
    [1e6, 2e7],      # eta.null
]


OUTPUT_DIR = "./dsnow_sensitivity_output"
PLOT_FILENAME = "morris_mu_star_barplot.png"
SAVE_PLOT = True

# ============================================================
# IMPORTS
# ============================================================

import os
import subprocess
import tempfile
from io import StringIO
from pathlib import Path
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze


# ============================================================
# LOAD OBSERVATIONS
# ============================================================

def load_seasonal_data():
    """
    Load and combine seasonal HS + SWE observations for all stations.

    For each station and season, this function:
      - Reads the corresponding CSV from DATA_DIR/year/station.csv
      - Checks that a SWE column exists
      - Renames HS and SWE to HS_obs and SWE_obs
      - Adds station and season metadata

    Returns
    -------
    pandas.DataFrame
        Long-format DataFrame with columns:
        ['station', 'season', 'date', 'HS_obs', 'SWE_obs']
    """
    rows = []

    for station in STATIONS:
        for year in SEASONS:
            csv_path = Path(DATA_DIR) / str(year) / f"{station}.csv"

            if not csv_path.exists():
                print(f"[WARN] Missing {csv_path}")
                continue

            df = pd.read_csv(csv_path)

            # SWE is required for RMSE calculation; bail out early if missing
            if "SWE" not in df.columns:
                raise ValueError(f"Missing SWE column in {csv_path}")

            df["station"] = station
            df["season"] = year

            df_sel = df[["station", "season", "date", "HS", "SWE"]].copy()
            df_sel.rename(
                columns={"HS": "HS_obs", "SWE": "SWE_obs"},
                inplace=True,
            )

            rows.append(df_sel)

    # Concatenate all station-season tables into a single DataFrame
    df_all = pd.concat(rows, ignore_index=True)
    df_all["date"] = pd.to_datetime(df_all["date"])

    return df_all


# ============================================================
# R RUNNER CALL
# ============================================================

def run_r_model(df, params):
    """
    Run the Δ-SNOW R model for a single station-season subset.

    The function:
      - Prepares a temporary CSV with 'date' and 'hs' (HS time series)
      - Enforces hs[0] = 0 to avoid unrealistic initial snowpack
      - Calls the R runner with parameter flags (e.g. --rho.max 400)
      - Parses the R stdout back into a DataFrame with model outputs

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame with at least ['date', 'HS_obs'] columns.
    params : dict
        Mapping of parameter name -> value, consistent with PARAM_NAMES.

    Returns
    -------
    pandas.DataFrame or None
        DataFrame with model outputs columns:
        ['date', 'HS_mod', 'SWE_mod']
        or None if the R model failed or returned invalid output.
    """
    # Prepare input DataFrame for R runner: date + HS (as "hs")
    df_r = df[["date", "HS_obs"]].copy()
    df_r.rename(columns={"HS_obs": "hs"}, inplace=True)

    # Enforce first HS = 0 and fill missing values with 0
    df_r["hs"] = df_r["hs"].fillna(0)
    df_r.at[df_r.index[0], "hs"] = 0

    # Write temporary CSV for R to consume
    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp:
        tmp_path = Path(tmp.name)
        df_r.to_csv(tmp_path, index=False)

    # Build R command: Rscript minimal_delta_snow_runner.R --in tmp.csv --param value ...
    cmd = [R_BIN, R_RUNNER, "--in", str(tmp_path)]
    for name, value in params.items():
        cmd += [f"--{name}", str(value)]

    # Execute R runner and capture stdout/stderr
    proc = subprocess.run(cmd, capture_output=True, text=True)

    # Clean up temporary file
    os.remove(tmp_path)

    # Non-zero return code → model failed
    if proc.returncode != 0:
        print("[WARN] R runner returned non-zero exit code.")
        print("stderr:", proc.stderr)
        return None

    # Try to parse R output as CSV from stdout
    try:
        df_out = pd.read_csv(StringIO(proc.stdout))
        df_out["date"] = pd.to_datetime(df_out["date"])
    except Exception as exc:
        print("[WARN] Failed to parse R output:", exc)
        return None

    # Standardize column names expected downstream
    df_out.rename(
        columns={"hs": "HS_mod", "swe_mod": "SWE_mod"},
        inplace=True,
    )

    return df_out[["date", "HS_mod", "SWE_mod"]]


# ============================================================
# RUN MODEL FOR ALL STATIONS × SEASONS
# ============================================================

def run_all_delta_snow(df_seasonal, params):
    """
    Run Δ-SNOW for all available station-season combinations.

    For each station and each season:
      - Extract subset of observations
      - Call the R model with current parameter set
      - Merge model output with observations on 'date'

    Parameters
    ----------
    df_seasonal : pandas.DataFrame
        Long-format observations from load_seasonal_data().
    params : dict
        Parameter dictionary used to configure Δ-SNOW.

    Returns
    -------
    pandas.DataFrame
        Combined DataFrame with observed and modelled time series:
        ['station', 'season', 'date', 'HS_obs', 'SWE_obs', 'HS_mod', 'SWE_mod']
        Empty DataFrame if all R calls fail.
    """
    results = []

    # Loop over stations
    for station in df_seasonal["station"].unique():
        df_station = df_seasonal[df_seasonal["station"] == station]

        # Loop over seasons
        for year in SEASONS:
            sub = df_station[df_station["season"] == year]

            # Skip if we have no data for this station-season
            if len(sub) == 0:
                continue

            df_r = run_r_model(sub, params)
            if df_r is None:
                # Model failed for this station-season; skip
                continue

            # Merge model output back onto observations by date
            merged = sub.merge(df_r, on="date", how="left")
            results.append(merged)

    if len(results) == 0:
        # No successful model runs → return empty DataFrame
        return pd.DataFrame()

    # Combine all station-season results
    return pd.concat(results, ignore_index=True)


# ============================================================
# RMSE ON SWE
# ============================================================

def compute_rmse(df):
    """
    Compute RMSE between observed and modelled SWE.

    The function returns a large penalty (1e9) if:
      - The DataFrame is empty, or
      - There are no overlapping non-NaN SWE_obs / SWE_mod values,
      - The resulting RMSE is NaN.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with 'SWE_obs' and 'SWE_mod' columns.

    Returns
    -------
    float
        Root Mean Squared Error (RMSE) or 1e9 as a failure penalty.
    """
    if df.empty:
        return 1e9

    mask = df["SWE_obs"].notna() & df["SWE_mod"].notna()
    if not np.any(mask):
        return 1e9

    err = df.loc[mask, "SWE_mod"] - df.loc[mask, "SWE_obs"]
    rmse = np.sqrt(np.mean(err**2))

    if np.isnan(rmse):
        return 1e9

    return float(rmse)


# ============================================================
# MORRIS WORKER
# ============================================================

def evaluate_theta(args):
    """
    Worker function to evaluate a single parameter vector (theta).

      - Converts the theta vector into a parameter dict
      - Runs Δ-SNOW for all stations × seasons
      - Computes RMSE on SWE

    Parameters
    ----------
    args : tuple
        (theta, df_seasonal), where:
          - theta : array-like of parameter values
          - df_seasonal : pandas.DataFrame of observations

    Returns
    -------
    float
        RMSE value for this theta (or 1e9 if anything fails).
    """
    theta, df_seasonal = args

    # Map parameter names to numeric values (ensure float for safety)
    params = {name: float(value) for name, value in zip(PARAM_NAMES, theta)}

    try:
        df_out = run_all_delta_snow(df_seasonal, params)
        return compute_rmse(df_out)
    except Exception as exc:
        print("Error for theta", theta, ":", exc)
        return 1e9


# ============================================================
# RUN MORRIS SAMPLING + PARALLEL MODEL EVALUATION
# ============================================================

def run_morris(df_seasonal):
    """
    Execute the Morris global sensitivity analysis for the Δ-SNOW model.

    Steps:
      1) Define SALib problem structure (parameter names, bounds)
      2) Generate Morris trajectories (sampled parameter vectors)
      3) Evaluate Δ-SNOW for each vector (parallel)
      4) Compute RMSE for each vector
      5) Analyze results with SALib.morris.analyze

    Parameters
    ----------
    df_seasonal : pandas.DataFrame
        Long-format seasonal HS + SWE observations.

    Returns
    -------
    dict
        Morris sensitivity indices with keys like:
        ['mu', 'mu_star', 'sigma', 'mu_star_conf', ...]
    """
    # 1) Define SALib problem dictionary
    problem = {
        "num_vars": len(PARAM_NAMES),
        "names": PARAM_NAMES,
        "bounds": BOUNDS,
    }

    # 2) Generate Morris samples
    print("Generating Morris samples…")
    X = morris_sample.sample(
        problem,
        N=MORRIS_N,
        num_levels=MORRIS_LEVELS,
    )
    n = len(X)

    # 3) Determine number of worker processes

    workers = max(5, cpu_count() - 1)

    print(f"Evaluating {n} parameter sets using {workers} worker(s)…")

    # 4) Build list of tasks: (theta, df_seasonal) for each sample row
    tasks = [(theta, df_seasonal) for theta in X]

    Y = []

    # 5) Evaluate all parameter sets (in parallel where possible)
    with Pool(workers) as pool:
        # pool.imap preserves the input order of X in the results
        for i, result in enumerate(pool.imap(evaluate_theta, tasks), start=1):
            print(f"[{i}/{n}] done")
            Y.append(result)

    # Convert to NumPy array for SALib
    Y = np.array(Y)

    # 6) Analyze Morris results
    Si = morris_analyze.analyze(
        problem,
        X,
        Y,
        print_to_console=False,
    )

    return Si


# ============================================================
# PLOT RESULTS
# ============================================================

def plot_morris(Si):
    """
    Plot Morris μ* importance with σ as error bars.

    Parameters
    ----------
    Si : dict
        Sensitivity indices as returned by run_morris().
        Must contain 'mu_star' and 'sigma'.
    """
    mu_star = Si["mu_star"]
    sigma = Si["sigma"]
    x = np.arange(len(PARAM_NAMES))

    plt.figure(figsize=(9, 5))
    plt.bar(x, mu_star, yerr=sigma, capsize=6)
    plt.xticks(x, PARAM_NAMES, rotation=45)
    plt.ylabel("Morris μ* (importance)")
    plt.title("Δ-SNOW Sensitivity (SWE RMSE)")
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    # Ensure output directory exists
    Path(OUTPUT_DIR).mkdir(exist_ok=True)
    outpath = Path(OUTPUT_DIR) / PLOT_FILENAME

    if SAVE_PLOT:
        plt.savefig(outpath, dpi=300)
        print(f"Saved plot → {outpath}")

    plt.show()


# ============================================================
# PLOT Morrison Scatter RESULTS
# ============================================================

def plot_morris_scatter(Si):
    """
    Plot the classical Morris sensitivity scatter plot:
    x-axis = μ* (mean absolute effect)
    y-axis = σ  (standard deviation of effects)

    Parameters
    ----------
    Si : dict
        Sensitivity index dictionary from SALib.morris.analyze()
        Must contain keys 'mu_star' and 'sigma'.
    """

    mu_star = Si["mu_star"]
    sigma = Si["sigma"]

    plt.figure(figsize=(8, 6))
    plt.scatter(mu_star, sigma, s=100)

    # Label each point with parameter name
    for name, x, y in zip(PARAM_NAMES, mu_star, sigma):
        plt.text(x, y, name, fontsize=10, ha="left", va="bottom")

    plt.xlabel("μ* (overall importance)")
    plt.ylabel("σ (interaction / nonlinearity)")
    plt.title("Morris Sensitivity Scatter Plot (Δ-SNOW)")
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # Save figure
    Path(OUTPUT_DIR).mkdir(exist_ok=True)
    outpath = Path(OUTPUT_DIR) / "morris_mu_star_sigma_scatter.png"
    plt.savefig(outpath, dpi=300)
    print(f"Saved scatter plot → {outpath}")

    plt.show()


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    
    print("Loading HS + SWE seasonal data…")
    df_all = load_seasonal_data()
    print(f"Loaded {len(df_all)} rows. SWE preview:")
    print(df_all["SWE_obs"].describe())

    print("\nRunning Morris sensitivity…")
    si = run_morris(df_all)

    print("\nMorris μ* values:")
    for name, mu_val in zip(PARAM_NAMES, si["mu_star"]):
        print(f"  {name:10s} : {mu_val:.6f}")

    plot_morris(si)

    plot_morris_scatter(si)

    print("\nDONE.")