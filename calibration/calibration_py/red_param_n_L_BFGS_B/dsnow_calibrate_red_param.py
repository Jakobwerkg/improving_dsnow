#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Δ-SNOW Calibration Pipeline (Python, L-BFGS-B)
----------------------------------------------

This script:

  - Loads HS + SWE data from a NetCDF file
  - Splits data into calibration (fit) and validation (val) winters
  - Runs the Δ-SNOW R model in parallel for many parameter sets
  - Calibrates ONLY a subset of parameters (defined in config)
    while keeping other parameters FIXED to given values
  - Uses an internal scaling for optimization stability
  - Logs RMSE, Bias and parameter values per iteration (CSV)
  - Exports final results as:
        - NetCDF with SWE_obs / SWE_mod / HS
        - RMSE/Bias time series plot
  - Validates calibrated parameters on held-out winters

Configuration is fully controlled via a YAML file
(e.g. config_calib.yml).
"""

# ============================================================
# IMPORTS
# ============================================================

import os
import csv
import yaml
import subprocess
import tempfile
from io import StringIO
from pathlib import Path
from datetime import timedelta
from multiprocessing import Pool, cpu_count, set_start_method

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.optimize import minimize


# Enable fast multiprocessing on macOS (and generally safe on POSIX)
try:
    set_start_method("fork")
except RuntimeError:
    # Already set in this interpreter session → ignore
    pass


# ============================================================
# CONFIG LOADER
# ============================================================

def load_config(config_path: str = "config_calib.yml") -> dict:
    """
    Load YAML configuration file for the calibration.

    Parameters
    ----------
    config_path : str
        Relative path to the configuration YAML.

    Returns
    -------
    dict
        Parsed configuration dictionary.
    """
    try:
        base_path = Path(__file__).parent
    except NameError:
        # Fallback for interactive environments (e.g. Jupyter)
        base_path = Path(os.getcwd())

    config_path = base_path / config_path
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    return cfg


# ============================================================
# LOGGING UTILITIES
# ============================================================

def init_logging(log_file: Path, param_names: list[str]) -> None:
    """
    Initialize CSV log file for calibration iterations.

    Parameters
    ----------
    log_file : Path
        Path to the CSV file to create/overwrite.
    param_names : list of str
        Names of calibrated parameters (columns).
    """
    with open(log_file, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["iteration", "rmse", "bias"] + param_names
        writer.writerow(header)


def append_logging(
    log_file: Path,
    iteration: int,
    rmse: float,
    bias: float,
    params_unscaled: dict,
    param_names: list[str],
) -> None:
    """
    Append a single row (iteration) to the calibration log CSV.

    Parameters
    ----------
    log_file : Path
        Path to the CSV log file.
    iteration : int
        Current iteration number (1-based).
    rmse : float
        RMSE for this parameter set.
    bias : float
        Mean bias for this parameter set.
    params_unscaled : dict
        Dictionary of all parameters (possibly including fixed).
    param_names : list of str
        List of calibrated parameter names to be logged.
    """
    with open(log_file, "a", newline="") as f:
        writer = csv.writer(f)
        row = [iteration, rmse, bias] + [params_unscaled[p] for p in param_names]
        writer.writerow(row)


def export_netcdf(merged_df: pd.DataFrame, outfile: Path) -> None:
    """
    Export final merged calibration results to a NetCDF file.

    Parameters
    ----------
    merged_df : pandas.DataFrame
        DataFrame containing at least ['date', 'hs', 'swe_obs', 'swe_mod'].
    outfile : Path
        Path to the NetCDF file to write.
    """
    ds = xr.Dataset(
        {
            "hs": (("time",), merged_df["hs"].values),
            "swe_obs": (("time",), merged_df["swe_obs"].values),
            "swe_mod": (("time",), merged_df["swe_mod"].values),
        },
        coords={"time": merged_df["date"].values},
    )
    ds.to_netcdf(outfile)
    print(f"Saved NetCDF → {outfile}")


def plot_rmse_bias(log_file: Path, outfile: Path) -> None:
    """
    Plot RMSE and Bias vs. iteration from the calibration log.

    Parameters
    ----------
    log_file : Path
        Path to the CSV log file created during calibration.
    outfile : Path
        Path to the PNG file to write.
    """
    df = pd.read_csv(log_file)
    plt.figure(figsize=(8, 5))
    plt.plot(df["iteration"], df["rmse"], label="RMSE")
    plt.plot(df["iteration"], df["bias"], label="Bias")
    plt.xlabel("Iteration")
    plt.ylabel("Value")
    plt.title("RMSE & Bias during Optimization")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    print(f"Saved plot → {outfile}")


# ============================================================
# DATA PREPARATION
# ============================================================

def build_tibble(d_obs_fit: dict, start_month: int = 8) -> pd.DataFrame:
    """
    Build a long-format 'tibble'-style DataFrame from station-wise data.

    Each station dataframe is assumed to have a DateTimeIndex and columns:
    - 'hs'       : snow depth
    - 'swe_obs'  : observed snow water equivalent

    Parameters
    ----------
    d_obs_fit : dict
        Mapping station_name -> pandas.DataFrame with index as datetime.
    start_month : int
        Month where a hydrological season starts (e.g. 8 for August).

    Returns
    -------
    pandas.DataFrame
        Columns: ['date', 'name', 'hs', 'swe_obs', 'block']
        where 'block' is an integer representing the winter season.
    """
    rows = []
    for name, df in d_obs_fit.items():
        df = df.copy()
        # Block (winter) index: year of the season start
        blocks = [
            (d.year - 1 if d.month < start_month else d.year)
            for d in df.index
        ]

        tib = pd.DataFrame(
            {
                "date": df.index,
                "name": name,
                "hs": df["hs"].values,
                "swe_obs": df["swe_obs"].values,
                "block": blocks,
            }
        )
        rows.append(tib)

    return pd.concat(rows, ignore_index=True)


def split_into_seasons(
    station_dfs: dict,
    season_start: str = "-08-01",
    season_end: str = "-07-31",
) -> tuple[dict, dict]:
    """
    Split station time series into winters and separate fit vs validation.

    Winters are split by hydrological years:
      - For year y, the season runs from y-08-01 to (y+1)-07-31.
      - Short winters (<150 days) are discarded.
      - Winters with starting year even go to the fit set,
        winters with starting year odd go to the validation set.

    Parameters
    ----------
    station_dfs : dict
        Mapping station_name -> DataFrame with columns ['date', 'hs', 'swe_obs'].
    season_start : str
        Seasonal start date suffix, e.g. "-08-01".
    season_end : str
        Seasonal end date suffix, e.g. "-07-31".

    Returns
    -------
    (dict, dict)
        d_obs_fit, d_obs_val
        Each mapping station_name -> concatenated DataFrame of winters.
    """
    d_obs_fit, d_obs_val = {}, {}

    for station, df in station_dfs.items():
        df = df.copy().set_index("date")
        years = df.index.year.unique()

        winters_fit, winters_val = [], []

        for y in years:
            start = pd.to_datetime(f"{y}{season_start}")
            end = pd.to_datetime(f"{y + 1}{season_end}") + timedelta(days=1)
            w = df[(df.index >= start) & (df.index < end)]

            # Skip too short winters
            if len(w) < 150:
                continue

            # Even year -> fit, Odd year -> validation
            if w.index[0].year % 2 == 0:
                winters_fit.append(w)
            else:
                winters_val.append(w)

        if winters_fit:
            d_obs_fit[station] = pd.concat(winters_fit)
        if winters_val:
            d_obs_val[station] = pd.concat(winters_val)

    return d_obs_fit, d_obs_val


# ============================================================
# Δ-SNOW PARALLEL RUNNER
# ============================================================

def run_single_block(args):
    """
    Helper function for multiprocessing: run Δ-SNOW for one station+block.

    Parameters
    ----------
    args : tuple
        (station, block, df_sub, r_script_path, params_unscaled, temp_dir)

    Returns
    -------
    pandas.DataFrame or None
        DataFrame with columns at least ['date', 'swe_mod', 'station', 'block'],
        or None if the R runner fails.
    """
    station, block, df_sub, r_script_path, params, temp_dir = args

    # Ensure a continuous daily time axis with hs filled as 0 where missing
    df_sub = df_sub.sort_values("date")
    full_range = pd.date_range(
        df_sub["date"].min(),
        df_sub["date"].max(),
        freq="D",
    )
    df_sub = (
        df_sub.set_index("date")
        .reindex(full_range)
        .rename_axis("date")
        .reset_index()
    )
    df_sub["hs"] = df_sub["hs"].fillna(0)

    # Write input CSV for R runner
    input_csv = Path(temp_dir) / f"hs_input_{station}_{block}.csv"
    df_sub[["date", "hs"]].to_csv(input_csv, index=False)

    # Build R command
    cmd = ["Rscript", r_script_path, "--in", str(input_csv)]
    for k, v in params.items():
        cmd += [f"--{k}", str(v)]

    p = subprocess.run(cmd, capture_output=True, text=True)

    if p.returncode != 0:
        print(f"ERROR {station} block {block}: {p.stderr}")
        return None

    try:
        out = pd.read_csv(StringIO(p.stdout), parse_dates=["date"])
        out["station"] = station
        out["block"] = block
        return out
    except Exception:
        return None


def evaluate_model_parallel(
    df_all: pd.DataFrame,
    r_script_path: str,
    params_unscaled: dict,
    workers: int,
    temp_dir: str = "/tmp",
) -> tuple[pd.DataFrame, dict]:
    """
    Evaluate the Δ-SNOW model on all station-block combinations in parallel.

    Parameters
    ----------
    df_all : pandas.DataFrame
        Long-format tibble with columns:
        ['date', 'name', 'hs', 'swe_obs', 'block'].
    r_script_path : str
        Path to the Δ-SNOW R runner script.
    params_unscaled : dict
        Dictionary of ALL parameters to pass to R (calibrated + fixed).
    workers : int
        Number of worker processes to use.
    temp_dir : str
        Directory for temporary CSV files.

    Returns
    -------
    (pandas.DataFrame, dict)
        merged : merged calibration dataframe with swe_mod attached.
        metrics : {'rmse': float, 'bias': float}
    """
    stations = df_all["name"].unique()
    blocks = df_all["block"].unique()

    jobs = []
    for st in stations:
        for blk in blocks:
            df_sub = df_all[(df_all["name"] == st) & (df_all["block"] == blk)]
            if len(df_sub) < 100:
                continue
            jobs.append(
                (st, blk, df_sub.copy(), r_script_path, params_unscaled, temp_dir)
            )

    # Parallel execution of all station-block runs
    with Pool(workers) as pool:
        results = list(
            tqdm(
                pool.imap_unordered(run_single_block, jobs),
                total=len(jobs),
                desc="Δ-SNOW blocks",
            )
        )

    # Keep only successful runs
    results = [r for r in results if r is not None]
    if not results:
        raise RuntimeError("No Δ-SNOW output from any block.")

    out = pd.concat(results, ignore_index=True)

    # Keep only necessary model outputs to avoid naming conflicts
    out = out[["date", "station", "block", "swe_mod"]]

    # Merge observations with model results by (date, name, block)
    merged = df_all.merge(
        out.rename(columns={"station": "name"}),
        on=["date", "name", "block"],
        how="inner",
    )

    # Filter valid SWE pairs and compute metrics
    valid = merged.dropna(subset=["swe_obs", "swe_mod"])
    rmse = np.sqrt(np.mean((valid["swe_mod"] - valid["swe_obs"]) ** 2))
    bias = np.mean(valid["swe_mod"] - valid["swe_obs"])

    return merged, {"rmse": rmse, "bias": bias}


# ============================================================
# PENALTY FUNCTION
# ============================================================

def exp_penalty(value: float, lower: float, upper: float) -> float:
    """
    Exponential penalty when a parameter is outside the central 50% of its range.

    The idea:
      - Middle 50% of [lower, upper] → safe zone, no penalty
      - Outside this zone → penalty grows exponentially towards the bounds

    Parameters
    ----------
    value : float
        Parameter value.
    lower : float
        Lower bound of the parameter.
    upper : float
        Upper bound of the parameter.

    Returns
    -------
    float
        Penalty contribution ≥ 0.0.
    """
    mid = 0.5 * (lower + upper)
    half_range = 0.25 * (upper - lower)  # central 50% safe zone

    # Inside central safe region → no penalty
    if (value >= lower + half_range) and (value <= upper - half_range):
        return 0.0

    # Normalized distance into penalty zone (0 → 1)
    if value < mid:
        distance = (lower + half_range - value) / half_range
    else:
        distance = (value - (upper - half_range)) / half_range

    # Exponential penalty that explodes near bounds
    return float(np.exp(5 * distance) - 1.0)


# ============================================================
# OPTIMIZATION — INTERNAL SCALING
# ============================================================

def optimize_params(
    df_fit: pd.DataFrame,
    r_script_path: str,
    param_names: list[str],
    initial_unscaled: list[float],
    bounds_unscaled: list[list[float]],
    scale: list[float],
    workers: int,
    log_file: Path,
    penalty_weight: float,
    maxiter: int,
    fixed_params: dict,
):
    """
    Optimize (calibrate) a subset of Δ-SNOW parameters using L-BFGS-B.

    Only parameters in `param_names` are calibrated. Additional
    parameters from `fixed_params` are kept constant but always passed
    to the R runner.

    Internal scaling is used for numerical stability:
      x_scaled = x_unscaled / scale

    Parameters
    ----------
    df_fit : pandas.DataFrame
        Tibble-style calibration data from build_tibble().
    r_script_path : str
        Path to the Δ-SNOW R runner.
    param_names : list of str
        Names of parameters to calibrate (subset of all params).
    initial_unscaled : list of float
        Initial guesses in original units.
    bounds_unscaled : list of [float, float]
        Lower/upper bounds for each calibrated parameter.
    scale : list of float
        Scaling factors (same length as param_names).
    workers : int
        Number of parallel workers.
    log_file : Path
        Path to the CSV log file.
    penalty_weight : float
        Weight factor multiplying the total penalty term.
    maxiter : int
        Maximum number of L-BFGS-B iterations.
    fixed_params : dict
        Parameter name → value for parameters that should NOT be calibrated.

    Returns
    -------
    (OptimizeResult, dict)
        res : scipy.optimize.OptimizeResult
        final_unscaled : dict of ALL parameters (calibrated + fixed).
    """
    initial_scaled = [v / s for v, s in zip(initial_unscaled, scale)]
    bounds_scaled = [(lo / s, hi / s) for (lo, hi), s in zip(bounds_unscaled, scale)]

    iteration = {"i": 0}

    # Parameters that should NOT receive penalties (optional)
    # You can adjust this set if desired.
    penalty_exempt = {"tau", "c.ov", "k.ov"}

    def objective(x_scaled: np.ndarray) -> float:
        """
        Objective function for the optimizer: RMSE + penalty.

        x_scaled are the scaled values of calibrated parameters.
        """
        iteration["i"] += 1

        # Unscale calibrated parameters
        params_unscaled = {
            p: x_scaled[i] * scale[i]
            for i, p in enumerate(param_names)
        }

        # Add fixed parameters
        params_unscaled.update(fixed_params)

        print(f"\n--- Iteration {iteration['i']} ---")
        print("Scaled params:", x_scaled)
        print("Unscaled calibrated params:", {p: params_unscaled[p] for p in param_names})

        # Run Δ-SNOW for this parameter set
        merged, metrics = evaluate_model_parallel(
            df_fit,
            r_script_path=r_script_path,
            params_unscaled=params_unscaled,
            workers=workers,
        )

        rmse, bias = metrics["rmse"], metrics["bias"]

        # Log iteration
        append_logging(
            log_file=log_file,
            iteration=iteration["i"],
            rmse=rmse,
            bias=bias,
            params_unscaled=params_unscaled,
            param_names=param_names,
        )

        # Penalty term
        penalty = 0.0
        penalty_components = {}

        for i, p in enumerate(param_names):
            if p in penalty_exempt:
                p_pen = 0.0
            else:
                v = params_unscaled[p]
                lo, hi = bounds_unscaled[i]
                p_pen = exp_penalty(v, lo, hi)

            penalty += p_pen
            penalty_components[p] = p_pen

        print("Penalty components:", penalty_components)

        rmse_penalized = rmse + penalty_weight * penalty

        print(
            f"Total penalty={penalty:.3f}, "
            f"RMSE={rmse:.3f}, "
            f"RMSE penalized={rmse_penalized:.3f}"
        )

        return rmse_penalized

    # L-BFGS-B optimization in scaled space
    res = minimize(
        objective,
        x0=np.array(initial_scaled, dtype=float),
        method="L-BFGS-B",
        bounds=bounds_scaled,
        options={"maxiter": maxiter},
    )

    # Recover unscaled calibrated parameters
    final_unscaled_calib = {
        p: res.x[i] * scale[i]
        for i, p in enumerate(param_names)
    }

    # Add fixed parameters to form the full parameter set
    final_unscaled = final_unscaled_calib.copy()
    final_unscaled.update(fixed_params)

    return res, final_unscaled


# ============================================================
# VALIDATION
# ============================================================

def run_validation(
    final_params: dict,
    d_val: dict,
    r_script: str,
    workers: int = 1,
) -> tuple[pd.DataFrame, dict]:
    """
    Run Δ-SNOW validation on held-out winters.

    Parameters
    ----------
    final_params : dict
        Dictionary of final UN-SCALED parameters (calibrated + fixed).
    d_val : dict
        Mapping station_name -> validation DataFrame(s).
    r_script : str
        Path to the Δ-SNOW R runner.
    workers : int
        Number of parallel workers.

    Returns
    -------
    (pandas.DataFrame, dict)
        merged_val : merged validation DataFrame with swe_mod
        metrics_val : {'rmse': float, 'bias': float}
    """
    print("\n=== Running Validation on Held-Out Winters ===")

    # Build tibble for validation winters
    df_val = build_tibble(d_val)

    # Run Δ-SNOW in parallel on the validation dataset
    merged_val, metrics_val = evaluate_model_parallel(
        df_all=df_val,
        r_script_path=r_script,
        params_unscaled=final_params,
        workers=workers,
    )

    rmse_val = metrics_val["rmse"]
    bias_val = metrics_val["bias"]

    print("\n--- VALIDATION RESULTS ---")
    print(f"Validation RMSE: {rmse_val:.3f}")
    print(f"Validation Bias: {bias_val:.3f}\n")

    return merged_val, metrics_val


# ============================================================
# MAIN
# ============================================================

def main():
    """
    Main entry point for the Δ-SNOW calibration pipeline.

    Steps:
      1) Load configuration
      2) Load NetCDF dataset and preprocess
      3) Split into fit and validation winters
      4) Run optimization (L-BFGS-B)
      5) Run validation on held-out winters
      6) Export NetCDF + RMSE/Bias plot + console summary
    """
    cfg = load_config("config_calib.yml")

    # Paths from config
    output_dir = Path(cfg["paths"]["output_dir"])
    nc_file = cfg["paths"]["nc_file"]
    r_script = cfg["paths"]["r_script"]

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Optimization settings
    maxiter = int(cfg["optimization"]["maxiter"])
    penalty_weight = float(cfg["optimization"]["PENALTY_WEIGHT"])

    # Parameter configuration
    param_names = cfg["parameters"]["names"]        # calibrated parameters
    initial_unscaled = [float(v) for v in cfg["parameters"]["initial"]]
    bounds_unscaled = [
        [float(b[0]), float(b[1])]
        for b in cfg["parameters"]["bounds"]
    ]
    scale = [float(s) for s in cfg["parameters"]["scale"]]
    fixed_params = cfg["parameters"]["fixed"]       # dict of fixed parameters

    # Logging
    log_file = output_dir / "calibration_log.csv"
    init_logging(log_file, param_names)

    # Number of workers (use all cores)
    workers = cpu_count()
    print(f"Using {workers} workers.")

    # --------------------------------------------------------
    # LOAD DATASET
    # --------------------------------------------------------
    print("\nLoading dataset...")
    ds = xr.open_dataset(nc_file, engine="netcdf4").rename(
        {"SWE": "swe_obs", "HS": "hs"}
    )
    ds_red = ds[["hs", "swe_obs"]]

    # Convert dataset to station-wise DataFrames
    station_dfs = {
        st: ds_red.sel(station=st)
        .to_dataframe()
        .reset_index()
        .rename(columns={"time": "date"})
        for st in ds_red["station"].values
    }

    print("Splitting dataset into winters (fit vs validation)...")
    d_fit, d_val = split_into_seasons(station_dfs)
    df_fit = build_tibble(d_fit)

    # --------------------------------------------------------
    # OPTIMIZATION
    # --------------------------------------------------------
    print("\nStarting optimization...\n")
    res, final_params = optimize_params(
        df_fit=df_fit,
        r_script_path=r_script,
        param_names=param_names,
        initial_unscaled=initial_unscaled,
        bounds_unscaled=bounds_unscaled,
        scale=scale,
        workers=workers,
        log_file=log_file,
        penalty_weight=penalty_weight,
        maxiter=maxiter,
        fixed_params=fixed_params,
    )

    print("\n=== Optimization complete ===")
    print("Optimizer status:", res.message)
    print("Final parameters (unscaled):")
    for k, v in final_params.items():
        print(f"  {k:10s} = {v:.6g}")

    # --------------------------------------------------------
    # VALIDATION
    # --------------------------------------------------------
    merged_val, metrics_val = run_validation(
        final_params=final_params,
        d_val=d_val,
        r_script=r_script,
        workers=workers,
    )

    # --------------------------------------------------------
    # FINAL MODEL RUN ON ALL FIT WINTERS (for export)
    # --------------------------------------------------------
    merged_final, _ = evaluate_model_parallel(
        df_all=df_fit,
        r_script_path=r_script,
        params_unscaled=final_params,
        workers=workers,
    )

    # --------------------------------------------------------
    # EXPORTS
    # --------------------------------------------------------
    export_netcdf(merged_final, output_dir / "results.nc")
    plot_rmse_bias(log_file, output_dir / "rmse_bias_plot.png")

    print("\nAll outputs saved to:", output_dir)
    print("Done.")


if __name__ == "__main__":
    main()