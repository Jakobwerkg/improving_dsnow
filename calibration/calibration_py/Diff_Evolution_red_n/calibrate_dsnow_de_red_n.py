#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Δ-SNOW calibration pipeline (Differential Evolution)
----------------------------------------------------

Features:
  - Reads dataset + model + parameter setup from a YAML config
  - Uses SciPy's differential_evolution in physical parameter space
  - Calibrates ONLY parameters listed under parameters.names
  - Keeps parameters listed under parameters.fixed constant
  - Evaluates Δ-SNOW (R script) in parallel over stations × winters
  - Logs each evaluation to CSV (iteration, rmse, bias, parameters)
  - Exports final NetCDF and RMSE/Bias plot
  - Optional quick test mode with very small maxiter (--test)
"""

import argparse
import csv
import os
import subprocess
from datetime import timedelta
from io import StringIO
from multiprocessing import Pool, cpu_count, set_start_method
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from scipy.optimize import differential_evolution
from tqdm import tqdm

# ---------------------------------------------------------
# Multiprocessing setup (macOS safe)
# ---------------------------------------------------------

try:
    set_start_method("fork")
except RuntimeError:
    # Already set, ignore
    pass


# ---------------------------------------------------------
# CONFIG LOADER
# ---------------------------------------------------------

def load_config(config_path: str = "config_calib.yml") -> dict:
    """
    Load YAML configuration.

    Parameters
    ----------
    config_path : str
        Name or path of the YAML config file, relative to script directory.

    Returns
    -------
    dict
        Parsed configuration dictionary.
    """
    try:
        base_path = Path(__file__).parent
    except NameError:
        # e.g. interactive
        base_path = Path(os.getcwd())

    cfg_path = base_path / config_path
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config file not found: {cfg_path}")

    with open(cfg_path, "r") as f:
        cfg = yaml.safe_load(f)

    return cfg


# ---------------------------------------------------------
# LOGGING UTILITIES
# ---------------------------------------------------------

def init_logging(log_file: Path, param_names: list[str]) -> None:
    """
    Create a fresh CSV log file for calibration.

    Columns: iteration, rmse, bias, <param1>, <param2>, ...
    """
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with open(log_file, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["iteration", "rmse", "bias"] + list(param_names)
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
    Append one row to the calibration log.

    Only calibrated parameter names are stored (not the fixed ones).
    """
    with open(log_file, "a", newline="") as f:
        writer = csv.writer(f)
        row = [iteration, rmse, bias] + [params_unscaled[p] for p in param_names]
        writer.writerow(row)


def export_netcdf(merged_df: pd.DataFrame, outfile: Path) -> None:
    """
    Export final calibration results (HS, SWE_obs, SWE_mod) as NetCDF.
    """
    ds = xr.Dataset(
        {
            "hs": (("time",), merged_df["hs"].values),
            "swe_obs": (("time",), merged_df["swe_obs"].values),
            "swe_mod": (("time",), merged_df["swe_mod"].values),
        },
        coords={"time": merged_df["date"].values},
    )
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(outfile)
    print(f"Saved NetCDF → {outfile}")


def plot_rmse_bias(log_file: Path, outfile: Path) -> None:
    """
    Plot RMSE and Bias versus iteration from the calibration log.
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

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outfile, dpi=200)
    print(f"Saved plot → {outfile}")


# ---------------------------------------------------------
# DATA PREPARATION
# ---------------------------------------------------------

def build_tibble(d_obs_fit: dict, start_month: int = 8) -> pd.DataFrame:
    """
    Build a long-format DataFrame ("tibble") from station-wise data.

    Each entry in d_obs_fit is a DataFrame indexed by date with columns:
      - 'hs'       : snow depth [m or cm]
      - 'swe_obs'  : observed snow water equivalent

    Returns
    -------
    DataFrame with columns:
      ['date', 'name', 'hs', 'swe_obs', 'block']
    where 'block' is an integer representing hydrological winter.
    """
    rows = []
    for name, df in d_obs_fit.items():
        df = df.copy()
        blocks = [(d.year - 1 if d.month < start_month else d.year) for d in df.index]

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
    station_dfs: dict, season_start: str = "-08-01", season_end: str = "-07-31"
) -> tuple[dict, dict]:
    """
    Split each station's time series into winters and partition the data
    into fit vs validation winters.

    Rules:
      - Winters <150 days are discarded
      - Even starting year → fit
      - Odd  starting year → validation
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

            if len(w) < 150:
                continue

            if w.index[0].year % 2 == 0:
                winters_fit.append(w)
            else:
                winters_val.append(w)

        if winters_fit:
            d_obs_fit[station] = pd.concat(winters_fit)
        if winters_val:
            d_obs_val[station] = pd.concat(winters_val)

    return d_obs_fit, d_obs_val


# ---------------------------------------------------------
# Δ-SNOW PARALLEL RUNNER
# ---------------------------------------------------------

def run_single_block(args):
    """
    Worker for one station+block run of Δ-SNOW (R script).

    args = (station, block, df_sub, r_script_path, params, temp_dir)
    """
    station, block, df_sub, r_script_path, params, temp_dir = args

    # Ensure continuous daily input with hs=0 where missing
    df_sub = df_sub.sort_values("date")
    full_range = pd.date_range(df_sub["date"].min(), df_sub["date"].max(), freq="D")
    df_sub = (
        df_sub.set_index("date")
        .reindex(full_range)
        .rename_axis("date")
        .reset_index()
    )
    df_sub["hs"] = df_sub["hs"].fillna(0)

    temp_dir = Path(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)

    input_csv = temp_dir / f"hs_input_{station}_{block}.csv"
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
    except Exception as e:
        print(f"Error reading Δ-SNOW output for {station} block {block}: {e}")
        return None


def evaluate_model_parallel(
    df_all: pd.DataFrame,
    r_script_path: str,
    params_unscaled: dict,
    workers: int,
    temp_dir: str = "/tmp",
) -> tuple[pd.DataFrame, dict]:
    """
    Evaluate Δ-SNOW for all station×block in df_all using multiprocessing.

    Parameters
    ----------
    df_all : DataFrame
        Tibble with ['date', 'name', 'hs', 'swe_obs', 'block'].
    r_script_path : str
        Path to Δ-SNOW R runner.
    params_unscaled : dict
        Parameter dictionary passed to R (includes calibrated + fixed).
    workers : int
        Number of worker processes.
    temp_dir : str
        Directory for temporary CSVs.

    Returns
    -------
    merged : DataFrame
        Merged observations + swe_mod.
    metrics : dict
        {'rmse': float, 'bias': float}
    """
    stations = df_all["name"].unique()
    blocks = df_all["block"].unique()

    jobs = []
    for st in stations:
        for blk in blocks:
            df_sub = df_all[(df_all["name"] == st) & (df_all["block"] == blk)]
            if len(df_sub) < 100:
                continue
            jobs.append((st, blk, df_sub.copy(), r_script_path, params_unscaled, temp_dir))

    if not jobs:
        raise RuntimeError("No jobs constructed for Δ-SNOW evaluation.")

    with Pool(workers) as pool:
        results = list(
            tqdm(
                pool.imap_unordered(run_single_block, jobs),
                total=len(jobs),
                desc="Δ-SNOW blocks",
            )
        )

    results = [r for r in results if r is not None]
    if not results:
        raise RuntimeError("No Δ-SNOW output from any block.")

    out = pd.concat(results, ignore_index=True)

    out = out[["date", "station", "block", "swe_mod"]]

    merged = df_all.merge(
        out.rename(columns={"station": "name"}),
        on=["date", "name", "block"],
        how="inner",
    )
    valid = merged.dropna(subset=["swe_obs", "swe_mod"])

    if valid.empty:
        raise RuntimeError("No valid swe_obs/swe_mod pairs after merge.")

    rmse = float(np.sqrt(np.mean((valid["swe_mod"] - valid["swe_obs"]) ** 2)))
    bias = float(np.mean(valid["swe_mod"] - valid["swe_obs"]))

    return merged, {"rmse": rmse, "bias": bias}


# ---------------------------------------------------------
# OPTIMIZATION — DIFFERENTIAL EVOLUTION (SIMPLE)
# ---------------------------------------------------------

def optimize_params_de(
    df_fit: pd.DataFrame,
    r_script_path: str,
    param_names: list[str],
    bounds_unscaled: list[list[float]],
    workers_model: int,
    log_file: Path,
    maxiter: int,
    test_mode: bool = False,
    fixed_params: dict | None = None,
) -> tuple[object, dict]:
    """
    Optimize parameters using SciPy's differential_evolution.

    - Works directly in physical units (unscaled bounds).
    - Only parameters in param_names are optimized.
    - fixed_params are always added and NOT optimized.
    """
    if fixed_params is None:
        fixed_params = {}

    iteration = {"i": 0}

    def objective(x: np.ndarray) -> float:
        """
        Objective function for differential_evolution.

        x is a 1D vector of length len(param_names).
        """
        iteration["i"] += 1

        # Combine free + fixed parameters
        params_unscaled = {p: float(x[i]) for i, p in enumerate(param_names)}
        params_unscaled.update(fixed_params)

        print(f"\n--- Evaluation {iteration['i']} ---")
        print("Parameters:", params_unscaled)

        merged, metrics = evaluate_model_parallel(
            df_fit,
            r_script_path=r_script_path,
            params_unscaled=params_unscaled,
            workers=workers_model,
        )

        rmse = metrics["rmse"]
        bias = metrics["bias"]

        append_logging(
            log_file=log_file,
            iteration=iteration["i"],
            rmse=rmse,
            bias=bias,
            params_unscaled=params_unscaled,
            param_names=param_names,
        )

        print(f"RMSE = {rmse:.3f}, Bias = {bias:.3f}")

        # Objective: pure RMSE
        return rmse

    print("\nUsing Differential Evolution optimizer.")
    print(f"Max generations (maxiter) = {maxiter} | Test mode = {test_mode}")

    res = differential_evolution(
        objective,
        bounds=bounds_unscaled,
        strategy="best1bin",
        maxiter=maxiter,
        disp=True,
        workers=1,  # keep DE itself serial; model eval is parallel
        polish=True,
        seed=42,
    )

    final_unscaled = {p: float(res.x[i]) for i, p in enumerate(param_names)}
    final_unscaled.update(fixed_params)

    return res, final_unscaled


# ---------------------------------------------------------
# VALIDATION
# ---------------------------------------------------------

def run_validation(
    final_params: dict,
    d_val: dict,
    r_script: str,
    workers: int = 1,
) -> tuple[pd.DataFrame, dict]:
    """
    Run Δ-SNOW on validation winters using the final parameter set.
    """
    print("\n=== Running Validation on Held-Out Winters ===")

    df_val = build_tibble(d_val)

    merged_val, metrics_val = evaluate_model_parallel(
        df_all=df_val,
        r_script_path=r_script,
        params_unscaled=final_params,
        workers=workers,
    )

    print("\n--- VALIDATION RESULTS ---")
    print(f"Validation RMSE: {metrics_val['rmse']:.3f}")
    print(f"Validation Bias: {metrics_val['bias']:.3f}\n")

    return merged_val, metrics_val


# ---------------------------------------------------------
# MAIN
# ---------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Δ-SNOW calibration using Differential Evolution."
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config_calib.yml",
        help="Path to YAML configuration file.",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Run a quick test with very few optimization iterations.",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)

    # Paths
    output_dir = Path(cfg["paths"]["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    log_file = output_dir / "calibration_log.csv"
    nc_file = cfg["paths"]["nc_file"]
    r_script = cfg["paths"]["r_script"]

    # Optimization setup
    maxiter_cfg = int(cfg["optimization"]["maxiter"])
    maxiter_test = int(cfg["optimization"].get("maxiter_test", 2))
    maxiter = maxiter_test if args.test else maxiter_cfg

    # Parameter setup
    param_names = list(cfg["parameters"]["names"])
    bounds_unscaled = [
        [float(b[0]), float(b[1])] for b in cfg["parameters"]["bounds"]
    ]
    fixed_params = cfg["parameters"].get("fixed", {})

    print("\nFree parameters:", param_names)
    print("Fixed parameters:", fixed_params)

    # Initialize logging
    init_logging(log_file, param_names)

    # Workers for model evaluation
    workers_model = cpu_count()
    print(f"\nUsing {workers_model} workers for Δ-SNOW evaluation.")

    # Load dataset
    print("\nLoading dataset...")
    ds = xr.open_dataset(nc_file, engine="netcdf4").rename(
        {"SWE": "swe_obs", "HS": "hs"}
    )
    ds_red = ds[["hs", "swe_obs"]]

    # Convert to station-wise DataFrames
    station_dfs = {
        st: ds_red.sel(station=st)
        .to_dataframe()
        .reset_index()
        .rename(columns={"time": "date"})
        for st in ds_red["station"].values
    }

    print("Splitting dataset into winters (fit / validation)...")
    d_fit, d_val = split_into_seasons(station_dfs)
    df_fit = build_tibble(d_fit)

    # Optimization
    print("\nStarting optimization with Differential Evolution...\n")
    res, final_params = optimize_params_de(
        df_fit=df_fit,
        r_script_path=r_script,
        param_names=param_names,
        bounds_unscaled=bounds_unscaled,
        workers_model=workers_model,
        log_file=log_file,
        maxiter=maxiter,
        test_mode=args.test,
        fixed_params=fixed_params,
    )

    print("\n=== Optimization complete ===")
    print("Success:", res.success)
    print("Message:", res.message)
    print("\nFinal unscaled parameters:")
    for k, v in final_params.items():
        print(f"  {k:10s} = {v:.6g}")

    # Validation
    print("\n=== Running validation on held-out winters ===")
    merged_val, metrics_val = run_validation(
        final_params=final_params,
        d_val=d_val,
        r_script=r_script,
        workers=workers_model,
    )

    # Re-run model on fit data with final parameters for export
    merged_final, _ = evaluate_model_parallel(
        df_all=df_fit,
        r_script_path=r_script,
        params_unscaled=final_params,
        workers=workers_model,
    )

    # Exports
    export_netcdf(merged_final, output_dir / "results.nc")
    plot_rmse_bias(log_file, output_dir / "rmse_bias_plot.png")

    print("\nAll outputs saved to:", output_dir)


if __name__ == "__main__":
    main()