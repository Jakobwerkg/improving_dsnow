#!/usr/bin/env python3
"""
HS2SWE calibration pipeline

Features:
 - YAML-based configuration
 - Parameter scaling for the optimizer
 - RMSE + bias objective with optional penalties
 - Training / validation split by winters
 - Multi-station support
 - NetCDF export of final fit
 - CSV log of iterations
 - RMSE/Bias figure export
"""

import os
import csv
from pathlib import Path
from datetime import timedelta
from io import StringIO

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import yaml
from scipy.optimize import minimize
from multiprocessing import Pool, cpu_count, set_start_method
from tqdm import tqdm

# HS2SWE model import (adjust if module name differs)
from swe_mod_slf import HS2SWE


# Try to enable fast multiprocessing on macOS
try:
    set_start_method("fork")
except RuntimeError:
    pass


# -------------------------------------------------------------------
# CONFIG
# -------------------------------------------------------------------

def load_config(config_path: str = "hs2swe_calib.yml") -> dict:
    base_path = Path(__file__).parent if "__file__" in globals() else Path.cwd()
    cfg_path = base_path / config_path
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config file not found: {cfg_path}")
    with open(cfg_path, "r") as f:
        return yaml.safe_load(f)


# -------------------------------------------------------------------
# LOGGING
# -------------------------------------------------------------------

def init_logging(log_file: Path, param_names):
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("w", newline="") as f:
        writer = csv.writer(f)
        header = ["iteration", "rmse", "bias"] + param_names
        writer.writerow(header)


def append_logging(log_file: Path, iteration, rmse, bias, params_unscaled, param_names):
    with log_file.open("a", newline="") as f:
        writer = csv.writer(f)
        row = [iteration, rmse, bias] + [params_unscaled[p] for p in param_names]
        writer.writerow(row)


def export_netcdf(merged_df: pd.DataFrame, outfile: Path):
    ds = xr.Dataset(
        data_vars={
            "hs": ("time", merged_df["hs"].values),
            "swe_obs": ("time", merged_df["swe_obs"].values),
            "swe_mod": ("time", merged_df["swe_mod"].values),
        },
        coords={"time": merged_df["date"].values},
    )
    outfile.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(outfile)
    print(f"Saved NetCDF → {outfile}")


def plot_rmse_bias(log_file: Path, outfile: Path):
    df = pd.read_csv(log_file)
    plt.figure(figsize=(8, 5))
    plt.plot(df["iteration"], df["rmse"], label="RMSE")
    plt.plot(df["iteration"], df["bias"], label="Bias")
    plt.xlabel("Iteration")
    plt.ylabel("Value")
    plt.title("RMSE & Bias during HS2SWE Calibration")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    outfile.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outfile, dpi=200)
    print(f"Saved plot → {outfile}")


# -------------------------------------------------------------------
# DATA PREPARATION
# -------------------------------------------------------------------

def load_station_csvs(data_dir: Path) -> dict:
    """
    Load all station CSVs from data_dir.
    Each file must contain: date, hs, swe_obs
    Returns dict: {station_name: DataFrame}
    """
    station_dfs = {}
    for f in sorted(data_dir.glob("*.csv")):
        df = pd.read_csv(f, parse_dates=["date"])
        if "hs" not in df.columns:
            raise ValueError(f"{f} missing 'hs' column")
        if "swe_obs" not in df.columns:
            # create empty observed SWE if not present
            df["swe_obs"] = np.nan
        station_name = f.stem  # e.g., Adelboden
        station_dfs[station_name] = df
    return station_dfs


def split_into_seasons(
    station_dfs: dict,
    season_start: str = "-08-01",
    season_end: str = "-07-31",
    min_days: int = 150,
):
    """
    Split each station into winters and alternate:
    - even start year -> fit
    - odd start year  -> validation
    """
    d_fit, d_val = {}, {}

    for station, df in station_dfs.items():
        df = df.copy().set_index("date").sort_index()
        years = df.index.year.unique()

        winters_fit, winters_val = [], []

        for y in years:
            start = pd.to_datetime(f"{y}{season_start}")
            end = pd.to_datetime(f"{y+1}{season_end}") + timedelta(days=1)
            w = df[(df.index >= start) & (df.index < end)]
            if len(w) < min_days:
                continue

            if start.year % 2 == 0:
                winters_fit.append(w)
            else:
                winters_val.append(w)

        if winters_fit:
            d_fit[station] = pd.concat(winters_fit)
        if winters_val:
            d_val[station] = pd.concat(winters_val)

    return d_fit, d_val


def build_tibble(d_obs_fit: dict, start_month: int = 8) -> pd.DataFrame:
    """
    Build a long-format tibble from dict of dataframes.
    """
    rows = []
    for name, df in d_obs_fit.items():
        df = df.copy().sort_index()
        blocks = [
            d.year - 1 if d.month < start_month else d.year
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


# -------------------------------------------------------------------
# HS2SWE RUNNER (INSTEAD OF R)
# -------------------------------------------------------------------

def hs2swe_from_block(station: str, block: int, df_sub: pd.DataFrame, params_unscaled: dict):
    """
    Run HS2SWE for one station/block combo and return a dataframe with swe_mod.
    """
    df_sub = df_sub.sort_values("date").copy()

    # Ensure full daily range
    full_range = pd.date_range(df_sub["date"].min(), df_sub["date"].max(), freq="D")
    df_sub = df_sub.set_index("date").reindex(full_range).rename_axis("date").reset_index()
    df_sub["hs"] = df_sub["hs"].fillna(0)

    # Prepare input array for HS2SWE: shape (T, 1)
    idata = np.transpose(np.array(df_sub["hs"].values, ndmin=2))

    # Map config parameters onto HS2SWE arguments
    hs2swe_kwargs = dict(
        RhoMax=params_unscaled["rho.max"],
        RhoNew=params_unscaled["rho.null"],
        c1=params_unscaled["c.ov"],
        c2=params_unscaled["k.ov"],
        c3=params_unscaled["k.exp"],
        c4=params_unscaled["tau"],
        Visc=params_unscaled["eta.null"],
    )

    swe_sim = HS2SWE(idata, **hs2swe_kwargs)
    swe_sim = swe_sim.squeeze()  # (T, 1) -> (T,)

    out = pd.DataFrame(
        {
            "date": df_sub["date"].values,
            "name": station,
            "block": block,
            "swe_mod": swe_sim,
        }
    )
    return out


def run_single_block(args):
    station, block, df_sub, params_unscaled = args
    try:
        return hs2swe_from_block(station, block, df_sub, params_unscaled)
    except Exception as e:
        print(f"Error in station={station}, block={block}: {e}")
        return None


def evaluate_model_parallel(
    df_all: pd.DataFrame,
    params_unscaled: dict,
    workers: int,
) -> tuple[pd.DataFrame, dict]:
    """
    Evaluate HS2SWE for all station/block combinations in parallel.
    """
    stations = df_all["name"].unique()
    blocks = df_all["block"].unique()

    jobs = []
    for st in stations:
        for blk in blocks:
            df_sub = df_all[(df_all["name"] == st) & (df_all["block"] == blk)]
            if len(df_sub) < 100:
                continue
            jobs.append((st, blk, df_sub.copy(), params_unscaled))

    if not jobs:
        raise RuntimeError("No jobs to evaluate. Check data and block setup.")

    with Pool(processes=workers) as pool:
        results = list(
            tqdm(
                pool.imap_unordered(run_single_block, jobs),
                total=len(jobs),
                desc="HS2SWE blocks",
            )
        )

    results = [r for r in results if r is not None]
    if not results:
        raise RuntimeError("No HS2SWE output produced.")

    out = pd.concat(results, ignore_index=True)

    merged = df_all.merge(out, on=["date", "name", "block"], how="inner")
    valid = merged.dropna(subset=["swe_obs", "swe_mod"])

    rmse = np.sqrt(np.mean((valid["swe_mod"] - valid["swe_obs"]) ** 2))
    bias = np.mean(valid["swe_mod"] - valid["swe_obs"])

    metrics = {"rmse": rmse, "bias": bias}
    return merged, metrics


# -------------------------------------------------------------------
# PENALTY
# -------------------------------------------------------------------

def exp_penalty(value: float, lower: float, upper: float) -> float:
    """
    Exponential penalty when parameter is outside the central 50% of [lower, upper].
    """
    mid = 0.5 * (lower + upper)
    half_range = 0.25 * (upper - lower)

    if (value >= lower + half_range) and (value <= upper - half_range):
        return 0.0

    if value < mid:
        distance = (lower + half_range - value) / half_range
    else:
        distance = (value - (upper - half_range)) / half_range

    return np.exp(5 * distance) - 1.0


# -------------------------------------------------------------------
# OPTIMIZATION
# -------------------------------------------------------------------

def optimize_params(
    df_fit: pd.DataFrame,
    param_names,
    initial_unscaled,
    bounds_unscaled,
    scale,
    workers: int,
    log_file: Path,
    penalty_weight: float,
    maxiter: int,
):
    initial_scaled = [v / s for v, s in zip(initial_unscaled, scale)]
    bounds_scaled = [(lo / s, hi / s) for (lo, hi), s in zip(bounds_unscaled, scale)]

    iteration = {"i": 0}
    penalty_exempt = set()  # or e.g. {"tau"} if you want to exempt some

    def objective(x_scaled):
        iteration["i"] += 1

        params_unscaled = {
            p: x_scaled[i] * scale[i] for i, p in enumerate(param_names)
        }

        print(f"\n--- Iteration {iteration['i']} ---")
        print("Scaled params:", np.round(x_scaled, 6))
        print("Unscaled params:", {k: round(v, 6) for k, v in params_unscaled.items()})

        merged, metrics = evaluate_model_parallel(
            df_all=df_fit,
            params_unscaled=params_unscaled,
            workers=workers,
        )

        rmse, bias = metrics["rmse"], metrics["bias"]

        append_logging(
            log_file,
            iteration["i"],
            rmse,
            bias,
            params_unscaled,
            param_names,
        )

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

        print("Penalty components:", {k: round(v, 4) for k, v in penalty_components.items()})

        rmse_penalized = rmse + penalty_weight * penalty
        print(
            f"Total penalty={penalty:.3f}, RMSE={rmse:.3f}, RMSE penalized={rmse_penalized:.3f}"
        )

        return rmse_penalized

    res = minimize(
        objective,
        initial_scaled,
        method="L-BFGS-B",
        bounds=bounds_scaled,
        options={"maxiter": maxiter},
    )

    final_unscaled = {
        p: res.x[i] * scale[i] for i, p in enumerate(param_names)
    }

    return res, final_unscaled


def run_validation(final_params: dict, d_val: dict, workers: int):
    print("\n=== Running validation on held-out winters ===")
    if not d_val:
        print("No validation data found (d_val is empty). Skipping validation.")
        return None, {"rmse": np.nan, "bias": np.nan}

    df_val = build_tibble(d_val)
    merged_val, metrics_val = evaluate_model_parallel(
        df_all=df_val,
        params_unscaled=final_params,
        workers=workers,
    )

    print("\n--- VALIDATION RESULTS ---")
    print(f"Validation RMSE: {metrics_val['rmse']:.3f}")
    print(f"Validation Bias: {metrics_val['bias']:.3f}")
    return merged_val, metrics_val


# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------

def main():
    cfg = load_config("hs2swe_calib.yml")

    data_dir = Path(cfg["paths"]["data_dir"])
    output_dir = Path(cfg["paths"]["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    maxiter = int(cfg["optimization"]["maxiter"])
    penalty_weight = float(cfg["optimization"]["PENALTY_WEIGHT"])

    param_names = cfg["parameters"]["names"]
    initial_unscaled = [float(v) for v in cfg["parameters"]["initial"]]
    bounds_unscaled = [[float(b[0]), float(b[1])] for b in cfg["parameters"]["bounds"]]
    scale = [float(s) for s in cfg["parameters"]["scale"]]

    log_file = output_dir / "calibration_log.csv"
    init_logging(log_file, param_names)

    workers = cpu_count()

    print("\nLoading station CSVs...")
    station_dfs = load_station_csvs(data_dir)

    print("Splitting dataset into winters (fit/validation)...")
    d_fit, d_val = split_into_seasons(station_dfs)
    df_fit = build_tibble(d_fit)

    print("\nStarting HS2SWE parameter optimization...\n")
    res, final_params = optimize_params(
        df_fit=df_fit,
        param_names=param_names,
        initial_unscaled=initial_unscaled,
        bounds_unscaled=bounds_unscaled,
        scale=scale,
        workers=workers,
        log_file=log_file,
        penalty_weight=penalty_weight,
        maxiter=maxiter,
    )

    print("\n=== Optimization complete ===")
    for k, v in final_params.items():
        print(f"{k}: {v}")

    merged_val, metrics_val = run_validation(
        final_params=final_params,
        d_val=d_val,
        workers=workers,
    )

    # Re-run on fit data to export full series with swe_mod
    merged_final, _ = evaluate_model_parallel(
        df_all=df_fit,
        params_unscaled=final_params,
        workers=workers,
    )

    export_netcdf(merged_final, output_dir / "results_fit.nc")
    plot_rmse_bias(log_file, output_dir / "rmse_bias.png")

    print("\nAll outputs saved to:", output_dir)


if __name__ == "__main__":
    main()