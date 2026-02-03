#!/usr/bin/env python3

"""
Δ-SNOW calibration pipeline with:
 - unscaled parameters from config
 - internal optimizer scaling
 - per-iteration parameter printing
 - robust parallelization
 - NetCDF export
 - CSV RMSE log
 - RMSE/Bias figure export
 - clean logging directory
"""

import argparse
import yaml
import numpy as np
import pandas as pd
import xarray as xr
from datetime import timedelta
from pathlib import Path
import subprocess
from io import StringIO
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import matplotlib.pyplot as plt
import csv
from scipy.optimize import minimize
import os

# Enable fast multiprocessing on macOS
from multiprocessing import set_start_method
try:
    set_start_method("fork")
except RuntimeError:
    pass

# ---------------------------------------------------------
# CONFIG LOADER
# ---------------------------------------------------------

def load_config(config_path="config_calib.yml"):
    try:
        base_path = Path(__file__).parent
    except NameError:
        base_path = Path(os.getcwd())

    config_path = base_path / config_path
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)
    return cfg


# ---------------------------------------------------------
# LOGGING UTILITIES
# ---------------------------------------------------------

def init_logging(log_file, param_names):
    with open(log_file, "w") as f:
        writer = csv.writer(f)
        header = ["iteration", "rmse", "bias"] + param_names
        writer.writerow(header)

def append_logging(log_file, iteration, rmse, bias, params_unscaled, param_names):
    with open(log_file, "a") as f:
        writer = csv.writer(f)
        row = [iteration, rmse, bias] + [params_unscaled[p] for p in param_names]
        writer.writerow(row)


def export_netcdf(merged_df, outfile):
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


def plot_rmse_bias(log_file, outfile):
    df = pd.read_csv(log_file)
    plt.figure(figsize=(8,5))
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


# ---------------------------------------------------------
# DATA PREPARATION
# ---------------------------------------------------------

def build_tibble(d_obs_fit, start_month=8):
    rows = []
    for name, df in d_obs_fit.items():
        df = df.copy()
        blocks = [(d.year - 1 if d.month < start_month else d.year) for d in df.index]

        tib = pd.DataFrame({
            "date": df.index,
            "name": name,
            "hs": df["hs"].values,
            "swe_obs": df["swe_obs"].values,
            "block": blocks,
        })
        rows.append(tib)

    return pd.concat(rows, ignore_index=True)


def split_into_seasons(station_dfs, season_start="-08-01", season_end="-07-31"):
    d_obs_fit, d_obs_val = {}, {}

    for station, df in station_dfs.items():
        df = df.copy().set_index("date")
        years = df.index.year.unique()

        winters_fit, winters_val = [], []

        for y in years:
            start = pd.to_datetime(f"{y}{season_start}")
            end = pd.to_datetime(f"{y+1}{season_end}") + timedelta(days=1)
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
    station, block, df_sub, r_script_path, params, temp_dir = args

    df_sub = df_sub.sort_values("date")
    full_range = pd.date_range(df_sub["date"].min(), df_sub["date"].max(), freq="D")
    df_sub = df_sub.set_index("date").reindex(full_range).rename_axis("date").reset_index()
    df_sub["hs"] = df_sub["hs"].fillna(0)

    input_csv = Path(temp_dir) / f"hs_input_{station}_{block}.csv"
    df_sub[["date", "hs"]].to_csv(input_csv, index=False)

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
    except:
        return None


def evaluate_model_parallel(df_all, r_script_path, params_unscaled, workers, temp_dir="/tmp"):

    stations = df_all["name"].unique()
    blocks = df_all["block"].unique()

    jobs = []
    for st in stations:
        for blk in blocks:
            df_sub = df_all[(df_all["name"] == st) & (df_all["block"] == blk)]
            if len(df_sub) < 100:
                continue
            jobs.append((st, blk, df_sub.copy(), r_script_path, params_unscaled, temp_dir))

    # --- PARALLEL Δ-SNOW EXECUTION (FASTER VERSION) ---
    with Pool(workers) as pool:
        results = list(
            tqdm(
                pool.imap_unordered(run_single_block, jobs),
                total=len(jobs),
                desc="Δ-SNOW blocks"
            )
        )

    # Filter None
    results = [r for r in results if r is not None]
    if not results:
        raise RuntimeError("No Δ-SNOW output.")

    out = pd.concat(results, ignore_index=True)

    # --- Keep only useful outputs to avoid hs_x / hs_y conflicts ---
    out = out[["date", "station", "block", "swe_mod"]]

    # Merge using date + name + block (safe, unique)
    merged = df_all.merge(
        out.rename(columns={"station": "name"}),
        on=["date", "name", "block"],
        how="inner"
    )
    valid = merged.dropna(subset=["swe_obs", "swe_mod"])

    rmse = np.sqrt(np.mean((valid["swe_mod"] - valid["swe_obs"]) ** 2))
    bias = np.mean(valid["swe_mod"] - valid["swe_obs"])

    return merged, {"rmse": rmse, "bias": bias}




# ---------------------------------------------------------
# Penalty — outside mid 50 percentile add an error
# ---------------------------------------------------------

def exp_penalty(value, lower, upper):
    """
    Exponential penalty when parameter is outside the middle 50% of its range.
    """
    mid = 0.5 * (lower + upper)
    half_range = 0.25 * (upper - lower)  # central 50% safe zone

    # Inside safe central region → no penalty
    if (value >= lower + half_range) and (value <= upper - half_range):
        return 0.0

    # Normalized distance into penalty zone (0 → 1)
    if value < mid:
        distance = (lower + half_range - value) / half_range
    else:
        distance = (value - (upper - half_range)) / half_range

    # Exponential penalty that explodes near bounds
    return np.exp(5 * distance) - 1


# ---------------------------------------------------------
# OPTIMIZATION — INTERNAL SCALING
# ---------------------------------------------------------

def optimize_params(df_fit, r_script_path, param_names,
                    initial_unscaled, bounds_unscaled, scale,
                    workers, log_file, PENALTY_WEIGHT, maxiter):

    initial_scaled = [v / s for v, s in zip(initial_unscaled, scale)]
    bounds_scaled  = [(lo / s, hi / s) for (lo, hi), s in zip(bounds_unscaled, scale)]

    iteration = {"i": 0}


        # Parameters that should NOT receive penalties
    penalty_exempt = {"tau", "c.ov", "k.ov"}

    # Parameters that should NOT receive penalties
    penalty_exempt = {"tau", "c.ov", "k.ov"}

    def objective(x_scaled):
        iteration["i"] += 1

        # Unscale parameters for R
        params_unscaled = {p: x_scaled[i] * scale[i] for i, p in enumerate(param_names)}

        print(f"\n--- Iteration {iteration['i']} ---")
        print("Scaled params:", x_scaled)
        print("Unscaled params:", params_unscaled)

        merged, metrics = evaluate_model_parallel(df_fit, r_script_path,
                                                params_unscaled, workers)

        rmse, bias = metrics["rmse"], metrics["bias"]
        append_logging(
            log_file,
            iteration["i"],
            rmse,
            bias,
            params_unscaled,
            param_names
        )

        # ============================
        # Penalty (now selective)
        # ============================
        penalty = 0.0
        penalty_components = {}

        for i, p in enumerate(param_names):
            if p in penalty_exempt:
                p_pen = 0.0   # no penalty for tau, c.ov, k.ov
            else:
                v = params_unscaled[p]
                lo, hi = bounds_unscaled[i]
                p_pen = exp_penalty(v, lo, hi)

            penalty += p_pen
            penalty_components[p] = p_pen

        print("Penalty components:", penalty_components)

        rmse_penalized = rmse + PENALTY_WEIGHT * penalty

        print(f"Total penalty={penalty:.3f}, RMSE={rmse:.3f}, RMSE penalized={rmse_penalized:.3f}")

        return rmse_penalized

        

    # Run optimizer
    res = minimize(
        objective,
        initial_scaled,
        method="L-BFGS-B",
        bounds=bounds_scaled,
        options={"maxiter": maxiter}
    )

    final_unscaled = {p: res.x[i] * scale[i] for i, p in enumerate(param_names)}
    return res, final_unscaled


def run_validation(final_params, d_val, r_script, workers=1):
    """
    Run Δ-SNOW validation using the OUT-OF-SAMPLE validation dataset.

    Parameters
    ----------
    final_params : dict
        Dictionary of UN-SCALED calibrated parameters.
    d_val : dict
        Dictionary of validation dataframes (same structure as d_obs_val).
    r_script : str
        Path to the Δ-SNOW R runner.
    workers : int
        Number of parallel workers to use.
    """

    print("\n=== Running Validation on Held-Out Winters ===")

    # Build tibble for validation winters
    df_val = build_tibble(d_val)

    # Run Δ-SNOW in parallel on the validation dataset
    merged_val, metrics_val = evaluate_model_parallel(
        df_val,
        r_script_path=r_script,
        params_unscaled=final_params,
        workers=workers
    )

    rmse_val = metrics_val["rmse"]
    bias_val = metrics_val["bias"]

    print("\n--- VALIDATION RESULTS ---")
    print(f"Validation RMSE: {rmse_val:.3f}")
    print(f"Validation Bias: {bias_val:.3f}\n")

    return merged_val, metrics_val



# ---------------------------------------------------------
# MAIN
# ---------------------------------------------------------

def main():

    cfg = load_config("config_calib.yml")

    # Logging directory
    output_dir = Path(cfg['paths']['output_dir'])


    log_file = output_dir / "calibration_log.csv"
    workers = cpu_count() 

    # --- LOAD PARAMETERS BEFORE LOGGING ---
    param_names      = cfg["parameters"]["names"]
    initial_unscaled = [float(v) for v in cfg["parameters"]["initial"]]
    bounds_unscaled  = [[float(b[0]), float(b[1])] for b in cfg["parameters"]["bounds"]]
    scale            = [float(s) for s in cfg["parameters"]["scale"]]

    # NOW param_names is defined → safe to call
    init_logging(log_file, param_names)



    log_file = output_dir / "calibration_log.csv"
    init_logging(log_file, param_names)

    nc_file = cfg["paths"]["nc_file"]
    r_script = cfg["paths"]["r_script"]
    maxiter = cfg['optimization']['maxiter']
    PENALTY_WEIGHT = cfg['optimization']['PENALTY_WEIGHT']

        # use ALL cores


    print("\nLoading dataset...")
    ds = xr.open_dataset(nc_file, engine="netcdf4").rename({"SWE": "swe_obs",
                                                            "HS": "hs"})
    ds_red = ds[["hs", "swe_obs"]]

    station_dfs = {
        st: ds_red.sel(station=st).to_dataframe()
                  .reset_index()
                  .rename(columns={"time": "date"})
        for st in ds_red["station"].values
    }

    print("Splitting dataset into winters...")
    d_fit, d_val = split_into_seasons(station_dfs)
    df_fit = build_tibble(d_fit)


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
        PENALTY_WEIGHT=PENALTY_WEIGHT,
        maxiter=maxiter          # ← add this
    )

    print("\n=== Optimization complete ===")
    for k, v in final_params.items():
        print(f"{k}: {v}")

    print("\n=== Running validation on held-out winters ===\n")
    merged_val, metrics_val = run_validation(
        final_params=final_params,
        d_val=d_val,
        r_script=r_script,
        workers=workers
    )

    merged_final, _ = evaluate_model_parallel(df_fit, r_script,
                                              final_params, workers)

    export_netcdf(merged_final, output_dir / "results.nc")
    plot_rmse_bias(log_file, output_dir / "rmse_bias_plot.png")

    print("\nAll outputs saved to:", output_dir)


if __name__ == "__main__":
    main()