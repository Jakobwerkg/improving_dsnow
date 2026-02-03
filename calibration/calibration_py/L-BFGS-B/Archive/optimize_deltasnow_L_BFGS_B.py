#!/usr/bin/env python3
"""
Simple L-BFGS-B calibration for Δ-SNOW using an external R runner.

Usage:
    python calibrate_delta_snow.py --config config.yaml [--disp] [--n-jobs N]

Config YAML needs at least:
    base_dir, data_dir, log_dir, r_runner,
    season_start, season_end,
    par_guess, par_scale, lower, upper
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess as sp
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import yaml
from concurrent.futures import ThreadPoolExecutor, as_completed
from scipy.optimize import minimize

BIG_PENALTY = 1e12
MIN_LEN_WINTER = 200  # minimal number of days per winter


# -------------------------- Parameter container -------------------------- #

@dataclass(frozen=True)
class Params:
    rho_max: float
    rho_null: float
    c_ov: float
    k_ov: float
    k: float
    tau: float
    eta_null: float


# ------------------------------ Data loading ----------------------------- #

def load_station_data(data_dir: Path) -> Dict[str, pd.DataFrame]:
    """
    Load all CSVs in data_dir.

    Each file must have: date, hs, swe_obs
    Station name = filename stem (optionally without '_hs_swe_obs').
    """
    files = sorted(data_dir.glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No CSV files found in: {data_dir}")

    stations: Dict[str, pd.DataFrame] = {}

    for f in files:
        df = pd.read_csv(f)
        if not {"date", "hs", "swe_obs"}.issubset(df.columns):
            print(f"[WARN] Skipping {f.name}: missing columns.")
            continue

        df = df[["date", "hs", "swe_obs"]].copy()
        df["date"] = pd.to_datetime(df["date"], errors="coerce")
        df["hs"] = pd.to_numeric(df["hs"], errors="coerce")
        df["swe_obs"] = pd.to_numeric(df["swe_obs"], errors="coerce")

        df = df.dropna(subset=["date"]).sort_values("date").drop_duplicates("date")
        if df.empty:
            print(f"[WARN] Skipping {f.name}: empty after cleaning.")
            continue

        st_name = f.stem
        st_name = re.sub(r"_hs_swe_obs$", "", st_name)
        stations[st_name] = df.reset_index(drop=True)

    if not stations:
        raise RuntimeError("No usable station CSVs after cleaning.")
    return stations


# ---------------------------- Winter splitting --------------------------- #

def season_slice(year: int, start: str, end: str) -> Tuple[pd.Timestamp, pd.Timestamp]:
    """
    Return [left, right) for hydrological season `year`.

    Example: start = "-08-01", end = "-07-31"
      season year = [year-08-01, (year+1)-07-31 + 1 day)
    """
    left = pd.Timestamp(f"{year}{start}")
    right = pd.Timestamp(f"{year + 1}{end}") + pd.Timedelta(days=1)
    return left, right


def split_into_winters(
    df: pd.DataFrame,
    season_start: str,
    season_end: str,
    min_len: int = MIN_LEN_WINTER,
) -> List[pd.DataFrame]:
    """Split a station time series into hydrological winters."""
    if df.empty:
        return []

    years = sorted(df["date"].dt.year.unique())
    winters: List[pd.DataFrame] = []

    for y in years[:-1]:
        left, right = season_slice(y, season_start, season_end)
        win = df[(df["date"] >= left) & (df["date"] < right)].copy()
        if len(win) >= min_len:
            winters.append(win.reset_index(drop=True))

    return winters


def alt_split_fit_val(
    stations: Dict[str, pd.DataFrame],
    season_start: str,
    season_end: str,
) -> Tuple[List[pd.DataFrame], List[pd.DataFrame]]:
    """
    Split all station winters into fit / validation (ALT-winter split).

    Even winters  -> fit
    Odd winters   -> validation
    """
    fit_blocks: List[pd.DataFrame] = []
    val_blocks: List[pd.DataFrame] = []

    for name, df in stations.items():
        winters = split_into_winters(df, season_start, season_end)
        if not winters:
            continue

        for i, w in enumerate(winters, start=1):
            if i % 2 == 0:
                fit_blocks.append(w)
            else:
                val_blocks.append(w)

    return fit_blocks, val_blocks


# ------------------------------ R bridge --------------------------------- #

def run_r_delta(
    hs_df: pd.DataFrame,
    par: Params,
    r_runner: Path,
    tz: str,
) -> pd.DataFrame:
    """
    Run external R script to compute modeled SWE from snow depth.

    Returns a DataFrame with columns: date, swe_mod
    """
    if not r_runner.exists():
        raise FileNotFoundError(f"R runner not found: {r_runner}")

    # --- Preprocess hs_df to daily, gap-free, non-negative series ---
    hs = hs_df.copy()
    hs["date"] = pd.to_datetime(hs["date"], errors="coerce").dt.normalize()
    hs["hs"] = pd.to_numeric(hs["hs"], errors="coerce")

    hs = hs.dropna(subset=["date"])
    if hs.empty:
        raise RuntimeError("No valid dates in hs_df.")

    # Full daily range
    date_range = pd.date_range(hs["date"].min(), hs["date"].max(), freq="D")
    full = pd.DataFrame({"date": date_range})
    hs = full.merge(hs[["date", "hs"]], on="date", how="left")

    hs["hs"] = pd.to_numeric(hs["hs"], errors="coerce").fillna(0.0)
    hs.loc[hs["hs"] < 0, "hs"] = 0.0

    # Skip almost flat winters
    if hs["hs"].max() - hs["hs"].min() < 0.01:  # < 1 cm variation
        raise RuntimeError("Block has almost no snow (flat hs); skipping.")

    with tempfile.TemporaryDirectory() as tmpdir_str:
        tmpdir = Path(tmpdir_str)
        in_csv = tmpdir / "hs.csv"
        out_csv = tmpdir / "swe.csv"

        hs[["date", "hs"]].to_csv(in_csv, index=False)

        cmd = [
            "Rscript",
            str(r_runner),
            "--in", str(in_csv),
            "--out", str(out_csv),
            "--tz", tz,
            "--rho.max", str(par.rho_max),
            "--rho.null", str(par.rho_null),
            "--c.ov", str(par.c_ov),
            "--k.ov", str(par.k_ov),
            "--k", str(par.k),
            "--tau", str(par.tau),
            "--eta.null", str(par.eta_null),
        ]

        cp = sp.run(cmd, capture_output=True, text=True)
        if cp.returncode != 0:
            raise RuntimeError(
                f"R runner failed (code {cp.returncode}).\n"
                f"STDOUT:\n{cp.stdout}\nSTDERR:\n{cp.stderr}"
            )

        if not out_csv.exists():
            raise RuntimeError("R runner did not produce output CSV.")

        res = pd.read_csv(out_csv)
        res.columns = [c.strip() for c in res.columns]

        if "date" not in res.columns or not any(c.lower() == "swe" for c in res.columns):
            raise RuntimeError("R output must contain 'date' and 'swe' columns.")

        swe_col = [c for c in res.columns if c.lower() == "swe"][0]
        res["date"] = pd.to_datetime(res["date"], errors="coerce")
        res = res.dropna(subset=["date"])
        res.rename(columns={swe_col: "swe_mod"}, inplace=True)

        return res[["date", "swe_mod"]]


# --------------------------- Metrics / evaluation ------------------------ #

def rmse_bias(joined: pd.DataFrame) -> Tuple[float, float]:
    """Compute RMSE and bias from a DataFrame with swe_obs, swe_mod."""
    obs = joined["swe_obs"].to_numpy()
    mod = joined["swe_mod"].to_numpy()
    mask = np.isfinite(obs) & np.isfinite(mod)
    if not np.any(mask):
        return BIG_PENALTY, np.nan

    err = mod[mask] - obs[mask]
    rmse = float(np.sqrt(np.mean(err ** 2)))
    bias = float(np.mean(err))
    return rmse, bias


def eval_blocks(
    blocks: List[pd.DataFrame],
    par: Params,
    r_runner: Path,
    tz: str,
    n_jobs: int,
) -> Tuple[float, float]:
    """Evaluate modeled SWE vs observed SWE over all blocks."""
    if not blocks:
        return BIG_PENALTY, np.nan

    def _eval_one(block: pd.DataFrame) -> Optional[pd.DataFrame]:
        try:
            mod = run_r_delta(block, par, r_runner, tz)
            left = block.copy()
            left["date"] = pd.to_datetime(left["date"], errors="coerce")
            mod["date"] = pd.to_datetime(mod["date"], errors="coerce")

            joined = pd.merge(left[["date", "swe_obs"]], mod, on="date", how="inner")
            joined = joined.dropna(subset=["swe_obs", "swe_mod"])
            if joined.empty:
                return None
            return joined
        except Exception as e:
            print(f"[WARN] Block evaluation failed: {e}")
            return None

    joined_all: List[pd.DataFrame] = []

    if n_jobs > 1:
        with ThreadPoolExecutor(max_workers=n_jobs) as ex:
            futures = [ex.submit(_eval_one, b) for b in blocks]
            for fut in as_completed(futures):
                res = fut.result()
                if res is not None and not res.empty:
                    joined_all.append(res)
    else:
        for b in blocks:
            res = _eval_one(b)
            if res is not None and not res.empty:
                joined_all.append(res)

    if not joined_all:
        return BIG_PENALTY, np.nan

    df_all = pd.concat(joined_all, ignore_index=True)
    return rmse_bias(df_all)


# ------------------------------- Logging --------------------------------- #

def write_log(
    log_dir: Path,
    comment: str,
    method: str,
    best_par: dict,
    rmse_fit: float,
    bias_fit: float,
    rmse_val: float,
    bias_val: float,
    n_fit_rows: int,
    station_names: List[str],
) -> Path:
    """Write a simple text log of the calibration."""
    log_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y_%m_%d_%H%M")
    safe_comment = re.sub(r"[^A-Za-z0-9_\-]+", "_", comment or "")
    fname = f"{ts}_calib_log_{safe_comment}.txt"
    log_file = log_dir / fname

    lines = []
    lines.append("===== Calibration Log =====")
    lines.append(f"Timestamp: {ts}")
    lines.append("")
    lines.append("---- Data / Stations ----")
    lines.append("Stations: " + ", ".join(sorted(station_names)))
    lines.append(f"Fit rows (approx.): {n_fit_rows}")
    lines.append("")
    lines.append("---- Optimization Result ----")
    lines.append(f"Method: {method}")
    lines.append(f"RMSE (fit) = {rmse_fit:.6g}")
    lines.append(f"Bias (fit) = {bias_fit:.6g}")
    lines.append(f"RMSE (val) = {rmse_val:.6g}")
    lines.append(f"Bias (val) = {bias_val:.6g}")
    lines.append("")
    lines.append("Best parameters (unscaled):")
    for k, v in best_par.items():
        lines.append(f"  {k} = {v:.10g}")
    lines.append("")
    lines.append("Comment:")
    lines.append(comment or "")
    lines.append("===== End of Log =====")

    log_file.write_text("\n".join(lines))
    print(f"[LOG] Written to {log_file}")
    return log_file


# ----------------------------- Script entry ------------------------------ #

if __name__ == "__main__":
    # 1) Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Δ-SNOW calibration (L-BFGS-B, ALT-winter)."
    )
    parser.add_argument("--config", required=True, help="Path to YAML config file.")
    parser.add_argument("--n-jobs", type=int, default=None, help="Number of parallel workers.")
    parser.add_argument("--maxiter", type=int, default=80, help="Max L-BFGS-B iterations.")
    parser.add_argument("--disp", action="store_true", help="Print RMSE/Bias during optimization.")
    args = parser.parse_args()

    # 2) Read YAML config directly here (no extra function)
    cfg_path = Path(args.config)
    with cfg_path.open("r") as f:
        cfg = yaml.safe_load(f)

    base_dir = Path(cfg["base_dir"])
    data_dir = Path(cfg["data_dir"])
    log_dir = Path(cfg["log_dir"])
    r_runner = Path(cfg["r_runner"])

    tz = cfg.get("timezone", "Europe/Vienna")
    comment = cfg.get("calib_comment", "")

    season_start = cfg.get("season_start", "-08-01")
    season_end = cfg.get("season_end", "-07-31")

    par_guess = cfg["par_guess"]
    par_scale = cfg["par_scale"]
    lower = cfg["lower"]
    upper = cfg["upper"]

    param_names = list(par_guess.keys())
    guess_vec = np.array([par_guess[n] for n in param_names], dtype=float)
    scale_vec = np.array([par_scale[n] for n in param_names], dtype=float)
    lower_vec = np.array([lower[n] for n in param_names], dtype=float)
    upper_vec = np.array([upper[n] for n in param_names], dtype=float)

    x0 = guess_vec / scale_vec
    bounds = list(zip(lower_vec / scale_vec, upper_vec / scale_vec))

    # 3) Number of parallel workers
    if args.n_jobs is None:
        n_jobs = max(1, (os.cpu_count() or 2) - 1)
    else:
        n_jobs = args.n_jobs
    print(f"[INFO] Using n_jobs = {n_jobs}")

    # 4) Work in base_dir
    os.chdir(base_dir)

    # 5) Load data and split into fit / validation winters
    stations = load_station_data(data_dir)
    fit_blocks, val_blocks = alt_split_fit_val(stations, season_start, season_end)

    # Debug: quick check of hs
    any_station = next(iter(stations.keys()))
    print("[DEBUG] Example station:", any_station)
    print(stations[any_station][["date", "hs"]].head())
    print(stations[any_station]["hs"].describe())

    if not fit_blocks:
        raise RuntimeError("No fit blocks found. Check season windows and winter length.")

    n_fit_rows = sum(len(b) for b in fit_blocks)
    print(f"[INFO] Number of fit blocks: {len(fit_blocks)}, rows: {n_fit_rows}")
    print(f"[INFO] Number of val blocks: {len(val_blocks)}")

    # 6) Define objective in scaled space (fit blocks only)
    def objective(x_scaled: np.ndarray) -> float:
        theta = x_scaled * scale_vec
        par = Params(*theta.tolist())
        rmse, bias = eval_blocks(fit_blocks, par, r_runner, tz, n_jobs=n_jobs)
        if args.disp:
            print(f"  RMSE(fit) = {rmse:.6g}, Bias(fit) = {bias:.6g}")
        return rmse if np.isfinite(rmse) else BIG_PENALTY

    # 7) Run optimization
    print("[INFO] Starting L-BFGS-B optimization...")
    res = minimize(
        objective,
        x0=x0,
        method="L-BFGS-B",
        bounds=bounds,
        options=dict(maxiter=args.maxiter, disp=args.disp),
    )

    if not res.success:
        print(f"[WARN] Optimizer did not converge cleanly: {res.message}")

    best_scaled = res.x
    best_real_vec = best_scaled * scale_vec
    best_par = {n: float(v) for n, v in zip(param_names, best_real_vec)}
    par_best = Params(*best_real_vec.tolist())

    # 8) Final metrics on fit and validation
    rmse_fit, bias_fit = eval_blocks(fit_blocks, par_best, r_runner, tz, n_jobs=n_jobs)
    rmse_val, bias_val = eval_blocks(val_blocks, par_best, r_runner, tz, n_jobs=n_jobs)

    # 9) Logging + printing
    write_log(
        log_dir=log_dir,
        comment=comment,
        method="L-BFGS-B",
        best_par=best_par,
        rmse_fit=rmse_fit,
        bias_fit=bias_fit,
        rmse_val=rmse_val,
        bias_val=bias_val,
        n_fit_rows=n_fit_rows,
        station_names=list(stations.keys()),
    )

    print("\nBest parameters (unscaled):")
    for k in param_names:
        print(f"  {k:10s} = {best_par[k]:.8g}")
    print(f"\nRMSE (fit) = {rmse_fit:.6g}")
    print(f"Bias (fit) = {bias_fit:.6g}")
    print(f"RMSE (val) = {rmse_val:.6g}")
    print(f"Bias (val) = {bias_val:.6g}")