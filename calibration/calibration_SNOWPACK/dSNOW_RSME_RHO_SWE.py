##############################################################################
# DELTASNOW PARAMETER OPTIMIZATION FROM SMET FILES
##############################################################################

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import minimize
from joblib import Parallel, delayed
from pydeltasnow import swe_deltasnow

##############################################################################
# PATHS AND SETTINGS
##############################################################################

INPUT_DIR = Path("/Users/jakobwerkgarner/code/mt_dsnow/par_sens/SNOWPACK_data")

SEASON_START_MONTH = 8
MIN_WINTER_LENGTH = 200
EPS = 1e-6
N_JOBS = -1

##############################################################################
# VARIABLES
##############################################################################

HS_COL = "HS_meas"
SWE_COL = "SWE"

##############################################################################
# SMET FUNCTIONS
##############################################################################

def read_smet(file_path):
    header = []
    fields = None
    data_rows = []
    data_start = False

    with open(file_path, "r") as f:
        for line in f:
            line = line.rstrip()

            if line.startswith("fields"):
                fields = line.split("=", 1)[1].strip().split()

            elif line.strip() == "[DATA]":
                data_start = True
                continue

            elif not data_start:
                header.append(line)

            else:
                if line:
                    data_rows.append(line.split())

    df = pd.DataFrame(data_rows, columns=fields)
    df["timestamp"] = pd.to_datetime(df["timestamp"])

    return header, df


def extract_station_id(header):
    for line in header:
        if line.startswith("station_id"):
            return line.split("=", 1)[1].strip()
    raise ValueError("station_id not found in SMET header")


##############################################################################
# GENERAL HELPER FUNCTIONS
##############################################################################

def rmse(obs, mod):
    obs = np.asarray(obs, dtype=float)
    mod = np.asarray(mod, dtype=float)

    ok = np.isfinite(obs) & np.isfinite(mod)
    if not np.any(ok):
        return np.nan

    return np.sqrt(np.mean((obs[ok] - mod[ok]) ** 2))


def combined_rmse(df, eps=1e-6):
    rmse_swe = rmse(df["swe_obs"], df["swe_mod"])

    rho_obs = np.where(df["hs"] > eps, df["swe_obs"] / df["hs"], np.nan)
    rho_mod = np.where(df["hs"] > eps, df["swe_mod"] / df["hs"], np.nan)

    rmse_rho = rmse(rho_obs, rho_mod)

    if not np.isfinite(rmse_swe) or not np.isfinite(rmse_rho):
        return np.inf

    return 0.5 * rmse_swe + 0.5 * rmse_rho


def is_valid_winter(df):
    if len(df) < MIN_WINTER_LENGTH:
        return False

    if len(df) < 365:
        first_hs = float(df["hs"].iloc[0])
        last_hs = float(df["hs"].iloc[-1])
        if first_hs > 0.05 or last_hs > 0.05:
            return False

    return True


##############################################################################
# LOAD ALL SMET FILES
##############################################################################

def load_all_smet_data(input_dir):
    all_data = []

    for smet_file in input_dir.glob("*.smet"):
        print(f"Reading {smet_file.name}")

        header, df = read_smet(smet_file)
        station_id = extract_station_id(header)

        needed = ["timestamp", HS_COL, SWE_COL]
        missing = [col for col in needed if col not in df.columns]
        if missing:
            raise ValueError(f"{smet_file.name} is missing columns: {missing}")

        df = df[needed].copy()

        df[[HS_COL, SWE_COL]] = (
            df[[HS_COL, SWE_COL]]
            .apply(pd.to_numeric, errors="coerce")
            .replace(-999, np.nan)
        )

        df_daily = (
            df.set_index("timestamp")[[HS_COL, SWE_COL]]
            .resample("D")
            .mean()
            .reset_index()
        )

        df_daily["name"] = station_id
        df_daily = df_daily.rename(
            columns={
                "timestamp": "date",
                HS_COL: "hs",
                SWE_COL: "swe_obs"
            }
        )

        all_data.append(df_daily)

    if len(all_data) == 0:
        raise ValueError("No SMET files found.")

    out = pd.concat(all_data, ignore_index=True)
    out["date"] = pd.to_datetime(out["date"])
    out = out.sort_values(["name", "date"]).reset_index(drop=True)

    return out


##############################################################################
# PREPARE FIT DATA
##############################################################################

def prepare_fit_data(df_all, start_month=8):
    fit_parts = []

    for station, df_station in df_all.groupby("name"):
        df_station = df_station.sort_values("date").copy()
        years = np.sort(df_station["date"].dt.year.unique())

        print(f"Processing station: {station}")

        for y in years[:-1]:
            start = pd.Timestamp(year=int(y), month=8, day=1)
            end = pd.Timestamp(year=int(y) + 1, month=7, day=31)

            winter = df_station[
                (df_station["date"] >= start) &
                (df_station["date"] <= end)
            ].copy()

            if winter.empty:
                continue

            if not is_valid_winter(winter):
                print(f"  skipping winter {y}/{y+1}")
                continue

            winter["block"] = int(y)
            fit_parts.append(winter)

    if len(fit_parts) == 0:
        raise ValueError("No valid winters found for fitting.")

    fit_data = pd.concat(fit_parts, ignore_index=True)
    return fit_data


##############################################################################
# RUN DELTASNOW FOR ONE STATION/BLOCK
##############################################################################

def run_one_block(df_block, par_real):
    df_block = df_block.sort_values("date").copy()

    full_dates = pd.date_range(
        start=df_block["date"].min(),
        end=df_block["date"].max(),
        freq="D"
    )

    joined = (
        df_block[["date", "hs", "swe_obs"]]
        .set_index("date")
        .reindex(full_dates)
        .rename_axis("date")
        .reset_index()
    )

    hs_series = pd.Series(
        joined["hs"].to_numpy(dtype=float),
        index=pd.DatetimeIndex(joined["date"]),
        name="hs"
    )

    try:
        swe_mod = swe_deltasnow(
            hs_series,
            rho_max=par_real[0],
            rho_null=par_real[1],
            c_ov=par_real[2],
            k_ov=par_real[3],
            k=par_real[4],
            tau=par_real[5],
            eta_null=par_real[6],
            hs_input_unit="m",
            swe_output_unit="mm"
        )

        joined["swe_mod"] = np.asarray(swe_mod, dtype=float)

    except Exception:
        joined["swe_mod"] = np.nan

    joined = joined.dropna(subset=["hs", "swe_obs", "swe_mod"]).copy()
    return joined


##############################################################################
# OBJECTIVE FUNCTION
##############################################################################

def objective(par_scaled, data, par_scale, verbose=False):
    par_real = np.asarray(par_scaled) * np.asarray(par_scale)

    groups = [g for _, g in data.groupby(["name", "block"], sort=False)]

    results = Parallel(n_jobs=N_JOBS)(
        delayed(run_one_block)(g, par_real) for g in groups
    )

    results = [r for r in results if r is not None and not r.empty]

    if len(results) == 0:
        return 1e12

    dff = pd.concat(results, ignore_index=True)

    score = combined_rmse(dff, eps=EPS)

    if not np.isfinite(score):
        return 1e12

    if verbose:
        rmse_swe = rmse(dff["swe_obs"], dff["swe_mod"])
        rho_obs = np.where(dff["hs"] > EPS, dff["swe_obs"] / dff["hs"], np.nan)
        rho_mod = np.where(dff["hs"] > EPS, dff["swe_mod"] / dff["hs"], np.nan)
        rmse_rho = rmse(rho_obs, rho_mod)
        bias_swe = np.nanmean(dff["swe_mod"] - dff["swe_obs"])

        print(
            f"| bias_swe = {bias_swe:.4f} "
            f"| rmse_swe = {rmse_swe:.4f} "
            f"| rmse_density = {rmse_rho:.4f} "
            f"| final_score = {score:.4f}"
        )

    return float(score)


##############################################################################
# START VALUES
##############################################################################

par_delta = np.array([
    396.9521,      # rho_max
    78.88324,      # rho_null
    0.0004986025,  # c_ov
    0.227786,      # k_ov
    0.0290256,     # k
    0.02489556,    # tau
    8792253        # eta_null
], dtype=float)

par_scale = np.array([
    1000,
    1000,
    0.001,
    1,
    0.1,
    0.1,
    1e7
], dtype=float)

##############################################################################
# BOUNDS (UNSCALED)
##############################################################################

lower_unscaled = np.array([
    400,    # rho_max
    100,    # rho_null
    0.0,    # c_ov
    0.05,   # k_ov
    0.01,   # k
    0.01,   # tau
    5e6     # eta_null
], dtype=float)

upper_unscaled = np.array([
    550,    # rho_max
    120,    # rho_null
    0.001,  # c_ov
    2.0,    # k_ov
    0.08,   # k
    0.05,   # tau
    9e6     # eta_null
], dtype=float)

##############################################################################
# SCALE START VALUES AND BOUNDS
##############################################################################

lower_scaled = lower_unscaled / par_scale
upper_scaled = upper_unscaled / par_scale

# force initial parameters into bounds
par_start = np.clip(par_delta, lower_unscaled, upper_unscaled) / par_scale

bounds = list(zip(lower_scaled, upper_scaled))

print("\nScaled start parameters:")
print(par_start)

print("\nUnscaled start parameters:")
print(par_start * par_scale)

print("\nBounds:")
for name, lo, hi in zip(
    ["rho_max", "rho_null", "c_ov", "k_ov", "k", "tau", "eta_null"],
    lower_unscaled,
    upper_unscaled
):
    print(f"{name:10s}: [{lo}, {hi}]")

##############################################################################
# CALLBACK TO PRINT OPTIMIZATION PROGRESS
##############################################################################

iteration = 0

def print_progress(xk):
    global iteration
    iteration += 1

    par_real = xk * par_scale
    score = objective(xk, fit_data, par_scale, verbose=False)

    print(
        f"\nITER {iteration}"
        f"\nparams = {np.round(par_real, 6)}"
        f"\nscore  = {score:.6f}"
    )

##############################################################################
# MAIN
##############################################################################

df_all = load_all_smet_data(INPUT_DIR)
fit_data = prepare_fit_data(df_all, start_month=SEASON_START_MONTH)

print("\nTesting objective function ...")
test_score = objective(par_start, fit_data, par_scale, verbose=True)
print("Initial score:", test_score)

print("\nStarting optimization ...")

res = minimize(
    fun=objective,
    x0=par_start,
    args=(fit_data, par_scale, False),
    method="L-BFGS-B",
    bounds=bounds,
    callback=print_progress,
    options={
        "maxiter": 500,
        "disp": True
    }
)

##############################################################################
# RESULTS
##############################################################################

best_par_scaled = res.x
best_par_real = best_par_scaled * par_scale

param_names = ["rho_max", "rho_null", "c_ov", "k_ov", "k", "tau", "eta_null"]
best_params = pd.Series(best_par_real, index=param_names, name="best_value")

print("\nBest parameters:")
print(best_params)

print("\nFinal score:")
print(res.fun)

print("\nOptimizer message:")
print(res.message)