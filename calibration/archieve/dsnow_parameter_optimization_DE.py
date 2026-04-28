##############################################################################
# DELTASNOW PARAMETER OPTIMIZATION FROM SMET FILES - DIFFERENTIAL EVOLUTION
# with FutureWarning suppression inside parallel workers
##############################################################################

import sys
import contextlib
import io
import warnings

# Suppress FutureWarnings globally (for the main process)
warnings.filterwarnings("ignore", category=FutureWarning)

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import differential_evolution
from joblib import Parallel, delayed

sys.path.insert(0, "/Users/jakobwerkgarner/code/mt_dsnow/snow_to_swe_master")
from main import SnowToSwe

##############################################################################
# PATHS AND SETTINGS
##############################################################################

INPUT_DIR = Path("/Users/jakobwerkgarner/code/mt_dsnow/par_sens/SNOWPACK_data")

SEASON_START_MONTH = 8
MIN_WINTER_LENGTH = 200
EPS = 1e-6
N_JOBS = -1

##############################################################################
# OBJECTIVE WEIGHTS  (should sum to 1)
##############################################################################

WEIGHT_SWE_NRMSE = 0.5   # normalized RMSE of SWE
WEIGHT_RHO_NRMSE = 0.5   # normalized RMSE of bulk density
WEIGHT_SWE_NBIAS = 0.0   # normalized absolute bias of SWE

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


def combined_score(df, eps=1e-6):
    swe_obs = df["swe_obs"].to_numpy()
    swe_mod = df["swe_mod"].to_numpy()
    hs      = df["hs"].to_numpy()

    # Bulk density (only where snow is present)
    rho_obs = np.where(hs > eps, swe_obs / hs, np.nan)
    rho_mod = np.where(hs > eps, swe_mod / hs, np.nan)

    rmse_swe = rmse(swe_obs, swe_mod)
    rmse_rho = rmse(rho_obs, rho_mod)

    mean_swe = np.nanmean(swe_obs)
    mean_rho = np.nanmean(rho_obs)

    ok = np.isfinite(swe_obs) & np.isfinite(swe_mod)
    bias_swe = np.mean(swe_mod[ok] - swe_obs[ok]) if np.any(ok) else np.nan

    nrmse_swe = rmse_swe / mean_swe  if (np.isfinite(rmse_swe) and mean_swe > 0) else np.inf
    nrmse_rho = rmse_rho / mean_rho  if (np.isfinite(rmse_rho) and mean_rho > 0) else np.inf
    nbias_swe = abs(bias_swe) / mean_swe if (np.isfinite(bias_swe) and mean_swe > 0) else np.inf

    score = (
        WEIGHT_SWE_NRMSE * nrmse_swe +
        WEIGHT_RHO_NRMSE * nrmse_rho +
        WEIGHT_SWE_NBIAS * nbias_swe
    )
    return score if np.isfinite(score) else np.inf


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
# RUN DELTASNOW FOR ONE STATION/BLOCK (with warning suppression)
##############################################################################

def run_one_block(df_block, par_real):
    # Suppress FutureWarnings inside this worker process
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)

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

        # SnowToSwe requires no NaNs and first value == 0
        hs_series = (
            joined["hs"]
            .interpolate(method="linear", limit_direction="both")
            .clip(lower=0.0)
        )
        if hs_series.isna().any():
            return pd.DataFrame()  # skip block if interpolation still leaves NaN
        hs_arr = hs_series.to_numpy(dtype=float)
        hs_arr[0] = 0.0

        try:
            with contextlib.redirect_stdout(io.StringIO()):
                model = SnowToSwe(
                    rho_max=par_real[0],
                    rho_null=par_real[1],
                    c_ov=par_real[2],
                    k_ov=par_real[3],
                    k=par_real[4],
                    tau=par_real[5],
                    eta_null=par_real[6],
                )
            swe_mod = model.convert_list(hs_arr.tolist(), timestep=24)
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

    score = combined_score(dff, eps=EPS)

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
# START VALUES AND BOUNDS (PHYSICAL, THEN SCALED)
##############################################################################

par_delta = np.array([
    401,
    81,
    5.1e-4,
    0.38,
    0.03,
    0.024,
    8.5e6,      
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

# Physical (unscaled) lower and upper bounds
lower_unscaled = np.array([200,  20,   1e-6, 0,   0,   0,   1e5])
upper_unscaled = np.array([600, 200,   0.01, 1,   0.5, 1,   1e8])

# Scale bounds
lower_scaled = lower_unscaled / par_scale
upper_scaled = upper_unscaled / par_scale

# Scaled start values
par_start = par_delta / par_scale

print("\nScaled start parameters:")
print(par_start)

print("\nUnscaled start parameters:")
print(par_start * par_scale)

print("\nScaled bounds (lower, upper):")
for i, (l, u) in enumerate(zip(lower_scaled, upper_scaled)):
    print(f"  param_{i}: {l:.6f} -> {u:.6f}")

##############################################################################
# MAIN
##############################################################################

df_all = load_all_smet_data(INPUT_DIR)
fit_data = prepare_fit_data(df_all, start_month=SEASON_START_MONTH)

print("\nTesting objective function ...")
test_score = objective(par_start, fit_data, par_scale, verbose=True)
print("Initial score:", test_score)

print("\nStarting Differential Evolution optimization ...")

# DE control parameters
popsize = 10 * len(par_start)   # population size
maxiter = 200                   # maximum number of generations
mutation = (0.5, 1.0)           # mutation factor (dithering)
recombination = 0.7             # crossover probability

# Run differential evolution
result = differential_evolution(
    func=objective,
    bounds=list(zip(lower_scaled, upper_scaled)),
    args=(fit_data, par_scale, False),
    maxiter=maxiter,
    popsize=popsize,
    mutation=mutation,
    recombination=recombination,
    seed=123,                   # for reproducibility
    disp=True,                  # print progress every generation
    workers=1                   # DE internal parallelisation is off (we use joblib inside objective)
)

##############################################################################
# RESULTS
##############################################################################

best_par_scaled = result.x
best_par_real = best_par_scaled * par_scale

param_names = ["rho_max", "rho_null", "c_ov", "k_ov", "k", "tau", "eta_null"]
best_params = pd.Series(best_par_real, index=param_names, name="best_value")

print("\n" + "="*60)
print("OPTIMIZATION FINISHED")
print("="*60)
print("\nBest parameters (unscaled):")
print(best_params)

print("\nFinal score (objective value):")
print(result.fun)

print("\nOptimizer message:")
print(result.message)

print("\nNumber of function evaluations:")
print(result.nfev)

# Optional: print a Python call example with the optimized parameters
fmt = lambda x: f"{x:.10f}"
print("\n# Example call to SnowToSwe with optimized parameters:")
print(f"model = SnowToSwe(")
print(f"    rho_max={fmt(best_par_real[0])},")
print(f"    rho_null={fmt(best_par_real[1])},")
print(f"    c_ov={fmt(best_par_real[2])},")
print(f"    k_ov={fmt(best_par_real[3])},")
print(f"    k={fmt(best_par_real[4])},")
print(f"    tau={fmt(best_par_real[5])},")
print(f"    eta_null={fmt(best_par_real[6])}")
print(f")")