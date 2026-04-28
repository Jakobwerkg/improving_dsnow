import sys
import os
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# ── Project root ─────────────────────────────────────────────────────────────
base_dir = "/Users/jakobwerkgarner/code/mt_dsnow"
os.chdir(base_dir)

# ── Local imports ────────────────────────────────────────────────────────────
# Validation helpers (density-scatter plots + metrics)
import HNW_validation_helper as val_helper

# SnowToSwe – Python port of the R nixmass::swe.delta.snow (Winkler et al. 2021)
sys.path.insert(0, os.path.join(base_dir, "snow_to_swe_master"))
from main import SnowToSwe

# SnowToSwe uses a deprecated pandas idiom that emits FutureWarnings.
# Values are correct; we just silence the noise.
warnings.filterwarnings("ignore", category=FutureWarning, module=r"main")

# ── Run configuration ────────────────────────────────────────────────────────
infile        = os.path.join(base_dir, "HNW_validation/validation_input_Mag25/Mag25_all.nc")
save_dir      = "HNW_validation/dSnow/validation_plots"
model_source  = "SnowToSwe"
calib_comment = "SnowToSwe_default_Winkler2021"
save_data     = False

# ── Winkler et al. 2021 original parameters (SnowToSwe defaults) ─────────────
snow_to_swe = SnowToSwe(
    rho_max  = 401,
    rho_null = 81,
    c_ov     = 5.1e-4,
    k_ov     = 0.38,
    k        = 0.03,
    tau      = 0.024,
    eta_null = 8.5e6,
)

print(f"model_source : {model_source}")
print(f"infile       : {infile}")
print(f"calib_comment: {calib_comment}")

from joblib import Parallel, delayed
from tqdm.auto import tqdm

# ── Load Mag25 multi-station NetCDF (FULL YEAR — do NOT pre-filter) ──────────
# Running SnowToSwe per-winter requires each call to end on a snow-free day.
# If we cut at Apr 30, stations with snow still on the ground that day crash
# SnowToSwe's internal loop. Using hydrological years Sep 1 -> Aug 31 lets
# the summer naturally bring HS to 0 at the boundaries.
Mag25_data = xr.open_dataset(infile)

Mag25_data_with_SWE = Mag25_data.copy()
Mag25_data_with_SWE["SWE_mod"] = xr.full_like(Mag25_data_with_SWE["HS"], np.nan)

# ── Group the time axis into hydrological years ──────────────────────────────
# Hydrological year Y = Sep 1 (Y) through Aug 31 (Y+1).
# Label each day: Sep–Dec -> Y, Jan–Aug -> Y-1.
times        = pd.to_datetime(Mag25_data["time"].values)
hyd_year     = np.where(times.month >= 9, times.year, times.year - 1)
winter_years = np.unique(hyd_year)


# ── Per-(station × hydrological year) worker ─────────────────────────────────
def process_station_winter(station_name, winter_year,
                           time_idx, hs_values, snow_to_swe):
    """Run SnowToSwe on ONE station × ONE hydrological year (Sep–Aug)."""
    import warnings
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    try:
        hs = pd.Series(hs_values).fillna(0).clip(lower=0).astype(float)
        if len(hs) == 0:
            return None
        # SnowToSwe preconditions: series must start at 0
        if hs.iloc[0] != 0:
            hs.iloc[0] = 0.0

        swe_list = snow_to_swe.convert_list(hs.tolist(),
                                            timestep=24, verbose=False)
        if swe_list is None:
            return None

        swe_arr = np.asarray(swe_list, dtype=float)
        assert swe_arr.shape[0] == len(hs)
        return (station_name, winter_year, time_idx, swe_arr)
    except Exception as e:
        return ("__ERROR__", station_name,
                f"winter {winter_year}/{winter_year + 1}: {e}")


# ── Build (station × winter) task list ───────────────────────────────────────
station_list = Mag25_data["station"].values
hs_by_stn    = {s: Mag25_data["HS"].sel(station=s).values for s in station_list}

tasks = []
for stn in station_list:
    hs_full = hs_by_stn[stn]
    for y in winter_years:
        mask     = hyd_year == y
        time_idx = times[mask].values
        hs_vals  = hs_full[mask]
        tasks.append((stn, int(y), time_idx, hs_vals))

n_jobs = -1
print(f"Running SnowToSwe on {len(station_list)} stations × "
      f"{len(winter_years)} hydrological years = {len(tasks)} tasks "
      f"(parallel, n_jobs={n_jobs})…")

results = []
with tqdm(total=len(tasks), desc="Station×Winter", unit="task") as pbar:
    for r in Parallel(n_jobs=n_jobs, return_as="generator")(
        delayed(process_station_winter)(stn, y, tidx, hsv, snow_to_swe)
        for (stn, y, tidx, hsv) in tasks
    ):
        results.append(r)
        pbar.update(1)


# ── Collect results back into SWE_mod ────────────────────────────────────────
errors      = []
ok_tasks    = 0
ok_stations = set()

for result in results:
    if result is None:
        continue
    if isinstance(result, tuple) and result and result[0] == "__ERROR__":
        errors.append((result[1], result[2]))
        continue

    station_name, wyr, time_idx, swe_arr = result
    Mag25_data_with_SWE["SWE_mod"].loc[
        dict(station=station_name, time=time_idx)
    ] = swe_arr
    ok_tasks += 1
    ok_stations.add(station_name)


# ── Observation counts per station ───────────────────────────────────────────
obs_counts_df = pd.DataFrame([
    {"station": s,
     "n_obs": int(Mag25_data["SWE"].sel(station=s).notnull().sum().item())}
    for s in station_list
])

print(f"\nDone. Tasks processed: {ok_tasks} / {len(tasks)}")
print(f"Stations with ≥1 successful winter: "
      f"{len(ok_stations)} / {len(station_list)}")
print(f"Hydrological years: {list(winter_years)} ({len(winter_years)} seasons)")
print(f"Total time steps (full year): {Mag25_data.sizes['time']}")

if errors:
    print(f"\nErrors on {len(errors)} (station, winter) pair(s):")
    for stn, msg in errors:
        print(f"  {stn}: {msg}")
# ── Derive HNW_mod from day-to-day SWE_mod increments ────────────────────────
# Negative diffs represent melt, not new snow, so clip them to 0.
# .reindex() restores the full time axis (diff drops the first timestep).

HNW_mod = Mag25_data_with_SWE["SWE_mod"].diff(dim="time").clip(min=0)
HNW_mod = HNW_mod.reindex(time=Mag25_data_with_SWE["time"])
Mag25_data_with_SWE["HNW_mod"] = HNW_mod

# ── Build one long-format frame with both pairs (obs & mod) ─────────────────
all_df = (Mag25_data_with_SWE[["HNW", "HNW_mod", "SWE", "SWE_mod"]]
          .to_dataframe()
          .reset_index()
          .rename(columns={"HNW": "HNW_obs", "SWE": "SWE_obs"}))

# ── Keep only days with a valid observed SWE ────────────────────────────────
# Mag25 SWE is biweekly snow-course data (NaN on most days). Validating on the
# observed days only yields statistically meaningful RMSE / R², and puts the
# HNW and SWE comparisons onto the same subset of rows.

all_df_valid_SWE_obs = all_df#[all_df["SWE_obs"].notna()].copy()
#all_df_valid_SWE_obs.index = pd.to_datetime(all_df_valid_SWE_obs["time"]).values


all_df_valid_HNW_obs = all_df[all_df["HNW_obs"].notna()].copy()
all_df_valid_HNW_obs.index = pd.to_datetime(all_df_valid_HNW_obs["time"]).values



print(f"Rows total      : {len(all_df):>7}")
print(f"Rows w/ SWE_obs : {len(all_df_valid_SWE_obs):>7}")
print(f"Rows w/ HNW_obs : {len(all_df_valid_HNW_obs):>7}")


## Validation — HNW (daily new snow water equivalent)


val_helper.validate_hnw_mag25(
    all_df_valid_HNW_obs,
    model_name = "dSnow (SnowToSwe)",
    obs_col    = "HNW_obs",
    mod_col    = "HNW_mod",
    save_dir   = save_dir,
    filename   = f"SnowToSwe_hnw_validation_{calib_comment}.png",
)

## Validation — SWE (biweekly snow-course observations)


val_helper.validate_swe_mag25(
    all_df_valid_SWE_obs,
    model_name = "dSnow (SnowToSwe)",
    obs_col    = "SWE_obs",
    mod_col    = "SWE_mod",
    save_dir   = save_dir,
    filename   = f"SnowToSwe_SWE_validation_{calib_comment}.png",
)
