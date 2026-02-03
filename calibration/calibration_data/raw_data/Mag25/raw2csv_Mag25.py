#!/usr/bin/env python3
"""
raw2csv_Mag25.py

Prepare SLF Mag25 dataset:
- Load SLF text files (HN, HNW, HS, SWE)
- Convert them to aligned xarray DataArrays
- Apply canonical station renaming
- Save full dataset to NetCDF (Mag25_all.nc)
- Export per-station CSV files for calibration with columns:
        date, hs, swe_obs

IMPORTANT:
SWE is kept in mm
HS is converted from cm → m

All logic follows exactly the original notebook.

Author: Jakob Werkgarner
"""

# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
import numpy as np
import pandas as pd
import os
import xarray as xr


# ---------------------------------------------------------
# Helper function
# ---------------------------------------------------------
def df_to_da(df, var_name, station_names):
    """
    Convert SLF wide-format dataframe to xarray DataArray aligned by station.
    """
    df = df.set_index("Time")
    df.index.name = "time"
    df = df.reindex(columns=station_names, fill_value=np.nan)

    arr = df.to_numpy(dtype=np.float64)
    time = pd.to_datetime(df.index)

    return xr.DataArray(
        arr,
        dims=("time", "station"),
        coords={"time": time, "station": station_names},
        name=var_name,
    )


# ---------------------------------------------------------
# Main processing script
# ---------------------------------------------------------
def main():

    print("\n=== Preparing Mag25 dataset ===\n")

    # -----------------------------------------------------
    # Paths
    # -----------------------------------------------------
    base_dir = "/Users/jakobwerkgarner/code/mt_dsnow"
    input_data_dir = "calibration/calibration_data/raw_data/Mag25/SLF_dataset"

    os.chdir(base_dir)

    # -----------------------------------------------------
    # Import SLF helper
    # -----------------------------------------------------
    import calibration.calibration_data.raw_data.Mag25.helpers_MAG25 as hMag25

    # -----------------------------------------------------
    # Load SLF data
    # -----------------------------------------------------
    print("Loading raw SLF data...")

    HN = hMag25.load_SLF_data(os.path.join(input_data_dir, "OBS-HN.txt"))
    HNW = hMag25.load_SLF_data(os.path.join(input_data_dir, "OBS-HNW.txt"))
    HS_man = hMag25.load_SLF_data(os.path.join(input_data_dir, "OBS-HS-PROFILE.txt"))
    HS_sta = hMag25.load_SLF_data(os.path.join(input_data_dir, "OBS-HS-STAKE.txt"))
    SWE = hMag25.load_SLF_data(os.path.join(input_data_dir, "OBS-SWE-PROFILE.txt"))

    station_list = pd.read_csv(
        os.path.join(input_data_dir, "meta/STATION-LIST.txt"), sep=","
    )

    station_names = station_list["Station"].values

    # -----------------------------------------------------
    # Coordinates
    # -----------------------------------------------------
    northing = xr.DataArray(
        station_list["Northing (m)"].values,
        dims="station",
        coords={"station": station_names},
        name="northing",
    )

    easting = xr.DataArray(
        station_list["Easting (m)"].values,
        dims="station",
        coords={"station": station_names},
        name="easting",
    )

    altitude = xr.DataArray(
        station_list["Altitude (m.a.s.l.)"].values,
        dims="station",
        coords={"station": station_names},
        name="altitude",
    )

    # -----------------------------------------------------
    # Convert data tables to DataArrays
    # -----------------------------------------------------
    HN_da = df_to_da(HN, "HN", station_names)
    HNW_da = df_to_da(HNW, "HNW", station_names)
    HS_da = df_to_da(HS_sta, "HS", station_names)
    SWE_da = df_to_da(SWE, "SWE", station_names)

    # Convert cm → m
    HN_da_m = HN_da / 100.0
    HS_da_m = HS_da / 100.0

    # -----------------------------------------------------
    # Build unified dataset
    # -----------------------------------------------------
    ds = xr.Dataset(
        {
            "HN": HN_da_m,
            "HNW": HNW_da,
            "HS": HS_da_m,
            "SWE": SWE_da,
            "northing": northing,
            "easting": easting,
            "altitude": altitude,
        }
    )

    # -----------------------------------------------------
    # Canonical station renaming
    # -----------------------------------------------------
    print("Applying canonical station names...")

    station_rename = {
        'SLF.1AD': 'Adelboden',
        'SLF.1GA': 'Gadmen',
        'SLF.1GB': 'Grindelwald_Bort',
        'SLF.1GS': 'Gsteig',
        'SLF.1GT': 'Gantrisch',
        'SLF.1LS': 'Leysin',
        'SLF.1MR': 'Muerren',
        'SLF.1SM': 'Saanenmoeser',
        'SLF.1WE': 'Wengen',
        'SLF.2SO': 'Srenberg',
        'SLF.2ST': 'Stoos',
        'SLF.3BR': 'Braunwald',
        'SLF.3MB': 'Malbun',
        'SLF.3MG': 'St_Margrethenberg',
        'SLF.4BN': 'Binn',
        'SLF.4BP': 'Bourg_St_Pierre',
        'SLF.4FY': 'Fionnay',
        'SLF.4GR': 'Grimentz',
        'SLF.4LA': 'Lauchernalp',
        'SLF.4MO': 'Montana',
        'SLF.4MS': 'Muenster',
        'SLF.4SF': 'Saas_Fee',
        'SLF.4SM': 'Simplon_Dorf',
        'SLF.4UL': 'Ulrichen',
        'SLF.4WI': 'Wiler',
        'SLF.5BI': 'Bivio',
        'SLF.5DF': 'Davos_Flueelastr',
        'SLF.5JU': 'Juf',
        'SLF.5OB': 'Obersaxen',
        'SLF.5PU': 'Pusserein',
        'SLF.5SA': 'St_Antoenien',
        'SLF.5SE': 'Sedrun',
        'SLF.5SP': 'Spluegen',
        'SLF.5VA': 'Vals',
        'SLF.5WJ': 'Weisfluh_Joch',
        'SLF.6BG': 'Bosco_Gurin',
        'SLF.6SB': 'San_Bernadino',
        'SLF.7MA': 'Maloja',
        'SLF.7MZ': 'Sankt_Moritz',
        'SLF.7SN': 'Samnaun',
        'SLF.7ZU': 'Zuoz'
    }

    labels = ds["station"].values
    missing = sorted(set(labels) - set(station_rename.keys()))
    if missing:
        raise ValueError(f"Unmapped station codes: {missing}")

    mapped = np.array([station_rename[k] for k in labels], dtype=object)
    ds = ds.assign_coords(station=("station", mapped))

    # -----------------------------------------------------
    # Add metadata
    # -----------------------------------------------------
    ds["HN"].attrs = {"units": "m"}
    ds["HNW"].attrs = {"units": "mm"}
    ds["HS"].attrs = {"units": "m"}
    ds["SWE"].attrs = {"units": "mm"}

    # -----------------------------------------------------
    # Save NetCDF
    # -----------------------------------------------------
    out_nc = os.path.join(input_data_dir, "Mag25_all.nc")
    ds.to_netcdf(out_nc)
    print(f"\nSaved NetCDF → {out_nc}\n")

    # -----------------------------------------------------
    # Export per-station CSV for calibration
    # -----------------------------------------------------
    print("Exporting per-station CSV files...")

    ds_reduced = ds.drop_vars(["HN", "HNW", "northing", "easting", "altitude"])
    out_dir = "calibration/calibration_data/output/HS_SWE_by_station"
    os.makedirs(out_dir, exist_ok=True)

    for station in ds_reduced["station"].values:

        df = (
            ds_reduced[["HS", "SWE"]]
            .sel(station=station)
            .to_dataframe()
            .reset_index()
        )

        df = df.rename(columns={
            "HS": "hs",
            "SWE": "swe_obs",
            "time": "date",
        })

        df["date"] = df["date"].dt.date
        df = df.drop(columns=["station"])

        safe = "".join(c if c.isalnum() else "_" for c in station)
        out_file = os.path.join(out_dir, f"{safe}_hs_swe_obs.csv")

        df.to_csv(out_file, index=False)
        print(f"Saved: {out_file}")

    print("\n=== DONE: Mag25 dataset ready. ===\n")


# ---------------------------------------------------------
if __name__ == "__main__":
    main()