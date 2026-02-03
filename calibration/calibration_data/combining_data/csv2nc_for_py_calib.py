def load_stations_to_nc(input_dir, output_nc):

    """
    Load all CSV station files and build a NetCDF with:
        dims: station x time
        vars: hs, swe_obs
    Handles missing SWE/HS automatically by aligning all data.
    """
    import pandas as pd
    import xarray as xr
    from pathlib import Path

    input_dir = Path(input_dir)
    csv_files = sorted(input_dir.glob("*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found in {input_dir}")

    df_list = []

    # ---------------------------------------------------------
    # LOAD ALL CSV FILES
    # ---------------------------------------------------------
    for f in csv_files:
        df = pd.read_csv(f)

        # station name from filename
        station = f.stem.split("_")[0]
        df["station"] = station

        df["date"] = pd.to_datetime(df["date"])
        df_list.append(df)

    df_all = pd.concat(df_list, ignore_index=True)
    print(f"✓ Loaded {len(csv_files)} station files ({len(df_all)} rows).")

    # ---------------------------------------------------------
    # UNIFIED AXES
    # ---------------------------------------------------------
    station_axis = sorted(df_all["station"].unique())
    time_axis = sorted(df_all["date"].unique())

    # Create empty dataset with unified axes
    ds = xr.Dataset(
        coords={
            "station": station_axis,
            "time": time_axis,
        }
    )

    # ---------------------------------------------------------
    # FILL VARIABLES (HS, SWE)
    # ---------------------------------------------------------
    for var in ["hs", "swe_obs"]:
        if var in df_all.columns:

            # Start with all-NaN matrix
            data = pd.DataFrame(
                index=station_axis,
                columns=time_axis,
                dtype=float
            )

            # Fill values
            for st in station_axis:
                sub = df_all[df_all["station"] == st]

                if var not in sub.columns:
                    continue

                for _, row in sub.iterrows():
                    data.loc[st, row["date"]] = row[var]

            ds[var] = (("station", "time"), data.values)

    # ---------------------------------------------------------
    # SAVE THE DATASET
    # ---------------------------------------------------------
    ds.to_netcdf(output_nc)
    print(f"Saved NetCDF → {output_nc}")


# ---------------------------------------------------------
# MAIN EXECUTION
# ---------------------------------------------------------
if __name__ == "__main__":

    ds = load_stations_to_nc(
    input_dir="/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_data/output/HS_SWE_by_station",
    output_nc="/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_data/output/merged_nc_files/HS_SWE_dataset.nc"
    )   