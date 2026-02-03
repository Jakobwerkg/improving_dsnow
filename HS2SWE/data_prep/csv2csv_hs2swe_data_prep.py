#!/usr/bin/env python3
"""
save_data_for_HS2SWE.py

Convert calibration CSV files into the format required by HS2SWE:

Required output format in HS2SWE/data_prep/:
    date,hs,swe_obs

This script:
 - Reads all CSV files from the SLF MAG25 output directory
 - Ensures correct column names and ordering
 - Most imortatly converts from m to cm
 - Writes cleaned files into HS2SWE/data_prep/

Author: Jakob Werkgarner
"""

import os
import pandas as pd

# ----------------------------
# Paths
# ----------------------------
input_dir = "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_data/output/HS_SWE_by_station"
output_dir = "/Users/jakobwerkgarner/code/mt_dsnow/HS2SWE/data_prep/hs2swe_input"
os.makedirs(output_dir, exist_ok=True)

# ----------------------------
# Process all station CSVs
# ----------------------------
for file in os.listdir(input_dir):
    if not file.endswith(".csv"):
        continue

    path_in = os.path.join(input_dir, file)
    df = pd.read_csv(path_in)

    # Expected columns from previous script:
    # date, hs, swe_obs
    # If station column exists, drop it
    df = df.rename(columns={
        "HS": "hs_obs",
        "SWE": "swe_obs",
        "HS_sta": "hs",
        "SWE_obs": "swe_obs",
        "time": "date"
    })

    # Enforce required columns
    required = ["date", "hs", "swe_obs"]
    for col in required:
        if col not in df.columns:
            df[col] = None

    df = df[required]


    #Converte to cm since the hs2swe requiers swe in cm
    df['hs'] = df['hs'] * 100


    # Save cleaned file
    station_name = file.replace("_hs_swe_obs", "").replace(".csv", "")
    out_path = os.path.join(output_dir, f"{station_name}.csv")
    df.to_csv(out_path, index=False)


    print(f"Saved: {out_path}")

print("\n=== All stations prepared for HS2SWE ===")