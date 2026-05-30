# DeltaSnow Parameter Calibration

Calibrates the 7 free parameters of the **deltasnow** model (`nixmass::swe.delta.snow`) against observed HS and SWE data from two independent datasets. The optimizer minimises a weighted combination of normalised SWE error, density error, and SWE bias.

---

## Pipeline Overview

```
Raw data                Data prep              Calibration              Results
─────────────────────   ────────────────────   ──────────────────────   ─────────────────────
SNOWPACK .smet files ──▶ prepare_SNOWPACK_      run_all_phases.sh
                         data.R                  └─ run_calibration.sh
                          │                           ├─ SNOWPACK NM ──▶ R_opt_logs/
                          ▼                           ├─ SNOWPACK DE ──▶ R_opt_logs_DE/
                     d_obs_SNOWPACK.rda               ├─ Win21 NM    ──▶ R_opt_logs/
                                                      └─ Win21 DE    ──▶ R_opt_logs_DE/
SLF Mag25 txt files ───▶ raw2csv_Mag25.py                                      │
Win21 H_SWE_obs.Rda ───▶ raw2csv_Win21.R       ◀────────────────────────────────┘
                          │                    collect_opt_results.R
                          ▼                      └─▶ opt_results_summary_grid_search.csv
                  HS_SWE_by_station/*.csv
                          │
                          ▼
                  csv2rda_calib.R
                          │
                          ▼
                   d_obs_WIN_MAG.rda
                   (used by Win21 scripts)
```

---

## File Descriptions

### Shell scripts (entry points)

| File | Purpose |
|------|---------|
| `run_all_phases.sh` | Top-level script. Runs a sequence of weight combinations (Phases 1–5). Each phase calls `run_calibration.sh`. Comment/uncomment phases to control what runs. |
| `run_calibration.sh` | Runs all 4 optimisation scripts sequentially for one weight triple `(w_swe, w_rho, w_bias)`. Usage: `./run_calibration.sh 0.4 0.6 0.0` |

### Data preparation

| File | Input | Output |
|------|-------|--------|
| `calibration_SNOWPACK/prepare_SNOWPACK_data.R` | `.smet` files from SNOWPACK simulation output (`par_sens/SNOWPACK_data/data_rain_gauge/raw_alpsolut/`) | `calibration_SNOWPACK/data/d_obs_SNOWPACK.rda` — named list of zoo objects, one per station, with columns `Hobs` (m) and `SWEobs` (mm) |
| `calibration_data/raw_data/Mag25/raw2csv_Mag25.py` | SLF raw text files: `OBS-HN.txt`, `OBS-HNW.txt`, `OBS-HS-STAKE.txt`, `OBS-SWE-PROFILE.txt`, station list | Per-station CSVs in `calibration_data/output/HS_SWE_by_station/` with columns `date, hs [m], swe_obs [mm]`. Also saves `Mag25_all.nc`. |
| `calibration_data/raw_data/dsnow/raw2csv_Win21.R` | `calibration_Win21/data/H_SWE_obs.Rda` (Win et al. 2021 dataset) | Per-station CSVs in `calibration_data/output/HS_SWE_by_station/` (same format as above) |
| `calibration_data/combining_data/csv2rda_calib.R` | Per-station CSVs from `HS_SWE_by_station/` | `calibration_data/output/calibration_rda_files/d_obs_WIN_MAG.rda` — named list of zoo objects used by the Win21 calibration scripts |

### Calibration scripts (4 optimisers)

All four scripts share the same logic — they differ only in dataset and optimisation algorithm.

| File | Dataset | Algorithm | Output dir |
|------|---------|-----------|-----------|
| `calibration_SNOWPACK/dsnow_parameter_optimization.R` | `d_obs_SNOWPACK.rda` | Nelder-Mead (`optimx`) | `calibration_SNOWPACK/data/R_opt_logs/` |
| `calibration_SNOWPACK/dsnow_parameter_optimization_DE.R` | `d_obs_SNOWPACK.rda` | Differential Evolution (`DEoptim`) | `calibration_SNOWPACK/data/R_opt_logs_DE/` |
| `calibration_Win21/WIN21_dsnow_paprameter_optimization.R` | `H_SWE_obs.Rda` | Nelder-Mead (`optimx`) | `calibration_Win21/data/R_opt_logs/` |
| `calibration_Win21/WIN21_dsnow_parameter_opitmization_DE.R` | `H_SWE_obs.Rda` | Differential Evolution (`DEoptim`) | `calibration_Win21/data/R_opt_logs_DE/` |

Each script saves one `.rds` file per run, named:
```
opt_results__SWE_NRMSE_<w>__RHO_NRMSE_<w>__SWE_NBIAS_<w>.rds
```

### Result collection

| File | Purpose |
|------|---------|
| `collect_opt_results.R` | Scans all `opt_results*.rds` files, extracts best parameters and scores, writes `opt_results_summary_grid_search.csv` and `.rds`. Run from the repo root: `Rscript calibration/collect_opt_results.R` |

---

## Objective Function

All scripts minimise the same weighted score:

```
score = w1 · NRMSE_SWE + w2 · NRMSE_rho + w3 · NBIAS_SWE

NRMSE_SWE = RMSE(SWE_mod, SWE_obs) / mean(SWE_obs)
NRMSE_rho = RMSE(rho_mod, rho_obs) / mean(rho_obs)     # bulk density = SWE / HS
NBIAS_SWE = |mean(SWE_mod - SWE_obs)| / mean(SWE_obs)
```

Weights are passed via command line and must sum to 1. The filename encodes them (e.g. `SWE_NRMSE_0p4__RHO_NRMSE_0p6__SWE_NBIAS_0p0`).

---

## Parameters Being Optimised

| Parameter | Physical meaning | Bounds |
|-----------|-----------------|--------|
| `rho.max` | Maximum bulk snow density [kg/m³] | 300 – 600 |
| `rho.null` | Fresh snow density [kg/m³] | 60 – 150 |
| `c.ov` | Overburden compaction coefficient | 1e-6 – 0.01 |
| `k.ov` | Overburden compaction exponent | 0.01 – 1.0 |
| `k` | Settling rate | 0.001 – 0.1 |
| `tau` | Time constant | 0.001 – 0.1 |
| `eta.null` | Snow viscosity [Pa·s] | 1e5 – 1e8 |

---

## Train / Validation Split

Per station, winters are split **alternating by season**:
- Even winters (counter 0, 2, 4, …) → **fit set**
- Odd winters (counter 1, 3, 5, …) → **validation set**

Winters shorter than 200 days or with snow remaining at the season boundary (incomplete record) are skipped entirely. The hydrological year runs **1 Aug → 31 Jul**.

---

## Running the Pipeline

### Full Phase 5 run (currently active)

```bash
cd ~/code/mt_dsnow/calibration
./run_all_phases.sh
```

This runs the 4 weight combinations defined in Phase 5 of `run_all_phases.sh`:

| Phase | w_SWE | w_RHO | w_BIAS |
|-------|-------|-------|--------|
| 5A | 0.8 | 0.1 | 0.1 |
| 5B | 0.1 | 0.8 | 0.1 |
| 5C | 0.4 | 0.4 | 0.2 |
| 5D | 0.25 | 0.25 | 0.5 |

Each combination runs all 4 optimisers → **4 × 4 = 16 `.rds` files** produced.

### Single weight combination

```bash
./run_calibration.sh 0.5 0.5 0.0
```

### Collect and compare all results

```bash
Rscript calibration/collect_opt_results.R
# → opt_results_summary_grid_search.csv
```

---

## Dependencies (R)

```r
install.packages(c("optimx", "DEoptim", "zoo", "foreach",
                   "doParallel", "lubridate", "nixmass", "tidyverse"))
```

The `nixmass` package provides `swe.delta.snow()`, the deltasnow forward model used inside the objective function.
