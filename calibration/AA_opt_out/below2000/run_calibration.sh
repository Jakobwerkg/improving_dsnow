#!/usr/bin/env bash
# Usage: ./run_calibration.sh <w_swe> <w_rho> <w_bias>
# Runs all 4 calibration scripts sequentially with the given objective weights.
# Results are saved automatically by each script (tagged by weight combination).
#
# Example:
#   ./run_calibration.sh 1.0 0.0 0.0   # SWE-only baseline
#   ./run_calibration.sh 0.5 0.5 0.0   # equal SWE + density
#   ./run_calibration.sh 0.3 0.7 0.0   # density-dominant

set -e

W1=${1:?Error: provide w_swe  (e.g. 1.0)}
W2=${2:?Error: provide w_rho  (e.g. 0.0)}
W3=${3:?Error: provide w_bias (e.g. 0.0)}

BASE="$(cd "$(dirname "$0")" && pwd)"

echo "========================================================"
echo "  DeltaSnow calibration  |  SWE=$W1  RHO=$W2  BIAS=$W3"
echo "========================================================"

echo ""
echo "[1/2] SNOWPACK — Nelder-Mead"
Rscript "$BASE/calibration_SNOWPACK/dsnow_parameter_optimization.R" "$W1" "$W2" "$W3"

echo ""
echo "[2/2] SNOWPACK — Differential Evolution"
Rscript "$BASE/calibration_SNOWPACK/dsnow_parameter_optimization_DE.R" "$W1" "$W2" "$W3"

# echo ""
# echo "[3/4] Win21 — Nelder-Mead"
# Rscript "$BASE/calibration_Win21/WIN21_dsnow_paprameter_optimization.R" "$W1" "$W2" "$W3"

# echo ""
# echo "[4/4] Win21 — Differential Evolution"
# Rscript "$BASE/calibration_Win21/WIN21_dsnow_parameter_opitmization_DE.R" "$W1" "$W2" "$W3"

echo ""
echo "========================================================"
echo "  All 2 calibrations complete."
echo "  Results saved to:"
echo "    calibration_SNOWPACK/data/R_opt_logs/"
echo "    calibration_SNOWPACK/data/R_opt_logs_DE/"2
echo "========================================================"