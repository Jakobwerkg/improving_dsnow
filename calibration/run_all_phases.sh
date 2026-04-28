#!/usr/bin/env bash
# Runs all calibration phases (1–4) sequentially.
# Each phase calls run_calibration.sh, which runs all 4 R scripts
# (SNOWPACK NM, SNOWPACK DE, Win21 NM, Win21 DE) for that weight combo.
#
# Usage: ./run_all_phases.sh
# Estimated runtime: ~several hours (each phase ≈ 4 × optimizer time)

set -e

BASE="$(cd "$(dirname "$0")" && pwd)"
RUN="$BASE/run_calibration.sh"

run_phase() {
  local label=$1 w1=$2 w2=$3 w3=$4
  echo ""
  echo "###################################################"
  echo "#  $label"
  echo "#  SWE=$w1  RHO=$w2  BIAS=$w3"
  echo "###################################################"
  "$RUN" "$w1" "$w2" "$w3"
}

# ------------------------------------------------------------
# PHASE 1 — SWE-only baseline (w = 1 / 0 / 0)
# ------------------------------------------------------------
run_phase "Phase 1 — SWE-only baseline"          1.0 0.0 0.0

# ------------------------------------------------------------
# PHASE 2 — Density-only extreme (w = 0 / 1 / 0)
# ------------------------------------------------------------
run_phase "Phase 2 — Density-only extreme"        0.0 1.0 0.0

# ------------------------------------------------------------
# PHASE 3 — Balanced weight sweep (no bias)
# ------------------------------------------------------------
run_phase "Phase 3A — SWE-dominant  (0.7 / 0.3)" 0.7 0.3 0.0
run_phase "Phase 3B — Equal weight  (0.5 / 0.5)" 0.5 0.5 0.0
run_phase "Phase 3C — Density-dom.  (0.3 / 0.7)" 0.3 0.7 0.0

# ------------------------------------------------------------
# PHASE 4 — Introduce bias penalty
# ------------------------------------------------------------
run_phase "Phase 4A — Balanced + bias    (0.6 / 0.2 / 0.2)" 0.6 0.2 0.2
run_phase "Phase 4B — SWE + bias only   (0.7 / 0.0 / 0.3)" 0.7 0.0 0.3
run_phase "Phase 4C — Density + bias    (0.3 / 0.5 / 0.2)" 0.3 0.5 0.2

echo ""
echo "###################################################"
echo "#  All phases complete."
echo "###################################################"
