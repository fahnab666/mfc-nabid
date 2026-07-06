#!/usr/bin/env bash
# MFC smoke driver: build -> validate -> run a 1D Sod shock tube -> assert output.
# Usage: .claude/skills/run-mfc/smoke.sh [case_dir] [ranks]
# Exit 0 = the full pre_process+simulation pipeline produced time-step output.
set -euo pipefail

cd "$(dirname "$0")/../../.."   # repo root (skill lives in .claude/skills/run-mfc)

CASE_DIR="${1:-examples/1D_sodshocktube}"
RANKS="${2:-2}"
CASE="$CASE_DIR/case.py"

echo "== MFC smoke: $CASE with $RANKS ranks =="

# 1. Build (no-op if up to date; first build takes ~5-15 min)
./mfc.sh build -t pre_process simulation -j "$(getconf _NPROCESSORS_ONLN)"

# 2. Static validation of the case file (catches parameter errors in seconds)
./mfc.sh validate "$CASE"

# 3. Run pre_process + simulation
rm -rf "$CASE_DIR/restart_data" "$CASE_DIR/D" "$CASE_DIR/p_all"
./mfc.sh run "$CASE" -t pre_process simulation -n "$RANKS"

# 4. Assert real output. With parallel_io (most cases) conserved fields land in
# restart_data/lustre_<step>.dat; serial-IO cases write D/*.dat instead.
n_out=$( (ls "$CASE_DIR"/restart_data/lustre_*.dat "$CASE_DIR"/D/*.dat 2>/dev/null || true) | wc -l | tr -d ' ')
if [ "$n_out" -lt 2 ]; then
    echo "SMOKE FAIL: expected >=2 output files (initial + later step), found $n_out"
    exit 1
fi
echo "SMOKE PASS: $n_out field-data output files written"
