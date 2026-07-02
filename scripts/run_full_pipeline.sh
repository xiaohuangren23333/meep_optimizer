#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
source /home/ubuntu/miniforge3/etc/profile.d/conda.sh
conda activate phc-meep
export PYTHONPATH=src

LOG="optimization_pipeline.log"
exec >> "$LOG" 2>&1

echo "===== PHC pipeline resume $(date) ====="

if [ ! -f results/best_bandgap_2d.json ]; then
  echo ">>> Stage 1: 2D bandgap scan (phase 1)"
  python src/scan_bandgap_2d.py --phase 1
else
  echo ">>> Stage 1: skip (results/best_bandgap_2d.json exists)"
fi

echo ">>> Stage 2: 2D cavity optimization (Q>=1000 target)"
rm -f results/cavity_optimization_all.csv
python src/optimize_cavity_2d.py

echo ">>> Stage 3: fast 3D bandgap verification"
python src/scan_bandgap_3d.py --nfreq 300 --resolution 14 --a-points 2 --rx-points 2 --ry-points 2 --a-span 0.01 --rx-span 0.01 --ry-span 0.01 --n-delta 2

echo ">>> Stage 4: 3D cavity verification"
python src/run_3d_cavity.py

echo "===== Pipeline finished $(date) ====="

python3 - <<'PY'
import json, os
from config import min_feature_nm, check_periodic_geometry, check_cavity_geometry

def load(path):
    return json.load(open(path)) if os.path.exists(path) else {}

bandgap = load("results_3d/best_bandgap_3d.json") or load("results/best_bandgap_params.json")
cavity = load("results/best_cavity_design.json")
cavity3d = load("results/cavity_3d/cavity_3d_result.json")

print("\n===== TARGET CHECK =====")
print(f"min_feature_nm={min_feature_nm}")
if bandgap:
    ok, v = check_periodic_geometry(bandgap["a"], bandgap["rx"], bandgap["ry"])
    print(f"Mirror: a={bandgap['a']} rx={bandgap['rx']} ry={bandgap['ry']} geom_ok={ok}")
    print(f"  center={bandgap.get('center_nm','?')} T_min={bandgap.get('T_min','?')}")
if cavity:
    ok, v = check_cavity_geometry(cavity['a_m'], cavity['rx_m'], cavity['ry_m'], cavity['a_c'], cavity['rx_c'], cavity['ry_c'], cavity['N_taper'])
    q, lam, tp = cavity.get('Q'), cavity.get('lambda0_nm'), cavity.get('T_peak')
    print(f"2D cavity: Q={q} lambda0={lam} T_peak={tp} geom_ok={ok}")
    met = (q or 0) >= 1000 and abs((lam or 0)-1550) <= 20 and (tp or 0) >= 0.9
    print(f"  2D targets met: {met}")
if cavity3d.get('3D_result'):
    r = cavity3d['3D_result']
    print(f"3D cavity: Q={r.get('Q')} lambda0={r.get('lambda0_nm')} T_peak={r.get('T_peak')}")
PY
