#!/usr/bin/env bash
# Full 3D-targeted optimization pipeline
set -euo pipefail
cd "$(dirname "$0")/.."
source /home/ubuntu/miniforge3/etc/profile.d/conda.sh
conda activate phc-meep
export PYTHONPATH=src

LOG="optimization_3d_pipeline.log"
exec > >(tee -a "$LOG") 2>&1

echo "===== 3D optimization pipeline $(date) ====="

# Stage 1: 3D bandgap verification around best 2D mirror (200nm-compliant)
if [ ! -f results_3d/best_bandgap_3d.json ]; then
  echo ">>> Stage 1: 3D bandgap scan"
  python src/scan_bandgap_3d.py \
    --best-2d results/best_bandgap_2d.json \
    --nfreq 400 --resolution 14 \
    --a-span 0.02 --rx-span 0.015 --ry-span 0.015 \
    --a-points 3 --rx-points 3 --ry-points 3 --n-delta 2
else
  echo ">>> Stage 1: skip (best_bandgap_3d.json exists)"
fi

# Stage 2: 3D cavity coarse scan
echo ">>> Stage 2: 3D cavity optimization (phase refine)"
python src/optimize_cavity_3d.py --phase refine --nfreq 600 --resolution 14 --max-runs 40

# Stage 3: 3D cavity fine refine if targets not met
if ! python3 - <<'PY'
import json, os
p = "results_3d/cavity/best_cavity_3d.json"
if not os.path.exists(p):
    raise SystemExit(1)
d = json.load(open(p))
if not d.get("targets_met"):
    raise SystemExit(1)
print("targets met")
PY
then
  echo ">>> Stage 3: 3D cavity fine refine"
  python src/optimize_cavity_3d.py --phase 2 --nfreq 900 --resolution 16 --max-runs 25
fi

# Stage 4: high-resolution verify best 3D design
echo ">>> Stage 4: 3D cavity high-resolution verify"
python src/run_3d_cavity.py

python3 - <<'PY'
import json, os
from config import Q_target, min_feature_nm, check_cavity_geometry

print("\n===== 3D TARGET CHECK =====")
bg = json.load(open("results_3d/best_bandgap_3d.json")) if os.path.exists("results_3d/best_bandgap_3d.json") else {}
cv = json.load(open("results_3d/cavity/best_cavity_3d.json")) if os.path.exists("results_3d/cavity/best_cavity_3d.json") else {}
vr = json.load(open("results/cavity_3d/cavity_3d_result.json")).get("3D_result") if os.path.exists("results/cavity_3d/cavity_3d_result.json") else None

if bg:
    print(f"3D bandgap: center={bg.get('center_nm')} T_min={bg.get('T_min')}")
if cv:
    ok, v = check_cavity_geometry(cv['a_m'], cv['rx_m'], cv['ry_m'], cv['a_c'], cv['rx_c'], cv['ry_c'], cv['N_taper'])
    print(f"3D scan best: Q={cv.get('Q')} λ={cv.get('lambda0_nm')} T={cv.get('T_peak')} geom_ok={ok}")
if vr:
    print(f"3D verify: Q={vr.get('Q')} λ={vr.get('lambda0_nm')} T={vr.get('T_peak')}")
print(f"min_feature_nm={min_feature_nm}, Q_target={Q_target}")
PY

echo "===== 3D pipeline finished $(date) ====="
