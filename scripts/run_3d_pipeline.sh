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

# Mirror params from 2D scan (200nm-compliant); 3D bandgap verified inside cavity optimizer
if [ -f results/best_bandgap_params.json ]; then
  cp -f results/best_bandgap_params.json results_3d/best_bandgap_3d.json
  echo ">>> Stage 1: use existing mirror params (skip lengthy 3D bandgap grid)"
else
  echo ">>> Stage 1: quick 3D bandgap verify"
  python src/scan_bandgap_3d.py --nfreq 300 --resolution 14 \
    --a-points 2 --rx-points 2 --ry-points 2 --a-span 0.01 --rx-span 0.01 --ry-span 0.01 --n-delta 0
fi

# Stage 2: 3D cavity coarse + refine
echo ">>> Stage 2: 3D cavity optimization (refine)"
python src/optimize_cavity_3d.py --phase refine --nfreq 400 --resolution 12 --fixed-until 120 --max-runs 24

# Stage 3: fine refine if 3D targets not met
if ! python3 - <<'PY'
import json, os
p = "results_3d/cavity/best_cavity_3d.json"
if not os.path.exists(p):
    raise SystemExit(1)
d = json.load(open(p))
if not d.get("targets_met"):
    raise SystemExit(1)
PY
then
  echo ">>> Stage 3: 3D cavity fine refine (higher resolution)"
  python src/optimize_cavity_3d.py --phase 2 --nfreq 800 --resolution 14 --fixed-until 0 --max-runs 12
fi

# Stage 4: final high-res verify
echo ">>> Stage 4: 3D cavity verification (resolution 18)"
python src/run_3d_cavity.py

python3 - <<'PY'
import json, os
from config import Q_target, min_feature_nm, check_cavity_geometry

print("\n===== 3D TARGET CHECK =====")
for label, path in [
    ("bandgap", "results_3d/best_bandgap_3d.json"),
    ("cavity_scan", "results_3d/cavity/best_cavity_3d.json"),
]:
    if os.path.exists(path):
        d = json.load(open(path))
        print(f"{label}: {d}")
if os.path.exists("results/cavity_3d/cavity_3d_result.json"):
    vr = json.load(open("results/cavity_3d/cavity_3d_result.json")).get("3D_result")
    print("cavity_verify:", vr)
print(f"min_feature_nm={min_feature_nm}, Q_target={Q_target}")
PY

echo "===== 3D pipeline finished $(date) ====="
