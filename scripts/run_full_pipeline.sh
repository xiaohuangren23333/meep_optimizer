#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
source /home/ubuntu/miniforge3/etc/profile.d/conda.sh
conda activate phc-meep
export PYTHONPATH=src

LOG="optimization_pipeline.log"
exec > >(tee -a "$LOG") 2>&1

echo "===== PHC full optimization pipeline started $(date) ====="

echo ">>> Stage 1: 2D bandgap scan (phase 1)"
python src/scan_bandgap_2d.py --phase 1

echo ">>> Stage 2: 3D bandgap verification scan"
python src/scan_bandgap_3d.py --nfreq 400 --resolution 16

echo ">>> Stage 3: 2D cavity optimization"
python src/optimize_cavity_2d.py

echo ">>> Stage 4: 3D cavity verification"
python src/run_3d_cavity.py

echo "===== Pipeline finished $(date) ====="

python3 - <<'PY'
import json, os, csv

def load_json(path):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}

targets = {"Q": 1000, "lambda0_nm": 1550, "T_peak": 0.9, "T_min": 0.1}
bandgap = load_json("results_3d/best_bandgap_3d.json") or load_json("results/best_bandgap_params.json")
cavity = load_json("results/best_cavity_design.json")
cavity3d = load_json("results/cavity_3d/cavity_3d_result.json")

print("\n===== TARGET CHECK =====")
if bandgap:
    print(f"Bandgap center: {bandgap.get('center_nm', 'n/a')} nm, T_min={bandgap.get('T_min', 'n/a')}")
if cavity:
    print(f"2D cavity: Q={cavity.get('Q')}, lambda0={cavity.get('lambda0_nm')} nm, T_peak={cavity.get('T_peak')}")
if cavity3d and cavity3d.get("3D_result"):
    r = cavity3d["3D_result"]
    print(f"3D cavity: Q={r.get('Q')}, lambda0={r.get('lambda0_nm')} nm, T_peak={r.get('T_peak')}")
PY
