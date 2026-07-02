#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
source /home/ubuntu/miniforge3/etc/profile.d/conda.sh
conda activate phc-meep
export PYTHONPATH=src
python src/run_3d_cavity.py
