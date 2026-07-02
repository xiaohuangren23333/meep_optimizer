#!/usr/bin/env bash
# Source in Cloud Agent shells: eval "$(bash scripts/cloud_env.sh)"
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MAMBA_ROOT_PREFIX="${HOME}/.micromamba"
export PATH="${HOME}/.local/bin:${PATH}"

if [ -x "${HOME}/.local/bin/micromamba" ]; then
  # shellcheck disable=SC1090
  eval "$("${HOME}/.local/bin/micromamba" shell hook -s bash)"
  micromamba activate phc-meep
fi

export PYTHONPATH="${ROOT}/src:${PYTHONPATH:-}"
echo "phc-meep ready (PYTHONPATH=${ROOT}/src)"
