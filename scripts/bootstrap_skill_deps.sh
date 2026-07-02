#!/usr/bin/env bash
# Extra Python/system deps for bundled .cursor/skills (idempotent)
set -euo pipefail

SENTINEL="${HOME}/.cache/phc-meep/.skill_deps_installed"
MAMBA_BIN="${HOME}/.local/bin/micromamba"

if [ -f "${SENTINEL}" ]; then
  echo "[skill-deps] already installed, skip."
  exit 0
fi

export MAMBA_ROOT_PREFIX="${HOME}/.micromamba"
if [ -x "${MAMBA_BIN}" ]; then
  # shellcheck disable=SC1090
  eval "$("${MAMBA_BIN}" shell hook -s bash)"
  micromamba activate phc-meep 2>/dev/null || true
fi

echo "[skill-deps] installing Python packages..."
python -m pip install --quiet --upgrade pip
python -m pip install --quiet \
  openpyxl \
  pypdf \
  pdfplumber \
  python-docx \
  python-pptx \
  pillow \
  chardet

# draw.io CLI for drawio-skill (headless export)
if ! command -v draw.io >/dev/null 2>&1 && ! command -v drawio >/dev/null 2>&1; then
  if command -v sudo >/dev/null 2>&1; then
    echo "[skill-deps] installing draw.io desktop CLI..."
    DRAWIO_DEB="/tmp/drawio-amd64.deb"
    curl -fsSL -o "${DRAWIO_DEB}" \
      "https://github.com/jgraph/drawio-desktop/releases/download/v26.0.4/drawio-amd64-26.0.4.deb"
    sudo apt-get update -qq
    sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
      "${DRAWIO_DEB}" \
      libgtk-3-0 libnotify4 libnss3 libxss1 libxtst6 xdg-utils \
      libatspi2.0-0 libdrm2 libgbm1 libasound2t64 \
      || sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq "${DRAWIO_DEB}" || true
    rm -f "${DRAWIO_DEB}"
  else
    echo "[skill-deps] draw.io not installed (no sudo); .drawio XML still works."
  fi
fi

python - <<'PY'
import importlib
for m in ("openpyxl", "pypdf", "pdfplumber", "docx", "pptx", "PIL"):
    importlib.import_module(m)
print("skill python deps ok")
PY

mkdir -p "$(dirname "${SENTINEL}")"
echo ok > "${SENTINEL}"
echo "[skill-deps] done."
