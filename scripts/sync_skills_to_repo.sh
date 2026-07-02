#!/usr/bin/env bash
# Sync selected skills from local Cursor skills dir into repo .cursor/skills/
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${SKILLS_SRC:-${HOME}/.cursor/skills}"
# Windows Cursor default (WSL)
if [ ! -d "${SRC}/optical-resonance-q" ] && [ -d "/mnt/c/Users/dell/.cursor/skills/optical-resonance-q" ]; then
  SRC="/mnt/c/Users/dell/.cursor/skills"
fi

DEST="${ROOT}/.cursor/skills"
mkdir -p "${DEST}"

SKILLS=(
  optical-resonance-q fdtdx-mode-source-overlap datafile-quicklook
  xlsx pdf python-debugger code-review git-commit-helper
  gs-researcher gs-search gs-advanced-search gs-cited-by gs-navigate-pages gs-fulltext gs-export
  ieee-researcher ieee-search ieee-advanced-search ieee-parse-results ieee-navigate-pages
  ieee-paper-detail ieee-journal-browse ieee-download ieee-export ieee-standards-search
  wos-researcher wos-search wos-parse-results wos-navigate-pages wos-paper-detail
  wos-download wos-export wos-dom
  cnki-researcher cnki-search cnki-parse-results cnki-paper-detail cnki-navigate-pages
  cnki-journal-toc cnki-journal-search cnki-journal-index cnki-export cnki-download cnki-advanced-search
  sd-researcher sd-search sd-paper-detail sd-parse-results sd-navigate-pages sd-journal-browse
  sd-export sd-download sd-advanced-search
  drawio-skill canvas-design docx pptx
)

echo "[sync] source: ${SRC}"
echo "[sync] dest:   ${DEST}"

for name in "${SKILLS[@]}"; do
  if [ ! -d "${SRC}/${name}" ]; then
    echo "[sync] SKIP missing: ${name}" >&2
    continue
  fi
  rsync -a --delete \
    --exclude '.git' \
    --exclude '__pycache__' \
    --exclude '*.pyc' \
    "${SRC}/${name}/" "${DEST}/${name}/"
  echo "[sync] OK ${name}"
done

# Rewrite global skill paths → repo-relative for Cloud
find "${DEST}" -type f \( -name '*.md' -o -name '*.py' -o -name '*.sh' \) -print0 \
  | while IFS= read -r -d '' f; do
      sed -i \
        -e 's|~/.cursor/skills/|.cursor/skills/|g' \
        -e 's|\$HOME/.cursor/skills/|.cursor/skills/|g' \
        -e 's|C:\\\\Users\\\\dell\\\\.cursor\\\\skills\\\\|.cursor/skills/|g' \
        -e 's|C:/Users/dell/.cursor/skills/|.cursor/skills/|g' \
        "$f" 2>/dev/null || true
    done

echo "[sync] done: $(find "${DEST}" -mindepth 1 -maxdepth 1 -type d | wc -l) skills"
