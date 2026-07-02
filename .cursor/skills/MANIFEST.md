# Project skills (Cursor Cloud)

Bundled under `.cursor/skills/` for Cloud Agents. Read `SKILL.md` in each folder when the task matches the skill description.

## Meep / FDTD

| Skill | Use |
|-------|-----|
| **`meep-phc-optimizer`** | **主 Skill**：二次渐变椭圆柱 PHC、Phase1/2、Meep 归一化与优化流程 |
| `optical-resonance-q` | Q from transmission spectra (FDTD: `peak_bracketed`) |
| `fdtdx-mode-source-overlap` | Mode source / overlap / transmission concepts |

## Data & code

| Skill | Use |
|-------|-----|
| `datafile-quicklook` | Inspect CSV/TSV/Excel logs |
| `xlsx` | Spreadsheet deliverables |
| `pdf` | Read/merge/process PDFs |
| `python-debugger` | Debug Meep scripts |
| `code-review` | Review simulation code |
| `git-commit-helper` | Commit messages |

## Literature (browser MCP on Cloud may need network + login)

| Orchestrator | Domain |
|--------------|--------|
| `gs-researcher` | Google Scholar |
| `ieee-researcher` | IEEE Xplore |
| `wos-researcher` | Web of Science |
| `cnki-researcher` | CNKI |
| `sd-researcher` | ScienceDirect |

Companion step skills: `gs-*`, `ieee-*`, `wos-*`, `cnki-*`, `sd-*` in this directory.

## Visualization / documents

| Skill | Use |
|-------|-----|
| `drawio-skill` | Structure / flow diagrams (needs `drawio` CLI) |
| `canvas-design` | Visual design artifacts |
| `docx` | Word reports |
| `pptx` | Slides |

## Scripts

- Q peak: `.cursor/skills/optical-resonance-q/scripts/q_peak.py`
- Re-sync from local machine: `bash scripts/sync_skills_to_repo.sh`
