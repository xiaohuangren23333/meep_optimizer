# Agent instructions

See `AGENT_INSTRUCTION.txt` for simulation protocol (GaussianSource, reference normalization, no per-hole optimization, 2D ≠ 3D).

## Cursor Cloud specific instructions

This repo ships a prebuilt Cloud Agent environment via `.cursor/environment.json` and `.cursor/Dockerfile`:

- Conda env: `phc-meep` (Python 3.11, Meep, numpy, scipy, matplotlib, pandas, h5py, mpi4py)
- Bootstrap (idempotent): `bash scripts/bootstrap_cloud.sh`
- Skill deps (idempotent): `bash scripts/bootstrap_skill_deps.sh`
- Activate in a shell: `source scripts/cloud_env.sh`

### Bundled Agent Skills (`.cursor/skills/`)

Skills are **in-repo** so Cloud Agents can read them without your local `C:\Users\dell\.cursor\skills`. Index: `.cursor/skills/MANIFEST.md`.

| Category | Skills |
|----------|--------|
| Meep / Q | `optical-resonance-q`, `fdtdx-mode-source-overlap` |
| Data / code | `datafile-quicklook`, `xlsx`, `pdf`, `python-debugger`, `code-review`, `git-commit-helper` |
| Literature | `gs-researcher`, `ieee-researcher`, `wos-researcher`, `cnki-researcher`, `sd-researcher` (+ companion `*-search`, `*-export`, …) |
| Docs / viz | `drawio-skill`, `canvas-design`, `docx`, `pptx` |

When a task matches a skill description, read `.cursor/skills/<name>/SKILL.md` before acting.

**Q from spectra (FDTD):**
```bash
source scripts/cloud_env.sh
python .cursor/skills/optical-resonance-q/scripts/q_peak.py \
  --method peak_bracketed --wl ... --T ...
```

**Literature skills** need browser MCP + network; CNKI/IEEE may require login on your side.

### Run simulations

Always use the `phc-meep` interpreter (via `run.py` after sourcing `cloud_env.sh`, or `scripts/run_cloud_pipeline.sh`):

```bash
source scripts/cloud_env.sh
python run.py mirror-v2       # Phase 1: 3D mirror bandgap
python run.py cavity-3d-v2    # Phase 2: 3D cavity (needs Phase 1 output)
python run.py phase1-2        # Phase 1 then Phase 2
bash scripts/run_cloud_pipeline.sh   # bootstrap + Phase 1 + Phase 2
```

Logs: `logs/mirror_phase1.log`, `logs/cavity_phase2.log`  
Outputs: `results_3d/best_mirror_phase1.json`, `results/best_cavity_design.json`

### Notes

- 3D Meep runs are CPU-heavy; expect hours for full Phase 1 + Phase 2.
- Do not reinstall Meep unless bootstrap fails; fix `.cursor/Dockerfile` instead.
- After changing `.cursor/environment.json`, `.cursor/Dockerfile`, or `.cursor/skills/`, push to GitHub and start a Cloud Agent from that branch to rebuild.
- Re-sync skills from your Windows machine: `bash scripts/sync_skills_to_repo.sh` then commit.
