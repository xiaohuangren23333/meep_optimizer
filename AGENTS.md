# AGENTS.md

## Cursor Cloud specific instructions

### What this project is
Meep (MIT FDTD) + Python tool for 1D photonic-crystal (elliptical-hole ridge
waveguide) bandgap scanning and defect-cavity Q optimization. Pure Python
scripts, no server/frontend. Core dependency is the real MIT **Meep**
(`import meep as mp`), which is distributed **only via conda-forge as
`pymeep`** — the `meep` package on PyPI is an unrelated placeholder, do not use it.

### Environment / how to run
- Meep lives in a conda env named `phc-meep` (Python 3.11: `pymeep`, `numpy`,
  `scipy`, `matplotlib`) under `~/miniforge3`. `conda` is not on `PATH` by
  default in non-interactive shells.
- Activate before running anything:
  `source $HOME/miniforge3/bin/activate phc-meep`
- Entry point is `run.py` (`bandgap-3d`, `cavity-2d`, `cavity-3d`, `all`), see
  `README.md`. Scripts read/write under `results/` and `results_3d/`.

### Non-obvious runtime caveats
- The `run.py` stages and the 3D scripts are **long-running (hours)**: they run
  full-resolution 3D FDTD twice per parameter point (structure + no-hole
  reference) with a slow field-decay stop condition. Do not expect a single 3D
  simulation (`src/optimize_3d_ridge.py:run_sim`) to finish within a few
  minutes even at low hole counts.
- For a fast end-to-end smoke test of the core FDTD pipeline, use the 2D bandgap
  Phase-0 scan (~15 s, 9 parameter points):
  `python src/scan_bandgap_2d.py --phase 0`
  It writes spectra CSVs to `results/spectra/`, figures to `results/figures/`,
  and candidates to `results/bandgap_candidates_2d.json`.
- Matplotlib is forced to the `Agg` backend in the scripts, so figures render
  headless without a display.
- `T_max > 1` in 2D output is a known EigenModeSource normalization artifact
  (see `PROJECT_SUMMARY.md`); the bandgap depth `T_min` is still meaningful.
- There are no automated tests, no linter config, and no build step in this
  repo; "running the app" means executing the scripts above.
