#!/usr/bin/env python3
"""
3D FDTD defect cavity optimizer.

Scans cavity parameters directly in 3D (not 2D proxy) until targets are met:
  Q >= 1000, lambda0 ~ 1550 nm, T_peak high, 200 nm min feature.
"""
import argparse
import csv
import json
import os
import time

import numpy as np

from config import (
    Q_target,
    check_cavity_geometry,
    check_periodic_geometry,
    min_feature_nm,
    nfreq_cavity_3d,
    resolution_3d,
)
from phc_3d import (
    analyze_bandgap_3d,
    build_cavity_geom_3d,
    build_periodic_geom_3d,
    build_ref_geom_3d,
    fit_cavity_peak,
    meets_3d_targets,
    normalized_spectrum,
    run_flux_sim_3d,
    score_cavity_3d,
)

RESULTS_DIR = "results_3d"
CAVITY_DIR = os.path.join(RESULTS_DIR, "cavity")
SPECTRA_DIR = os.path.join(CAVITY_DIR, "spectra")
CSV_PATH = os.path.join(CAVITY_DIR, "cavity_3d_scan.csv")
BEST_PATH = os.path.join(CAVITY_DIR, "best_cavity_3d.json")
COMPAT_PATH = "results/best_cavity_design.json"
BANDGAP_PATH = os.path.join(RESULTS_DIR, "best_bandgap_3d.json")
MIRROR_FALLBACK = "results/best_bandgap_params.json"

CSV_FIELDS = [
    "run_id", "label", "a_m", "rx_m", "ry_m", "a_c", "rx_c", "ry_c",
    "N_taper", "N_mirror", "lambda0_nm", "FWHM_nm", "Q", "T_peak", "fit_r2",
    "score", "targets_met", "elapsed_s", "spectrum_csv", "error",
]


def ensure_dirs():
    os.makedirs(SPECTRA_DIR, exist_ok=True)
    os.makedirs(CAVITY_DIR, exist_ok=True)


def load_mirror():
    for path in (BANDGAP_PATH, MIRROR_FALLBACK):
        if os.path.exists(path):
            with open(path) as f:
                raw = json.load(f)
            mirror = {
                "a_m": float(raw.get("a_m") or raw["a"]),
                "rx_m": float(raw.get("rx_m") or raw["rx"]),
                "ry_m": float(raw.get("ry_m") or raw["ry"]),
                "bandgap_start_nm": float(raw.get("gap_start_nm", 1480)),
                "bandgap_end_nm": float(raw.get("gap_end_nm", 1620)),
                "gap_width_nm": float(raw.get("gap_width_nm", 0)),
                "T_min_gap": float(raw.get("T_min") or raw.get("T_min_gap", 1)),
            }
            ok, viol = check_periodic_geometry(mirror["a_m"], mirror["rx_m"], mirror["ry_m"])
            if not ok:
                raise RuntimeError(f"Mirror violates {min_feature_nm}nm: {viol}")
            return mirror
    raise FileNotFoundError("Run bandgap-3d first to create best_bandgap_3d.json")


def load_seeds():
    seeds = []
    for path in (
        "results/best_cavity_design_highQ.json",
        "results/best_cavity_design_balanced.json",
        COMPAT_PATH,
    ):
        if not os.path.exists(path):
            continue
        with open(path) as f:
            d = json.load(f)
        seeds.append({
            "name": os.path.basename(path),
            "a_c": float(d["a_c"]),
            "rx_c": float(d.get("rx_c", d.get("rx_m", 0.22))),
            "ry_c": float(d.get("ry_c", d.get("ry_m", 0.20))),
            "N_taper": int(d.get("N_taper", 3)),
            "N_mirror": int(d.get("N_mirror", 16)),
        })
    if not seeds:
        seeds.append({
            "name": "default",
            "a_c": 0.785, "rx_c": 0.230, "ry_c": 0.207,
            "N_taper": 3, "N_mirror": 16,
        })
    return seeds


def run_periodic_bandgap_3d(a, rx, ry, n_period, nfreq, resolution):
    geom, cell = build_periodic_geom_3d(a, rx, ry, n_period)
    fh, fr = run_flux_sim_3d(geom, cell, nfreq, resolution)
    rr = run_flux_sim_3d(build_ref_geom_3d(), cell, nfreq, resolution)
    wl, tr = normalized_spectrum(fh, fr, rr[0], rr[1])
    return analyze_bandgap_3d(wl, tr)


def evaluate_cavity_3d(run_id, mirror, cavity, nfreq, resolution, bandgap=None,
                         fixed_until=120, high_accuracy=False):
    label = (f"ac{cavity['a_c']:.4f}_rxc{cavity['rx_c']:.3f}_"
             f"Nm{cavity['N_mirror']}_Nt{cavity['N_taper']}")
    ok, viol = check_cavity_geometry(
        mirror["a_m"], mirror["rx_m"], mirror["ry_m"],
        cavity["a_c"], cavity["rx_c"], cavity["ry_c"], cavity["N_taper"],
    )
    if not ok:
        return {"run_id": run_id, "label": label, "error": "; ".join(viol),
                "score": -1e6, "targets_met": False}

    print(f"[{run_id:03d}] 3D {label}", end="", flush=True)
    t0 = time.time()
    try:
        geom, sx, sy, sz = build_cavity_geom_3d(mirror, cavity)
        cell = (sx, sy, sz)
        kw = {} if high_accuracy else {"fixed_until": fixed_until}
        freqs_h, flux_h = run_flux_sim_3d(geom, cell, nfreq, resolution, **kw)
        freqs_r, flux_r = run_flux_sim_3d(build_ref_geom_3d(), cell, nfreq, resolution, **kw)
        wl, tr = normalized_spectrum(freqs_h, flux_h, freqs_r, flux_r)
        peak = fit_cavity_peak(wl, tr, bandgap=mirror)
        elapsed = time.time() - t0
        spectrum_path = os.path.join(SPECTRA_DIR, f"cavity3d_{run_id:04d}_{label}.csv")
        np.savetxt(spectrum_path, np.column_stack([wl, tr]),
                   delimiter=",", header="wavelength_nm,transmission", comments="")

        score = score_cavity_3d(peak)
        met, reasons = meets_3d_targets(peak, bandgap=bandgap)
        row = {
            "run_id": run_id,
            "label": label,
            "a_m": mirror["a_m"], "rx_m": mirror["rx_m"], "ry_m": mirror["ry_m"],
            "a_c": cavity["a_c"], "rx_c": cavity["rx_c"], "ry_c": cavity["ry_c"],
            "N_taper": cavity["N_taper"], "N_mirror": cavity["N_mirror"],
            "lambda0_nm": peak["lambda0_nm"] if peak else 0,
            "FWHM_nm": peak["FWHM_nm"] if peak else 0,
            "Q": peak["Q"] if peak else 0,
            "T_peak": peak["T_peak"] if peak else 0,
            "fit_r2": peak["fit_r2"] if peak else 0,
            "score": round(score, 2),
            "targets_met": met,
            "elapsed_s": round(elapsed, 1),
            "spectrum_csv": spectrum_path,
            "error": "" if peak else "no peak",
        }
        if peak:
            flag = "🎯" if met else ""
            print(f" | Q={peak['Q']:.0f} λ={peak['lambda0_nm']:.1f} "
                  f"T={peak['T_peak']:.3f} [{elapsed:.0f}s] {flag}")
        else:
            print(f" | no peak [{elapsed:.0f}s]")
        return row
    except Exception as exc:
        print(f" | ERROR {exc}")
        return {"run_id": run_id, "label": label, "error": str(exc),
                "score": -1e6, "targets_met": False, "elapsed_s": round(time.time() - t0, 1)}


def append_csv(row):
    exists = os.path.isfile(CSV_PATH)
    with open(CSV_PATH, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        if not exists:
            w.writeheader()
        w.writerow(row)


def save_best(row, mirror):
    payload = {
        "a_m": mirror["a_m"], "rx_m": mirror["rx_m"], "ry_m": mirror["ry_m"],
        "a_c": row["a_c"], "rx_c": row["rx_c"], "ry_c": row["ry_c"],
        "N_taper": int(row["N_taper"]), "N_mirror": int(row["N_mirror"]),
        "lambda0_nm": row["lambda0_nm"], "FWHM_nm": row["FWHM_nm"],
        "Q": row["Q"], "T_peak": row["T_peak"], "fit_r2": row["fit_r2"],
        "score": row["score"], "targets_met": row["targets_met"],
        "source": "3D FDTD cavity scan",
        "scan_csv": CSV_PATH,
    }
    with open(BEST_PATH, "w") as f:
        json.dump(payload, f, indent=2)
    with open(COMPAT_PATH, "w") as f:
        json.dump(payload, f, indent=2)


def build_phase1_grid(seeds, mirror):
    grid = []
    for seed in seeds:
        ac0 = seed["a_c"]
        nm0 = min(seed["N_mirror"], 18)
        for nm in sorted({10, 14, nm0}):
            for ac in np.round(np.linspace(ac0 - 0.015, ac0 + 0.015, 4), 4):
                cavity = {
                    "a_c": float(ac),
                    "rx_c": seed["rx_c"],
                    "ry_c": seed["ry_c"],
                    "N_taper": seed["N_taper"],
                    "N_mirror": int(nm),
                }
                ok, _ = check_cavity_geometry(
                    mirror["a_m"], mirror["rx_m"], mirror["ry_m"],
                    cavity["a_c"], cavity["rx_c"], cavity["ry_c"], cavity["N_taper"],
                )
                if ok:
                    grid.append(cavity)
    # dedupe
    seen = set()
    unique = []
    for c in grid:
        key = (c["a_c"], c["rx_c"], c["ry_c"], c["N_taper"], c["N_mirror"])
        if key not in seen:
            seen.add(key)
            unique.append(c)
    return unique


def build_refine_grid(best_row, mirror):
    grid = []
    ac0 = float(best_row["a_c"])
    nm0 = int(best_row["N_mirror"])
    for nm in range(max(nm0 - 2, 12), min(nm0 + 8, 36), 2):
        for ac in np.round(np.linspace(ac0 - 0.025, ac0 + 0.015, 9), 4):
            cavity = {
                "a_c": float(ac),
                "rx_c": float(best_row["rx_c"]),
                "ry_c": float(best_row["ry_c"]),
                "N_taper": int(best_row["N_taper"]),
                "N_mirror": int(nm),
            }
            ok, _ = check_cavity_geometry(
                mirror["a_m"], mirror["rx_m"], mirror["ry_m"],
                cavity["a_c"], cavity["rx_c"], cavity["ry_c"], cavity["N_taper"],
            )
            if ok:
                grid.append(cavity)
    return grid


def _as_bool(val):
    if isinstance(val, bool):
        return val
    return str(val).lower() in ("1", "true", "yes")


def pick_best(rows):
    valid = [r for r in rows if float(r.get("Q", 0) or 0) > 0]
    if not valid:
        return None
    hits = [r for r in valid if _as_bool(r.get("targets_met"))]
    if hits:
        return max(hits, key=lambda r: float(r["Q"]))
    return max(valid, key=lambda r: (float(r["score"]), float(r["Q"])))


def verify_mirror_3d(mirror, nfreq, resolution, fixed_until=80):
    print("Verifying mirror bandgap in 3D...")
    geom, cell = build_periodic_geom_3d(mirror["a_m"], mirror["rx_m"], mirror["ry_m"], 20)
    fh, _ = run_flux_sim_3d(geom, cell, nfreq, resolution, fixed_until=fixed_until)
    rr = run_flux_sim_3d(build_ref_geom_3d(), cell, nfreq, resolution, fixed_until=fixed_until)
    wl, tr = normalized_spectrum(fh, _, rr[0], rr[1])
    bg = analyze_bandgap_3d(wl, tr)
    print(f"  3D bandgap: {bg.get('gap_start_nm',0):.0f}-{bg.get('gap_end_nm',0):.0f} nm "
          f"T_min_gap={bg.get('T_min_gap',1):.4f}")
    mirror["bandgap_start_nm"] = bg.get("gap_start_nm", mirror["bandgap_start_nm"])
    mirror["bandgap_end_nm"] = bg.get("gap_end_nm", mirror["bandgap_end_nm"])
    mirror["T_min_gap"] = bg.get("T_min_gap", mirror["T_min_gap"])
    return bg


def main():
    parser = argparse.ArgumentParser(description="3D defect cavity optimizer")
    parser.add_argument("--phase", default="1", choices=["1", "2", "refine"],
                        help="1=coarse grid, 2=refine best, refine=coarse+refine loop")
    parser.add_argument("--nfreq", type=int, default=400)
    parser.add_argument("--resolution", type=int, default=12)
    parser.add_argument("--fixed-until", type=int, default=120,
                        help="Fixed FDTD time steps for fast scan (0=use field decay)")
    parser.add_argument("--max-runs", type=int, default=24)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    ensure_dirs()
    mirror = load_mirror()

    if args.phase == "refine":
        grid = build_phase1_grid(load_seeds(), mirror)
    elif args.phase == "1":
        grid = build_phase1_grid(load_seeds(), mirror)
    elif args.phase == "2":
        rows = list(csv.DictReader(open(CSV_PATH))) if os.path.exists(CSV_PATH) else []
        best = pick_best(rows)
        if not best:
            print("No prior results; running phase 1 grid")
            grid = build_phase1_grid(load_seeds(), mirror)
        else:
            grid = build_refine_grid(best, mirror)

    print(f"3D cavity scan: {len(grid)} points, nfreq={args.nfreq}, res={args.resolution}")
    print(f"Targets: Q>={Q_target}, λ≈1550nm, min feature {min_feature_nm}nm")
    if args.dry_run:
        for c in grid[:15]:
            print(f"  a_c={c['a_c']:.4f} Nm={c['N_mirror']}")
        return

    bandgap = verify_mirror_3d(mirror, min(args.nfreq, 400), args.resolution)

    run_id = 0
    if os.path.exists(CSV_PATH):
        run_id = max(int(r["run_id"]) for r in csv.DictReader(open(CSV_PATH)))

    all_rows = []
    best_row = None
    target_hit = None

    fixed_until = args.fixed_until if args.fixed_until > 0 else None

    for cavity in grid[: args.max_runs]:
        run_id += 1
        row = evaluate_cavity_3d(run_id, mirror, cavity, args.nfreq, args.resolution,
                                 bandgap, fixed_until=fixed_until)
        append_csv(row)
        if float(row.get("Q", 0) or 0) > 0:
            all_rows.append(row)
            if best_row is None or row["score"] > best_row["score"]:
                best_row = row
            if row.get("targets_met"):
                target_hit = row
                save_best(row, mirror)
                print(f"\n🎯 3D TARGETS MET: Q={row['Q']:.0f} λ={row['lambda0_nm']:.1f} "
                      f"T={row['T_peak']:.3f}")
                break

    if target_hit is None and best_row:
        save_best(best_row, mirror)
        print(f"\nBest 3D so far: Q={best_row['Q']:.0f} λ={best_row['lambda0_nm']:.1f} "
              f"T={best_row['T_peak']:.3f} -> {BEST_PATH}")

    if args.phase == "refine" and target_hit is None and best_row:
        print("\nRefining around best 3D design...")
        for cavity in build_refine_grid(best_row, mirror)[:30]:
            run_id += 1
            row = evaluate_cavity_3d(run_id, mirror, cavity, args.nfreq + 200,
                                     min(args.resolution + 2, 16), bandgap,
                                     fixed_until=None, high_accuracy=True)
            append_csv(row)
            if row.get("targets_met"):
                save_best(row, mirror)
                print(f"\n🎯 3D TARGETS MET on refine: Q={row['Q']:.0f}")
                return
            if float(row.get("Q", 0) or 0) > float(best_row.get("Q", 0)):
                best_row = row
                save_best(row, mirror)


if __name__ == "__main__":
    main()
