#!/usr/bin/env python3
"""
3D FDTD graded (mirror-free) defect cavity optimizer.

Design philosophy (user request): remove the abrupt uniform mirror block and
use a long, smooth taper of holes. The lattice PERIOD is graded quadratically
from a small center value ``a_center`` (which forms the defect) out to the full
mirror period ``a_end`` (which reflects); the taper end itself acts as the
mirror. Hole size is kept constant to preserve the 3D bandgap.

Scans (a_center, N_taper) directly in 3D FDTD until targets are met:
  Q >= 1000, lambda0 ~ 1550 nm, high T_peak, 200 nm min feature.
"""
import argparse
import csv
import json
import os
import time

import numpy as np

from config import (
    Q_target,
    check_hole_geometry,
    min_feature_nm,
    min_feature_um,
    w_wg,
)
from phc_3d import (
    analyze_bandgap_3d,
    build_graded_cavity_geom_3d,
    build_periodic_geom_3d,
    build_ref_geom_3d,
    fit_cavity_peak,
    graded_hole_layout,
    meets_3d_targets,
    normalized_spectrum,
    run_flux_sim_3d,
    score_cavity_3d,
)

RESULTS_DIR = "results_3d"
CAVITY_DIR = os.path.join(RESULTS_DIR, "cavity")
SPECTRA_DIR = os.path.join(CAVITY_DIR, "spectra")
CSV_PATH = os.path.join(CAVITY_DIR, "cavity_3d_graded_scan.csv")
BEST_PATH = os.path.join(CAVITY_DIR, "best_cavity_3d.json")
COMPAT_PATH = "results/best_cavity_design.json"
BANDGAP_PATH = os.path.join(RESULTS_DIR, "best_bandgap_3d.json")
MIRROR_FALLBACK = "results/best_bandgap_params.json"

CSV_FIELDS = [
    "run_id", "label", "a_end", "a_center", "rx", "ry", "N_taper", "N_mirror",
    "total_holes", "device_um", "lambda0_nm", "FWHM_nm", "Q", "T_peak", "fit_r2",
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
                "a_end": float(raw.get("a_m") or raw["a"]),
                "rx": float(raw.get("rx_m") or raw["rx"]),
                "ry": float(raw.get("ry_m") or raw["ry"]),
                "bandgap_start_nm": float(raw.get("gap_start_nm", 1480)),
                "bandgap_end_nm": float(raw.get("gap_end_nm", 1620)),
                "T_min_gap": float(raw.get("T_min") or raw.get("T_min_gap", 1)),
            }
            ok, viol = check_hole_geometry(mirror["a_end"], mirror["rx"], mirror["ry"])
            if not ok:
                raise RuntimeError(f"Mirror violates {min_feature_nm}nm: {viol}")
            return mirror
    raise FileNotFoundError("Run bandgap first to create mirror params JSON")


def geometry_ok(a_end, a_center, rx, ry, N_taper, N_mirror):
    """All graded holes must respect the 200 nm min feature (smallest period is a_center)."""
    ok, _ = check_hole_geometry(a_center, rx, ry, w_wg=w_wg)
    if not ok:
        return False
    ok, _ = check_hole_geometry(a_end, rx, ry, w_wg=w_wg)
    return ok


def evaluate(run_id, mirror, a_center, N_taper, N_mirror, nfreq, resolution,
             bandgap=None, fixed_until=None):
    a_end = mirror["a_end"]
    rx = mirror["rx"]
    ry = mirror["ry"]
    label = f"ac{a_center:.3f}_Nt{N_taper}_Nm{N_mirror}"
    if not geometry_ok(a_end, a_center, rx, ry, N_taper, N_mirror):
        return {"run_id": run_id, "label": label, "error": "min feature",
                "score": -1e6, "targets_met": False}

    positions, _ = graded_hole_layout(a_end, a_center, N_taper, N_mirror)
    total_holes = 2 * len(positions)
    device_um = 2 * (positions[-1] + rx)

    print(f"[{run_id:03d}] {label} holes={total_holes} L={device_um:.1f}um",
          end="", flush=True)
    t0 = time.time()
    try:
        geom, sx, sy, sz = build_graded_cavity_geom_3d(
            a_end, a_center, rx, ry, N_taper, N_mirror)
        cell = (sx, sy, sz)
        kw = {"fixed_until": fixed_until} if fixed_until else {}
        fh, flux_h = run_flux_sim_3d(geom, cell, nfreq, resolution, **kw)
        fr, flux_r = run_flux_sim_3d(build_ref_geom_3d(), cell, nfreq, resolution, **kw)
        wl, tr = normalized_spectrum(fh, flux_h, fr, flux_r)
        peak = fit_cavity_peak(wl, tr, bandgap=mirror)
        elapsed = time.time() - t0
        spectrum_path = os.path.join(SPECTRA_DIR, f"graded_{run_id:04d}_{label}.csv")
        np.savetxt(spectrum_path, np.column_stack([wl, tr]),
                   delimiter=",", header="wavelength_nm,transmission", comments="")

        score = score_cavity_3d(peak)
        met, _ = meets_3d_targets(peak, bandgap=bandgap)
        row = {
            "run_id": run_id, "label": label,
            "a_end": a_end, "a_center": a_center, "rx": rx, "ry": ry,
            "N_taper": N_taper, "N_mirror": N_mirror,
            "total_holes": total_holes, "device_um": round(device_um, 2),
            "lambda0_nm": peak["lambda0_nm"] if peak else 0,
            "FWHM_nm": peak["FWHM_nm"] if peak else 0,
            "Q": peak["Q"] if peak else 0,
            "T_peak": peak["T_peak"] if peak else 0,
            "fit_r2": peak["fit_r2"] if peak else 0,
            "score": round(score, 2), "targets_met": met,
            "elapsed_s": round(elapsed, 1), "spectrum_csv": spectrum_path,
            "error": "" if peak else "no peak",
        }
        if peak:
            flag = "TARGET" if met else ""
            print(f" | Q={peak['Q']:.0f} lam={peak['lambda0_nm']:.1f} "
                  f"T={peak['T_peak']:.3f} [{elapsed:.0f}s] {flag}")
        else:
            print(f" | no peak [{elapsed:.0f}s]")
        return row
    except Exception as exc:
        print(f" | ERROR {exc}")
        return {"run_id": run_id, "label": label, "error": str(exc),
                "score": -1e6, "targets_met": False,
                "elapsed_s": round(time.time() - t0, 1)}


def append_csv(row):
    exists = os.path.isfile(CSV_PATH)
    with open(CSV_PATH, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        if not exists:
            w.writeheader()
        w.writerow(row)


def save_best(row, mirror):
    payload = {
        "design": "graded_mirror_free",
        "a_m": mirror["a_end"], "rx_m": mirror["rx"], "ry_m": mirror["ry"],
        "a_end": row["a_end"], "a_center": row["a_center"],
        "rx_c": row["rx"], "ry_c": row["ry"],
        "N_taper": int(row["N_taper"]), "N_mirror": int(row["N_mirror"]),
        "total_holes": int(row["total_holes"]), "device_um": row["device_um"],
        "lambda0_nm": row["lambda0_nm"], "FWHM_nm": row["FWHM_nm"],
        "Q": row["Q"], "T_peak": row["T_peak"], "fit_r2": row["fit_r2"],
        "score": row["score"], "targets_met": row["targets_met"],
        "source": "3D FDTD graded cavity scan", "scan_csv": CSV_PATH,
        # legacy keys so run_3d_cavity can still read it
        "a_c": row["a_center"],
    }
    with open(BEST_PATH, "w") as f:
        json.dump(payload, f, indent=2)
    with open(COMPAT_PATH, "w") as f:
        json.dump(payload, f, indent=2)


def _as_bool(val):
    return val is True or str(val).lower() in ("1", "true", "yes")


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
    geom, cell = build_periodic_geom_3d(mirror["a_end"], mirror["rx"], mirror["ry"], 20)
    fh, flux_h = run_flux_sim_3d(geom, cell, nfreq, resolution, fixed_until=fixed_until)
    fr, flux_r = run_flux_sim_3d(build_ref_geom_3d(), cell, nfreq, resolution,
                                 fixed_until=fixed_until)
    wl, tr = normalized_spectrum(fh, flux_h, fr, flux_r)
    bg = analyze_bandgap_3d(wl, tr)
    print(f"  3D bandgap: {bg.get('gap_start_nm',0):.0f}-{bg.get('gap_end_nm',0):.0f} nm "
          f"T_min_gap={bg.get('T_min_gap',1):.4f}")
    mirror["bandgap_start_nm"] = bg.get("gap_start_nm", mirror["bandgap_start_nm"])
    mirror["bandgap_end_nm"] = bg.get("gap_end_nm", mirror["bandgap_end_nm"])
    mirror["T_min_gap"] = bg.get("T_min_gap", mirror["T_min_gap"])
    return bg


def coarse_grid(a_end):
    """Center-period defect depths and taper lengths (mirror-free).

    Center period spans ~0.84-0.97 of a_end (shallow to moderate defect).
    """
    lo = round(0.84 * a_end, 3)
    hi = round(0.96 * a_end, 3)
    centers = [round(c, 3) for c in np.arange(lo, hi + 1e-9, 0.02)]
    tapers = [12, 18, 24]
    return centers, tapers


def main():
    parser = argparse.ArgumentParser(description="3D graded defect cavity optimizer")
    parser.add_argument("--phase", default="coarse", choices=["coarse", "refine"])
    parser.add_argument("--nfreq", type=int, default=500)
    parser.add_argument("--resolution", type=int, default=12)
    parser.add_argument("--n-mirror", type=int, default=0,
                        help="extra uniform holes after taper (0 = fully mirror-free)")
    parser.add_argument("--max-runs", type=int, default=15)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    ensure_dirs()
    mirror = load_mirror()
    a_end = mirror["a_end"]
    centers, tapers = coarse_grid(a_end)

    if args.phase == "coarse":
        grid = [(c, nt) for nt in tapers for c in centers]
    else:
        rows = list(csv.DictReader(open(CSV_PATH))) if os.path.exists(CSV_PATH) else []
        best = pick_best(rows)
        if not best:
            grid = [(c, nt) for nt in tapers for c in centers]
        else:
            c0 = float(best["a_center"])
            nt0 = int(best["N_taper"])
            cs = [round(c, 3) for c in np.arange(c0 - 0.02, c0 + 0.021, 0.01)
                  if 0.30 <= c <= a_end]
            nts = sorted({max(8, nt0 - 6), nt0, nt0 + 6, nt0 + 12})
            grid = [(c, nt) for nt in nts for c in cs]

    print(f"3D graded cavity scan: {len(grid)} points, N_mirror={args.n_mirror}, "
          f"nfreq={args.nfreq}, res={args.resolution}")
    print(f"Targets: Q>={Q_target}, lambda~1550nm, min feature {min_feature_nm}nm")
    if args.dry_run:
        for c, nt in grid[:20]:
            pos, _ = graded_hole_layout(a_end, c, nt, args.n_mirror)
            print(f"  a_center={c:.3f} N_taper={nt} holes={2*len(pos)} "
                  f"L={2*(pos[-1]+mirror['rx']):.1f}um")
        return

    bandgap = verify_mirror_3d(mirror, min(args.nfreq, 400), args.resolution)

    run_id = 0
    if os.path.exists(CSV_PATH):
        run_id = max((int(r["run_id"]) for r in csv.DictReader(open(CSV_PATH))),
                     default=0)

    best_row = None
    for c, nt in grid[: args.max_runs]:
        run_id += 1
        row = evaluate(run_id, mirror, c, nt, args.n_mirror,
                       args.nfreq, args.resolution, bandgap)
        append_csv(row)
        if float(row.get("Q", 0) or 0) > 0:
            if best_row is None or row["score"] > best_row["score"]:
                best_row = row
                save_best(best_row, mirror)
            if row.get("targets_met"):
                print(f"\n*** 3D TARGETS MET: Q={row['Q']:.0f} "
                      f"lam={row['lambda0_nm']:.1f} T={row['T_peak']:.3f} ***")
                save_best(row, mirror)
                return

    if best_row:
        print(f"\nBest 3D graded so far: Q={best_row['Q']:.0f} "
              f"lam={best_row['lambda0_nm']:.1f} T={best_row['T_peak']:.3f} "
              f"Nt={best_row['N_taper']} a_c={best_row['a_center']}")
        print(f"Saved -> {BEST_PATH}")


if __name__ == "__main__":
    main()
