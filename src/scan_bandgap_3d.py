#!/usr/bin/env python3
"""
3D periodic bandgap scan seeded by best 2D candidate.

Workflow mirrors Meep waveguide-cavity tutorials:
1) broadband excitation
2) transmission monitor
3) reference run without holes
4) normalized transmission spectrum T = flux_holes / flux_ref
"""
import argparse
import csv
import itertools
import json
import os
import time

import meep as mp
import numpy as np

from config import (
    df,
    dpml,
    fcen,
    h_ridge,
    h_slab,
    h_total,
    lambda_max,
    lambda_min,
    n_air,
    n_sub,
    n_wg,
    nfreq_bandgap_3d,
    pad,
    resolution_3d,
    w_wg,
)

RESULTS_DIR = "results_3d"
SPECTRA_DIR = os.path.join(RESULTS_DIR, "spectra")
CSV_PATH = os.path.join(RESULTS_DIR, "bandgap_3d_scan.csv")
BEST_PATH = os.path.join(RESULTS_DIR, "best_bandgap_3d.json")
COMPAT_PATH = "results/best_bandgap_params.json"
BEST_2D_PATH = "results/best_bandgap_2d.json"


def build_geom_3d(a, rx, ry, n_period):
    sx = 2 * dpml + n_period * a + 2 * pad
    sy = 2 * dpml + w_wg + 2 * pad
    sz = 2 * dpml + h_total + 1.5

    sub = mp.Block(
        material=mp.Medium(index=n_sub),
        center=mp.Vector3(0, 0, -0.5),
        size=mp.Vector3(mp.inf, mp.inf, 1.0),
    )
    slab = mp.Block(
        material=mp.Medium(index=n_wg),
        center=mp.Vector3(0, 0, h_slab / 2),
        size=mp.Vector3(mp.inf, mp.inf, h_slab),
    )
    ridge = mp.Block(
        material=mp.Medium(index=n_wg),
        center=mp.Vector3(0, 0, h_slab + h_ridge / 2),
        size=mp.Vector3(mp.inf, w_wg, h_ridge),
    )
    geom = [sub, slab, ridge]

    x0 = -(n_period - 1) * a / 2.0
    for i in range(n_period):
        geom.append(
            mp.Ellipsoid(
                material=mp.Medium(index=n_air),
                center=mp.Vector3(x0 + i * a, 0, h_total / 2),
                size=mp.Vector3(2 * rx, 2 * ry, h_total),
            )
        )
    return geom, mp.Vector3(sx, sy, sz)


def build_ref_3d():
    return [
        mp.Block(
            material=mp.Medium(index=n_sub),
            center=mp.Vector3(0, 0, -0.5),
            size=mp.Vector3(mp.inf, mp.inf, 1.0),
        ),
        mp.Block(
            material=mp.Medium(index=n_wg),
            center=mp.Vector3(0, 0, h_slab / 2),
            size=mp.Vector3(mp.inf, mp.inf, h_slab),
        ),
        mp.Block(
            material=mp.Medium(index=n_wg),
            center=mp.Vector3(0, 0, h_slab + h_ridge / 2),
            size=mp.Vector3(mp.inf, w_wg, h_ridge),
        ),
    ]


def run_transmission(geom, cell, nfreq):
    sx = cell.x
    src_x = -sx / 2 + dpml + 0.5
    mon_x = sx / 2 - dpml - 0.5
    mon_point = mp.Vector3(mon_x, 0, h_total / 2)

    sources = [
        mp.Source(
            mp.GaussianSource(fcen, fwidth=df),
            component=mp.Ey,
            center=mp.Vector3(src_x, 0, h_total / 2),
            size=mp.Vector3(0, w_wg, h_total),
        )
    ]
    sim = mp.Simulation(
        cell_size=cell,
        resolution=resolution_3d,
        geometry=geom,
        sources=sources,
        boundary_layers=[mp.PML(dpml)],
    )
    fr = mp.FluxRegion(center=mon_point, size=mp.Vector3(0, 2 * w_wg, 2 * h_total))
    tr = sim.add_flux(fcen, df, nfreq, fr)
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ey, mon_point, 1e-4))
    freqs = np.array(mp.get_flux_freqs(tr))
    flux = np.array(mp.get_fluxes(tr))
    sim.reset_meep()
    return freqs, flux


def analyze_bandgap(wl_nm, transmission):
    mask = (wl_nm >= 1450.0) & (wl_nm <= 1650.0)
    wl = wl_nm[mask]
    tr = transmission[mask]
    if wl.size < 8:
        return {"has_gap": False, "gap_start_nm": 0.0, "gap_end_nm": 0.0, "gap_width_nm": 0.0, "T_min_gap": float(np.min(transmission))}

    best = {
        "has_gap": False,
        "gap_start_nm": 0.0,
        "gap_end_nm": 0.0,
        "gap_width_nm": 0.0,
        "T_min_gap": float(np.min(tr)),
    }
    for th in np.arange(0.03, 0.25, 0.01):
        below = tr < th
        edges = np.diff(np.concatenate(([False], below, [False])).astype(int))
        starts = np.where(edges == 1)[0]
        ends = np.where(edges == -1)[0]
        if starts.size == 0:
            continue
        lengths = ends - starts
        pick = int(np.argmax(lengths))
        if lengths[pick] < 4:
            continue
        s_idx = starts[pick]
        e_idx = ends[pick] - 1
        width = float(abs(wl[e_idx] - wl[s_idx]))
        if width <= best["gap_width_nm"]:
            continue
        s_nm = float(wl[s_idx])
        e_nm = float(wl[e_idx])
        in_gap = (wl_nm >= s_nm) & (wl_nm <= e_nm)
        best = {
            "has_gap": True,
            "gap_start_nm": s_nm,
            "gap_end_nm": e_nm,
            "gap_width_nm": width,
            "T_min_gap": float(np.min(transmission[in_gap])) if np.any(in_gap) else float(np.min(tr)),
        }
    return best


def score_result(analysis):
    # wide gap + low in-gap transmission
    return 2.0 * analysis["gap_width_nm"] - 220.0 * analysis["T_min_gap"]


def load_best_2d(path):
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found. Run `python run.py bandgap-2d` first to generate best_bandgap_2d.json."
        )
    with open(path) as f:
        payload = json.load(f)
    return {
        "a": float(payload["a"]),
        "rx": float(payload["rx"]),
        "ry": float(payload["ry"]),
        "N": int(payload["N"]),
    }


def build_scan_space(seed, args):
    a_vals = np.round(np.linspace(seed["a"] - args.a_span, seed["a"] + args.a_span, args.a_points), 4)
    rx_vals = np.round(np.linspace(seed["rx"] - args.rx_span, seed["rx"] + args.rx_span, args.rx_points), 4)
    ry_vals = np.round(np.linspace(seed["ry"] - args.ry_span, seed["ry"] + args.ry_span, args.ry_points), 4)
    n_vals = sorted(
        {
            max(8, seed["N"] - args.n_delta),
            seed["N"],
            seed["N"] + args.n_delta,
        }
    )
    return list(itertools.product(a_vals, rx_vals, ry_vals, n_vals))


def ensure_dirs():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(SPECTRA_DIR, exist_ok=True)
    os.makedirs("results", exist_ok=True)


def main():
    parser = argparse.ArgumentParser(description="3D periodic bandgap verification scan")
    parser.add_argument("--best-2d", default=BEST_2D_PATH, help="Best 2D parameter JSON")
    parser.add_argument("--a-span", type=float, default=0.01)
    parser.add_argument("--rx-span", type=float, default=0.02)
    parser.add_argument("--ry-span", type=float, default=0.03)
    parser.add_argument("--a-points", type=int, default=3)
    parser.add_argument("--rx-points", type=int, default=3)
    parser.add_argument("--ry-points", type=int, default=3)
    parser.add_argument("--n-delta", type=int, default=2)
    parser.add_argument("--nfreq", type=int, default=nfreq_bandgap_3d)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    ensure_dirs()
    seed = load_best_2d(args.best_2d)
    params = build_scan_space(seed, args)
    print(f"Seed from 2D best: a={seed['a']:.4f}, rx={seed['rx']:.4f}, ry={seed['ry']:.4f}, N={seed['N']}")
    print(f"3D scan points: {len(params)}")
    if args.dry_run:
        for item in params[:20]:
            print(f"  a={item[0]:.4f}, rx={item[1]:.4f}, ry={item[2]:.4f}, N={item[3]}")
        if len(params) > 20:
            print("  ...")
        return

    fieldnames = [
        "label",
        "a",
        "rx",
        "ry",
        "N",
        "gap_width_nm",
        "gap_start_nm",
        "gap_end_nm",
        "T_min_gap",
        "T_max",
        "score",
        "elapsed_s",
        "spectrum_csv",
    ]

    best_row = None
    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for idx, (a, rx, ry, n_period) in enumerate(params, 1):
            label = f"a{a:.4f}_rx{rx:.4f}_ry{ry:.4f}_N{n_period}"
            print(f"[{idx:03d}/{len(params)}] {label}", end="", flush=True)
            t0 = time.time()
            geom, cell = build_geom_3d(a, rx, ry, n_period)
            freqs_h, flux_h = run_transmission(geom, cell, args.nfreq)
            freqs_r, flux_r = run_transmission(build_ref_3d(), cell, args.nfreq)
            tr = np.divide(flux_h, flux_r, out=np.zeros_like(flux_h), where=flux_r > 1e-15)
            freqs = freqs_r if freqs_r.size == flux_r.size else freqs_h
            wl_nm = 1000.0 / freqs
            order = np.argsort(wl_nm)
            wl_nm = wl_nm[order]
            tr = tr[order]
            analysis = analyze_bandgap(wl_nm, tr)
            score = score_result(analysis)
            elapsed = time.time() - t0
            spectrum_path = os.path.join(SPECTRA_DIR, f"{label}.csv")
            np.savetxt(
                spectrum_path,
                np.column_stack([wl_nm, tr]),
                delimiter=",",
                header="wavelength_nm,transmission",
                comments="",
            )
            row = {
                "label": label,
                "a": a,
                "rx": rx,
                "ry": ry,
                "N": n_period,
                "gap_width_nm": round(analysis["gap_width_nm"], 3),
                "gap_start_nm": round(analysis["gap_start_nm"], 3),
                "gap_end_nm": round(analysis["gap_end_nm"], 3),
                "T_min_gap": round(analysis["T_min_gap"], 6),
                "T_max": round(float(np.max(tr)), 6),
                "score": round(score, 4),
                "elapsed_s": round(elapsed, 1),
                "spectrum_csv": spectrum_path,
            }
            writer.writerow(row)
            if best_row is None or row["score"] > best_row["score"]:
                best_row = row
            print(f" | gap={row['gap_width_nm']:.1f}nm Tmin={row['T_min_gap']:.4f} score={row['score']:.2f} [{elapsed:.1f}s]")

    if best_row is None:
        raise RuntimeError("No valid 3D scan results generated")

    best_payload = {
        "a": best_row["a"],
        "rx": best_row["rx"],
        "ry": best_row["ry"],
        "N": int(best_row["N"]),
        "gap_start_nm": best_row["gap_start_nm"],
        "gap_end_nm": best_row["gap_end_nm"],
        "gap_width_nm": best_row["gap_width_nm"],
        "center_nm": round((best_row["gap_start_nm"] + best_row["gap_end_nm"]) / 2.0, 3),
        "T_min": best_row["T_min_gap"],
        "selection_score": best_row["score"],
        "source": "3D scan seeded from best 2D candidate",
        "scan_csv": CSV_PATH,
    }
    with open(BEST_PATH, "w") as f:
        json.dump(best_payload, f, indent=2)
    with open(COMPAT_PATH, "w") as f:
        json.dump(best_payload, f, indent=2)
    print(f"\nBest 3D bandgap saved: {BEST_PATH}")
    print(f"Compatibility best params: {COMPAT_PATH}")


if __name__ == "__main__":
    main()
