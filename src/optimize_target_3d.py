#!/usr/bin/env python3
"""
3D 目标驱动缺陷腔优化调度器
=========================
在 2D 缺陷腔扫描结果基础上，批量执行 3D FDTD 验证，并按目标函数排序：
  - Q >= 1000
  - lambda0 接近 1550nm
  - T_peak 接近 1
  - 谐振峰附近局部谷值接近 0
  - 镜区禁带宽度尽量大
"""
import argparse
import csv
import json
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

import run_3d_cavity as c3d


RESULT_DIR = "results/cavity_3d"
SUMMARY_JSON = os.path.join(RESULT_DIR, "target_optimization_summary.json")
BEST_JSON = os.path.join(RESULT_DIR, "best_target_design_3d.json")
SUMMARY_CSV = os.path.join(RESULT_DIR, "target_optimization_summary.csv")


def load_mirror() -> Dict:
    for p in ["results/best_bandgap_params.json", "results_3d/optimized_params.json"]:
        if os.path.isfile(p):
            with open(p) as f:
                m = json.load(f)
            gs = float(m.get("gap_start_nm", m.get("bandgap_start_nm", 0)) or 0)
            ge = float(m.get("gap_end_nm", m.get("bandgap_end_nm", 0)) or 0)
            gw = float(m.get("gap_width_nm", max(0.0, ge - gs)) or 0)
            return {
                "a_m": float(m.get("a_m", m.get("a", 0.44))),
                "rx_m": float(m.get("rx_m", m.get("rx", 0.12))),
                "ry_m": float(m.get("ry_m", m.get("ry", 0.22))),
                "gap_start_nm": gs,
                "gap_end_nm": ge,
                "gap_width_nm": gw,
                "source": p,
            }
    return {
        "a_m": 0.44,
        "rx_m": 0.12,
        "ry_m": 0.22,
        "gap_start_nm": 1529.0,
        "gap_end_nm": 1589.0,
        "gap_width_nm": 60.0,
        "source": "fallback-default",
    }


def _f(x: str, d: float = 0.0) -> float:
    try:
        return float(x)
    except Exception:
        return d


def _i(x: str, d: int = 0) -> int:
    try:
        return int(float(x))
    except Exception:
        return d


def load_2d_candidates(path: str, top_k: int) -> List[Dict]:
    out: List[Dict] = []
    if os.path.isfile(path):
        with open(path) as f:
            rd = csv.DictReader(f)
            for r in rd:
                out.append({
                    "source": "csv",
                    "run_id": r.get("run_id", ""),
                    "a_c": _f(r.get("a_c", ""), 0.44),
                    "rx_c": _f(r.get("rx_c", ""), 0.12),
                    "ry_c": _f(r.get("ry_c", ""), 0.22),
                    "N_taper": _i(r.get("N_taper", ""), 3),
                    "N_mirror": _i(r.get("N_mirror", ""), 10),
                    "q_2d": _f(r.get("Q", ""), 0),
                    "t_peak_2d": _f(r.get("T_peak", ""), 0),
                    "score_2d": _f(r.get("best_peak_score", ""), -100),
                    "t_floor_local_2d": _f(r.get("T_floor_local", ""), 1.0),
                })
    if not out and os.path.isfile("results/best_cavity_design.json"):
        with open("results/best_cavity_design.json") as f:
            b = json.load(f)
        out.append({
            "source": "best_cavity_design.json",
            "run_id": "best",
            "a_c": float(b.get("a_c", 0.44)),
            "rx_c": float(b.get("rx_c", 0.12)),
            "ry_c": float(b.get("ry_c", 0.22)),
            "N_taper": int(b.get("N_taper", 3)),
            "N_mirror": int(b.get("N_mirror", 10)),
            "q_2d": float(b.get("Q", 0)),
            "t_peak_2d": float(b.get("T_peak", 0)),
            "score_2d": float(b.get("target_score_2d", -100)),
            "t_floor_local_2d": float(b.get("T_floor_local", 1.0)),
        })
    out.sort(key=lambda x: (x["score_2d"], x["q_2d"], x["t_peak_2d"]), reverse=True)
    return out[:top_k]


def lorentzian(x, x0, gamma, A, offset):
    return offset + A * gamma**2 / ((x - x0)**2 + gamma**2)


def fit_peak(wl_seg: np.ndarray, t_seg: np.ndarray) -> Optional[Dict]:
    if len(wl_seg) < 7:
        return None
    try:
        x0_guess = float(wl_seg[np.argmax(t_seg)])
        A_guess = float(np.max(t_seg) - np.min(t_seg))
        gamma_guess = float((wl_seg[-1] - wl_seg[0]) / 4)
        offset_guess = float(np.min(t_seg))
        popt, _ = curve_fit(
            lorentzian, wl_seg, t_seg,
            p0=[x0_guess, gamma_guess, A_guess, offset_guess],
            maxfev=5000
        )
        x0, gamma, _, _ = popt
        fwhm = float(2 * gamma)
        q = float(x0 / fwhm) if fwhm > 0 else 0.0
        tpk = float(lorentzian(x0, *popt))
        res = t_seg - lorentzian(wl_seg, *popt)
        ss_res = float(np.sum(res**2))
        ss_tot = float(np.sum((t_seg - np.mean(t_seg))**2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        return {"lambda0_nm": float(x0), "FWHM_nm": fwhm, "Q": q, "T_peak": tpk, "fit_r2": float(r2)}
    except Exception:
        return None


def local_floor(wl: np.ndarray, t: np.ndarray, x0: float, fwhm: float) -> float:
    half_span_nm = 30.0
    exclude_nm = max(1.5 * fwhm, 0.5)
    m_local = (wl >= x0 - half_span_nm) & (wl <= x0 + half_span_nm)
    m_valley = m_local & ((wl < x0 - exclude_nm) | (wl > x0 + exclude_nm))
    if np.any(m_valley):
        return float(np.min(t[m_valley]))
    if np.any(m_local):
        return float(np.min(t[m_local]))
    return float(np.min(t))


def score_3d_target(peak: Dict, mirror_gap_width_nm: float) -> float:
    q = float(peak["Q"])
    lam = float(peak["lambda0_nm"])
    tpk = float(peak["T_peak"])
    tfl = float(peak["T_floor_local"])
    r2 = float(peak["fit_r2"])

    score = 0.0
    if q >= 1000:
        score += 150
    score += min(q / 8.0, 120)
    score += max(0.0, 80 - 2.0 * abs(lam - 1550.0))
    score += max(0.0, 60 - 180.0 * abs(tpk - 1.0))
    score += max(0.0, 80 - 300.0 * tfl)
    score += min(max(mirror_gap_width_nm, 0.0), 120) * 0.35
    score += 20.0 * r2
    if lam < 1500 or lam > 1600:
        score -= 80
    return float(score)


def target_hit(peak: Dict, mirror_gap_width_nm: float) -> bool:
    return (
        peak["Q"] >= 1000
        and abs(peak["lambda0_nm"] - 1550) <= 20
        and peak["T_peak"] >= 0.85
        and peak["T_floor_local"] <= 0.10
        and mirror_gap_width_nm >= 50
    )


def evaluate_design(mirror: Dict, cand: Dict, quick: bool = False) -> Dict:
    if quick:
        c3d.nfreq = 1200
        c3d.resolution = 16

    cavity = {
        "a_c": cand["a_c"],
        "rx_c": cand["rx_c"],
        "ry_c": cand["ry_c"],
        "N_taper": cand["N_taper"],
        "N_mirror": cand["N_mirror"],
    }

    geom, sx, sy, sz = c3d.build_cavity_3d(mirror, cavity)
    freqs_h, flux_h = c3d.run_sim(geom, sx, sy, sz)
    ref_geom = c3d.build_ref_geom()
    ref_sx = 2 * c3d.dpml + 4.0 + 2 * c3d.pad
    freqs_r, flux_r = c3d.run_sim(ref_geom, ref_sx, sy, sz)

    T = np.divide(flux_h, flux_r, out=np.zeros_like(flux_h), where=flux_r > 1e-15)
    wl = 1000.0 / freqs_r
    idx = np.argsort(wl)
    wl = wl[idx]
    T = T[idx]

    mask = (wl >= 1450) & (wl <= 1650)
    wl_f = wl[mask]
    T_f = T[mask]
    peaks, _ = find_peaks(T_f, height=0.01, prominence=0.005, width=3)
    peak_rows = []
    for ip in peaks:
        left = max(ip - 15, 0)
        right = min(ip + 15, len(wl_f) - 1)
        fit = fit_peak(wl_f[left:right + 1], T_f[left:right + 1])
        if fit is None:
            continue
        fit["T_floor_local"] = local_floor(wl_f, T_f, fit["lambda0_nm"], fit["FWHM_nm"])
        fit["score_3d"] = score_3d_target(fit, mirror["gap_width_nm"])
        fit["target_hit"] = target_hit(fit, mirror["gap_width_nm"])
        peak_rows.append(fit)

    best_peak = max(peak_rows, key=lambda x: x["score_3d"]) if peak_rows else None
    return {
        "candidate": cand,
        "n_peaks": len(peak_rows),
        "best_peak_3d": best_peak,
        "all_peaks_3d": peak_rows,
        "T_min_full": float(np.min(T)),
        "T_max_full": float(np.max(T)),
        "mirror_gap_width_nm": float(mirror["gap_width_nm"]),
    }


def write_summary_csv(rows: List[Dict]) -> None:
    with open(SUMMARY_CSV, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow([
            "idx", "source", "run_id", "a_c", "rx_c", "ry_c", "N_taper", "N_mirror",
            "Q_3d", "lambda0_nm_3d", "T_peak_3d", "T_floor_local_3d",
            "fit_r2_3d", "score_3d", "target_hit"
        ])
        for i, r in enumerate(rows, 1):
            c = r["candidate"]
            b = r.get("best_peak_3d") or {}
            wr.writerow([
                i, c.get("source", ""), c.get("run_id", ""),
                c.get("a_c", 0), c.get("rx_c", 0), c.get("ry_c", 0),
                c.get("N_taper", 0), c.get("N_mirror", 0),
                b.get("Q", 0), b.get("lambda0_nm", 0), b.get("T_peak", 0),
                b.get("T_floor_local", 0), b.get("fit_r2", 0),
                b.get("score_3d", -1e9), b.get("target_hit", False)
            ])


def main():
    parser = argparse.ArgumentParser(description="3D target-driven cavity optimization")
    parser.add_argument("--best-csv", default="results/cavity_optimization_best.csv")
    parser.add_argument("--top-k", type=int, default=3)
    parser.add_argument("--quick", action="store_true", help="lower resolution/nfreq for faster iteration")
    parser.add_argument("--stop-on-target", action="store_true", help="stop once target is met")
    parser.add_argument("--dry-run", action="store_true", help="show selected candidates only")
    args = parser.parse_args()

    os.makedirs(RESULT_DIR, exist_ok=True)
    mirror = load_mirror()
    candidates = load_2d_candidates(args.best_csv, args.top_k)
    if not candidates:
        raise RuntimeError("No 2D candidates found. Run optimize_cavity_2d.py first.")

    print("=" * 70)
    print("3D 目标驱动优化")
    print("=" * 70)
    print(f"mirror from: {mirror['source']}")
    print(f"mirror gap: {mirror['gap_start_nm']:.1f}-{mirror['gap_end_nm']:.1f} nm "
          f"(width={mirror['gap_width_nm']:.1f} nm)")
    print(f"selected candidates: {len(candidates)}")
    for i, c in enumerate(candidates, 1):
        print(f"  [{i}] a_c={c['a_c']:.4f} rx_c={c['rx_c']:.4f} ry_c={c['ry_c']:.4f} "
              f"Nt={c['N_taper']} Nm={c['N_mirror']} | 2D Q={c['q_2d']:.0f} score={c['score_2d']:.1f}")
    if args.dry_run:
        return

    results: List[Dict] = []
    for i, c in enumerate(candidates, 1):
        print(f"\n[{i}/{len(candidates)}] run 3D candidate...", flush=True)
        res = evaluate_design(mirror, c, quick=args.quick)
        results.append(res)
        bp = res.get("best_peak_3d")
        if bp:
            print(f"  -> Q={bp['Q']:.0f} Tpk={bp['T_peak']:.3f} Tfloor={bp['T_floor_local']:.3f} "
                  f"lambda={bp['lambda0_nm']:.1f}nm score={bp['score_3d']:.1f} "
                  f"hit={bp['target_hit']}")
        else:
            print("  -> no valid 3D peak")
        if args.stop_on_target and bp and bp["target_hit"]:
            print("target met, stop early.")
            break

    sortable = [r for r in results if r.get("best_peak_3d")]
    sortable.sort(key=lambda x: x["best_peak_3d"]["score_3d"], reverse=True)

    summary = {
        "mirror": mirror,
        "config": {
            "best_csv": args.best_csv,
            "top_k": args.top_k,
            "quick": args.quick,
            "stop_on_target": args.stop_on_target,
        },
        "n_evaluated": len(results),
        "n_with_peak": len(sortable),
        "best_result": sortable[0] if sortable else None,
        "results": results,
    }
    with open(SUMMARY_JSON, "w") as f:
        json.dump(summary, f, indent=2)
    if sortable:
        with open(BEST_JSON, "w") as f:
            json.dump(sortable[0], f, indent=2)
    write_summary_csv(results)

    print("\n" + "=" * 70)
    if sortable:
        b = sortable[0]["best_peak_3d"]
        c = sortable[0]["candidate"]
        print("best 3D candidate:")
        print(f"  a_c={c['a_c']:.4f} rx_c={c['rx_c']:.4f} ry_c={c['ry_c']:.4f} Nt={c['N_taper']} Nm={c['N_mirror']}")
        print(f"  Q={b['Q']:.0f} lambda={b['lambda0_nm']:.2f}nm Tpk={b['T_peak']:.3f} "
              f"Tfloor={b['T_floor_local']:.3f} score={b['score_3d']:.1f}")
        print(f"  target_hit={b['target_hit']}")
    else:
        print("no valid 3D peak found.")
    print(f"summary: {SUMMARY_JSON}")
    print(f"best: {BEST_JSON}")
    print("=" * 70)


if __name__ == "__main__":
    main()
