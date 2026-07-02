#!/usr/bin/env python3
"""
精细扫描 a_c，对准 1550 nm（使用可靠拟合模块）
================================================
在 2D 中快速验证，固定其余腔参数，只扫 a_c。
"""
import csv
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import optimize_cavity_2d as oc
from spectrum_analysis import extract_peaks, load_bandgap_from_json, save_spectrum_bundle, validate_normalized_spectrum

OUT_CSV = "results/cavity_ac_scan.csv"
OUT_JSON = "results/best_cavity_design.json"


def run_one(a_c, mirror, params):
    rx_c, ry_c, Nt, Nm = params
    geom, sx = oc.build_cavity_geom(
        mirror["a_m"], mirror["rx_m"], mirror["ry_m"],
        a_c, rx_c, ry_c, Nt, Nm,
    )
    wl, T = oc.run_cavity_sim(geom, sx)
    delta_wl = (oc.lambda_max - oc.lambda_min) * 1000 / oc.nfreq
    bg = load_bandgap_from_json()
    val = validate_normalized_spectrum(wl, T, label=f"a_c={a_c}")
    peaks = extract_peaks(
        wl, T,
        bandgap_start=bg["bandgap_start_nm"],
        bandgap_end=bg["bandgap_end_nm"],
        target_lambda=1550.0,
        delta_wl=delta_wl,
    )
    valid = [p for p in peaks if p["valid"]]
    best = valid[0] if valid else None
    return wl, T, val, peaks, best


def main():
    with open(OUT_JSON) as f:
        base = json.load(f)

    mirror = load_bandgap_from_json()
    mirror.update({
        "a_m": base["a_m"], "rx_m": base["rx_m"], "ry_m": base["ry_m"],
        "bandgap_start_nm": mirror["bandgap_start_nm"],
        "bandgap_end_nm": mirror["bandgap_end_nm"],
    })
    params = (base["rx_c"], base["ry_c"], base["N_taper"], base["N_mirror"])

    # 以当前最佳 a_c 为中心，向 1550 nm 方向扫描
    center = base["a_c"]
    a_vals = np.round(np.linspace(center - 0.025, center + 0.025, 11), 4)

    header = (
        "a_c,lambda0_nm,FWHM_nm,Q,T_peak,fit_r2,valid,reliable,"
        "T_min,T_max,issues,warnings,elapsed_s\n"
    )
    write_header = not os.path.isfile(OUT_CSV)
    rows = []

    print("=" * 60)
    print(f"a_c 精细扫描 (rx_c={params[0]} ry_c={params[1]} Nt={params[2]} Nm={params[3]})")
    print("=" * 60)

    for a_c in a_vals:
        t0 = time.time()
        print(f"\n▶ a_c={a_c:.4f}", end="", flush=True)
        wl, T, val, peaks, best = run_one(a_c, mirror, params)
        elapsed = time.time() - t0

        if best:
            print(
                f" | Q={best['Q']:.0f} T={best['T_peak']:.3f} "
                f"λ={best['lambda0_nm']:.1f} R²={best['fit_r2']:.3f} ({elapsed:.0f}s)"
            )
        else:
            print(f" | no valid peak ({elapsed:.0f}s)")

        os.makedirs("results/spectra", exist_ok=True)
        tag = f"ac{a_c:.4f}".replace(".", "p")
        np.savetxt(
            f"results/spectra/scan_{tag}.csv",
            np.column_stack([wl, T]),
            delimiter=",", header="wavelength_nm,transmission", comments="",
        )

        row = {
            "a_c": a_c,
            "best": best,
            "val": val,
            "elapsed": elapsed,
        }
        rows.append(row)

        with open(OUT_CSV, "a") as f:
            if write_header:
                f.write(header)
                write_header = False
            if best:
                f.write(
                    f"{a_c},{best['lambda0_nm']:.2f},{best['FWHM_nm']:.4f},"
                    f"{best['Q']:.0f},{best['T_peak']:.4f},{best['fit_r2']:.4f},"
                    f"1,{val['reliable']},{val['T_min']:.4f},{val['T_max']:.4f},"
                    f"\"{';'.join(val['issues'])}\",\"{';'.join(val['warnings'])}\","
                    f"{elapsed:.1f}\n"
                )
            else:
                f.write(
                    f"{a_c},0,0,0,0,0,0,{val['reliable']},"
                    f"{val['T_min']:.4f},{val['T_max']:.4f},"
                    f"\"{';'.join(val['issues'])}\",\"{';'.join(val['warnings'])}\","
                    f"{elapsed:.1f}\n"
                )

    scored = []
    for r in rows:
        b = r["best"]
        if not b:
            continue
        dl = abs(b["lambda0_nm"] - 1550)
        s = (
            -dl * 2
            + min(b["Q"] / 50, 30)
            + min(b["T_peak"] * 20, 20)
            + b["fit_r2"] * 10
            + (50 if b["Q"] >= 1000 else 0)
        )
        scored.append((s, r["a_c"], b))

    scored.sort(key=lambda x: -x[0])
    if scored:
        _, best_ac, b = scored[0]
        design = {**base, "a_c": float(best_ac)}
        for k in ("lambda0_nm", "FWHM_nm", "Q", "T_peak", "fit_r2"):
            design[k] = b[k]
        design["scan_source"] = "scan_cavity_ac.py"
        with open(OUT_JSON, "w") as f:
            json.dump(design, f, indent=2)
        print(f"\n🏆 更新 best_cavity_design.json: a_c={design['a_c']:.4f}")
        print(f"   Q={design['Q']:.0f} λ={design['lambda0_nm']:.1f}nm T={design['T_peak']:.3f}")

    print(f"\n结果: {OUT_CSV}")


if __name__ == "__main__":
    main()
