#!/usr/bin/env python3
"""
3D FDTD 缺陷腔验证（可靠归一化 + 约束 Lorentzian 拟合）
======================================================
- 有孔/参考使用相同 cell（与 optimize_3d_ridge.py 一致）
- 保存原始 flux + 校验元数据，便于独立 FDTD 复现
"""
import json
import os
import sys
import time

import meep as mp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spectrum_analysis import (
    extract_peaks,
    load_bandgap_from_json,
    normalize_transmission,
    save_spectrum_bundle,
    validate_normalized_spectrum,
    lorentzian,
)

# ============================================================================
# 结构参数（与 optimize_3d_ridge / optimize_cavity_2d 一致）
# ============================================================================
n_wg = 2.18
n_air = 1.0
n_sub = 1.44
w_wg = 1.5
h_slab = 0.2
h_ridge = 0.2
h_total = 0.4

lambda_min = 1.25
lambda_max = 1.75
fcen = 1.0 / 1.50
df = (1.0 / lambda_min - 1.0 / lambda_max) / 2.0
nfreq = 2000

resolution = 20
dpml = 1.0
pad = 2.0

OUT_DIR = "results/cavity_3d"
os.makedirs(OUT_DIR, exist_ok=True)


def build_cavity_3d(mirror, cavity):
    a_m, rx_m, ry_m = mirror["a_m"], mirror["rx_m"], mirror["ry_m"]
    a_c, rx_c, ry_c = cavity["a_c"], cavity["rx_c"], cavity["ry_c"]
    Nt = cavity.get("N_taper", 3)
    Nm = cavity.get("N_mirror", 5)

    total_holes = 2 * (Nm + Nt)
    sx = 2 * dpml + max(total_holes * a_m, 4.0) + 2 * pad
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

    for i in range(Nm):
        x = -(Nm - i + Nt) * a_m
        geom.append(
            mp.Ellipsoid(
                material=mp.Medium(index=n_air),
                center=mp.Vector3(x, 0, h_total / 2),
                size=mp.Vector3(2 * rx_m, 2 * ry_m, h_total),
            )
        )
    for i in range(1, Nt + 1):
        t = i / Nt
        ai = a_c + (a_m - a_c) * t**2
        rxi = rx_c + (rx_m - rx_c) * t**2
        ryi = ry_c + (ry_m - ry_c) * t**2
        x = -(Nt - i + 0.5) * ai
        geom.append(
            mp.Ellipsoid(
                material=mp.Medium(index=n_air),
                center=mp.Vector3(x, 0, h_total / 2),
                size=mp.Vector3(2 * rxi, 2 * ryi, h_total),
            )
        )
    for i in range(1, Nt + 1):
        t = i / Nt
        ai = a_c + (a_m - a_c) * t**2
        rxi = rx_c + (rx_m - rx_c) * t**2
        ryi = ry_c + (ry_m - ry_c) * t**2
        x = (Nt - i + 0.5) * ai
        geom.append(
            mp.Ellipsoid(
                material=mp.Medium(index=n_air),
                center=mp.Vector3(x, 0, h_total / 2),
                size=mp.Vector3(2 * rxi, 2 * ryi, h_total),
            )
        )
    for i in range(Nm):
        x = (Nm - i + Nt) * a_m
        geom.append(
            mp.Ellipsoid(
                material=mp.Medium(index=n_air),
                center=mp.Vector3(x, 0, h_total / 2),
                size=mp.Vector3(2 * rx_m, 2 * ry_m, h_total),
            )
        )
    return geom, sx, sy, sz


def build_ref_geom_3d():
    """无孔参考波导（与 optimize_3d_ridge build_geom(N=0) 相同）"""
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
    return [sub, slab, ridge]


def run_cavity_fdtd(geom_hole, geom_ref, sx, sy, sz, decay=1e-5):
    """
    有孔 + 参考 FDTD，相同 cell / 源 / 监视器。
    返回 freqs, flux_h, flux_r
    """
    cell = mp.Vector3(sx, sy, sz)
    src_x = -sx / 2 + dpml + 0.5
    mon_x = sx / 2 - dpml - 0.5
    sources = [
        mp.Source(
            mp.GaussianSource(fcen, fwidth=df),
            component=mp.Ey,
            center=mp.Vector3(src_x, 0, h_total / 2),
            size=mp.Vector3(0, w_wg, h_total),
        )
    ]
    freg = mp.FluxRegion(
        center=mp.Vector3(mon_x, 0, h_total / 2),
        size=mp.Vector3(0, 2 * w_wg, 2 * h_total),
    )
    pt = mp.Vector3(mon_x, 0, h_total / 2)

    sim = mp.Simulation(
        cell_size=cell,
        resolution=resolution,
        geometry=geom_hole,
        sources=sources,
        boundary_layers=[mp.PML(dpml)],
    )
    trans = sim.add_flux(fcen, df, nfreq, freg)
    t0 = time.time()
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ey, pt, decay))
    t_hole = time.time() - t0
    freqs = np.array(mp.get_flux_freqs(trans))
    flux_h = np.array(mp.get_fluxes(trans))
    sim.reset_meep()

    sim_r = mp.Simulation(
        cell_size=cell,
        resolution=resolution,
        geometry=geom_ref,
        sources=sources,
        boundary_layers=[mp.PML(dpml)],
    )
    trans_r = sim_r.add_flux(fcen, df, nfreq, freg)
    t0 = time.time()
    sim_r.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ey, pt, decay))
    t_ref = time.time() - t0
    flux_r = np.array(mp.get_fluxes(trans_r))
    sim_r.reset_meep()

    meta = {
        "cell_sx_um": sx,
        "cell_sy_um": sy,
        "cell_sz_um": sz,
        "resolution": resolution,
        "dpml": dpml,
        "fcen": fcen,
        "df": df,
        "nfreq": nfreq,
        "source": "GaussianSource(Ey)",
        "source_size": f"(0, {w_wg}, {h_total})",
        "flux_region": f"(0, {2*w_wg}, {2*h_total})",
        "decay_threshold": decay,
        "time_hole_s": t_hole,
        "time_ref_s": t_ref,
        "normalization": "T = flux_hole / flux_ref, same cell",
    }
    return freqs, flux_h, flux_r, meta


def plot_spectrum(wl, T, bg, best_peak, cavity, out_path):
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(wl, T, "b-", lw=1.5, label="3D FDTD (normalized)")
    ax.axvspan(bg["bandgap_start_nm"], bg["bandgap_end_nm"], alpha=0.1, color="red",
               label=f"bandgap ({bg['bandgap_start_nm']:.0f}-{bg['bandgap_end_nm']:.0f} nm)")
    ax.axvline(1550, color="gray", ls=":", lw=1, alpha=0.7, label="target 1550 nm")
    if best_peak and best_peak.get("valid"):
        x0, g = best_peak["lambda0_nm"], best_peak["gamma_nm"]
        wl_f = np.linspace(x0 - 8 * g, x0 + 8 * g, 200)
        T_f = lorentzian(wl_f, x0, g, best_peak["A"], best_peak["offset"])
        ax.plot(wl_f, T_f, "r-", lw=2, alpha=0.7,
                label=f"fit Q={best_peak['Q']:.0f} λ={x0:.1f}nm")
        ax.axvline(x0, color="r", ls="--", lw=1.5)
    ax.set_xlim(1450, 1650)
    ax.set_ylim(-0.02, min(1.05, T.max() * 1.15))
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Transmission")
    ax.set_title(
        f"3D Cavity: a_c={cavity['a_c']:.4f} Nt={cavity['N_taper']} Nm={cavity['N_mirror']}"
    )
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


if __name__ == "__main__":
    print("=" * 70)
    print("3D 缺陷腔 FDTD 验证 (fixed normalization + fitting)")
    print("=" * 70)

    with open("results/best_cavity_design.json") as f:
        design = json.load(f)

    mirror = {"a_m": design["a_m"], "rx_m": design["rx_m"], "ry_m": design["ry_m"]}
    cavity = {
        "a_c": design["a_c"],
        "rx_c": design["rx_c"],
        "ry_c": design["ry_c"],
        "N_taper": design["N_taper"],
        "N_mirror": design["N_mirror"],
    }
    bg = load_bandgap_from_json()

    print(f"\n2D 设计: Q={design['Q']:.0f} T={design['T_peak']:.3f} λ={design['lambda0_nm']:.1f}nm")
    print(f"禁带 ({bg['source']}): {bg['bandgap_start_nm']:.0f}-{bg['bandgap_end_nm']:.0f} nm")

    geom_h, sx, sy, sz = build_cavity_3d(mirror, cavity)
    geom_r = build_ref_geom_3d()
    print(f"Cell (same for hole/ref): {sx:.2f} x {sy:.2f} x {sz:.2f} μm")

    print("\n运行 3D FDTD (hole + ref)...", flush=True)
    freqs, flux_h, flux_r, sim_meta = run_cavity_fdtd(geom_h, geom_r, sx, sy, sz)
    wl, T, fh, fr = normalize_transmission(freqs, flux_h, freqs, flux_r)
    delta_wl = (lambda_max - lambda_min) * 1000 / nfreq

    validation = validate_normalized_spectrum(wl, T, fh, fr, label="3D_cavity")
    print(f"\n数据校验: reliable={validation['reliable']}")
    print(f"  T range: [{validation['T_min']:.4f}, {validation['T_max']:.4f}]")
    if validation["issues"]:
        print(f"  ❌ issues: {validation['issues']}")
    if validation["warnings"]:
        print(f"  ⚠ warnings: {validation['warnings']}")

    peaks = extract_peaks(
        wl, T,
        bandgap_start=bg["bandgap_start_nm"],
        bandgap_end=bg["bandgap_end_nm"],
        target_lambda=1550.0,
        delta_wl=delta_wl,
    )
    valid_peaks = [p for p in peaks if p["valid"]]
    best = valid_peaks[0] if valid_peaks else (peaks[0] if peaks else None)

    print(f"\n检测到 {len(peaks)} 个峰，有效 {len(valid_peaks)} 个")
    for i, p in enumerate(peaks[:5]):
        flag = "✅" if p["valid"] else "❌"
        print(
            f"  {flag} #{i+1}: λ={p['lambda0_nm']:.1f} Q={p['Q']:.0f} "
            f"T={p['T_peak']:.3f} FWHM={p['FWHM_nm']:.3f} R²={p['fit_r2']:.3f}"
            + (f" ({p['invalid_reason']})" if not p["valid"] else "")
        )

    save_spectrum_bundle(
        OUT_DIR, "cavity_3d",
        wl, T, fh, fr, validation, peaks, best, sim_meta,
    )
    # 兼容旧文件名
    np.savetxt(
        f"{OUT_DIR}/cavity_3d_spectrum.csv",
        np.column_stack([wl, T]),
        delimiter=",", header="wavelength_nm,transmission", comments="",
    )
    plot_spectrum(wl, T, bg, best, cavity, f"{OUT_DIR}/cavity_3d_spectrum.png")

    result = {
        "2D_design": design,
        "3D_result": best,
        "validation": validation,
        "bandgap": bg,
        "simulation": sim_meta,
        "n_peaks_found": len(peaks),
        "n_valid_peaks": len(valid_peaks),
        "T_range": [validation["T_min"], validation["T_max"]],
        "data_reliable": validation["reliable"] and best is not None and best.get("valid"),
    }
    with open(f"{OUT_DIR}/cavity_3d_result.json", "w") as f:
        json.dump(result, f, indent=2)

    print(f"\n{'='*70}")
    if best and best.get("valid"):
        print(
            f"3D 腔模: Q={best['Q']:.0f}, λ₀={best['lambda0_nm']:.1f}nm, "
            f"T_peak={best['T_peak']:.3f}, FWHM={best['FWHM_nm']:.3f}nm"
        )
    else:
        print("⚠ 未找到禁带内可靠腔模；请检查频谱或增大 N_mirror")
    print(f"原始 flux: {OUT_DIR}/cavity_3d_flux.csv")
    print(f"分析: {OUT_DIR}/cavity_3d_analysis.json")
    print(f"{'='*70}")
