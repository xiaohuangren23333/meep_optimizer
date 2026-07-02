#!/usr/bin/env python3
"""
3D FDTD 缺陷腔验证
=================
基于 2D 优化结果 (Q=113, T=0.81 在 2D)，在 3D 正确脊型波导中验证实际 Q 值。
"""
import meep as mp
import numpy as np
import os, time, json, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from config import min_feature_nm, check_cavity_geometry

# ============================================================================
# 结构参数
# ============================================================================
n_wg   = 2.18
n_air  = 1.0
n_sub  = 1.44
w_wg   = 1.5
h_slab = 0.2
h_ridge = 0.2
h_total = 0.4

lambda_min = 1.25
lambda_max = 1.75
fcen = 1.0 / 1.50
df   = (1.0/lambda_min - 1.0/lambda_max) / 2.0
nfreq = 2000

resolution = 20
dpml = 1.0
pad = 2.0

os.makedirs("results/cavity_3d", exist_ok=True)

# ============================================================================
# Lorentzian 拟合
# ============================================================================
def lorentzian(x, x0, gamma, A, offset):
    return offset + A * gamma**2 / ((x - x0)**2 + gamma**2)

# ============================================================================
# 3D 缺陷腔构建
# ============================================================================
def build_cavity_3d(mirror, cavity):
    """构建 3D 二次渐变缺陷腔"""
    a_m = mirror["a_m"]; rx_m = mirror["rx_m"]; ry_m = mirror["ry_m"]
    a_c = cavity["a_c"]; rx_c = cavity["rx_c"]; ry_c = cavity["ry_c"]
    Nt = cavity.get("N_taper", 3)
    Nm = cavity.get("N_mirror", 5)

    total_holes = 2 * (Nm + Nt)
    sx = 2*dpml + max(total_holes * a_m, 4.0) + 2*pad
    sy = 2*dpml + w_wg + 2*pad
    sz = 2*dpml + h_total + 1.5

    # 衬底
    sub = mp.Block(material=mp.Medium(index=n_sub),
                   center=mp.Vector3(0,0,-0.5),
                   size=mp.Vector3(mp.inf, mp.inf, 1.0))
    # 平板层
    slab = mp.Block(material=mp.Medium(index=n_wg),
                    center=mp.Vector3(0,0,h_slab/2),
                    size=mp.Vector3(mp.inf, mp.inf, h_slab))
    # 脊型层
    ridge = mp.Block(material=mp.Medium(index=n_wg),
                     center=mp.Vector3(0,0,h_slab+h_ridge/2),
                     size=mp.Vector3(mp.inf, w_wg, h_ridge))

    geom = [sub, slab, ridge]

    # 左侧: 镜区 → 渐变区
    for i in range(Nm):
        x = -(Nm - i + Nt) * a_m
        geom.append(mp.Ellipsoid(material=mp.Medium(index=n_air),
                                center=mp.Vector3(x, 0, h_total/2),
                                size=mp.Vector3(2*rx_m, 2*ry_m, h_total)))
    for i in range(1, Nt + 1):
        t = i / Nt
        ai = a_c + (a_m - a_c) * t**2
        rxi = rx_c + (rx_m - rx_c) * t**2
        ryi = ry_c + (ry_m - ry_c) * t**2
        x = -(Nt - i + 0.5) * ai
        geom.append(mp.Ellipsoid(material=mp.Medium(index=n_air),
                                center=mp.Vector3(x, 0, h_total/2),
                                size=mp.Vector3(2*rxi, 2*ryi, h_total)))

    # 右侧: 渐变区 → 镜区
    for i in range(1, Nt + 1):
        t = i / Nt
        ai = a_c + (a_m - a_c) * t**2
        rxi = rx_c + (rx_m - rx_c) * t**2
        ryi = ry_c + (ry_m - ry_c) * t**2
        x = (Nt - i + 0.5) * ai
        geom.append(mp.Ellipsoid(material=mp.Medium(index=n_air),
                                center=mp.Vector3(x, 0, h_total/2),
                                size=mp.Vector3(2*rxi, 2*ryi, h_total)))
    for i in range(Nm):
        x = (Nm - i + Nt) * a_m
        geom.append(mp.Ellipsoid(material=mp.Medium(index=n_air),
                                center=mp.Vector3(x, 0, h_total/2),
                                size=mp.Vector3(2*rx_m, 2*ry_m, h_total)))

    return geom, sx, sy, sz


def build_ref_geom():
    """参考波导（无孔）"""
    sub = mp.Block(material=mp.Medium(index=n_sub),
                   center=mp.Vector3(0,0,-0.5),
                   size=mp.Vector3(mp.inf, mp.inf, 1.0))
    slab = mp.Block(material=mp.Medium(index=n_wg),
                    center=mp.Vector3(0,0,h_slab/2),
                    size=mp.Vector3(mp.inf, mp.inf, h_slab))
    ridge = mp.Block(material=mp.Medium(index=n_wg),
                     center=mp.Vector3(0,0,h_slab+h_ridge/2),
                     size=mp.Vector3(mp.inf, w_wg, h_ridge))
    return [sub, slab, ridge]


# ============================================================================
# 仿真函数
# ============================================================================
def run_sim(geom, cell_sx, cell_sy, cell_sz):
    """运行单次 3D FDTD，返回波长和透射率"""
    cell = mp.Vector3(cell_sx, cell_sy, cell_sz)
    src_x = -cell_sx/2 + dpml + 0.5
    mon_x = cell_sx/2 - dpml - 0.5

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ey,
                         center=mp.Vector3(src_x, 0, h_total/2),
                         size=mp.Vector3(0, w_wg, h_total))]

    sim = mp.Simulation(cell_size=cell, resolution=resolution,
                        geometry=geom, sources=sources,
                        boundary_layers=[mp.PML(dpml)])

    freg = mp.FluxRegion(center=mp.Vector3(mon_x,0,h_total/2),
                          size=mp.Vector3(0, 2*w_wg, 2*h_total))
    trans = sim.add_flux(fcen, df, nfreq, freg)

    # 腔需要更长的运行时间
    sim.run(until_after_sources=mp.stop_when_fields_decayed(
        50, mp.Ey, mp.Vector3(mon_x,0,h_total/2), 1e-5))

    freqs = np.array(mp.get_flux_freqs(trans))
    flux = np.array(mp.get_fluxes(trans))
    sim.reset_meep()
    return freqs, flux


# ============================================================================
# 主程序
# ============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("3D 缺陷腔 FDTD 验证")
    print("=" * 70)

    # 读取 3D 优化最佳设计（优先）或 2D 回退
    design_path = "results_3d/cavity/best_cavity_3d.json"
    if not os.path.exists(design_path):
        design_path = "results/best_cavity_design.json"
    with open(design_path) as f:
        design = json.load(f)

    mirror = {
        "a_m": design["a_m"], "rx_m": design["rx_m"], "ry_m": design["ry_m"]
    }
    cavity = {
        "a_c": design["a_c"], "rx_c": design["rx_c"], "ry_c": design["ry_c"],
        "N_taper": design["N_taper"], "N_mirror": design["N_mirror"]
    }

    print(f"\n设计来源: {design_path}")
    print(f"  镜区: a_m={mirror['a_m']:.3f} rx_m={mirror['rx_m']:.3f} ry_m={mirror['ry_m']:.3f}")
    print(f"  缺陷: a_c={cavity['a_c']:.4f} rx_c={cavity['rx_c']:.3f} ry_c={cavity['ry_c']:.3f}")
    print(f"  N_taper={cavity['N_taper']} N_mirror={cavity['N_mirror']}")
    print(f"  3D/2D 预测: Q={design.get('Q',0):.0f} T_peak={design.get('T_peak',0):.3f} λ₀={design.get('lambda0_nm',0):.1f}nm")

    ok, violations = check_cavity_geometry(
        mirror["a_m"], mirror["rx_m"], mirror["ry_m"],
        cavity["a_c"], cavity["rx_c"], cavity["ry_c"],
        cavity["N_taper"],
    )
    if not ok:
        print(f"\n⚠️ 设计不满足最小特征尺寸 {min_feature_nm}nm:")
        for v in violations:
            print(f"   - {v}")
        print("   建议先重新运行 bandgap-3d / cavity-2d 优化。")
        sys.exit(1)

    # 构建几何
    geom, sx, sy, sz = build_cavity_3d(mirror, cavity)
    print(f"\n3D 计算区域: {sx:.1f} x {sy:.1f} x {sz:.1f}")

    # 有孔腔仿真
    t0 = time.time()
    print(f"\n运行 3D 腔仿真...", end="", flush=True)
    freqs_h, flux_h = run_sim(geom, sx, sy, sz)
    print(f" {time.time()-t0:.1f}s")

    # 参考波导 (无孔): 使用完全相同 cell，避免归一化偏差
    ref_sx = sx
    ref_geom = build_ref_geom()
    t0 = time.time()
    print(f"运行 3D 参考波导...", end="", flush=True)
    freqs_r, flux_r = run_sim(ref_geom, ref_sx, sy, sz)
    print(f" {time.time()-t0:.1f}s")

    # 归一化
    T = np.divide(flux_h, flux_r, out=np.zeros_like(flux_h), where=flux_r > 1e-15)
    wl = 1000.0 / freqs_r

    # 排序
    idx = np.argsort(wl)
    wl, T = wl[idx], T[idx]

    # 保存 CSV
    csv_path = "results/cavity_3d/cavity_3d_spectrum.csv"
    np.savetxt(csv_path, np.column_stack([wl, T]),
               delimiter=",", header="wavelength_nm,transmission", comments="")

    # 找峰
    from scipy.signal import find_peaks
    mask = (wl >= 1450) & (wl <= 1650)
    wl_focus, T_focus = wl[mask], T[mask]
    peaks, props = find_peaks(T_focus, height=0.01, prominence=0.005, width=3)

    print(f"\n{'='*50}")
    print(f"3D 缺陷腔结果")
    print(f"{'='*50}")
    print(f"  T 范围: [{T.min():.4f}, {T.max():.4f}]")
    print(f"  找到 {len(peaks)} 个峰")

    best_fit = None
    if len(peaks) > 0:
        for k, idx in enumerate(peaks):
            left = max(idx - 15, 0)
            right = min(idx + 15, len(wl_focus) - 1)
            wl_seg, T_seg = wl_focus[left:right+1], T_focus[left:right+1]
            if len(wl_seg) < 7:
                continue
            try:
                x0_g = wl_focus[idx]
                A_g = T_focus[idx] - T_focus[min(left + 3, len(wl_focus)-1)]
                gamma_g = 5.0
                popt, _ = curve_fit(lorentzian, wl_seg, T_seg,
                                    p0=[x0_g, gamma_g, A_g, min(T_seg)],
                                    maxfev=5000)
                x0, gamma, A, offset = popt
                Q = x0 / (2*gamma)
                T_peak = lorentzian(x0, *popt)
                residuals = T_seg - lorentzian(wl_seg, *popt)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((T_seg - np.mean(T_seg))**2)
                r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
                print(f"\n  峰 #{k+1}:")
                print(f"    λ₀ = {x0:.2f} nm")
                print(f"    FWHM = {2*gamma:.3f} nm")
                print(f"    Q = {Q:.0f}")
                print(f"    T_peak = {T_peak:.4f}")
                print(f"    R² = {r2:.4f}")
                best_fit = {"lambda0_nm": float(x0), "FWHM_nm": float(2*gamma),
                           "Q": float(Q), "T_peak": float(T_peak), "fit_r2": float(r2)}
            except Exception as e:
                print(f"  峰 #{k+1}: 拟合失败 ({e})")

    # 绘图
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(wl, T, 'b-', lw=1.5, label='3D FDTD')
    ax.axvspan(1530, 1575, alpha=0.08, color='red', label='mirror bandgap')
    if best_fit:
        ax.axvline(best_fit["lambda0_nm"], color='r', ls='--', lw=2,
                   label=f"Q={best_fit['Q']:.0f}, λ₀={best_fit['lambda0_nm']:.1f}nm")
        # 绘制拟合曲线
        wl_f = np.linspace(best_fit["lambda0_nm"] - 20, best_fit["lambda0_nm"] + 20, 200)
        T_f = lorentzian(wl_f, best_fit["lambda0_nm"], best_fit["FWHM_nm"]/2,
                        best_fit["T_peak"] - min(T_focus), min(T_focus))
        ax.plot(wl_f, T_f, 'r-', lw=3, alpha=0.4)
    ax.set_xlim(1450, 1650)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Transmission")
    ax.set_title(f"3D Cavity: a_c={cavity['a_c']:.4f} rx_c={cavity['rx_c']:.3f} ry_c={cavity['ry_c']:.3f}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    png_path = "results/cavity_3d/cavity_3d_spectrum.png"
    plt.savefig(png_path, dpi=150)
    plt.close()
    print(f"\n✅ 频谱图: {png_path}")

    # 保存结果 JSON
    result = {
        "2D_design": design,
        "3D_result": best_fit,
        "T_range": [float(T.min()), float(T.max())],
        "T_max_wavelength_nm": float(wl[T.argmax()]),
    }
    with open("results/cavity_3d/cavity_3d_result.json", "w") as f:
        json.dump(result, f, indent=2)

    print(f"\n{'='*70}")
    if best_fit:
        print(f"🏆 3D 缺陷腔验证: Q = {best_fit['Q']:.0f}, λ₀ = {best_fit['lambda0_nm']:.1f}nm, T_peak = {best_fit['T_peak']:.3f}")
    else:
        print(f"⚠️ 未检测到腔模")
    print(f"{'='*70}")