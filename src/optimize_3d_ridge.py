#!/usr/bin/env python3
"""
3D 正确脊型波导自动优化器
=========================
从最佳起点 (a=0.44, rx=0.12, ry=0.22, N=20) 出发，
使用坐标下降法自动搜索最优参数，最大化禁带宽度并对准 1550nm。

策略:
  1. 固定 N=20
  2. 扫描 a ∈ [0.43, 0.45] → 选择中心最接近 1550nm 且 gap 最宽的组合
  3. 固定最佳 a, 扫描 ry ∈ [0.20, 0.24] → 最大化宽度
  4. 固定最佳 a, ry, 扫描 rx ∈ [0.10, 0.14] → 微调宽度
  5. 迭代收敛

目标函数: f = gap_width * penalty_center
  penalty_center = exp(-((center - 1550)/20)^2)
"""
import sys, os, time, json, importlib.util

import meep as mp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from config import (
    n_wg, n_air, n_sub, w_wg, h_slab, h_ridge, h_total,
    min_feature_nm, check_periodic_geometry,
)

lambda_min = 1.25
lambda_max = 1.75
fcen = 1.0 / 1.50
df   = (1.0/lambda_min - 1.0/lambda_max) / 2.0
nfreq = 800

resolution = 20
dpml = 1.0
pad = 2.0

os.makedirs("results_3d/spectra", exist_ok=True)
os.makedirs("results_3d/figures", exist_ok=True)


# ============================================================================
# 几何构建
# ============================================================================
def build_geom(a, rx, ry, N):
    sx = 2*dpml + N*a + 2*pad
    sy = 2*dpml + w_wg + 2*pad
    sz = 2*dpml + h_total + 1.5

    sub = mp.Block(material=mp.Medium(index=n_sub),
                   center=mp.Vector3(0, 0, -0.5),
                   size=mp.Vector3(mp.inf, mp.inf, 1.0))
    slab = mp.Block(material=mp.Medium(index=n_wg),
                    center=mp.Vector3(0, 0, h_slab/2),
                    size=mp.Vector3(mp.inf, mp.inf, h_slab))
    ridge = mp.Block(material=mp.Medium(index=n_wg),
                     center=mp.Vector3(0, 0, h_slab + h_ridge/2),
                     size=mp.Vector3(mp.inf, w_wg, h_ridge))
    geom = [sub, slab, ridge]

    if N > 0:
        hole_start = -(N-1)*a / 2.0
        for i in range(N):
            geom.append(mp.Ellipsoid(
                material=mp.Medium(index=n_air),
                center=mp.Vector3(hole_start + i*a, 0, h_total/2),
                size=mp.Vector3(2*rx, 2*ry, h_total)))
    return geom, sx, sy, sz


# ============================================================================
# 单次 3D 仿真 (返回波长和透射率)
# ============================================================================
def run_sim(a, rx, ry, N, label=""):
    """运行一次仿真，返回 (wl_nm, transmission)"""
    geom_h, sx, sy, sz = build_geom(a, rx, ry, N)
    cell = mp.Vector3(sx, sy, sz)
    src_x = -sx/2 + dpml + 0.5
    mon_x = sx/2 - dpml - 0.5

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ey,
                         center=mp.Vector3(src_x, 0, h_total/2),
                         size=mp.Vector3(0, w_wg, h_total))]

    # ========== 有孔仿真 ==========
    sim = mp.Simulation(cell_size=cell, resolution=resolution,
                        geometry=geom_h, sources=sources,
                        boundary_layers=[mp.PML(dpml)])
    freg = mp.FluxRegion(center=mp.Vector3(mon_x,0,h_total/2),
                          size=mp.Vector3(0, 2*w_wg, 2*h_total))
    trans = sim.add_flux(fcen, df, nfreq, freg)
    sim.run(until_after_sources=mp.stop_when_fields_decayed(
        50, mp.Ey, mp.Vector3(mon_x,0,h_total/2), 1e-4))
    freqs_h = np.array(mp.get_flux_freqs(trans))
    flux_h = np.array(mp.get_fluxes(trans))
    sim.reset_meep()

    # ========== 参考仿真 (N=0, 使用相同 cell 尺寸!) ==========
    geom_r, _, _, _ = build_geom(a, rx, ry, 0)
    sim_r = mp.Simulation(cell_size=cell, resolution=resolution,
                          geometry=geom_r, sources=sources,
                          boundary_layers=[mp.PML(dpml)])
    trans_r = sim_r.add_flux(fcen, df, nfreq, freg)
    sim_r.run(until_after_sources=mp.stop_when_fields_decayed(
        50, mp.Ey, mp.Vector3(mon_x,0,h_total/2), 1e-4))
    flux_r = np.array(mp.get_fluxes(trans_r))

    T = np.divide(flux_h, flux_r, out=np.zeros_like(flux_h), where=flux_r > 1e-15)
    wl = 1000.0 / freqs_h
    idx = np.argsort(wl)
    wl = wl[idx]; T = T[idx]
    sim_r.reset_meep()
    return wl, T


# ============================================================================
# 禁带分析
# ============================================================================
def analyze_gap(wl, T, threshold=0.1, min_pts=5,
                target_center=1550.0, sigma=20.0):
    low = T < threshold
    runs = []; cur = 0
    for i, v in enumerate(low):
        if v:
            if cur == 0: si = i
            cur += 1
        else:
            if cur >= min_pts:
                runs.append((wl[si], wl[i-1], cur))
            cur = 0
    if cur >= min_pts:
        runs.append((wl[si], wl[-1], cur))

    result = {"has_gap": len(runs) > 0, "T_min": float(T.min()), "T_avg": float(T.mean())}

    if not runs:
        result.update({"gap_start": 0, "gap_end": 0, "gap_width": 0,
                       "center": 0, "score": -100.0})
        return result

    best = max(runs, key=lambda x: x[1] - x[0])
    gs, ge, npts = best
    gw = ge - gs
    center = (gs + ge) / 2.0
    penalty = np.exp(-((center - target_center) / sigma)**2)
    score = gw * penalty
    if T.min() > threshold * 0.5:
        score *= 0.5

    result.update({"gap_start": float(gs), "gap_end": float(ge),
                   "gap_width": float(gw), "center": float(center),
                   "npts": npts, "score": float(score), "penalty": float(penalty)})
    return result


# ============================================================================
# 保存谱图
# ============================================================================
def save_spectrum(wl, T, params, label, result):
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(wl, T, 'b-', lw=1.5, label='3D FDTD')
    ax.axhline(y=0.1, color='r', ls='--', lw=1, label='T=0.1')
    ax.axhline(y=1.0, color='gray', ls=':', lw=0.5)
    ax.set_xlim(1250, 1750); ax.set_ylim(-0.05, 1.15)
    if result.get("has_gap"):
        ax.axvspan(result["gap_start"], result["gap_end"],
                   alpha=0.15, color='green', label=f"gap={result['gap_width']:.0f}nm")
        ax.annotate(f"{result['gap_start']:.0f}-{result['gap_end']:.0f}nm",
                    xy=(result["center"], 0.5), ha='center', fontsize=11, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
    a, rx, ry, N = params
    ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Transmission")
    ax.set_title(f"Optimizer: a={a:.3f} rx={rx:.3f} ry={ry:.3f} N={N} | "
                 f"score={result['score']:.1f} width={result['gap_width']:.0f}nm")
    ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"results_3d/figures/{label}_spectrum.png", dpi=150)
    plt.close()


# ============================================================================
# 坐标下降优化器
# ============================================================================
def coordinate_descent(start_params, scan_ranges, fixed_N=20,
                        max_iterations=3, convergence_tol=0.5,
                        log_file="results_3d/optimization_log.csv"):
    best = dict(start_params)
    best["N"] = fixed_N
    best_result = None
    iteration = 0
    param_order = ["a", "ry", "rx"]

    ok, violations = check_periodic_geometry(
        best["a"], best["rx"], best["ry"], w_wg=w_wg
    )
    if not ok:
        print(f"⚠️ 起始参数不满足最小特征尺寸 {min_feature_nm}nm: {', '.join(violations)}")

    with open(log_file, "w") as f:
        f.write("iteration,param,a,rx,ry,N,T_min,T_avg,gap_start,gap_end,gap_width,center,score,penalty\n")

    while iteration < max_iterations:
        iteration += 1
        print(f"\n{'='*70}")
        print(f"迭代 #{iteration}")
        print(f"{'='*70}")
        any_improved = False

        for param in param_order:
            values = scan_ranges[param]
            print(f"\n--- 扫描 {param}: {values} ---")
            print(f"    当前最佳: {best}")

            candidates = []
            for val in values:
                if abs(val - best[param]) < 1e-6:
                    continue

                params = dict(best)
                params[param] = val
                a, rx, ry, N = params["a"], params["rx"], params["ry"], params["N"]
                ok, violations = check_periodic_geometry(a, rx, ry, w_wg=w_wg)
                if not ok:
                    print(f"\n  ⛔ SKIP a={a:.3f} rx={rx:.3f} ry={ry:.3f} | "
                          f"min feature {min_feature_nm}nm: {', '.join(violations)}")
                    continue

                label = f"opt_i{iteration}_{param}{val:.3f}"

                print(f"\n  ▶ a={a:.3f} rx={rx:.3f} ry={ry:.3f} N={N}", end="", flush=True)
                t0 = time.time()
                try:
                    wl, T = run_sim(a, rx, ry, N)
                    elapsed = time.time() - t0
                    result = analyze_gap(wl, T)
                    result["a"] = a; result["rx"] = rx; result["ry"] = ry; result["N"] = N
                    result["elapsed"] = elapsed
                    candidates.append(result)
                    save_spectrum(wl, T, (a, rx, ry, N), label, result)

                    status = "✅" if result["has_gap"] else "❌"
                    print(f" {elapsed:.0f}s | {status} score={result['score']:.1f} "
                          f"width={result['gap_width']:.0f}nm center={result['center']:.0f}nm "
                          f"T_min={result['T_min']:.4f}")

                    with open(log_file, "a") as f:
                        f.write(f"{iteration},{param},{a},{rx},{ry},{N},"
                                f"{result['T_min']:.6f},{result['T_avg']:.6f},"
                                f"{result['gap_start']:.1f},{result['gap_end']:.1f},"
                                f"{result['gap_width']:.1f},{result['center']:.1f},"
                                f"{result['score']:.2f},{result['penalty']:.4f}\n")
                except Exception as e:
                    print(f" ❌ ERROR: {e}")
                    continue

            if candidates:
                best_candidate = max(candidates, key=lambda x: x["score"])
                print(f"\n  → {param} 最佳值: {best_candidate[param]:.3f} "
                      f"(score={best_candidate['score']:.1f}, "
                      f"width={best_candidate['gap_width']:.0f}nm)")

                if best_result is None:
                    best_result = best_candidate
                    for k in ["a", "rx", "ry", "N"]:
                        best[k] = best_candidate[k]
                    any_improved = True
                    print(f"  → 设置基线")
                elif best_candidate["score"] > best_result["score"] + convergence_tol:
                    old_score = best_result["score"]
                    best_result = best_candidate
                    for k in ["a", "rx", "ry", "N"]:
                        best[k] = best_candidate[k]
                    any_improved = True
                    print(f"  → ✓ 改进: {best_result['score'] - old_score:.1f}")
                else:
                    print(f"  → 无显著改进 (≤{convergence_tol})")

        if not any_improved:
            print(f"\n✅ 收敛! 所有参数无进一步改进")
            break

    # 最终结果
    if best_result is None:
        print("\n❌ 未找到有效结果")
        return best

    print(f"\n{'='*70}")
    print("🏆 优化完成!")
    print(f"{'='*70}")
    print(f"  最佳参数: a={best_result['a']:.3f} rx={best_result['rx']:.3f} "
          f"ry={best_result['ry']:.3f} N={best_result['N']}")
    print(f"  禁带: {best_result['gap_start']:.0f}-{best_result['gap_end']:.0f}nm "
          f"(宽度 {best_result['gap_width']:.0f}nm)")
    print(f"  中心波长: {best_result['center']:.0f}nm")
    print(f"  T_min: {best_result['T_min']:.4f}")
    print(f"  综合得分: {best_result['score']:.1f}")

    best_params = {
        "a": best_result["a"], "rx": best_result["rx"], "ry": best_result["ry"],
        "N": best_result["N"],
        "gap_start_nm": best_result["gap_start"],
        "gap_end_nm": best_result["gap_end"],
        "gap_width_nm": best_result["gap_width"],
        "center_nm": best_result["center"],
        "T_min": best_result["T_min"],
        "score": best_result["score"],
    }
    with open("results_3d/optimized_params.json", "w") as f:
        json.dump(best_params, f, indent=2)
    print(f"  → 已保存 results_3d/optimized_params.json")
    print(f"  → 优化日志: {log_file}")
    return best


# ============================================================================
# 主程序
# ============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("3D 正确脊型波导自动优化器")
    print("=" * 70)

    start = {"a": 0.640, "rx": 0.200, "ry": 0.220}

    scan_ranges = {
        "a":  [0.600, 0.620, 0.640, 0.660, 0.680, 0.700],
        "rx": [0.200, 0.210, 0.220],
        "ry": [0.200, 0.210, 0.220, 0.230, 0.240],
    }

    print(f"\n📌 起始点: a={start['a']} rx={start['rx']} ry={start['ry']} N=20")
    print(f"   最小特征尺寸限制: {min_feature_nm}nm")
    print("   (c1_ridge, 已知 3D 禁带 31nm @1549nm)")
    print()

    coordinate_descent(
        start_params=start,
        scan_ranges=scan_ranges,
        fixed_N=20,
        max_iterations=3,
        convergence_tol=0.5,
        log_file="results_3d/optimization_log.csv"
    )

    print(f"\n{'='*70}")
    print("完成! 可查看:")
    print("  results_3d/optimized_params.json — 最佳参数")
    print("  results_3d/optimization_log.csv — 全部扫描记录")
    print("  results_3d/figures/opt_*.png — 透射谱图")
    print(f"{'='*70}")
