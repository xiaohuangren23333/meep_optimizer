#!/usr/bin/env python3
"""
2D FDTD 周期孔禁带大规模扫描
=============================
3D 结构的 2D 等效:
  - 波导: n=2.18 (LT), 宽度 w_wg=1.5um
  - 椭圆柱孔: n=1.0, 全刻蚀 (贯穿)
  - 2D 模拟: 在 x-y 平面, 忽略 z 方向

扫描策略: 四阶段
  Phase 0: 快速粗扫 (验证脚本正确性, 10组)
  Phase 1: rx × ry × a × N 大规模扫描
  Phase 2: 聚焦最佳区域精细扫描
  Phase 3: 输出候选参数 (3D验证用)

2D 每组 ~3-5 秒, 1000 组约 1-2 小时
"""
import meep as mp
import numpy as np
import os, time, json, csv, argparse, itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import traceback

from config import min_feature_nm, min_feature_um, w_wg, check_periodic_geometry

# ============================================================================
# 物理参数 (同3D)
# ============================================================================
n_wg = 2.18; n_clad = 1.0

lambda_min = 1.25; lambda_max = 1.75
fcen = 1.0 / 1.50
df = (1.0/lambda_min - 1.0/lambda_max) / 2.0
nfreq = 1000                   # 2D 可用较少频率点
resolution = 20
dpml = 1.0
pad = 2.0

# 存储目录
RESULTS_DIR = "results"
SPECTRA_DIR = f"{RESULTS_DIR}/spectra"
FIGURES_DIR = f"{RESULTS_DIR}/figures"
os.makedirs(SPECTRA_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

RESULTS_CSV = f"{RESULTS_DIR}/bandgap_2d_scan.csv"
CANDIDATES_JSON = f"{RESULTS_DIR}/bandgap_candidates_2d.json"
BEST_2D_JSON = f"{RESULTS_DIR}/best_bandgap_2d.json"

# ============================================================================
# 几何构建 (2D)
# ============================================================================
def build_geom_and_cell_2d(a, rx, ry, N):
    """2D 波导 + 椭圆柱孔"""
    sx = 2*dpml + N*a + 2*pad
    sy = 2*dpml + w_wg + 2*pad
    cell = mp.Vector3(sx, sy)

    # 波导 (2D slab)
    wg = mp.Block(material=mp.Medium(index=n_wg),
                  center=mp.Vector3(0, 0),
                  size=mp.Vector3(mp.inf, w_wg))

    geom = [wg]

    # 椭圆柱孔
    for i in range(N):
        x = (i - (N-1)/2) * a
        geom.append(mp.Ellipsoid(material=mp.Medium(index=n_clad),
                                center=mp.Vector3(x, 0),
                                size=mp.Vector3(2*rx, 2*ry, mp.inf)))
    return geom, cell


def build_ref_geom_2d(cell):
    """无孔参考波导"""
    wg = mp.Block(material=mp.Medium(index=n_wg),
                  center=mp.Vector3(0, 0),
                  size=mp.Vector3(mp.inf, w_wg))
    return [wg]


# ============================================================================
# 仿真 (2D EigenModeSource)
# ============================================================================
def run_sim_2d(geom, cell):
    """2D 仿真 (eig_parity=ODD_Z, 同 Meep 官方教程)"""
    sx = cell.x
    src_x = -sx/2 + dpml + 0.5
    mon_x = sx/2 - dpml - 0.5

    src_size = mp.Vector3(0, 3*w_wg)
    sources = [mp.EigenModeSource(
        mp.GaussianSource(fcen, fwidth=df),
        eig_band=1,
        eig_parity=mp.ODD_Z,  # 2D 标准 TE 模式
        center=mp.Vector3(src_x, 0),
        size=src_size,
    )]

    sim = mp.Simulation(cell_size=cell, resolution=resolution,
                        geometry=geom, sources=sources,
                        boundary_layers=[mp.PML(dpml)])

    mon_size = mp.Vector3(0, 3*w_wg)
    freg = mp.FluxRegion(center=mp.Vector3(mon_x, 0), size=mon_size)
    trans = sim.add_flux(fcen, df, nfreq, freg)

    t0 = time.time()
    sim.run(until_after_sources=mp.stop_when_fields_decayed(
        50, mp.Ey, mp.Vector3(mon_x, 0), 1e-5))
    elapsed = time.time() - t0

    freqs = np.array(mp.get_flux_freqs(trans))
    flux = np.array(mp.get_fluxes(trans))
    sim.reset_meep()
    return freqs, flux, elapsed


# ============================================================================
# 禁带分析
# ============================================================================
def analyze_bandgap(wl, T):
    """分析禁带: 在 1480-1650nm 范围内找最深禁带"""
    mask = (wl >= 1480) & (wl <= 1650)
    wl_f, T_f = wl[mask], T[mask]
    if len(wl_f) < 20:
        return {"T_min": float(T.min()), "gap_nm": 0, "has_gap": False, "T_min_gap": float(T.min())}

    Tmin = float(T_f.min())
    best = {"T_min": Tmin, "gap_nm": 0, "has_gap": False, "T_min_gap": Tmin}

    for th in np.arange(0.01, 0.90, 0.02):
        below = T_f < th
        if below.sum() < 3: continue
        d = np.diff(np.concatenate([[False], below, [False]]).astype(int))
        st = np.where(d == 1)[0]; en = np.where(d == -1)[0]
        lens = en - st
        if len(lens) == 0: continue
        bi = np.argmax(lens)
        gw = abs(wl_f[st[bi]] - wl_f[en[bi]-1])
        if gw > best["gap_nm"] and gw > 3:
            gs, ge = float(wl_f[st[bi]]), float(wl_f[en[bi]-1])
            gm = (wl >= gs) & (wl <= ge)
            gt = float(T[gm].min()) if gm.sum() > 0 else Tmin
            best = {"T_min": Tmin, "gap_nm": gw,
                    "T_min_gap": gt, "gap_start": gs, "gap_end": ge,
                    "threshold": th, "has_gap": True}
    return best


def candidate_score(row):
    """Score favors wide gap and low in-gap transmission."""
    if not row.get("has_gap"):
        return -1e9
    gap_nm = float(row.get("gap_nm", 0.0))
    t_gap = float(row.get("T_min_gap", 1.0))
    # Wider stopband is better; lower transmission is better.
    # The transmission term mirrors Meep tutorial logic: prioritize deep bandgaps.
    return gap_nm * 2.0 - 180.0 * t_gap


# ============================================================================
# 绘图 (仅对有禁带的)
# ============================================================================
def plot_spectrum(wl, T, params, result, fpath):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(wl, T, 'b-', lw=1.5, alpha=0.8)
    ax.axhline(1.0, color='gray', ls=':', lw=0.5, alpha=0.5)
    ax.axhline(0.1, color='red', ls='--', lw=1, alpha=0.5)

    if result.get("has_gap"):
        ax.axvspan(result["gap_start"], result["gap_end"],
                   alpha=0.12, color='red',
                   label=f"GAP {result['gap_nm']:.0f}nm Tmin={result['T_min_gap']:.4f}")

    title = (f"2D: a={params['a']:.3f} rx={params['rx']:.3f} ry={params['ry']:.3f} N={params['N']}"
             f" | Tmin={result['T_min']:.4f}")
    if result.get("has_gap"):
        title += f" GAP {result['gap_nm']:.0f}nm"

    ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Transmission")
    ax.set_title(title)
    ax.set_xlim(1400, 1700); ax.set_ylim(-0.02, 1.25)
    if result.get("has_gap"):
        ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout(); plt.savefig(fpath, dpi=150); plt.close()


# ============================================================================
# 单次运行
# ============================================================================
def run_parameter(a, rx, ry, N, run_id, total, plot_gap_only=True):
    label = f"a{a:.3f}_rx{rx:.3f}_ry{ry:.3f}_N{N}"
    csv_path = f"{SPECTRA_DIR}/2d_{label}.csv"
    png_path = f"{FIGURES_DIR}/2d_{label}.png"

    ok, violations = check_periodic_geometry(a, rx, ry, w_wg=w_wg)
    if not ok:
        print(f"[{run_id:02d}/{total}] ⛔ SKIP {label} | min feature {min_feature_nm}nm: {', '.join(violations)}")
        return None

    # 断点续传
    if os.path.isfile(csv_path):
        # 读已有结果
        d = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        if len(d) > 0:
            wl, T = d[:, 0], d[:, 1]
            result = analyze_bandgap(wl, T)
            T_max = float(T.max())
            print(f"[{run_id:02d}/{total}] ⏭ SKIP {label} | Tmin={result['T_min']:.4f}", end="")
            if result.get("has_gap"):
                print(f" GAP={result['gap_nm']:.0f}nm", end="")
            print()
            return {"label": label, "a": a, "rx": rx, "ry": ry, "N": N,
                    "T_min": result["T_min"], "T_min_gap": result.get("T_min_gap", result["T_min"]),
                    "gap_nm": result.get("gap_nm", 0), "has_gap": result.get("has_gap", False),
                    "T_max": T_max, "time_s": 0, "filename": f"2d_{label}.csv"}

    print(f"[{run_id:02d}/{total}] ▶ {label}", end="", flush=True)
    try:
        geom, cell = build_geom_and_cell_2d(a, rx, ry, N)
        freqs_h, flux_h, t_h = run_sim_2d(geom, cell)
        ref_geom = build_ref_geom_2d(cell)
        freqs_r, flux_r, t_r = run_sim_2d(ref_geom, cell)

        T = np.divide(flux_h, flux_r, out=np.zeros_like(flux_h), where=flux_r > 1e-15)
        wl = 1000.0 / freqs_r
        idx = np.argsort(wl); wl, T = wl[idx], T[idx]

        result = analyze_bandgap(wl, T)
        T_max = float(T.max())

        # 保存频谱
        np.savetxt(csv_path, np.column_stack([wl, T]),
                   delimiter=",", header="wavelength_nm,transmission", comments="")

        # 绘图 (只画有禁带的，除非 plot_gap_only=False)
        has_gap = result.get("has_gap", False) and result.get("gap_nm", 0) > 10
        if has_gap or not plot_gap_only:
            plot_spectrum(wl, T, {"a": a, "rx": rx, "ry": ry, "N": N}, result, png_path)

        status = "✅" if result.get("has_gap") else "  "
        score = ""
        if result.get("has_gap"):
            gt = result.get("T_min_gap", 1)
            if gt < 0.05: score = "🏆"
            elif gt < 0.1: score = "⭐"
            elif gt < 0.3: score = "👍"

        gap_info = f" | GAP {result['gap_nm']:.0f}nm Tgap={result.get('T_min_gap',0):.4f}" if result.get("has_gap") else f" | Tmin={result['T_min']:.4f}"
        print(f" | {status}{gap_info} Tmax={T_max:.3f} {score} [{t_h+t_r:.1f}s]")

        return {"label": label, "a": a, "rx": rx, "ry": ry, "N": N,
                "T_min": result["T_min"], "T_min_gap": result.get("T_min_gap", result["T_min"]),
                "gap_nm": result.get("gap_nm", 0), "has_gap": result.get("has_gap", False),
                "T_max": T_max, "time_s": round(t_h + t_r, 1),
                "filename": f"2d_{label}.csv"}
    except Exception as e:
        print(f" ❌ ERROR: {e}")
        traceback.print_exc()
        return None


# ============================================================================
# 主程序
# ============================================================================
def main():
    parser = argparse.ArgumentParser(description="2D 禁带大规模扫描")
    parser.add_argument("--phase", type=int, default=1, help="扫描阶段 (0=验证, 1=粗扫, 2=精细)")
    parser.add_argument("--dry-run", action="store_true", help="只打印参数列表，不运行")
    args = parser.parse_args()

    if args.phase == 0:
        # Phase 0: 快速验证 (10组)
        params = []
        for rx in [0.20, 0.21, 0.22]:
            for ry in [0.20, 0.30, 0.40]:
                params.append((0.64, rx, ry, 16))
        params = params[:10]  # 最多10组

    elif args.phase == 1:
        # Phase 1: 大规模粗扫 (最小特征尺寸 >= 200nm)
        # a 需满足 a >= 2*rx + 200nm，故 a 从 0.60 μm 起扫
        a_vals = np.round(np.arange(0.60, 0.78, 0.02), 3)
        rx_vals = np.round(np.arange(0.20, 0.26, 0.02), 3)
        ry_vals = np.round(np.arange(0.20, 0.47, 0.04), 3)
        N_vals = [12, 16, 20]

        params = [
            (a, rx, ry, N)
            for a, rx, ry, N in itertools.product(a_vals, rx_vals, ry_vals, N_vals)
            if check_periodic_geometry(a, rx, ry, w_wg=w_wg)[0]
        ]
        print(f"Phase 1: 有效组合 {len(params)} 组 (已过滤 <{min_feature_nm}nm 特征)")
        print(f"最小特征尺寸限制: {min_feature_nm}nm (rx, ry, 孔间距, 波导剩余宽度)")

    elif args.phase == 2:
        # Phase 2: 聚焦最佳区域精细扫描
        a_vals = np.round(np.arange(0.60, 0.72, 0.01), 3)
        rx_vals = np.round(np.arange(0.20, 0.25, 0.01), 3)
        ry_vals = np.round(np.arange(0.20, 0.45, 0.02), 3)
        N_vals = [16, 20, 24]

        params = [
            (a, rx, ry, N)
            for a, rx, ry, N in itertools.product(a_vals, rx_vals, ry_vals, N_vals)
            if check_periodic_geometry(a, rx, ry, w_wg=w_wg)[0]
        ]
        print(f"Phase 2: 有效组合 {len(params)} 组 (已过滤 <{min_feature_nm}nm 特征)")

    else:
        print(f"Unknown phase: {args.phase}")
        return

    total = len(params)
    print(f"预计时间: ~{total * 5 / 60:.1f} 分钟 ({total * 5 / 3600:.1f} 小时)")
    print(f"结果: {RESULTS_CSV}")
    print(f"频谱: {SPECTRA_DIR}/")
    print(f"图片: {FIGURES_DIR}/")
    print()

    if args.dry_run:
        for a, rx, ry, N in params[:20]:
            print(f"  a={a:.3f} rx={rx:.3f} ry={ry:.3f} N={N}")
        if len(params) > 20:
            print(f"  ... ({len(params)} total)")
        return

    results = []
    t_start = time.time()
    completed = 0

    for run_id, (a, rx, ry, N) in enumerate(params, 1):
        res = run_parameter(a, rx, ry, N, run_id, total, plot_gap_only=True)
        if res is not None:
            results.append(res)
            completed += 1
            # 实时写入 CSV
            fn = ["label","a","rx","ry","N","T_min","T_min_gap",
                  "gap_nm","has_gap","T_max","time_s","filename"]
            fe = os.path.isfile(RESULTS_CSV)
            with open(RESULTS_CSV, "a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fn)
                if not fe: w.writeheader()
                w.writerow({k: res.get(k, "") for k in fn})

        # 进度和剩余时间
        if run_id % 20 == 0 or run_id == total:
            elapsed = time.time() - t_start
            if run_id > 0:
                per_run = elapsed / run_id
                remaining = (total - run_id) * per_run
                print(f"      ⏱ 已过{elapsed/60:.1f}分, 预计剩余{remaining/60:.1f}分")
        print()

    elapsed = time.time() - t_start
    print("=" * 70)
    print(f"🎉 扫描完成! {completed}/{total} 组, {elapsed/60:.1f}分 ({elapsed/3600:.2f}小时)")

    # 分析结果
    if results:
        # 有禁带的 (gap > 10nm), 按“宽禁带 + 低透射”综合评分排序
        valid = [r for r in results if r.get("has_gap") and r.get("gap_nm", 0) > 10]
        valid.sort(key=candidate_score, reverse=True)

        print(f"\n有禁带 (>10nm): {len(valid)}/{len(results)}")
        if valid:
            print(f"\n🏆 Top 20 禁带候选 (宽禁带 + 低透射):")
            print(f"{'#':>3} {'a':>5} {'rx':>5} {'ry':>5} {'N':>3} {'gap_nm':>7} {'T_min':>7} {'score':>8}")
            print("-" * 58)
            for i, r in enumerate(valid[:20]):
                badge = "🏆" if r["T_min_gap"] < 0.05 else ("⭐" if r["T_min_gap"] < 0.1 else "👍" if r["T_min_gap"] < 0.3 else "")
                print(f"{i+1:>3} {r['a']:.3f} {r['rx']:.3f} {r['ry']:.3f} {r['N']:>3} "
                      f"{r['gap_nm']:>5.0f}nm {r['T_min_gap']:>7.4f} {candidate_score(r):>8.2f} {badge}")

        # 保存候选
        candidates = []
        for i, r in enumerate(valid[:20]):
            candidates.append({
                "rank": i+1, "a": r["a"], "rx": r["rx"], "ry": r["ry"],
                "N": r["N"], "gap_nm": r["gap_nm"], "T_min_gap": r["T_min_gap"],
                "T_max": r["T_max"]
            })

        # 也保存 all T_min < 0.3
        deep = [r for r in results if r.get("T_min", 1) < 0.3]
        deep_candidates = []
        for i, r in enumerate(sorted(deep, key=lambda x: x.get("T_min", 999))[:20]):
            deep_candidates.append({
                "rank": i+1, "a": r["a"], "rx": r["rx"], "ry": r["ry"],
                "N": r["N"], "gap_nm": r.get("gap_nm", 0), "T_min": r["T_min"],
                "T_max": r.get("T_max", 1)
            })

        with open(CANDIDATES_JSON, "w") as f:
            json.dump({"top_bandgap": candidates, "deep_Tmin": deep_candidates}, f, indent=2)
        print(f"\n📁 候选: {CANDIDATES_JSON}")

        if valid:
            best = valid[0]
            best_payload = {
                "a": best["a"],
                "rx": best["rx"],
                "ry": best["ry"],
                "N": best["N"],
                "gap_width_nm": best["gap_nm"],
                "T_min_gap": best["T_min_gap"],
                "T_max": best["T_max"],
                "selection_score": candidate_score(best),
                "selection_rule": "maximize (2*gap_nm - 180*T_min_gap) with gap_nm>10 and has_gap=True",
                "source_csv": RESULTS_CSV,
            }
            with open(BEST_2D_JSON, "w") as f:
                json.dump(best_payload, f, indent=2)
            print(f"📌 2D 最优参数: {BEST_2D_JSON}")

        # 检查 T_max
        maxT = max(r.get("T_max", 1) for r in results)
        print(f"全局 T_max = {maxT:.4f}" + (" ✅ 全部 T≤1" if maxT <= 1.05 else f" ⚠️ 有 T>{1}"))
    else:
        print("⚠️ 无结果数据")

    print("=" * 70)


if __name__ == "__main__":
    main()