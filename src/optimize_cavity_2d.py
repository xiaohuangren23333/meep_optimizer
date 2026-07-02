#!/usr/bin/env python3
"""
一维光子晶体二次渐变缺陷腔自动优化器 (2D)
=========================================
根据规范分 6 步优化，而非全维度扫描。

优化顺序:
  1. 扫描 a_c (固定 Nt=3, Nm=5, rx_c=rx_m, ry_c=ry_m)
  2. 扫描 rx_c (固定最佳 a_c)
  3. 扫描 ry_c (固定最佳 a_c, rx_c)
  4. 局部联合微调 a_c, rx_c, ry_c
  5. 扫描 N_taper
  6. 扫描 N_mirror (减少孔数)

数据格式严格遵循规范要求。
"""
import sys, os, time, json, warnings
import numpy as np
import meep as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from config import (
    n_wg, w_wg, min_feature_nm,
    check_cavity_geometry, check_periodic_geometry,
)

warnings.filterwarnings('ignore')

# ============================================================================
# 常量
# ============================================================================
n_air = 1.0
resolution = 30
dpml = 1.0
pad = 2.0

lambda_min = 1.25
lambda_max = 1.75
fcen = 1.0 / 1.50
df = (1.0/lambda_min - 1.0/lambda_max) / 2.0
nfreq = 2000  # 高分辨率

os.makedirs("results/spectra", exist_ok=True)
os.makedirs("results/figures", exist_ok=True)

# ============================================================================
# Mirror 候选读取
# ============================================================================
def load_mirror_candidate(path=None):
    """加载 mirror candidate JSON，自动映射键名"""
    raw = None
    if path and os.path.exists(path):
        with open(path) as f:
            raw = json.load(f)
    else:
        raw = {
            "a": 0.640, "rx": 0.200, "ry": 0.220, "N_period": 20,
            "gap_start_nm": 1529.0, "gap_end_nm": 1589.0,
            "gap_width_nm": 60.0, "T_min": 0.017, "T_avg": 0.5,
        }

    # 映射各种可能的键名 → 统一字段
    mirror = {
        "a_m": raw.get("a_m") or raw.get("a") or raw.get("a_m"),
        "rx_m": raw.get("rx_m") or raw.get("rx") or raw.get("rx_m"),
        "ry_m": raw.get("ry_m") or raw.get("ry") or raw.get("ry_m"),
        "N_period": raw.get("N_period") or raw.get("N") or 20,
        "bandgap_start_nm": raw.get("bandgap_start_nm") or raw.get("gap_start_nm", 1529),
        "bandgap_end_nm": raw.get("bandgap_end_nm") or raw.get("gap_end_nm", 1589),
        "gap_width_nm": raw.get("gap_width_nm") or raw.get("gap_width_nm", 60),
        "T_min_1500_1600": raw.get("T_min_1500_1600") or raw.get("T_min", 0.017),
        "T_avg_1500_1600": raw.get("T_avg_1500_1600") or raw.get("T_avg", 0.5),
    }
    ok, violations = check_periodic_geometry(
        mirror["a_m"], mirror["rx_m"], mirror["ry_m"], w_wg=w_wg
    )
    mirror["geometry_valid"] = ok
    mirror["geometry_violations"] = violations
    return mirror


# ============================================================================
# 几何构建
# ============================================================================
def build_cavity_geom(a_m, rx_m, ry_m, a_c, rx_c, ry_c, N_taper, N_mirror):
    """
    二次渐变缺陷腔:
      镜区: N_mirror 个周期孔 (a_m, rx_m, ry_m)
      渐变: N_taper 个孔按二次渐变过渡到缺陷
      中心: 空缺陷
      对称结构
    """
    total_holes = 2 * (N_mirror + N_taper)
    sx = 2*dpml + max(total_holes * a_m, 4.0) + 2*pad

    geom = [
        mp.Block(material=mp.Medium(index=1.0),
                 center=mp.Vector3(0, 0), size=mp.Vector3(mp.inf, mp.inf, mp.inf)),
        mp.Block(material=mp.Medium(index=n_wg),
                 center=mp.Vector3(0, 0), size=mp.Vector3(mp.inf, w_wg, 0)),
    ]

    # 左侧: 镜区(远) → 渐变区(近)
    # 注意: 2D 模拟中 z 方向必须用 mp.inf 而非 0
    for i in range(N_mirror):
        x = -(N_mirror - i + N_taper) * a_m
        geom.append(mp.Ellipsoid(material=mp.Medium(index=n_air),
                                  center=mp.Vector3(x, 0, 0),
                                  size=mp.Vector3(2*rx_m, 2*ry_m, mp.inf)))
    for i in range(1, N_taper + 1):
        t = i / N_taper
        ai = a_c + (a_m - a_c) * t**2
        rxi = rx_c + (rx_m - rx_c) * t**2
        ryi = ry_c + (ry_m - ry_c) * t**2
        x = -(N_taper - i + 0.5) * ai
        geom.append(mp.Ellipsoid(material=mp.Medium(index=n_air),
                                  center=mp.Vector3(x, 0, 0),
                                  size=mp.Vector3(2*rxi, 2*ryi, mp.inf)))

    # 右侧: 渐变区(近) → 镜区(远)
    for i in range(1, N_taper + 1):
        t = i / N_taper
        ai = a_c + (a_m - a_c) * t**2
        rxi = rx_c + (rx_m - rx_c) * t**2
        ryi = ry_c + (ry_m - ry_c) * t**2
        x = (N_taper - i + 0.5) * ai
        geom.append(mp.Ellipsoid(material=mp.Medium(index=n_air),
                                  center=mp.Vector3(x, 0, 0),
                                  size=mp.Vector3(2*rxi, 2*ryi, mp.inf)))
    for i in range(N_mirror):
        x = (N_mirror - i + N_taper) * a_m
        geom.append(mp.Ellipsoid(material=mp.Medium(index=n_air),
                                  center=mp.Vector3(x, 0, 0),
                                  size=mp.Vector3(2*rx_m, 2*ry_m, mp.inf)))
    return geom, sx


# ============================================================================
# 仿真
# ============================================================================
def run_cavity_sim(geom, cell_sx):
    """运行 2D FDTD 腔仿真 + 参考波导归一化"""
    cell = mp.Vector3(cell_sx, 6.0)
    # 源在左侧 PML 内侧 + 0.5um 缓冲
    src_x = -cell_sx/2 + dpml + 0.5
    # 通量监视器在右侧 PML 内侧 - 0.5um 缓冲
    mon_x = cell_sx/2 - dpml - 0.5

    # === 有孔结构 ===
    src1 = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ey,
                      center=mp.Vector3(src_x, 0),
                      size=mp.Vector3(0, w_wg))]
    sim = mp.Simulation(cell_size=cell, resolution=resolution,
                        geometry=geom, sources=src1,
                        boundary_layers=[mp.PML(dpml)])
    freg1 = mp.FluxRegion(center=mp.Vector3(mon_x, 0), size=mp.Vector3(0, 2*w_wg))
    trans1 = sim.add_flux(fcen, df, nfreq, freg1)
    sim.run(until_after_sources=mp.stop_when_fields_decayed(
        50, mp.Ey, mp.Vector3(mon_x, 0), 1e-4))
    freqs = np.array(mp.get_flux_freqs(trans1))
    flux_h = np.array(mp.get_fluxes(trans1))
    sim.reset_meep()

    # === 参考波导（无孔，其余完全相同）===
    ref_geom = [
        mp.Block(material=mp.Medium(index=1.0),
                 center=mp.Vector3(0,0), size=mp.Vector3(mp.inf, mp.inf, mp.inf)),
        mp.Block(material=mp.Medium(index=n_wg),
                 center=mp.Vector3(0,0), size=mp.Vector3(mp.inf, w_wg, 0)),
    ]
    src2 = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ey,
                      center=mp.Vector3(src_x, 0),
                      size=mp.Vector3(0, w_wg))]
    sim_r = mp.Simulation(cell_size=cell, resolution=resolution,
                          geometry=ref_geom, sources=src2,
                          boundary_layers=[mp.PML(dpml)])
    freg2 = mp.FluxRegion(center=mp.Vector3(mon_x, 0), size=mp.Vector3(0, 2*w_wg))
    trans2 = sim_r.add_flux(fcen, df, nfreq, freg2)
    sim_r.run(until_after_sources=mp.stop_when_fields_decayed(
        50, mp.Ey, mp.Vector3(mon_x, 0), 1e-4))
    flux_r = np.array(mp.get_fluxes(trans2))
    sim_r.reset_meep()

    # === 归一化 ===
    T = np.divide(flux_h, flux_r, out=np.zeros_like(flux_h), where=flux_r > 1e-15)
    wl = 1000.0 / freqs
    idx = np.argsort(wl)
    wl = wl[idx]; T = T[idx]

    # 快速诊断：检查是否有异常
    sum_h = float(np.sum(flux_h))
    sum_r = float(np.sum(flux_r))
    if sum_h >= sum_r * 0.99:  # 有孔和参考几乎一样，可能腔结构有问题
        print(f" ⚠️ flux_h={sum_h:.1f} ≈ flux_r={sum_r:.1f}", end="")
    return wl, T


# ============================================================================
# Lorentzian 拟合
# ============================================================================
def lorentzian(x, x0, gamma, A, offset):
    return offset + A * gamma**2 / ((x - x0)**2 + gamma**2)


def fit_peak(wl_seg, T_seg):
    if len(wl_seg) < 7:
        return None
    try:
        x0_guess = wl_seg[np.argmax(T_seg)]
        A_guess = max(T_seg) - min(T_seg)
        gamma_guess = (wl_seg[-1] - wl_seg[0]) / 4
        offset_guess = min(T_seg)
        popt, _ = curve_fit(lorentzian, wl_seg, T_seg,
                            p0=[x0_guess, gamma_guess, A_guess, offset_guess],
                            maxfev=5000)
        x0, gamma, A, offset = popt
        T_peak = lorentzian(x0, *popt)
        FWHM = 2 * gamma
        residuals = T_seg - lorentzian(wl_seg, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((T_seg - np.mean(T_seg))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        return {"lambda0_nm": float(x0), "FWHM_nm": float(FWHM),
                "Q": float(x0 / FWHM) if FWHM > 0 else 0,
                "T_peak": float(T_peak), "A": float(A), "gamma": float(gamma),
                "offset": float(offset), "fit_r2": float(r2)}
    except:
        return None


# ============================================================================
# 峰提取
# ============================================================================
def extract_peaks(wl, T, bg_start, bg_end, delta_wl):
    mask = (wl >= bg_start - 10) & (wl <= bg_end + 10)
    wl_bg, T_bg = wl[mask], T[mask]
    if len(wl_bg) < 20:
        return []

    peaks, props = find_peaks(T_bg, height=0.01, prominence=0.005,
                               width=3, distance=3, rel_height=0.5)
    if len(peaks) == 0:
        return []

    candidates = []
    for k, idx in enumerate(peaks):
        left = max(idx - 12, 0)
        right = min(idx + 12, len(wl_bg) - 1)

        fit = fit_peak(wl_bg[left:right+1], T_bg[left:right+1])
        if fit is None:
            continue

        x0, FWHM, Q, T_pk, r2 = fit["lambda0_nm"], fit["FWHM_nm"], fit["Q"], fit["T_peak"], fit["fit_r2"]

        prom = float(props["prominences"][k]) if k < len(props["prominences"]) else 0

        valid = True; reasons = []
        if x0 < 1500 or x0 > 1600: valid = False; reasons.append("lambda范围")
        if x0 < bg_start or x0 > bg_end: valid = False; reasons.append("不在禁带")
        if T_pk < 0.05: valid = False; reasons.append("T_peak<0.05")
        if FWHM < 3 * delta_wl: valid = False; reasons.append("FWHM过窄")
        if r2 < 0.95: valid = False; reasons.append(f"R²={r2:.3f}<0.95")

        candidates.append({
            "lambda0_nm": float(x0), "FWHM_nm": float(FWHM), "Q": float(Q),
            "T_peak": float(T_pk), "prominence": float(prom),
            "fit_r2": float(r2),
            "distance_to_gap_edge_nm": float(min(abs(x0-bg_start), abs(x0-bg_end))),
            "valid": valid, "invalid_reason": "; ".join(reasons) if reasons else "",
        })
    return candidates


def score_peak(peak):
    if not peak["valid"]:
        return -100
    s = 0
    dl = abs(peak["lambda0_nm"] - 1550)
    if dl > 40:
        s -= 120
    elif dl > 20:
        s -= 50
    elif dl > 10:
        s -= 15
    else:
        s += max(0, 60 - dl * 4)
    if peak["Q"] >= 1000:
        s += 200
    elif peak["Q"] >= 500:
        s += 80
    s += min(peak["Q"] / 50, 80)
    s += min(peak["T_peak"] * 40, 40)
    s += peak["fit_r2"] * 15
    return s


# ============================================================================
# 保存
# ============================================================================
CSV_HEADER = (
    "run_id,source_mirror_candidate_id,"
    "a_m,rx_m,ry_m,"
    "a_c,rx_c,ry_c,"
    "N_taper,N_mirror,defect_gap,N_total,"
    "valid_geometry,"
    "lambda0_nm,FWHM_nm,Q,T_peak,prominence,fit_r2,"
    "best_peak_score,highest_q_peak_Q,"
    "spectrum_csv,spectrum_png,fit_png"
    ",error_message\n"
)

def save_run(wl, T, run_id, params, peaks, mirror, elapsed=0):
    """
    保存所有结果: CSV, 谱图, 拟合图。
    返回 fill_data dict。
    """
    a_c, rx_c, ry_c, Nt, Nm = params
    a_m, rx_m, ry_m = mirror["a_m"], mirror["rx_m"], mirror["ry_m"]
    bg_s, bg_e = mirror["bandgap_start_nm"], mirror["bandgap_end_nm"]
    N_total = 2 * (Nm + Nt)
    delta_wl = (lambda_max - lambda_min) * 1000 / nfreq

    valid_peaks = [p for p in peaks if p["valid"]]
    best_peak = max(valid_peaks, key=score_peak) if valid_peaks else None
    high_q_peak = max(valid_peaks, key=lambda x: x["Q"]) if valid_peaks else None

    # 保存 CSV
    csv_path = f"results/spectra/cavity_run_{run_id:04d}.csv"
    np.savetxt(csv_path, np.column_stack([wl, T]),
               delimiter=",", header="wavelength_nm,transmission", comments="")

    # ===== 谱图 =====
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(wl, T, 'b-', lw=1.5, alpha=0.8)
    ax.axvspan(bg_s, bg_e, alpha=0.08, color='red', label='bandgap')
    for p in peaks:
        c = 'green' if p["valid"] else 'gray'
        ax.axvline(p["lambda0_nm"], color=c, ls='--', lw=1, alpha=0.5)
        ax.annotate(f"Q={p['Q']:.0f}" if p["valid"] else "inv",
                    xy=(p["lambda0_nm"], p["T_peak"]), fontsize=7, color=c)
    if best_peak:
        ax.axvline(best_peak["lambda0_nm"], color='red', lw=2,
                   label=f"best: λ={best_peak['lambda0_nm']:.1f} Q={best_peak['Q']:.0f}")
    ax.set_xlim(1400, 1700); ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Transmission")
    ax.set_title(f"Cavity: a_c={a_c:.3f} rx_c={rx_c:.3f} ry_c={ry_c:.3f} Nt={Nt} Nm={Nm}")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    png_path = f"results/figures/cavity_run_{run_id:04d}_spectrum.png"
    plt.savefig(png_path, dpi=150)
    plt.close()

    # ===== 拟合图 =====
    fit_png = ""
    if best_peak and best_peak["valid"]:
        fit_png = f"results/figures/cavity_run_{run_id:04d}_fit.png"
        x0, gamma = best_peak["lambda0_nm"], best_peak["FWHM_nm"] / 2
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(wl, T, 'b-', lw=1, alpha=0.3, label='full')
        mask = (wl >= bg_s - 10) & (wl <= bg_e + 10)
        ax.plot(wl[mask], T[mask], 'b-', lw=1.5, label='in-bandgap')
        wl_fit = np.linspace(x0 - 3*gamma, x0 + 3*gamma, 200)
        T_fit = lorentzian(wl_fit, x0, gamma, best_peak["T_peak"] - best_peak.get("offset", 0),
                           best_peak.get("offset", 0))
        ax.plot(wl_fit, T_fit, 'r-', lw=2,
                label=f"Lorentz: Q={best_peak['Q']:.0f} R²={best_peak['fit_r2']:.3f}")
        ax.set_xlabel("nm"); ax.set_ylabel("T")
        ax.set_title(f"Fit: λ₀={x0:.2f}nm FWHM={best_peak['FWHM_nm']:.3f}nm")
        ax.legend(); ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(fit_png, dpi=150)
        plt.close()

    # 构建 data row
    row = {
        "run_id": run_id, "source_mirror_candidate_id": "default",
        "a_m": a_m, "rx_m": rx_m, "ry_m": ry_m,
        "a_c": a_c, "rx_c": rx_c, "ry_c": ry_c,
        "N_taper": Nt, "N_mirror": Nm, "defect_gap": 0, "N_total": N_total,
        "valid_geometry": 1,
        "lambda0_nm": best_peak["lambda0_nm"] if best_peak else 0,
        "FWHM_nm": best_peak["FWHM_nm"] if best_peak else 0,
        "Q": best_peak["Q"] if best_peak else 0,
        "T_peak": best_peak["T_peak"] if best_peak else 0,
        "prominence": best_peak["prominence"] if best_peak else 0,
        "fit_r2": best_peak["fit_r2"] if best_peak else 0,
        "best_peak_score": score_peak(best_peak) if best_peak else -100,
        "highest_q_peak_Q": high_q_peak["Q"] if high_q_peak else 0,
        "spectrum_csv": csv_path, "spectrum_png": png_path, "fit_png": fit_png,
        "error_message": "",
    }
    return row


def append_csv_row(row, fpath="results/cavity_optimization_all.csv"):
    """追加一行到 CSV"""
    with open(fpath, "a") as f:
        f.write(
            f"{row['run_id']},{row['source_mirror_candidate_id']},"
            f"{row['a_m']},{row['rx_m']},{row['ry_m']},"
            f"{row['a_c']},{row['rx_c']},{row['ry_c']},"
            f"{row['N_taper']},{row['N_mirror']},{row['defect_gap']},{row['N_total']},"
            f"{row['valid_geometry']},"
            f"{row['lambda0_nm']:.2f},{row['FWHM_nm']:.4f},"
            f"{row['Q']:.0f},{row['T_peak']:.4f},"
            f"{row['prominence']:.4f},{row['fit_r2']:.4f},"
            f"{row['best_peak_score']:.1f},{row['highest_q_peak_Q']:.0f},"
            f"{row['spectrum_csv']},{row['spectrum_png']},{row['fit_png']},"
            f"{row['error_message']}\n"
        )


def invalid_geometry_row(run_id, mirror, a_c, rx_c, ry_c, Nt, Nm, reason):
    """Return a skipped row for designs that violate minimum feature size."""
    return {
        "run_id": run_id, "source_mirror_candidate_id": "default",
        "a_m": mirror["a_m"], "rx_m": mirror["rx_m"], "ry_m": mirror["ry_m"],
        "a_c": a_c, "rx_c": rx_c, "ry_c": ry_c,
        "N_taper": Nt, "N_mirror": Nm, "defect_gap": 0, "N_total": 2*(Nm+Nt),
        "valid_geometry": 0,
        "lambda0_nm": 0, "FWHM_nm": 0, "Q": 0, "T_peak": 0,
        "prominence": 0, "fit_r2": 0,
        "best_peak_score": -100, "highest_q_peak_Q": 0,
        "spectrum_csv": "", "spectrum_png": "", "fit_png": "",
        "error_message": reason,
    }


def filter_cavity_scan_values(values, mirror, param, best_ac, best_rx_c, best_ry_c, Nt):
    """Keep cavity scan values that satisfy the 200nm minimum feature constraint."""
    filtered = []
    for val in values:
        ac = val if param == "a_c" else best_ac
        rxc = val if param == "rx_c" else best_rx_c
        ryc = val if param == "ry_c" else best_ry_c
        ok, _ = check_cavity_geometry(
            mirror["a_m"], mirror["rx_m"], mirror["ry_m"],
            ac, rxc, ryc, Nt, w_wg=w_wg,
        )
        if ok:
            filtered.append(val)
    return filtered


# ============================================================================
# 单次评估
# ============================================================================
def evaluate(run_id, a_c, rx_c, ry_c, Nt, Nm, mirror):
    """运行一次仿真 + 分析 + 保存，返回 row"""
    print(f"[{run_id}] a_c={a_c:.3f} rx_c={rx_c:.3f} ry_c={ry_c:.3f} Nt={Nt} Nm={Nm}", end="", flush=True)
    ok, violations = check_cavity_geometry(
        mirror["a_m"], mirror["rx_m"], mirror["ry_m"],
        a_c, rx_c, ry_c, Nt, w_wg=w_wg,
    )
    if not ok:
        reason = f"min feature {min_feature_nm}nm: {', '.join(violations)}"
        print(f" ⛔ SKIP | {reason}")
        row = invalid_geometry_row(run_id, mirror, a_c, rx_c, ry_c, Nt, Nm, reason)
        append_csv_row(row)
        return row

    t0 = time.time()
    try:
        geom, sx = build_cavity_geom(mirror["a_m"], mirror["rx_m"], mirror["ry_m"],
                                      a_c, rx_c, ry_c, Nt, Nm)
        wl, T = run_cavity_sim(geom, sx)
        elapsed = time.time() - t0
        delta_wl = (lambda_max - lambda_min) * 1000 / nfreq
        peaks = extract_peaks(wl, T, mirror["bandgap_start_nm"],
                              mirror["bandgap_end_nm"], delta_wl)
        row = save_run(wl, T, run_id, (a_c, rx_c, ry_c, Nt, Nm), peaks, mirror, elapsed)
        row["valid_geometry"] = 1
        bp = max([p for p in peaks if p["valid"]], key=score_peak) if any(p["valid"] for p in peaks) else None
        if bp:
            print(f" | Q={bp['Q']:.0f} T={bp['T_peak']:.2f} λ={bp['lambda0_nm']:.1f} R²={bp['fit_r2']:.3f}")
        else:
            print(" | no valid peak")
        row["error_message"] = ""
        append_csv_row(row)
        return row
    except Exception as e:
        elapsed = time.time() - t0
        print(f" ❌ {e}")
        row = {
            "run_id": run_id, "source_mirror_candidate_id": "default",
            "a_m": mirror["a_m"], "rx_m": mirror["rx_m"], "ry_m": mirror["ry_m"],
            "a_c": a_c, "rx_c": rx_c, "ry_c": ry_c,
            "N_taper": Nt, "N_mirror": Nm, "defect_gap": 0, "N_total": 2*(Nm+Nt),
            "valid_geometry": 0,
            "lambda0_nm": 0, "FWHM_nm": 0, "Q": 0, "T_peak": 0,
            "prominence": 0, "fit_r2": 0,
            "best_peak_score": -100, "highest_q_peak_Q": 0,
            "spectrum_csv": "", "spectrum_png": "", "fit_png": "",
            "error_message": str(e),
        }
        append_csv_row(row)
        return row


# ============================================================================
# 主优化
# ============================================================================
def optimize(mirror):
    a_m, rx_m, ry_m = mirror["a_m"], mirror["rx_m"], mirror["ry_m"]

    # 初始化 CSV
    if not os.path.isfile("results/cavity_optimization_all.csv"):
        with open("results/cavity_optimization_all.csv", "w") as f:
            f.write(CSV_HEADER)

    run_id = sum(1 for _ in open("results/cavity_optimization_all.csv")) - 1
    if run_id < 0: run_id = 0

    best_overall = None  # (row, score)
    base_nm = 10  # 高 Q 需要更多镜区孔

    # ================================================================
    # Step 1: 扫描 a_c
    # ================================================================
    print(f"\n{'='*60}")
    print(f"Step 1: 扫描 a_c (固定 Nt=3 Nm={base_nm} rx_c=rx_m ry_c=ry_m)")
    print(f"{'='*60}")
    a_c_vals = filter_cavity_scan_values(
        np.round(np.linspace(0.80*a_m, 1.10*a_m, 13), 4),
        mirror, "a_c", a_m, rx_m, ry_m, 3,
    )
    step1_best = None
    for a_c in a_c_vals:
        run_id += 1
        row = evaluate(run_id, a_c, rx_m, ry_m, Nt=3, Nm=base_nm, mirror=mirror)
        if best_overall is None or row["best_peak_score"] > best_overall[1]:
            best_overall = (row, row["best_peak_score"])
        if step1_best is None or row["best_peak_score"] > step1_best[1]:
            step1_best = (row, row["best_peak_score"])
    best_ac = step1_best[0]["a_c"] if step1_best else a_m
    print(f"  → 最佳 a_c = {best_ac:.3f}")

    # ================================================================
    # Step 2: 扫描 rx_c
    # ================================================================
    print(f"\n{'='*60}")
    print(f"Step 2: 扫描 rx_c (固定 a_c={best_ac:.3f} Nt=3 Nm={base_nm} ry_c=ry_m)")
    print(f"{'='*60}")
    rx_c_vals = filter_cavity_scan_values(
        np.round(np.linspace(0.70*rx_m, 1.20*rx_m, 7), 4),
        mirror, "rx_c", best_ac, rx_m, ry_m, 3,
    )
    step2_best = None
    for rx_c in rx_c_vals:
        run_id += 1
        row = evaluate(run_id, best_ac, rx_c, ry_m, Nt=3, Nm=base_nm, mirror=mirror)
        if best_overall is None or row["best_peak_score"] > best_overall[1]:
            best_overall = (row, row["best_peak_score"])
        if step2_best is None or row["best_peak_score"] > step2_best[1]:
            step2_best = (row, row["best_peak_score"])
    best_rx_c = step2_best[0]["rx_c"] if step2_best else rx_m
    print(f"  → 最佳 rx_c = {best_rx_c:.3f}")

    # ================================================================
    # Step 3: 扫描 ry_c
    # ================================================================
    print(f"\n{'='*60}")
    print(f"Step 3: 扫描 ry_c (固定 a_c={best_ac:.3f} rx_c={best_rx_c:.3f} Nt=3 Nm={base_nm})")
    print(f"{'='*60}")
    ry_c_vals = filter_cavity_scan_values(
        np.round(np.linspace(0.70*ry_m, 1.20*ry_m, 7), 4),
        mirror, "ry_c", best_ac, best_rx_c, ry_m, 3,
    )
    step3_best = None
    for ry_c in ry_c_vals:
        run_id += 1
        row = evaluate(run_id, best_ac, best_rx_c, ry_c, Nt=3, Nm=base_nm, mirror=mirror)
        if best_overall is None or row["best_peak_score"] > best_overall[1]:
            best_overall = (row, row["best_peak_score"])
        if step3_best is None or row["best_peak_score"] > step3_best[1]:
            step3_best = (row, row["best_peak_score"])
    best_ry_c = step3_best[0]["ry_c"] if step3_best else ry_m
    print(f"  → 最佳 ry_c = {best_ry_c:.3f}")

    # ================================================================
    # Step 4: 局部联合微调
    # ================================================================
    print(f"\n{'='*60}")
    print(f"Step 4: 局部微调 a_c, rx_c, ry_c")
    print(f"{'='*60}")
    for a_step in [0.005, 0.003]:
        for param, vals in [("a_c", np.round(np.linspace(best_ac-2*a_step, best_ac+2*a_step, 5), 4)),
                            ("rx_c", np.round(np.linspace(best_rx_c-2*a_step, best_rx_c+2*a_step, 5), 4)),
                            ("ry_c", np.round(np.linspace(best_ry_c-2*a_step, best_ry_c+2*a_step, 5), 4))]:
            vals = filter_cavity_scan_values(vals, mirror, param, best_ac, best_rx_c, best_ry_c, 3)
            for val in vals:
                ac = val if param == "a_c" else best_ac
                rxc = val if param == "rx_c" else best_rx_c
                ryc = val if param == "ry_c" else best_ry_c
                run_id += 1
                row = evaluate(run_id, ac, rxc, ryc, Nt=3, Nm=base_nm, mirror=mirror)
                if best_overall is None or row["best_peak_score"] > best_overall[1]:
                    best_overall = (row, row["best_peak_score"])

    # 更新最佳参数
    if best_overall:
        best_ac, best_rx_c, best_ry_c = best_overall[0]["a_c"], best_overall[0]["rx_c"], best_overall[0]["ry_c"]

    # ================================================================
    # Step 5: 扫描 N_taper
    # ================================================================
    print(f"\n{'='*60}")
    print(f"Step 5: 扫描 N_taper (a_c={best_ac:.3f} rx_c={best_rx_c:.3f} ry_c={best_ry_c:.3f})")
    print(f"{'='*60}")
    for Nt in [2, 3, 4, 5, 6]:
        run_id += 1
        row = evaluate(run_id, best_ac, best_rx_c, best_ry_c, Nt=Nt, Nm=base_nm, mirror=mirror)
        if best_overall is None or row["best_peak_score"] > best_overall[1]:
            best_overall = (row, row["best_peak_score"])
    best_Nt = best_overall[0]["N_taper"] if best_overall else 3

    # ================================================================
    # Step 6: 扫描 N_mirror (增加镜区以提升 Q)
    # ================================================================
    print(f"\n{'='*60}")
    print(f"Step 6: 扫描 N_mirror (增加镜区孔数)")
    print(f"{'='*60}")
    best_Nm = base_nm
    for Nm in [10, 14, 18, 22, 26, 30, 34]:
        run_id += 1
        row = evaluate(run_id, best_ac, best_rx_c, best_ry_c, Nt=best_Nt, Nm=Nm, mirror=mirror)
        if best_overall is None or row["best_peak_score"] > best_overall[1]:
            best_overall = (row, row["best_peak_score"])
            best_Nm = Nm

    # ================================================================
    # Step 7: 对准 1550 nm 精细扫描 a_c
    # ================================================================
    print(f"\n{'='*60}")
    print(f"Step 7: 精细扫描 a_c 对准 1550 nm (Nm={best_Nm})")
    print(f"{'='*60}")
    fine_ac = filter_cavity_scan_values(
        np.round(np.linspace(best_ac - 0.04, best_ac + 0.04, 17), 4),
        mirror, "a_c", best_ac, best_rx_c, best_ry_c, best_Nt,
    )
    for a_c in fine_ac:
        run_id += 1
        row = evaluate(run_id, a_c, best_rx_c, best_ry_c, Nt=best_Nt, Nm=best_Nm, mirror=mirror)
        if best_overall is None or row["best_peak_score"] > best_overall[1]:
            best_overall = (row, row["best_peak_score"])
            best_ac = a_c

    # ================================================================
    # Step 8: 高 Q 搜索 — 增加镜区 + 超精细 a_c
    # ================================================================
    print(f"\n{'='*60}")
    print(f"Step 8: 高 Q 搜索 (Nm={best_Nm}, 目标 Q>=1000, λ≈1550nm)")
    print(f"{'='*60}")
    if best_overall:
        best_ac = best_overall[0]["a_c"]
        best_rx_c = best_overall[0]["rx_c"]
        best_ry_c = best_overall[0]["ry_c"]
        best_Nt = best_overall[0]["N_taper"]
        best_Nm = best_overall[0]["N_mirror"]

    for Nm in range(max(best_Nm, 24), 41, 2):
        run_id += 1
        row = evaluate(run_id, best_ac, best_rx_c, best_ry_c, Nt=best_Nt, Nm=Nm, mirror=mirror)
        if best_overall is None or row["best_peak_score"] > best_overall[1]:
            best_overall = (row, row["best_peak_score"])
            best_Nm = Nm
        if float(row.get("Q", 0) or 0) >= 1000 and abs(float(row.get("lambda0_nm", 0) or 0) - 1550) <= 15:
            print(f"  🎯 达标: Q={row['Q']:.0f} λ={row['lambda0_nm']:.1f}nm")
            break

    ultra_ac = filter_cavity_scan_values(
        np.round(np.linspace(best_ac - 0.015, best_ac + 0.015, 13), 4),
        mirror, "a_c", best_ac, best_rx_c, best_ry_c, best_Nt,
    )
    for a_c in ultra_ac:
        run_id += 1
        row = evaluate(run_id, a_c, best_rx_c, best_ry_c, Nt=best_Nt, Nm=best_Nm, mirror=mirror)
        if best_overall is None or row["best_peak_score"] > best_overall[1]:
            best_overall = (row, row["best_peak_score"])
            best_ac = a_c
        if float(row.get("Q", 0) or 0) >= 1000 and abs(float(row.get("lambda0_nm", 0) or 0) - 1550) <= 10:
            print(f"  🎯 达标: Q={row['Q']:.0f} λ={row['lambda0_nm']:.1f}nm")
            break

    # ================================================================
    # 输出最佳结果
    # ================================================================
    # 从 CSV 读取所有结果
    import csv
    all_rows = []
    with open("results/cavity_optimization_all.csv") as f:
        reader = csv.DictReader(f)
        for r in reader:
            q = float(r.get("Q", 0) or 0)
            r2 = float(r.get("fit_r2", 0) or 0)
            if q > 0 and r2 > 0.95:
                all_rows.append(r)

    all_rows.sort(key=lambda x: (
        -float(x.get("best_peak_score", -100) or -100),
        -float(x.get("Q", 0) or 0),
        -float(x.get("T_peak", 0) or 0),
        int(x.get("N_total", 999) or 999)
    ))

    # 保存 top 20
    with open("results/cavity_optimization_best.csv", "w") as f:
        f.write(CSV_HEADER)
        seen = set()
        written = 0
        for r in all_rows:
            rid = r["run_id"]
            if rid in seen:
                continue
            seen.add(rid)
            written += 1
            if written > 20:
                break
            f.write(
                f"{r['run_id']},{r.get('source_mirror_candidate_id','default')},"
                f"{r['a_m']},{r['rx_m']},{r['ry_m']},"
                f"{r['a_c']},{r['rx_c']},{r['ry_c']},"
                f"{r['N_taper']},{r['N_mirror']},"
                f"{r.get('defect_gap','0')},{r['N_total']},"
                f"{r['valid_geometry']},"
                f"{r['lambda0_nm']},{r['FWHM_nm']},{r['Q']},{r['T_peak']},"
                f"{r.get('prominence','0')},{r['fit_r2']},"
                f"{r['best_peak_score']},{r.get('highest_q_peak_Q','0')},"
                f"{r['spectrum_csv']},{r['spectrum_png']},{r['fit_png']},"
                f"{r.get('error_message','')}\n"
            )

    # 保存 best JSON
    if all_rows:
        best = all_rows[0]
        design = {
            "a_m": float(best["a_m"]), "rx_m": float(best["rx_m"]), "ry_m": float(best["ry_m"]),
            "a_c": float(best["a_c"]), "rx_c": float(best["rx_c"]), "ry_c": float(best["ry_c"]),
            "N_taper": int(best["N_taper"]), "N_mirror": int(best["N_mirror"]),
            "defect_gap": int(best.get("defect_gap", 0)),
            "N_total": int(best["N_total"]),
            "lambda0_nm": float(best["lambda0_nm"]),
            "FWHM_nm": float(best["FWHM_nm"]),
            "Q": float(best["Q"]),
            "T_peak": float(best["T_peak"]),
            "fit_r2": float(best["fit_r2"]),
        }
        with open("results/best_cavity_design.json", "w") as f:
            json.dump(design, f, indent=2)

    return all_rows


# ============================================================================
# 主程序
# ============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("一维光子晶体二次渐变缺陷腔自动优化器 (2D)")
    print("=" * 70)

    mc_path = "results_3d/best_bandgap_3d.json"
    if not os.path.exists(mc_path):
        mc_path = "results/best_bandgap_params.json"
    if not os.path.exists(mc_path):
        mc_path = "results_3d/optimized_params.json"
    mirror = load_mirror_candidate(mc_path)
    print(f"镜区: a_m={mirror['a_m']:.3f} rx_m={mirror['rx_m']:.3f} "
          f"ry_m={mirror['ry_m']:.3f}")
    print(f"禁带: {mirror.get('bandgap_start_nm',1529):.0f}-"
          f"{mirror.get('bandgap_end_nm',1589):.0f}nm")
    print(f"最小特征尺寸限制: {min_feature_nm}nm")
    if not mirror.get("geometry_valid", True):
        print(f"⚠️ 镜区参数不满足最小尺寸约束: {', '.join(mirror.get('geometry_violations', []))}")
        print("   请先重新运行 bandgap-2d / bandgap-3d 或更新 best_bandgap_params.json")
        sys.exit(1)
    print()

    results = optimize(mirror)

    print(f"\n{'='*70}")
    print("🏆 缺陷腔优化完成!")
    print(f"{'='*70}")
    q1k = sum(1 for r in results if float(r.get("Q",0) or 0) > 1000)
    print(f"  总仿真数: {len(results)}")
    print(f"  Q > 1000: {q1k}")
    if results:
        b = results[0]
        print(f"\n  最佳设计:")
        print(f"    a_c={b['a_c']} rx_c={b['rx_c']} ry_c={b['ry_c']}")
        print(f"    N_taper={b['N_taper']} N_mirror={b['N_mirror']} N_total={b['N_total']}")
        print(f"    λ₀={float(b['lambda0_nm']):.1f}nm  Q={float(b['Q']):.0f}")
        print(f"    T_peak={float(b['T_peak']):.3f}  R²={float(b['fit_r2']):.4f}")
    print(f"\n  输出文件:")
    print(f"    results/cavity_optimization_all.csv")
    print(f"    results/cavity_optimization_best.csv")
    print(f"    results/best_cavity_design.json")
    print(f"    results/figures/cavity_run_*.png")
    print(f"    results/spectra/cavity_run_*.csv")
    print(f"{'='*70}")
