"""
透射谱归一化校验 + Lorentzian 峰拟合（2D/3D 共用）

可靠性要求（AGENT_INSTRUCTION）:
  - 必须用无孔参考仿真归一化，且参考 cell 与有孔结构相同
  - 拟合参数带物理约束（gamma>0, FWHM 合理）
  - 输出原始 flux 与校验元数据，便于 FDTD 独立验证
"""
import json
import os
from typing import Optional

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


def lorentzian(x, x0, gamma, A, offset):
    return offset + A * gamma**2 / ((x - x0) ** 2 + gamma**2)


def normalize_transmission(freqs_h, flux_h, freqs_r, flux_r):
    """同频率网格归一化；freq 数组应来自相同 fcen/df/nfreq 设置。"""
    freqs_h = np.asarray(freqs_h, dtype=float)
    flux_h = np.asarray(flux_h, dtype=float)
    freqs_r = np.asarray(freqs_r, dtype=float)
    flux_r = np.asarray(flux_r, dtype=float)

    if len(freqs_h) != len(freqs_r) or np.max(np.abs(freqs_h - freqs_r)) > 1e-12:
        raise ValueError("有孔/参考仿真频率网格不一致，归一化无效")

    T = np.divide(flux_h, flux_r, out=np.zeros_like(flux_h), where=flux_r > 1e-15)
    wl = 1000.0 / freqs_h
    idx = np.argsort(wl)
    return wl[idx], T[idx], flux_h[idx], flux_r[idx]


def validate_normalized_spectrum(wl, T, flux_h=None, flux_r=None, label=""):
    """返回校验 dict；issues 非空表示数据可能不可靠。"""
    wl = np.asarray(wl)
    T = np.asarray(T)
    issues = []
    warnings = []

    if len(wl) < 50:
        issues.append("频率点过少")
    if np.any(T < -0.02):
        issues.append(f"T 出现负值 (min={T.min():.4f})")
    if np.max(T) > 1.05:
        warnings.append(f"T_max={T.max():.4f}>1.05，检查参考仿真或 PML")
    if np.max(T) > 1.2:
        issues.append(f"T_max={T.max():.4f} 严重超 1，归一化错误")

    if flux_h is not None and flux_r is not None:
        ratio = np.sum(flux_h) / max(np.sum(flux_r), 1e-30)
        if ratio > 0.98:
            warnings.append(f"sum(flux_h)/sum(flux_r)={ratio:.3f}≈1，结构可能未形成有效腔")

    band = (wl >= 1450) & (wl <= 1650)
    if band.sum() > 0:
        t_band = T[band]
        if t_band.max() < 0.02:
            issues.append("1450-1650nm 无可见透射峰")
        if t_band.min() > 0.5:
            warnings.append("禁带区 T_min 偏高，镜区反射可能不足")

    return {
        "label": label,
        "n_freq": int(len(wl)),
        "wl_min_nm": float(wl.min()),
        "wl_max_nm": float(wl.max()),
        "T_min": float(T.min()),
        "T_max": float(T.max()),
        "T_max_wl_nm": float(wl[np.argmax(T)]),
        "issues": issues,
        "warnings": warnings,
        "reliable": len(issues) == 0,
    }


def fit_peak(wl_seg, T_seg, delta_wl=None):
    """
    Lorentzian 拟合，带 bounds 防止负 gamma / 负 FWHM。
    Q = lambda0 / FWHM,  FWHM = 2*gamma (nm)
    """
    wl_seg = np.asarray(wl_seg, dtype=float)
    T_seg = np.asarray(T_seg, dtype=float)
    if len(wl_seg) < 9:
        return None

    if delta_wl is None:
        delta_wl = (wl_seg[-1] - wl_seg[0]) / max(len(wl_seg) - 1, 1)

    x0_g = wl_seg[np.argmax(T_seg)]
    offset_g = float(np.min(T_seg))
    A_g = max(float(np.max(T_seg) - offset_g), 0.01)
    span = max(wl_seg[-1] - wl_seg[0], 3 * delta_wl)
    gamma_g = max(span / 8, delta_wl)

    try:
        popt, _ = curve_fit(
            lorentzian,
            wl_seg,
            T_seg,
            p0=[x0_g, gamma_g, A_g, offset_g],
            bounds=(
                [wl_seg.min(), delta_wl * 0.5, 0.0, -0.05],
                [wl_seg.max(), span, 2.0, 1.0],
            ),
            maxfev=8000,
        )
    except Exception:
        return None

    x0, gamma, A, offset = popt
    if gamma <= 0 or A < 0:
        return None

    FWHM = 2.0 * gamma
    T_peak = float(lorentzian(x0, *popt))
    residuals = T_seg - lorentzian(wl_seg, *popt)
    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((T_seg - np.mean(T_seg)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return {
        "lambda0_nm": float(x0),
        "FWHM_nm": float(FWHM),
        "gamma_nm": float(gamma),
        "Q": float(x0 / FWHM) if FWHM > 0 else 0.0,
        "T_peak": T_peak,
        "A": float(A),
        "offset": float(offset),
        "fit_r2": float(r2),
        "fit_window_nm": float(span),
    }


def extract_peaks(
    wl,
    T,
    bandgap_start=1530.0,
    bandgap_end=1565.0,
    target_lambda=1550.0,
    delta_wl=None,
    min_T_peak=0.05,
    min_r2=0.90,
    min_fwhm_factor=2.0,
):
    """在禁带区内找峰并拟合；返回按综合评分排序的候选列表。"""
    wl = np.asarray(wl)
    T = np.asarray(T)
    if delta_wl is None:
        delta_wl = (wl[-1] - wl[0]) / max(len(wl) - 1, 1)

    mask = (wl >= bandgap_start - 15) & (wl <= bandgap_end + 15)
    wl_bg, T_bg = wl[mask], T[mask]
    if len(wl_bg) < 30:
        return []

    peaks, props = find_peaks(
        T_bg, height=0.03, prominence=0.015, width=2, distance=5
    )
    if len(peaks) == 0:
        return []

    candidates = []
    for k, pidx in enumerate(peaks):
        # 自适应窗口：至少 ±15 点，或 ±3*估计半宽
        half_pts = max(15, int(3.0 / delta_wl) if delta_wl > 0 else 15)
        left = max(pidx - half_pts, 0)
        right = min(pidx + half_pts, len(wl_bg) - 1)
        fit = fit_peak(wl_bg[left : right + 1], T_bg[left : right + 1], delta_wl)
        if fit is None:
            continue

        x0 = fit["lambda0_nm"]
        valid = True
        reasons = []
        if x0 < bandgap_start or x0 > bandgap_end:
            valid = False
            reasons.append("不在镜区禁带内")
        if fit["T_peak"] < min_T_peak:
            valid = False
            reasons.append(f"T_peak<{min_T_peak}")
        if fit["FWHM_nm"] < min_fwhm_factor * delta_wl:
            valid = False
            reasons.append("FWHM 窄于频率分辨率")
        if fit["fit_r2"] < min_r2:
            valid = False
            reasons.append(f"R²={fit['fit_r2']:.3f}<{min_r2}")
        if fit["Q"] < 1 or fit["FWHM_nm"] <= 0:
            valid = False
            reasons.append("拟合参数非物理")

        prom = float(props["prominences"][k]) if k < len(props["prominences"]) else 0.0
        dl = abs(x0 - target_lambda)
        score = (
            (100 if valid else -50)
            + min(fit["Q"] / 100, 40)
            + min(fit["T_peak"] * 25, 25)
            + fit["fit_r2"] * 10
            + max(0, 20 - dl / 2)
        )

        candidates.append({
            **fit,
            "prominence": prom,
            "valid": valid,
            "invalid_reason": "; ".join(reasons),
            "score": score,
            "in_bandgap": bandgap_start <= x0 <= bandgap_end,
        })

    candidates.sort(key=lambda p: (-p["score"], -p["Q"], -p["T_peak"]))
    return candidates


def load_bandgap_from_json(path=None):
    if path is None:
        for p in (
            "results/best_bandgap_params.json",
            "results_3d/bandgap_verify.json",
            "results_3d/optimized_params.json",
        ):
            if os.path.isfile(p):
                path = p
                break
    if not path or not os.path.isfile(path):
        return {"bandgap_start_nm": 1530.0, "bandgap_end_nm": 1565.0, "source": "default"}
    with open(path) as f:
        raw = json.load(f)
    return {
        "bandgap_start_nm": float(
            raw.get("gap_start_nm") or raw.get("bandgap_start_nm") or 1530
        ),
        "bandgap_end_nm": float(
            raw.get("gap_end_nm") or raw.get("bandgap_end_nm") or 1565
        ),
        "T_min_mirror": raw.get("T_min"),
        "source": path,
    }


def save_spectrum_bundle(
    out_dir,
    prefix,
    wl,
    T,
    flux_h=None,
    flux_r=None,
    validation=None,
    peaks=None,
    best_peak=None,
    meta=None,
):
    os.makedirs(out_dir, exist_ok=True)
    np.savetxt(
        f"{out_dir}/{prefix}_T.csv",
        np.column_stack([wl, T]),
        delimiter=",",
        header="wavelength_nm,transmission", comments="",
    )
    if flux_h is not None and flux_r is not None:
        np.savetxt(
            f"{out_dir}/{prefix}_flux.csv",
            np.column_stack([wl, flux_h, flux_r, T]),
            delimiter=",",
            header="wavelength_nm,flux_hole,flux_ref,transmission", comments="",
        )
    bundle = {
        "validation": validation,
        "best_peak": best_peak,
        "all_peaks": peaks,
        "simulation": meta or {},
    }
    with open(f"{out_dir}/{prefix}_analysis.json", "w") as f:
        json.dump(bundle, f, indent=2)
    return bundle
