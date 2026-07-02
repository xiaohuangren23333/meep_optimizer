#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
谐振峰 Q = λ_peak / FWHM（可复用、无项目依赖）。

方法:
  - peak_bracketed: FDTD/原始谱，峰附近肩台基线 + 包住峰的半高交点（推荐）
  - transmission_peak: TMM/PPT，透射峰局部谷底基线 + 0.5*(T_peak-T_base)
  - dip_interval: 两谷底之间全区间（易在宽禁带误算，慎用）
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Literal

import numpy as np

Method = Literal["peak_bracketed", "transmission_peak", "dip_interval"]


def _all_crossings(
    wl: np.ndarray,
    T: np.ndarray,
    i_left: int,
    i_right: int,
    level: float,
) -> list[float]:
    wl = np.asarray(wl, dtype=np.float64)
    T = np.asarray(T, dtype=np.float64)
    i_left = int(max(0, i_left))
    i_right = int(min(T.size - 1, i_right))
    xs: list[float] = []
    for i in range(i_left, i_right):
        t0, t1 = float(T[i]), float(T[i + 1])
        if not (np.isfinite(t0) and np.isfinite(t1)):
            continue
        if (t0 - level) * (t1 - level) > 0:
            continue
        if t1 == t0:
            lam = 0.5 * (wl[i] + wl[i + 1])
        else:
            frac = (level - t0) / (t1 - t0)
            lam = float(wl[i] + frac * (wl[i + 1] - wl[i]))
        xs.append(lam)
    return xs


def _idx_near(wl: np.ndarray, lam: float) -> int:
    return int(np.argmin(np.abs(np.asarray(wl, dtype=np.float64) - float(lam))))


def find_peak_in_window(
    wl: np.ndarray,
    T: np.ndarray,
    *,
    lam_center: float,
    half_window: float,
    maximize: bool = True,
) -> dict:
    """在窗口内找峰/谷索引。"""
    wl = np.asarray(wl, dtype=np.float64).reshape(-1)
    T = np.asarray(T, dtype=np.float64).reshape(-1)
    m = (wl >= lam_center - half_window) & (wl <= lam_center + half_window)
    if not np.any(m):
        return {"ok": False, "reason": "empty_window"}
    widx = np.where(m)[0]
    i_local = int(np.argmax(T[widx]) if maximize else np.argmin(T[widx]))
    i = int(widx[i_local])
    return {
        "ok": True,
        "index": i,
        "lambda_um": float(wl[i]),
        "T": float(T[i]),
    }


def peak_fwhm_bracketed(
    wl_um: np.ndarray,
    T: np.ndarray,
    peak_lam_um: float | None = None,
    *,
    half_win_um: float = 0.006,
    shoulder_um: float = 0.003,
    lam_center_um: float = 1.55,
    search_half_window_um: float = 0.05,
) -> dict:
    """
    FDTD 原始谱推荐：峰 ±half_win 内，肩台基线，取包住峰的一对半高交点。
    wl_um 单位 µm。
    """
    fail = {"ok": False, "method": "peak_bracketed"}
    wl_um = np.asarray(wl_um, dtype=np.float64).reshape(-1)
    T = np.asarray(T, dtype=np.float64).reshape(-1)

    if peak_lam_um is None or not np.isfinite(peak_lam_um):
        pk = find_peak_in_window(
            wl_um, T, lam_center=lam_center_um, half_window=search_half_window_um, maximize=True
        )
        if not pk["ok"]:
            return {**fail, "reason": "peak_not_found"}
        peak_lam_um = pk["lambda_um"]

    ip = _idx_near(wl_um, peak_lam_um)
    T_peak = float(T[ip])
    lam_p = float(wl_um[ip])

    lmask = (wl_um >= lam_p - shoulder_um - half_win_um) & (wl_um <= lam_p - shoulder_um)
    rmask = (wl_um >= lam_p + shoulder_um) & (wl_um <= lam_p + shoulder_um + half_win_um)
    bases: list[float] = []
    if np.any(lmask):
        bases.append(float(np.min(T[lmask])))
    if np.any(rmask):
        bases.append(float(np.min(T[rmask])))
    T_base = min(bases) if bases else float(np.min(T))

    if T_peak <= T_base + 1e-12:
        return {**fail, "lambda_peak_um": lam_p, "T_peak": T_peak, "T_base": T_base}

    T_half = T_base + 0.5 * (T_peak - T_base)
    i_lo = int(np.where(wl_um >= lam_p - half_win_um)[0][0])
    i_hi = int(np.where(wl_um <= lam_p + half_win_um)[0][-1])
    xs = sorted(_all_crossings(wl_um, T, i_lo, i_hi, T_half))

    lam_l = lam_r = None
    for i in range(len(xs) - 1):
        if xs[i] <= lam_p <= xs[i + 1]:
            lam_l, lam_r = float(xs[i]), float(xs[i + 1])
            break

    if lam_l is None or lam_r is None:
        return {
            **fail,
            "lambda_peak_um": lam_p,
            "T_peak": T_peak,
            "T_base": T_base,
            "n_crossings": len(xs),
            "reason": "no_bracketing_crossings",
        }

    fwhm = lam_r - lam_l
    if fwhm <= 0:
        return {**fail, "lambda_peak_um": lam_p, "reason": "fwhm_nonpos"}

    dlam_nm = float(np.median(np.diff(np.sort(wl_um)))) * 1e3 if wl_um.size > 1 else np.nan
    return {
        "ok": True,
        "method": "peak_bracketed",
        "Q_peak": lam_p / fwhm,
        "fwhm_um": fwhm,
        "lambda_peak_um": lam_p,
        "T_peak": T_peak,
        "T_base": T_base,
        "T_half": T_half,
        "lambda_L_um": lam_l,
        "lambda_R_um": lam_r,
        "n_crossings": len(xs),
        "sampling_dlam_nm": dlam_nm,
        "fwhm_oversampled": bool(fwhm * 1e3 < 2.5 * dlam_nm) if np.isfinite(dlam_nm) else False,
    }


def transmission_peak_fwhm(
    wl: np.ndarray,
    T: np.ndarray,
    *,
    lam_guess: float | None = None,
    peak_search_half_window: float = 5e-9,
    local_half_window: float = 20e-9,
    wavelength_unit: Literal["m", "um"] = "m",
) -> dict:
    """
    TMM/PPT 透射峰：局部 min 为基线，半高 0.5*(T_peak-T_base)。
    wl 默认米；若 wavelength_unit='um' 则 wl 为 µm，输出 lambda_*_um。
    """
    fail = {"ok": False, "method": "transmission_peak"}
    wl = np.asarray(wl, dtype=np.float64).reshape(-1)
    T = np.asarray(T, dtype=np.float64).reshape(-1)
    scale = 1e6 if wavelength_unit == "m" else 1.0

    if lam_guess is None:
        ip = int(np.argmax(T))
    else:
        sw = peak_search_half_window if wavelength_unit == "m" else peak_search_half_window * 1e-6
        m = (wl > lam_guess - sw) & (wl < lam_guess + sw)
        if not np.any(m):
            ip = int(np.argmax(T))
        else:
            widx = np.where(m)[0]
            ip = int(widx[int(np.argmax(T[widx]))])

    lam_p = float(wl[ip])
    Tp = float(T[ip])
    lw = local_half_window if wavelength_unit == "m" else local_half_window * 1e-6
    m = (wl > lam_p - lw) & (wl < lam_p + lw)
    Tb = float(np.min(T[m])) if np.any(m) else float(np.min(T))
    half = Tb + 0.5 * (Tp - Tb)
    if Tp <= half + 1e-15:
        return {**fail, "lambda_peak_um": lam_p * scale, "T_peak": Tp}

    li = ip
    while li > 0 and T[li] > half:
        li -= 1
    ri = ip
    while ri < len(T) - 1 and T[ri] > half:
        ri += 1
    if li == 0 or ri >= len(T) - 1:
        return {**fail, "lambda_peak_um": lam_p * scale, "T_peak": Tp, "reason": "edge_hit"}

    lam_l = float(np.interp(half, [T[li], T[li + 1]], [wl[li], wl[li + 1]]))
    lam_r = float(np.interp(half, [T[ri - 1], T[ri]], [wl[ri - 1], wl[ri]]))
    fwhm = lam_r - lam_l
    if fwhm <= 0:
        return {**fail, "lambda_peak_um": lam_p * scale, "T_peak": Tp}
    return {
        "ok": True,
        "method": "transmission_peak",
        "Q_peak": lam_p / fwhm,
        "fwhm_um": fwhm * scale,
        "lambda_peak_um": lam_p * scale,
        "T_peak": Tp,
        "T_base": Tb,
        "lambda_L_um": lam_l * scale,
        "lambda_R_um": lam_r * scale,
    }


def dip_interval_fwhm(
    wl_um: np.ndarray,
    T: np.ndarray,
    peak_lam_um: float,
    dip1_lam_um: float,
    dip2_lam_um: float,
) -> dict:
    """两谷底区间内求半高（宽禁带时 FWHM 易偏大）。"""
    fail = {"ok": False, "method": "dip_interval"}
    if not all(np.isfinite(x) for x in (peak_lam_um, dip1_lam_um, dip2_lam_um)):
        return fail
    wl_um = np.asarray(wl_um, dtype=np.float64).reshape(-1)
    T = np.asarray(T, dtype=np.float64).reshape(-1)
    i_p = _idx_near(wl_um, peak_lam_um)
    i_1, i_2 = _idx_near(wl_um, dip1_lam_um), _idx_near(wl_um, dip2_lam_um)
    i_lo, i_hi = sorted((i_1, i_2))
    if i_p < i_lo or i_p > i_hi:
        i_p = int(i_lo + np.argmax(T[i_lo : i_hi + 1]))
    T_peak = float(T[i_p])
    lam_p = float(wl_um[i_p])
    T_base = float(np.min(T[i_lo : i_hi + 1]))
    T_half = T_base + 0.5 * (T_peak - T_base)
    xs = _all_crossings(wl_um, T, i_lo, i_hi, T_half)
    if len(xs) < 2:
        return {**fail, "lambda_peak_um": lam_p, "n_crossings": len(xs)}
    lam_l, lam_r = float(min(xs)), float(max(xs))
    fwhm = lam_r - lam_l
    if fwhm <= 0:
        return {**fail, "lambda_peak_um": lam_p}
    return {
        "ok": True,
        "method": "dip_interval",
        "Q_peak": lam_p / fwhm,
        "fwhm_um": fwhm,
        "lambda_peak_um": lam_p,
        "T_peak": T_peak,
        "T_base": T_base,
        "lambda_L_um": lam_l,
        "lambda_R_um": lam_r,
        "n_crossings": len(xs),
    }


def compute_q(
    wl: np.ndarray,
    y: np.ndarray,
    *,
    method: Method = "peak_bracketed",
    wavelength_unit: Literal["m", "um", "auto"] = "auto",
    **kwargs,
) -> dict:
    """统一入口。auto：|wl|>0.01 视为 µm，否则为 m。"""
    wl = np.asarray(wl, dtype=np.float64).reshape(-1)
    y = np.asarray(y, dtype=np.float64).reshape(-1)
    if wavelength_unit == "auto":
        wavelength_unit = "um" if float(np.nanmax(np.abs(wl))) > 0.01 else "m"

    if method == "peak_bracketed":
        if wavelength_unit == "m":
            return peak_fwhm_bracketed(wl * 1e6, y, **kwargs)
        return peak_fwhm_bracketed(wl, y, **kwargs)
    if method == "transmission_peak":
        return transmission_peak_fwhm(wl, y, wavelength_unit=wavelength_unit, **kwargs)
    if method == "dip_interval":
        if wavelength_unit == "m":
            wl_um, peak = wl * 1e6, kwargs.get("peak_lam_um")
            if peak is not None:
                kwargs = {**kwargs, "peak_lam_um": peak * 1e6}
            return dip_interval_fwhm(
                wl_um,
                y,
                kwargs["peak_lam_um"],
                kwargs["dip1_lam_um"] * 1e6,
                kwargs["dip2_lam_um"] * 1e6,
            )
        return dip_interval_fwhm(wl, y, **kwargs)
    raise ValueError(f"unknown method: {method}")


def _cli() -> int:
    ap = argparse.ArgumentParser(description="Compute resonance Q = lambda_peak / FWHM")
    ap.add_argument("--wl", nargs="+", type=float, help="wavelength samples (with --T)")
    ap.add_argument("--T", nargs="+", type=float, help="spectrum samples")
    ap.add_argument("--npz", type=str, help="npz with wavelength_um/wavelength_m + spectrum/T")
    ap.add_argument("--wl-key", default="")
    ap.add_argument("--y-key", default="")
    ap.add_argument(
        "--method",
        choices=("peak_bracketed", "transmission_peak", "dip_interval"),
        default="peak_bracketed",
    )
    ap.add_argument("--peak-um", type=float, default=None)
    ap.add_argument("--dip1-um", type=float, default=None)
    ap.add_argument("--dip2-um", type=float, default=None)
    ap.add_argument("--json-out", type=str, default="")
    args = ap.parse_args()

    if args.npz:
        z = np.load(args.npz)
        keys = list(z.files)
        wl_key = args.wl_key or next(
            (k for k in ("wavelength_um", "wavelength_m", "lambda", "wl") if k in keys), keys[0]
        )
        y_key = args.y_key or next(
            (k for k in ("spectrum", "spectrum_raw", "T", "T_defect") if k in keys), keys[min(1, len(keys) - 1)]
        )
        wl = np.asarray(z[wl_key]).reshape(-1)
        y = np.asarray(z[y_key]).reshape(-1)
        if wl_key.endswith("_m"):
            wunit = "m"
        elif wl_key.endswith("_um"):
            wunit = "um"
        else:
            wunit = "auto"
    elif args.wl and args.T:
        wl = np.array(args.wl, dtype=np.float64)
        y = np.array(args.T, dtype=np.float64)
        wunit = "auto"
    else:
        ap.error("provide --npz or --wl and --T")
        return 2

    kw: dict = {}
    if args.peak_um is not None:
        kw["peak_lam_um"] = args.peak_um
    if args.method == "dip_interval":
        if args.dip1_um is None or args.dip2_um is None:
            ap.error("dip_interval needs --dip1-um and --dip2-um")
            return 2
        kw["peak_lam_um"] = args.peak_um
        kw["dip1_lam_um"] = args.dip1_um
        kw["dip2_lam_um"] = args.dip2_um

    out = compute_q(wl, y, method=args.method, wavelength_unit=wunit, **kw)
    text = json.dumps({k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in out.items()}, indent=2)
    if args.json_out:
        Path(args.json_out).write_text(text, encoding="utf-8")
    print(text)
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(_cli())
