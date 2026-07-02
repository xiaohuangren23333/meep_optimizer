"""Shared 3D ridge waveguide geometry and FDTD helpers."""

import meep as mp
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from config import (
    Q_target,
    T_peak_target,
    T_min_target,
    df,
    dpml,
    fcen,
    h_ridge,
    h_slab,
    h_total,
    lambda_target,
    n_air,
    n_sub,
    n_wg,
    pad,
    w_wg,
)


def build_ref_geom_3d():
  return [
      mp.Block(material=mp.Medium(index=n_sub), center=mp.Vector3(0, 0, -0.5),
               size=mp.Vector3(mp.inf, mp.inf, 1.0)),
      mp.Block(material=mp.Medium(index=n_wg), center=mp.Vector3(0, 0, h_slab / 2),
               size=mp.Vector3(mp.inf, mp.inf, h_slab)),
      mp.Block(material=mp.Medium(index=n_wg),
               center=mp.Vector3(0, 0, h_slab + h_ridge / 2),
               size=mp.Vector3(mp.inf, w_wg, h_ridge)),
  ]


def build_periodic_geom_3d(a, rx, ry, n_period):
    sx = 2 * dpml + n_period * a + 2 * pad
    sy = 2 * dpml + w_wg + 2 * pad
    sz = 2 * dpml + h_total + 1.5
    geom = list(build_ref_geom_3d())
    x0 = -(n_period - 1) * a / 2.0
    for i in range(n_period):
        geom.append(mp.Ellipsoid(
            material=mp.Medium(index=n_air),
            center=mp.Vector3(x0 + i * a, 0, h_total / 2),
            size=mp.Vector3(2 * rx, 2 * ry, h_total),
        ))
    return geom, mp.Vector3(sx, sy, sz)


def build_cavity_geom_3d(mirror, cavity):
    a_m = mirror["a_m"]
    rx_m = mirror["rx_m"]
    ry_m = mirror["ry_m"]
    a_c = cavity["a_c"]
    rx_c = cavity["rx_c"]
    ry_c = cavity["ry_c"]
    nt = int(cavity.get("N_taper", 3))
    nm = int(cavity.get("N_mirror", 10))

    total_holes = 2 * (nm + nt)
    sx = 2 * dpml + max(total_holes * a_m, 4.0) + 2 * pad
    sy = 2 * dpml + w_wg + 2 * pad
    sz = 2 * dpml + h_total + 1.5
    geom = list(build_ref_geom_3d())

    for i in range(nm):
        x = -(nm - i + nt) * a_m
        geom.append(mp.Ellipsoid(material=mp.Medium(index=n_air),
                                center=mp.Vector3(x, 0, h_total / 2),
                                size=mp.Vector3(2 * rx_m, 2 * ry_m, h_total)))
    for i in range(1, nt + 1):
        t = i / nt
        ai = a_c + (a_m - a_c) * t**2
        rxi = rx_c + (rx_m - rx_c) * t**2
        ryi = ry_c + (ry_m - ry_c) * t**2
        x = -(nt - i + 0.5) * ai
        geom.append(mp.Ellipsoid(material=mp.Medium(index=n_air),
                                center=mp.Vector3(x, 0, h_total / 2),
                                size=mp.Vector3(2 * rxi, 2 * ryi, h_total)))
    for i in range(1, nt + 1):
        t = i / nt
        ai = a_c + (a_m - a_c) * t**2
        rxi = rx_c + (rx_m - rx_c) * t**2
        ryi = ry_c + (ry_m - ry_c) * t**2
        x = (nt - i + 0.5) * ai
        geom.append(mp.Ellipsoid(material=mp.Medium(index=n_air),
                                center=mp.Vector3(x, 0, h_total / 2),
                                size=mp.Vector3(2 * rxi, 2 * ryi, h_total)))
    for i in range(nm):
        x = (nm - i + nt) * a_m
        geom.append(mp.Ellipsoid(material=mp.Medium(index=n_air),
                                center=mp.Vector3(x, 0, h_total / 2),
                                size=mp.Vector3(2 * rx_m, 2 * ry_m, h_total)))
    return geom, sx, sy, sz


def run_flux_sim_3d(geom, cell, nfreq, resolution=16, decay=1e-3, decay_time=30,
                      fixed_until=None):
    sx = cell.x if hasattr(cell, "x") else cell[0]
    src_x = -sx / 2 + dpml + 0.5
    mon_x = sx / 2 - dpml - 0.5
    mon_point = mp.Vector3(mon_x, 0, h_total / 2)
    cell_v = mp.Vector3(sx, cell.y if hasattr(cell, "y") else cell[1],
                        cell.z if hasattr(cell, "z") else cell[2])

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ey,
                         center=mp.Vector3(src_x, 0, h_total / 2),
                         size=mp.Vector3(0, w_wg, h_total))]
    sim = mp.Simulation(cell_size=cell_v, resolution=resolution, geometry=geom,
                        sources=sources, boundary_layers=[mp.PML(dpml)])
    fr = mp.FluxRegion(center=mon_point, size=mp.Vector3(0, 2 * w_wg, 2 * h_total))
    tr = sim.add_flux(fcen, df, nfreq, fr)
    if fixed_until is not None:
        sim.run(until_after_sources=fixed_until)
    else:
        sim.run(until_after_sources=mp.stop_when_fields_decayed(
            decay_time, mp.Ey, mon_point, decay))
    freqs = np.array(mp.get_flux_freqs(tr))
    flux = np.array(mp.get_fluxes(tr))
    sim.reset_meep()
    return freqs, flux


def normalized_spectrum(freqs_h, flux_h, freqs_r, flux_r):
    tr = np.divide(flux_h, flux_r, out=np.zeros_like(flux_h), where=flux_r > 1e-15)
    wl = 1000.0 / freqs_r
    order = np.argsort(wl)
    return wl[order], tr[order]


def lorentzian(x, x0, gamma, a, offset):
    return offset + a * gamma**2 / ((x - x0)**2 + gamma**2)


def analyze_bandgap_3d(wl_nm, transmission):
    mask = (wl_nm >= 1450.0) & (wl_nm <= 1650.0)
    wl = wl_nm[mask]
    tr = transmission[mask]
    if wl.size < 8:
        return {"has_gap": False, "gap_start_nm": 0.0, "gap_end_nm": 0.0,
                "gap_width_nm": 0.0, "T_min_gap": float(np.min(transmission)),
                "center_nm": 0.0}

    best = {"has_gap": False, "gap_start_nm": 0.0, "gap_end_nm": 0.0,
            "gap_width_nm": 0.0, "T_min_gap": float(np.min(tr)),
            "center_nm": 0.0}
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
        tmin = float(np.min(transmission[in_gap])) if np.any(in_gap) else float(np.min(tr))
        best = {"has_gap": True, "gap_start_nm": s_nm, "gap_end_nm": e_nm,
                "gap_width_nm": width, "T_min_gap": tmin,
                "center_nm": 0.5 * (s_nm + e_nm)}
    return best


def fit_cavity_peak(wl, transmission, bandgap=None):
    if bandgap:
        lo = bandgap.get("gap_start_nm", 1480) - 20
        hi = bandgap.get("gap_end_nm", 1620) + 20
    else:
        lo, hi = 1480, 1620
    mask = (wl >= lo) & (wl <= hi)
    wl_f, t_f = wl[mask], transmission[mask]
    if wl_f.size < 20:
        return None

    peaks, _ = find_peaks(t_f, height=0.005, prominence=0.003, width=1, distance=2)
    if peaks.size == 0:
        return None

    best = None
    for idx in peaks:
        left = max(int(idx) - 12, 0)
        right = min(int(idx) + 12, len(wl_f) - 1)
        seg_wl, seg_t = wl_f[left:right + 1], t_f[left:right + 1]
        if seg_wl.size < 7:
            continue
        try:
            x0g = wl_f[idx]
            popt, _ = curve_fit(
                lorentzian, seg_wl, seg_t,
                p0=[x0g, 5.0, max(seg_t) - min(seg_t), min(seg_t)],
                maxfev=8000,
            )
            x0, gamma, amp, offset = popt
            if gamma <= 0:
                continue
            fwhm = 2 * gamma
            q = x0 / fwhm if fwhm > 0 else 0
            t_peak = float(lorentzian(x0, *popt))
            resid = seg_t - lorentzian(seg_wl, *popt)
            ss_res = float(np.sum(resid**2))
            ss_tot = float(np.sum((seg_t - np.mean(seg_t))**2))
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            cand = {"lambda0_nm": float(x0), "FWHM_nm": float(fwhm), "Q": float(q),
                    "T_peak": t_peak, "fit_r2": float(r2)}
            if best is None or cand["Q"] > best["Q"]:
                best = cand
        except Exception:
            continue
    return best


def score_cavity_3d(peak):
    if not peak:
        return -1e6
    s = 0.0
    dl = abs(peak["lambda0_nm"] - lambda_target * 1000)
    if dl <= 10:
        s += 120
    elif dl <= 20:
        s += 70
    elif dl <= 40:
        s += 20
    else:
        s -= dl
    if peak["Q"] >= Q_target:
        s += 250
    elif peak["Q"] >= 500:
        s += 100
    s += min(peak["Q"] / 40.0, 120)
    s += min(peak["T_peak"] * 80, 80)
    s += peak.get("fit_r2", 0) * 20
    return s


def meets_3d_targets(peak, bandgap=None, lambda_tol_nm=20, t_peak_min=None):
    if not peak:
        return False, ["no peak"]
    tpk = t_peak_min if t_peak_min is not None else T_peak_target * 0.5
    reasons = []
    if peak["Q"] < Q_target:
        reasons.append(f"Q={peak['Q']:.0f}<{Q_target}")
    if abs(peak["lambda0_nm"] - lambda_target * 1000) > lambda_tol_nm:
        reasons.append(f"lambda={peak['lambda0_nm']:.1f}nm")
    if peak["T_peak"] < tpk:
        reasons.append(f"T_peak={peak['T_peak']:.3f}<{tpk}")
    if bandgap and bandgap.get("has_gap") and bandgap.get("T_min_gap", 1) > T_min_target:
        reasons.append(f"T_min_gap={bandgap['T_min_gap']:.3f}>{T_min_target}")
    return len(reasons) == 0, reasons
