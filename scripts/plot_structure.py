#!/usr/bin/env python3
"""Plot the 3D ridge photonic-crystal defect cavity structure (top + side view)."""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from config import w_wg, h_slab, h_ridge, h_total, min_feature_um
from phc_3d import size_graded_layout


def load_design():
    mirror = json.load(open("results_3d/best_bandgap_3d.json"))
    a_m = mirror["a"]; rx_m = mirror["rx"]; ry_m = mirror["ry"]
    for p in ("results_3d/cavity/best_cavity_3d.json",
              "results/best_cavity_design.json"):
        if os.path.exists(p):
            d = json.load(open(p))
            return {
                "a": d.get("a", a_m), "rx": d.get("rx", rx_m),
                "ry_end": d.get("ry_end", ry_m),
                "ry_center": d.get("ry_center", ry_m * 0.6),
                "N_taper": int(d.get("N_taper", 12)),
                "N_mirror": int(d.get("N_mirror", 8)),
                "Q": d.get("Q"), "lam": d.get("lambda0_nm"), "src": p,
            }
    return {"a": a_m, "rx": rx_m, "ry_end": ry_m, "ry_center": ry_m * 0.6,
            "N_taper": 12, "N_mirror": 8, "Q": None, "lam": None, "src": "defaults"}


def hole_positions(design):
    """Size-graded layout: (x, rx, ry, kind)."""
    return size_graded_layout(design["a"], design["rx"], design["ry_end"],
                              design["ry_center"], design["N_taper"],
                              design["N_mirror"])


def main():
    design = load_design()
    a_m = design["a"]; rx_m = design["rx"]; ry_m = design["ry_end"]
    ry_c = design["ry_center"]; Nt = design["N_taper"]; Nm = design["N_mirror"]
    Q = design["Q"]; lam = design["lam"]; src = design["src"]
    holes = hole_positions(design)
    span = max(abs(h[0]) for h in holes) + a_m

    colors = {"mirror": "#c0392b", "taper": "#f39c12"}

    fig, (ax_top, ax_side) = plt.subplots(
        2, 1, figsize=(15, 6), gridspec_kw={"height_ratios": [3, 1]})

    # ---- Top view (x-y plane) ----
    ax_top.add_patch(Rectangle((-span, -w_wg / 2), 2 * span, w_wg,
                               facecolor="#2980b9", alpha=0.30,
                               edgecolor="#1b4f72", lw=1.2, label="ridge (n=2.18)"))
    for x, rx, ry, kind in holes:
        ax_top.add_patch(Ellipse((x, 0), 2 * rx, 2 * ry,
                                 facecolor="white", edgecolor=colors[kind], lw=1.3))
    # legend proxies
    ax_top.add_patch(Ellipse((0, 1e6), 1, 1, facecolor="white",
                             edgecolor=colors["mirror"], label="mirror holes"))
    ax_top.add_patch(Ellipse((0, 1e6), 1, 1, facecolor="white",
                             edgecolor=colors["taper"], label="taper holes"))
    ax_top.axvline(0, color="gray", ls="--", lw=0.8, alpha=0.6)
    ax_top.annotate("cavity center", (0, w_wg / 2 + 0.15), ha="center",
                    fontsize=9, color="gray")
    ax_top.set_xlim(-span, span)
    ax_top.set_ylim(-w_wg / 2 - 0.6, w_wg / 2 + 0.6)
    ax_top.set_aspect("equal")
    ax_top.set_xlabel("x (μm)")
    ax_top.set_ylabel("y (μm)")
    title = (f"3D Ridge PhC Size-Graded Defect Cavity (1st-order gap) — Top view\n"
             f"a={a_m:.3f} rx={rx_m:.3f} ry: {ry_c:.3f}(center)→{ry_m:.3f}(mirror) μm | "
             f"N_taper={Nt} N_mirror={Nm} | total holes={len(holes)}")
    if Q:
        title += f"\n3D: Q={Q:.0f}, λ₀={lam:.1f} nm  (source: {os.path.basename(src)})"
    ax_top.set_title(title, fontsize=10)
    ax_top.legend(loc="upper right", fontsize=8, ncol=3)

    # ---- Side view (x-z plane) ----
    ax_side.add_patch(Rectangle((-span, -1.0), 2 * span, 1.0,
                                facecolor="#7f8c8d", alpha=0.35,
                                edgecolor="none", label="substrate (n=1.44)"))
    ax_side.add_patch(Rectangle((-span, 0), 2 * span, h_slab,
                                facecolor="#2980b9", alpha=0.25,
                                edgecolor="none", label="slab"))
    ax_side.add_patch(Rectangle((-span, h_slab), 2 * span, h_ridge,
                                facecolor="#2980b9", alpha=0.45,
                                edgecolor="#1b4f72", lw=0.8, label="ridge"))
    for x, rx, ry, kind in holes:
        ax_side.add_patch(Rectangle((x - rx, 0), 2 * rx, h_total,
                                    facecolor="white", edgecolor=colors[kind], lw=0.8))
    ax_side.set_xlim(-span, span)
    ax_side.set_ylim(-1.0, h_total + 0.3)
    ax_side.set_aspect("equal")
    ax_side.set_xlabel("x (μm)")
    ax_side.set_ylabel("z (μm)")
    ax_side.set_title("Side view (x-z), holes fully etched through h=0.4 μm", fontsize=9)
    ax_side.legend(loc="upper right", fontsize=7, ncol=4)

    plt.tight_layout()
    os.makedirs("results_3d/figures", exist_ok=True)
    out = "results_3d/figures/cavity_structure.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)

    # ---- Zoomed center view ----
    fig2, ax = plt.subplots(figsize=(12, 3.2))
    ax.add_patch(Rectangle((-span, -w_wg / 2), 2 * span, w_wg,
                           facecolor="#2980b9", alpha=0.30, edgecolor="#1b4f72"))
    right = sorted([h for h in holes if h[0] > 0], key=lambda h: h[0])
    for j, (x, rx, ry, kind) in enumerate(holes):
        ax.add_patch(Ellipse((x, 0), 2 * rx, 2 * ry,
                             facecolor="white", edgecolor=colors[kind], lw=1.5))
    # annotate the graded hole y-diameter for the inner right-side holes
    for j in range(min(len(right), 7)):
        x, rx, ry, kind = right[j]
        ax.annotate(f"{2*ry*1000:.0f}", (x, ry + 0.08), ha="center", fontsize=7,
                    color="#555")
    zoom = min((Nt + 2) * a_m, span)
    ax.set_xlim(-zoom, zoom)
    ax.set_ylim(-w_wg / 2 - 0.3, w_wg / 2 + 0.45)
    ax.set_aspect("equal")
    ax.axvline(0, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("x (μm)")
    ax.set_ylabel("y (μm)")
    ax.set_title(f"Center zoom — hole ry graded {2*ry_c*1000:.0f}→{2*ry_m*1000:.0f} nm "
                 f"(numbers = hole y-diameter); fixed period {a_m*1000:.0f} nm; "
                 f"min feature {min_feature_um*1000:.0f} nm", fontsize=9)
    plt.tight_layout()
    out2 = "results_3d/figures/cavity_center_zoom.png"
    plt.savefig(out2, dpi=150, bbox_inches="tight")
    print("saved", out2)


if __name__ == "__main__":
    main()
