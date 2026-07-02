#!/usr/bin/env python3
"""Targeted scan: Q>=1000 with lambda near 1550 nm."""
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from optimize_cavity_2d import (
    evaluate,
    filter_cavity_scan_values,
    load_mirror_candidate,
    row_to_design,
    select_design_candidates,
)


def main():
    mirror = load_mirror_candidate("results/best_bandgap_params.json")
    rows = list(csv.DictReader(open("results/cavity_optimization_all.csv")))
    run_id = max(int(r["run_id"]) for r in rows)

    rxc, ryc, nt = 0.2303, 0.2067, 3
    best_hit = None

    print("Targeted scan: high N_mirror + a_c near 1550 nm")
    for nm in [40, 42, 44, 46, 48]:
        ac_vals = filter_cavity_scan_values(
            np.round(np.linspace(0.785, 0.815, 11), 4),
            mirror, "a_c", 0.80, rxc, ryc, nt,
        )
        for a_c in ac_vals:
            run_id += 1
            row = evaluate(run_id, a_c, rxc, ryc, Nt=nt, Nm=nm, mirror=mirror)
            q = float(row.get("Q", 0) or 0)
            lam = float(row.get("lambda0_nm", 0) or 0)
            tp = float(row.get("T_peak", 0) or 0)
            print(f"  -> Q={q:.0f} λ={lam:.1f} T={tp:.3f} Nm={nm} a_c={a_c}")
            if q >= 1000 and abs(lam - 1550) <= 25:
                best_hit = row
                print(f"🎯 TARGET HIT Q={q:.0f} λ={lam:.1f}nm T={tp:.3f}")
                break
        if best_hit:
            break

    rows = []
    with open("results/cavity_optimization_all.csv") as f:
        for r in csv.DictReader(f):
            if float(r.get("Q", 0) or 0) > 0 and float(r.get("fit_r2", 0) or 0) > 0.95:
                rows.append(r)
    best_target, best_q, best_score = select_design_candidates(rows)
    pick = best_hit or best_target or best_q
    for name, row in [
        ("best_cavity_design.json", pick),
        ("best_cavity_design_highQ.json", best_q),
        ("best_cavity_design_balanced.json", best_score),
    ]:
        with open(f"results/{name}", "w") as f:
            json.dump(row_to_design(row), f, indent=2)
        print(f"Saved {name}: Q={row['Q']} λ={row['lambda0_nm']} T={row['T_peak']}")


if __name__ == "__main__":
    main()
