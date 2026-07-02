#!/usr/bin/env python3
"""Resume high-Q refinement from existing cavity optimization CSV."""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from optimize_cavity_2d import (
    load_mirror_candidate,
    select_design_candidates,
    row_to_design,
    refine_high_q,
)


def main():
    csv_path = "results/cavity_optimization_all.csv"
    if not os.path.exists(csv_path):
        print("Missing", csv_path)
        sys.exit(1)

    mirror = load_mirror_candidate("results/best_bandgap_params.json")
    rows = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            if float(r.get("Q", 0) or 0) > 0 and float(r.get("fit_r2", 0) or 0) > 0.95:
                rows.append(r)

    _, best_q, best_score = select_design_candidates(rows)
    print(f"Current best Q: {best_q['Q']} @ {best_q['lambda0_nm']} nm (Nm={best_q['N_mirror']})")

    run_id = max(int(r["run_id"]) for r in rows)
    refine_high_q(mirror, best_q, run_id)

    rows = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            if float(r.get("Q", 0) or 0) > 0 and float(r.get("fit_r2", 0) or 0) > 0.95:
                rows.append(r)

    best_target, best_q, best_score = select_design_candidates(rows)
    for name, row in [
        ("best_cavity_design.json", best_target),
        ("best_cavity_design_highQ.json", best_q),
        ("best_cavity_design_balanced.json", best_score),
    ]:
        with open(f"results/{name}", "w") as f:
            json.dump(row_to_design(row), f, indent=2)
        print(f"Wrote results/{name}: Q={row['Q']} λ={row['lambda0_nm']} T={row['T_peak']}")


if __name__ == "__main__":
    main()
