#!/usr/bin/env python3
"""
MEEP 一维光子晶体优化 — 主入口

用法:
  python run.py bandgap-3d     # 3D 禁带坐标下降优化
  python run.py cavity-2d      # 2D 缺陷腔 6 步优化 (Q, T_peak)
  python run.py cavity-3d      # 3D 缺陷腔验证
  python run.py all            # 依次执行 bandgap → cavity-2d → cavity-3d
"""
import argparse
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(ROOT, "src")


def run_script(name):
    path = os.path.join(SRC, name)
    if not os.path.isfile(path):
        print(f"Error: {path} not found")
        sys.exit(1)
    env = os.environ.copy()
    env["PYTHONPATH"] = SRC + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.call([sys.executable, path], cwd=ROOT, env=env)


def main():
    parser = argparse.ArgumentParser(description="PHC optimizer pipeline")
    parser.add_argument(
        "command",
        choices=["bandgap-3d", "cavity-2d", "cavity-3d", "scan-ac", "all"],
        help="optimization stage to run",
    )
    args = parser.parse_args()

    stages = {
        "bandgap-3d": ["optimize_3d_ridge.py"],
        "cavity-2d": ["optimize_cavity_2d.py"],
        "cavity-3d": ["run_3d_cavity.py"],
        "scan-ac": ["scan_cavity_ac.py"],
    }

    if args.command == "all":
        for cmd in ["bandgap-3d", "cavity-2d", "cavity-3d"]:
            print(f"\n{'='*70}\n>>> {cmd}\n{'='*70}")
            for script in stages[cmd]:
                rc = run_script(script)
                if rc != 0:
                    sys.exit(rc)
    else:
        for script in stages[args.command]:
            rc = run_script(script)
            if rc != 0:
                sys.exit(rc)


if __name__ == "__main__":
    main()
