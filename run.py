#!/usr/bin/env python3
"""
MEEP 一维光子晶体优化 — 主入口

用法:
  python run.py bandgap-2d     # 2D 周期结构禁带扫描 (选最优参数)
  python run.py bandgap-3d     # 3D 周期结构验证扫描 (基于 2D 最优)
  python run.py cavity-2d      # 2D 缺陷腔优化 (可选，提供初值)
  python run.py cavity-3d-opt  # 3D 缺陷腔扫描优化 ⭐
  python run.py cavity-3d      # 3D 缺陷腔高分辨率验证
  python run.py all-3d         # 3D 目标流水线
  python run.py all            # 2D+3D 全流程
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
        choices=["bandgap-2d", "bandgap-3d", "cavity-2d", "cavity-3d-opt", "cavity-3d", "all", "all-3d"],
        help="optimization stage to run",
    )
    args = parser.parse_args()

    stages = {
        "bandgap-2d": ["scan_bandgap_2d.py"],
        "bandgap-3d": ["scan_bandgap_3d.py"],
        "cavity-2d": ["optimize_cavity_2d.py"],
        "cavity-3d-opt": ["optimize_cavity_3d.py"],
        "cavity-3d": ["run_3d_cavity.py"],
    }

    if args.command == "all":
        for cmd in ["bandgap-2d", "bandgap-3d", "cavity-2d", "cavity-3d"]:
            print(f"\n{'='*70}\n>>> {cmd}\n{'='*70}")
            for script in stages[cmd]:
                rc = run_script(script)
                if rc != 0:
                    sys.exit(rc)
    elif args.command == "all-3d":
        for cmd in ["bandgap-3d", "cavity-3d-opt", "cavity-3d"]:
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
