# MEEP 一维光子晶体优化器

基于 [Meep](https://meep.readthedocs.io/) 的一维光子晶体（椭圆孔脊型波导）禁带扫描与缺陷腔 Q 值优化。

## 目标

| 指标 | 目标 |
|------|------|
| Q 值 | ≥ 1000 |
| 谐振波长 λ₀ | ~1550 nm |
| 透射峰 T_peak | ~1 |
| 禁带底 T_min | ~0 |

## 环境

```bash
conda activate phc-meep
cd ~/projects/meep_optimizer
```

依赖：Python 3、Meep 1.33+、numpy、scipy、matplotlib

## 快速开始

```bash
# 1. 3D 禁带优化（坐标下降，约数小时）
python run.py bandgap-3d

# 2. 2D 缺陷腔优化（6 步扫描，约数小时）
python run.py cavity-2d

# 3. 3D 缺陷腔验证
python run.py cavity-3d

# 或一键全流程
python run.py all
```

## 项目结构

```
meep_optimizer/
├── run.py                    # 主入口
├── AGENT_INSTRUCTION.txt     # 仿真规范
├── PROJECT_SUMMARY.md        # 历史进展报告
├── src/
│   ├── config.py             # 统一物理参数
│   ├── optimize_3d_ridge.py    # 3D 禁带优化 ⭐
│   ├── optimize_cavity_2d.py   # 2D 缺陷腔 Q 优化 ⭐
│   ├── run_3d_cavity.py        # 3D 腔体验证
│   └── scan_bandgap_2d.py      # 2D 禁带粗扫
└── results/
    ├── best_bandgap_params.json   # 镜区参数
    └── best_cavity_design.json    # 最佳腔设计（优化后生成）
```

## 输出文件

- `results_3d/optimized_params.json` — 3D 禁带最佳参数
- `results/best_cavity_design.json` — 2D 腔最佳设计（Q、λ₀、T_peak）
- `results/cavity_optimization_all.csv` — 全部腔优化记录
- `results/cavity_3d/cavity_3d_result.json` — 3D 验证 Q 值

## 物理参数

- 波导芯层 n = 2.18，衬底 n = 1.44
- 脊型波导：w = 1.5 μm，h = 0.4 μm
- 扫描：周期 a、椭圆 rx/ry、孔数 N

## 规范要点

- 必须使用 GaussianSource + 无孔参考仿真归一化
- 禁止逐孔自由优化
- 2D 结果不代表真实 3D 器件，需 3D 验证

## 手机端操作

克隆仓库后，在 WSL/服务器上：

```bash
git pull
conda activate phc-meep
nohup python run.py cavity-2d > cavity.log 2>&1 &
tail -f cavity.log
```

可用 GitHub Actions / Codespaces 远程触发（需自行配置）。
