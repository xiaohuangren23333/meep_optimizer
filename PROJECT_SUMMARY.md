# MEEP Photonic Crystal Optimizer — 项目总结报告

## 项目概述

使用 Meep + Python 优化一维光子晶体椭圆孔缺陷腔，通过 FDTD 仿真寻找最优禁带参数和缺陷腔设计。

## 环境

- **WSL Ubuntu** + conda 环境 `phc-meep`
- **项目路径**: `~/projects/meep_optimizer`
- **Meep 版本**: 已安装，Python 接口可用

---

## 文件结构

```
meep_optimizer/
├── AGENT_INSTRUCTION.txt              # 项目规范和约束
├── PROJECT_SUMMARY.md                 # 本报告
│
├── src/                               # 源代码
│   ├── scan_bandgap_2d.py             # 2D 禁带大规模扫描 (Phase 0/1/2)
│   ├── optimize_3d_ridge.py           # 3D 脊型波导坐标下降法优化 ⭐最佳
│   ├── run_3d_periodic_scan.py        # 3D 周期孔扫描 (v8, EigenModeSource)
│   ├── run_3d_scan_overnight.py       # 3D 过夜扫描 (备用)
│   ├── run_3d_validate_candidates.py  # 3D 验证 2D 候选 (当前主用 ❌问题)
│   ├── run_3d_correct_ridge.py        # 3D 正确脊型波导 (旧版)
│   ├── run_3d_cavity.py               # 3D 缺陷腔仿真
│   ├── launch_cavity_after_3d.py      # 3D 缺陷腔启动脚本
│   ├── optimize_cavity_2d.py          # 2D 缺陷腔6步优化器
│   └── run_3d_after_2d_cavity.sh      # 2D→3D 腔 shell 脚本
│
├── results/                           # 2D 结果
│   ├── bandgap_2d_scan.csv            # 2D Phase 1 扫描结果
│   ├── bandgap_candidates_2d.json     # 2D 候选 (9组, Phase 0)
│   ├── bandgap_candidates_2d_test5.json  # 2D 候选 (测试)
│   ├── bandgap_candidates_2d_top30.json  # 2D Top 30 候选
│   ├── cavity_3d/                     # 3D 缺陷腔结果
│   ├── spectra/                       # 2D 频谱 CSV (9 组验证)
│   └── figures/                       # 2D 频谱 PNG (9 组验证)
│
├── results_3d/                        # 3D 结果
│   ├── periodic_3d_scan_results.csv   # 3D 周期扫描结果
│   ├── bandgap_candidates_3d.json     # 3D 候选 JSON
│   ├── validation_3d_results.csv      # 3D 验证 30 组结果 (EigenModeSource)
│   ├── spectra/                       # 3D 频谱 CSV
│   ├── figures/                       # 3D 频谱 PNG
│   ├── fields/                        # 3D 场分布
│   └── models/                        # 3D 模型
│
└── *.log                              # 运行日志 (~20个)
```

---

## 物理参数

### 材料
- **波导芯层**: n_wg = 2.18 (LiTaO₃ 等效)
- **包层/孔**: n_air = 1.0
- **衬底**: n_sub = 1.44 (SiO₂ 等效)

### 3D 结构 (脊型波导)
- 波导宽度: w_wg = 1.5 μm
- 平板层: h_slab = 0.2 μm
- 脊型层: h_ridge = 0.2 μm
- 总高度: h_total = 0.4 μm
- 椭圆柱孔: 贯穿孔 (全刻蚀 0.4 μm)

### 2D 等效结构
- 波导: n=2.18, 宽度 w_wg=1.5 μm
- 椭圆柱孔: n=1.0, 2D 截面 (x-y 平面)

### 仿真参数
- 频率范围: λ = 1250~1750 nm


### 待扫描参数
| 参数 | 含义 | 范围 |
|------|------|------|
| a | 周期 (μm) | 0.36 ~ 0.52 |
| rx | 椭圆 x 半轴 (μm) | 0.08 ~ 0.26 |
| ry | 椭圆 y 半轴 (μm) | 0.15 ~ 0.47 |
| N | 孔数 | 12 ~ 24 |

---

## 已完成工作

### 1. 2D 禁带大规模扫描 (`scan_bandgap_2d.py`)
- **Phase 0**: 9 组快速验证 ✅
  - a=0.44, rx∈[0.10,0.14,0.18], ry∈[0.20,0.30,0.40], N=16
  - 全部有禁带 (96~160 nm)
  - T_min_gap 最低 0.0017
- **Phase 1**: 1728 组粗扫 ✅
  - a∈[0.36,0.52], rx∈[0.08,0.22], ry∈[0.15,0.43], N∈[12,16,20]
  - 结果保存在 `results/bandgap_2d_scan.csv`
- **Phase 2**: 精细扫描 (未完整运行)

**已知问题**: 2D 中 T_max > 1 (EigenModeSource 归一化偏差)，但 T_min_gap 仍可靠。

### 2. 3D 脊型波导坐标下降法 (`optimize_3d_ridge.py`) ⭐最佳结果
- **方法**: 从 a=0.44, rx=0.12, ry=0.22, N=20 出发，坐标下降法搜索
- **光源**: `mp.Source(GaussianSource, Ey)`，源尺寸 (0, w_wg, h_total)
- **监视器**: FluxRegion (0, 2*w_wg, 2*h_total)
- **参考仿真**: 使用相同 cell，N=0 的同一个 `build_geom` 函数
- **nfreq**: 800
- **结果**: 
  - T_min ≈ 0.064, gap ≈ 161 nm
  - **T_max ≤ 1.0** ✅
  - 禁带中心对准 1550 nm

### 3. 3D 验证 2D Top 30 (`run_3d_validate_candidates.py`)
- **30 组全部跑完** (EigenModeSource 版本)
- **结果**: 
  - 全部有禁带 (33~161 nm)
  - T_min 范围: 0.266 ~ 0.876
  - 最佳: a=0.48, rx=0.22, ry=0.43, N=20 → T_min=0.266, T_max=0.694
  - **问题**: T_min 远不如旧优化脚本 (0.064)

### 4. 3D 周期扫描 (`run_3d_periodic_scan.py` / `run_3d_scan_overnight.py`)
- v8 版本使用 EigenModeSource，验证 T≤1 ✅
- 但结果中 T_min 普遍偏高

### 5. 2D 缺陷腔优化 (`optimize_cavity_2d.py`)
- **功能**: 6步自动优化缺陷腔参数
- **方法**: 二次渐变 + Lorentzian 拟合
- 700 行完整优化器，包含数据保存、绘图、评分

---

## 关键问题分析

### 问题 1: 脚本碎片化严重
- **7 个 src/*.py 文件**功能重叠（多个 3D 扫描脚本）
- 参数定义分散在各脚本中，修改需同步多处
- 无统一配置文件或 argparse

### 问题 2: optimize_3d_ridge.py vs run_3d_validate_candidates.py 差异

| 参数 | optimize_3d_ridge.py (最佳) | run_3d_validate_candidates.py (问题) |
|------|---------------------------|-----------------------------------|
| 光源 | GaussianSource(Ey) | GaussianSource(Ey) *[已修复]* |
| 源尺寸 | (0, w_wg, h_total) | (0, w_wg, h_total) *[已修复]* |
| 监视器 | (0, 2*w_wg, 2*h_total) | (0, 2*w_wg, 2*h_total) *[已修复]* |
| nfreq | 800 | 1500 |
| 参考仿真 | build_geom(N=0) 同函数 | build_ref_geom() 独立函数 |
| cell 高度 (sz) | 2*dpml + h_total + 1.5 | 2*dpml + h_total + 2.0 |
| 结果 T_min | 0.064 | 0.266 (第一批) |

**差异分析**:
1. **nfreq**: 800 vs 1500 — 更多频率点应给出更精确结果，不导致变差
2. **参考仿真构建**: `build_geom(N=0)` vs `build_ref_geom()` — 关键差异！
   - `build_geom(N=0)` 返回 `[sub, slab, ridge]`，cell 未使用但几何相同
   - `build_ref_geom(cell)` 也返回 `[sub, slab, ridge]`，但 cell 传入后未使用
   - **实际上两者应等价**，cell 在 run_sim 中根据 geom 和 cell 参数重新构建
3. **cell sz**: +1.5 vs +2.0 — 1.5 更合理，减少仿真体积
4. **参数集不同**: 
   - optimize_3d_ridge.py 扫描 a=0.43~0.45, rx=0.10~0.14, ry=0.20~0.24 (小范围精细)
   - validate_candidates 扫描 a=0.40~0.48, rx=0.08~0.22, ry=0.27~0.43 (2D 筛选的大范围)
   - **参数空间不同可能是 T_min 差异的主因** — 2D 最优参数 ≠ 3D 最优参数

### 问题 3: 2D → 3D 迁移效果
- 2D 预测 T_min_gap ~ 0.001 (极深禁带)
- 3D 实际 T_min ~ 0.27 (浅很多)
- **原因**: 2D 忽略 z 方向模式泄漏，禁带深度被大幅高估
- **结论**: 2D 粗筛有效但需要修正预期，3D 验证不可或缺

### 问题 4: 运行日志堆积
- ~20 个 `run_3d_*.log` 文件，名称混乱
- 难以追踪哪个对应哪个版本

---

## 成果亮点

### 1. 已确认 3D 禁带存在
- 30/30 组都有禁带 (EigenModeSource 版)
- 禁带宽度 33~161 nm


### 2. 最佳参数候选项
- **3D 坐标下降法**: a=0.44, rx=0.12, ry=0.22, N=20 → T_min=0.064, T≤1
- **2D 扫描**: a=0.44, rx=0.14, ry=0.40, N=16 → T_min_gap=0.0017 (2D 高估)

### 3. 完整工具链
- 2D 扫描 → 候选筛选 → 3D 验证 → 缺陷腔设计
- 所有脚本支持断点续传 (检查已有 CSV 跳过)
- 自动生成频谱图和拟合图

---

## 待解决任务

| 优先级 | 任务 | 说明 |
|--------|------|------|
| 🔴 P0 | 统一验证方法 | 确认 optimize_3d_ridge.py 的参数在当前框架是否能复现 T_min=0.064 |
| 🔴 P0 | 定位根因 | 确定 run_3d_validate_candidates.py 与 optimize_3d_ridge.py 结果差异的真正原因 |
| 🟡 P1 | 清理项目 | 删除废弃脚本和日志，合并重复功能 |
| 🟡 P1 | 参数集中化 | 创建统一配置文件，避免参数分散 |
| 🟢 P2 | 重新 3D 扫描 | 用正确方法扫描有前景的参数区域 (a=0.44~0.48) |
| 🟢 P2 | 缺陷腔 3D | 用最佳禁带参数构建 3D 缺陷腔并仿真 |

---
