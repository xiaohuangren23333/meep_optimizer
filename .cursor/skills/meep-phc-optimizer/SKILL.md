---
name: meep-phc-optimizer
description: >-
  Optimize 1D photonic crystal ridge waveguides with elliptical air holes and
  quadratic-taper defect cavities in Meep (Python). Use for meep_optimizer
  tasks: mirror bandgap scan, cavity Q/T/lambda optimization, 3D FDTD validation,
  GaussianSource + reference normalization, mp.Ellipsoid holes, coordinate descent.
  Triggers on Meep, PHC, 光子晶体, 椭圆柱, 二次渐变, 缺陷腔, Q factor, bandgap, T_peak.
---

# Meep PHC Optimizer（二次渐变椭圆柱一维光子晶体）

本 Skill 指导 Agent 在本仓库用 **Meep Python** 完成：**镜区禁带优化 → 二次渐变缺陷腔优化 → 3D 验证**。与官方教程 [Resonant Modes and Transmission in a Waveguide Cavity](https://meep.readthedocs.io/en/latest/Python_Tutorials/Resonant_Modes_and_Transmission_in_a_Waveguide_Cavity/) 及 [Basics (bend-flux normalization)](https://meep.readthedocs.io/en/latest/Python_Tutorials/Basics/) 一致。

## 必读项目文件（先于改代码）

| 文件 | 内容 |
|------|------|
| `AGENT_INSTRUCTION.txt` | 禁止项、官方教程链接 |
| `SIMULATION_PROTOCOL.md` | 归一化、校验、Q 拟合 |
| `src/config.py` | 共享物理/频谱参数 |
| `src/spectrum_analysis.py` | T 归一化、Lorentzian、评分 |
| `src/run_3d_cavity.py` | `build_cavity_3d()` 几何与 3D FDTD |
| `src/optimize_3d_ridge.py` | 镜区 `build_geom()` + `run_sim()` |

## 结构物理模型

**脊型波导 + 椭圆柱空气孔**（3D 用 `mp.Ellipsoid`，2D 用 `size=(2rx, 2ry, mp.inf)`）：

```
衬底 n_sub=1.44 | 平板 h_slab=0.2 | 脊 w_wg=1.5, h_ridge=0.2 | 孔贯穿 h_total=0.4 μm
芯层 n_wg=2.18
```

**镜区（Mirror）**：N 个等间距孔，周期 `a_m`，椭圆半轴 `rx_m`, `ry_m`。

**缺陷腔（Cavity）**：中心为空缺陷；两侧各 `N_mirror` 个镜区孔 + `N_taper` 个**二次渐变**孔：

```python
t = i / N_taper  # i = 1..N_taper，从缺陷侧→镜区侧
a_i  = a_c  + (a_m  - a_c)  * t**2
rx_i = rx_c + (rx_m - rx_c) * t**2
ry_i = ry_c + (ry_m - ry_c) * t**2
```

实现见 `run_3d_cavity.build_cavity_3d()` / `optimize_cavity_2d.build_cavity_geom()`。

## 硬性约束（违反即错误）

1. **GaussianSource** + **Ey**（TM-like，E 在平面内）
2. **必须先无孔参考仿真**，再跑有孔结构；`T = flux_hole / flux_ref`
3. **有孔与参考必须相同 `cell_size`**、相同源、相同 `FluxRegion`、相同 `fcen/df/nfreq`
4. **禁止逐孔自由优化**（只允许规范中的坐标下降维度）
5. **2D Q 不能当 3D 最终结果**；报告必须区分维度
6. **不要重装 Meep**；用 `source scripts/cloud_env.sh` 或 `conda activate phc-meep`

## Meep FDTD 标准流程（每次透射谱）

```python
import meep as mp
import numpy as np

# 1) 几何 + cell
geom_h, sx, sy, sz = build_cavity_3d(mirror, cavity)  # 或 build_geom(a,rx,ry,N)
cell = mp.Vector3(sx, sy, sz)
src_x = -sx/2 + dpml + 0.5
mon_x =  sx/2 - dpml - 0.5

sources = [mp.Source(
    mp.GaussianSource(fcen, fwidth=df),
    component=mp.Ey,
    center=mp.Vector3(src_x, 0, h_total/2),
    size=mp.Vector3(0, w_wg, h_total),  # 3D；2D 用 size=(0, w_wg)
)]
freg = mp.FluxRegion(
    center=mp.Vector3(mon_x, 0, h_total/2),
    size=mp.Vector3(0, 2*w_wg, 2*h_total),
)

# 2) 有孔仿真
sim = mp.Simulation(cell_size=cell, resolution=resolution,
                    geometry=geom_h, sources=sources,
                    boundary_layers=[mp.PML(dpml)])
trans = sim.add_flux(fcen, df, nfreq, freg)
sim.run(until_after_sources=mp.stop_when_fields_decayed(
    50, mp.Ey, mp.Vector3(mon_x,0,h_total/2), decay))
flux_h = np.array(mp.get_fluxes(trans))
freqs  = np.array(mp.get_flux_freqs(trans))

# 3) 参考仿真（相同 cell，无孔波导）
geom_ref = build_ref_geom_3d()  # 仅 sub+slab+ridge
sim_r = mp.Simulation(cell_size=cell, ...)  # 同上
trans_r = sim_r.add_flux(fcen, df, nfreq, freg)
sim_r.run(...)
flux_r = np.array(mp.get_fluxes(trans_r))

# 4) 归一化（用 spectrum_analysis）
from spectrum_analysis import normalize_transmission, validate_normalized_spectrum
wl, T, fh, fr = normalize_transmission(freqs, flux_h, freqs, flux_r)
val = validate_normalized_spectrum(wl, T, fh, fr)
```

参考：`optimize_3d_ridge.run_sim()`（镜区）、`run_3d_cavity.run_cavity_fdtd()`（腔体）。

## 优化流水线

| 阶段 | 命令 | 脚本 | 输出 |
|------|------|------|------|
| Phase 1 镜区 | `python run.py mirror-v2` | `optimize_mirror_3d_v2.py` | `results_3d/best_mirror_phase1.json` |
| Phase 2 腔体 | `python run.py cavity-3d-v2` | `optimize_cavity_3d_v2.py` | `results/best_cavity_design.json` |
| 一键 | `python run.py phase1-2` | 顺序执行 | 上两者 |
| 3D 验证 | `python run.py cavity-3d` | `run_3d_cavity.py` | `results/cavity_3d/` |
| 2D 腔（快筛） | `python run.py cavity-2d` | `optimize_cavity_2d.py` | 仅筛选，须 3D 确认 |

**Phase 2 六步坐标下降**（不可打乱顺序）：

1. `a_c`（Nt=3, Nm=base）
2. `rx_c`
3. `ry_c`
4. 局部微调 a_c, rx_c, ry_c
5. `N_taper`
6. `N_mirror`

## 目标与评分

| 指标 | 目标 | 说明 |
|------|------|------|
| Q | ≥ 1000 | Lorentzian 拟合，γ>0；禁带内取峰 |
| λ₀ | 1500–1600 nm 可接受 | 1550 nm 权重降低 |
| T_peak | 高（~1 理想） | 与 Q 存在耦合权衡 |
| T_min（镜区） | 低 | Phase 1 宽深禁带 |

- 腔体评分：`spectrum_analysis.score_cavity_peak()`
- FDTD 原始谱 Q 交叉验证：读 `.cursor/skills/optical-resonance-q/SKILL.md`，方法 **`peak_bracketed`**

## 频谱与分辨率默认值

来自 `src/config.py`：

- λ: 1.25–1.75 μm，`fcen=1/1.55`，`df=(1/λ_min-1/λ_max)/2`
- `resolution_3d=20`, `resolution_2d=30`, `dpml=1.0`, `pad=2.0`
- `nfreq_bandgap_3d=800`, `nfreq_cavity_3d=2000`

## 常见错误与修复

| 症状 | 原因 | 修复 |
|------|------|------|
| T>1 或 T 全异常 | 参考 cell 与有孔 cell 不一致 | 对齐 `sx,sy,sz`（已修复于 `run_3d_cavity.py`） |
| Q 假高/假低 | FDTD 宽禁带用错误 FWHM | 用 `peak_bracketed` 或约束 Lorentzian |
| 2D Q≫3D Q | z 向泄漏未建模 | 以 3D 为准 |
| 无腔模峰 | 禁带太窄/未对准 λ | 先重跑 Phase 1 扩禁带 |
| 仿真不收敛 | decay 太严或 cell 太小 | 检查 `stop_when_fields_decayed`；增大 cell |

## 2D vs 3D 几何差异

| | 2D | 3D |
|--|----|----|
| 波导 | `Block(..., size=(inf, w_wg, 0))` | sub + slab + ridge 三层 |
| 孔 | `Ellipsoid(..., size=(2rx,2ry, mp.inf))` | `Ellipsoid(..., size=(2rx,2ry, h_total))` |
| 源 z | 可省略 | `center z = h_total/2` |
| 可信度 | 快速扫描 | **最终指标** |

## Agent 工作 checklist

- [ ] 读 `best_mirror_phase1.json` / `best_cavity_design.json` 再改参
- [ ] 新脚本复用 `build_cavity_3d` / `run_cavity_fdtd` / `spectrum_analysis`
- [ ] 保存 `flux_h`, `flux_ref`, `T(λ)` 与 `validation.reliable`
- [ ] 长任务用 `nohup` 或 `scripts/run_cloud_pipeline.sh`
- [ ] 报告注明 2D 或 3D、归一化是否同 cell

## 延伸阅读

- 详细公式与 API：本目录 `reference.md`
- 几何示意与参数表：本目录 `geometry.md`
- Meep 示例：`https://github.com/NanoComp/meep/tree/master/python/examples`
