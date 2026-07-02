# 仿真数据可靠性说明

本文档说明本项目 FDTD 数据的生成方式，便于你在 Meep 或其他 FDTD 中独立验证。

## 归一化（必须满足）

```
T(λ) = flux_hole(λ) / flux_ref(λ)
```

| 要求 | 实现位置 |
|------|----------|
| 有孔与参考 **相同 cell 尺寸** | `optimize_3d_ridge.py`, `run_3d_cavity.py`, `optimize_cavity_2d.py` |
| 相同光源 (GaussianSource, Ey) | 所有脚本 |
| 相同 flux 监视器位置与尺寸 | 所有脚本 |
| 相同 fcen / df / nfreq | 见各脚本或 `config.py` |
| 参考 = 无孔波导 | 3D: `build_ref_geom_3d()`; 2D: 两 Block 波导 |

**已修复的问题**：旧版 `run_3d_cavity.py` 参考仿真用了更短的 cell（4 μm），导致归一化偏差。现已改为与腔体相同 cell。

## 输出文件（可复现验证）

每次 3D 腔验证生成：

| 文件 | 内容 |
|------|------|
| `results/cavity_3d/cavity_3d_flux.csv` | wavelength, flux_hole, flux_ref, T |
| `results/cavity_3d/cavity_3d_T.csv` | wavelength, transmission |
| `results/cavity_3d/cavity_3d_analysis.json` | 校验结果、峰拟合、仿真元数据 |
| `results/cavity_3d/cavity_3d_result.json` | 汇总 |

`cavity_3d_analysis.json` 中 `validation.reliable == true` 表示通过自动校验。

## 自动校验项

- T ≥ -0.02（无严重负透射）
- T_max ≤ 1.05（轻微超出仅警告）
- 1450–1650 nm 内有可见峰
- 有孔/参考频率网格一致

## Q 值提取

- Lorentzian 拟合，**gamma > 0 有界约束**
- Q = λ₀ / FWHM，FWHM = 2γ
- 只接受 **镜区禁带内** 的峰（见 `best_bandgap_params.json`）
- R² ≥ 0.90（3D）/ 0.95（2D 优化）

## 3D 结构参数

```
n_wg=2.18, n_sub=1.44, w_wg=1.5 μm
h_slab=0.2, h_ridge=0.2, h_total=0.4 μm
椭圆孔贯穿脊型层 (z 方向 h_total)
resolution=20, dpml=1.0
field decay: 1e-5 (3D cavity), 1e-4 (bandgap scan)
```

## 独立验证建议

1. 用 `best_cavity_design.json` 中的几何参数重建结构
2. 复现相同 cell、源、监视器、PML
3. 分别跑有孔/无孔，计算 T = flux_h/ref
4. 对比 `cavity_3d_flux.csv` 中的 flux 与 T 曲线

2D 与 3D Q 值会有差异（2D 忽略 z 泄漏），以 **3D 结果** 为准。
