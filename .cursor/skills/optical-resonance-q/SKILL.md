---
name: optical-resonance-q
description: >-
  Computes optical resonance Q factor (Q = lambda_peak / FWHM) from transmission
  spectra with method selection for FDTD raw spectra vs TMM/PPT peaks. Use when
  the user mentions Q factor, Q值, 品质因子, FWHM, 半高宽, resonance peak
  extraction, or validating FDTD/TMM spectra against dataset Q_peak_fwhm.
---

# 光学谐振 Q 因子计算

## 公式与单位

- **Q = λ_peak / FWHM**（无量纲）
- **λ_peak、FWHM** 必须与波长轴同单位；输出字段统一为 **`lambda_*_um`、`fwhm_um`（µm）**
- 波长轴为 **米** 时，调用 `transmission_peak` 或 `compute_q(..., wavelength_unit="m")`

## 方法选择（必做）

| 谱类型 | 方法 | 原因 |
|--------|------|------|
| **FDTD 原始透射**（阻带宽、多交点） | `peak_bracketed` | 峰肩台基线 + 只取**包住峰**的一对半高交点，避免禁带远处假 FWHM |
| **TMM / PPT 无损尖峰** | `transmission_peak` | 峰附近局部谷底 + 0.5(T_peak−T_base)，与 `tmm_lcavity_opt.extract_peak` 一致 |
| **已知 dip1/dip2 且禁带窄** | `dip_interval` | 谷底区间内全交点 min/max；**宽禁带勿用** |

**禁止**：对 FDTD 原始谱用 `0.5*T_peak` 且在全谱找交点 → FWHM 可达数十 nm，Q 假高/假低。

## 快速使用

### 1. 运行自带脚本（推荐）

脚本路径（全局 skill）：

`.cursor/skills/optical-resonance-q/scripts/q_peak.py`

```bash
# FDTD：npz 中 wavelength_um + spectrum
python .cursor/skills/optical-resonance-q/scripts/q_peak.py \
  --npz path/to/n1.550.npz --method peak_bracketed

# TMM：波长米 + 透射率
python .cursor/skills/optical-resonance-q/scripts/q_peak.py \
  --method transmission_peak --wl ... --T ...
```

### 2. 在 Python 中 import

将 `scripts` 目录加入 `sys.path` 后：

```python
from q_peak import compute_q, peak_fwhm_bracketed, transmission_peak_fwhm

# FDTD 原始谱（wl 单位 µm）
det = peak_fwhm_bracketed(wl_um, T_raw, peak_lam_um=1.5504)
# 或
det = compute_q(wl_um, T_raw, method="peak_bracketed", wavelength_unit="um")

# TMM（wl 单位 m）
det = compute_q(wl_m, T_tmm, method="transmission_peak", wavelength_unit="m", lam_guess=1550e-9)
```

成功时读取：`det["Q_peak"]`, `det["fwhm_um"]`, `det["lambda_peak_um"]`, `det.get("T_peak")`。

## 工作流

1. **确认波长单位**（µm 或 m）与谱类型（原始 / 尖峰 / 平滑）。
2. **选方法**（见上表）；不确定 FDTD → `peak_bracketed`。
3. **峰位**：未给 `peak_lam_um` 时，在 `lam_center±50 nm`（或 `±50e-9 m`）内取 `argmax(T)`。
4. **计算 Q**；检查 `ok`、`fwhm_oversampled`（采样过粗时 FWHM < 2.5×dλ）。
5. **对照**：数据集 `Q_peak_fwhm`、同 KCG 的 FDTD 报告；TMM 无损谱 Q 常低于 FDTD（T_peak≈1）。

## 输出字段（标准）

| 字段 | 含义 |
|------|------|
| `ok` | 是否成功 |
| `method` | 使用的方法名 |
| `Q_peak` | 品质因子 |
| `fwhm_um` | 半高全宽 (µm) |
| `lambda_peak_um` | 峰位 (µm) |
| `lambda_L_um` / `lambda_R_um` | 半高左右交点 |
| `T_peak` / `T_base` | 峰高 / 基线（透射法） |
| `reason` | 失败原因（`ok=False`） |

## 项目内映射（Example_10）

| 位置 | 实现 |
|------|------|
| FDTD 灵敏度 / 原始谱 | `11-GAN图形生成/tools/fdtd_sensitivity_analyze_raw.py` → `peak_fwhm_bracketed` |
| TMM Lcavity 优化 | `_tmm_bragg_review/tmm_lcavity_opt.py` → `extract_peak`（= `transmission_peak`） |
| 数据集批量 | `10-光谱预测/tools/compute_q_fwhm_spikeonly.py` → `dip_interval`（需 dip1/dip2） |

新代码**优先**调用本 skill 的 `q_peak.py`，避免复制三份逻辑。

## 失败与排查

| 现象 | 处理 |
|------|------|
| `no_bracketing_crossings` | 收窄/平谱；检查是否阻带内需先找峰；略增 `half_win_um` |
| `fwhm` 极大、Q 极小 | 可能用了 `dip_interval` 或全谱交点 → 改 `peak_bracketed` |
| Q 与 npz 差 >10% | 对照是否同一谱（raw vs clean）、同一峰位、采样点数 |
| TMM Q~500 vs FDTD Q~900 | 模型差异，非 Q 公式错误；见 TMM/FDTD 对比说明 |

## 报告模板

```text
λ_peak = xxx nm, FWHM = xxx nm, Q = xxx
方法: peak_bracketed | transmission_peak
T_peak = xxx, T_base = xxx（如有）
```

## 附加说明

- 高分辨率 TMM 峰应用 **≥20001 点/±2 nm** 再算 Q；稀疏网格（如 200 点）仅适合看峰位。
- 需要批量处理 npz 时，可写循环调用 `compute_q`，不要对整表重算 spike-clean 后再标 Q 除非用户明确要求。
