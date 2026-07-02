---
name: fdtdx-mode-source-overlap
description: Summarize and apply the fdtdx ModeSource and ModeOverlapDetector tutorial. Use when the user mentions fdtdx ModeSource, ModeOverlapDetector, mode overlap, TE-like/TM-like modes, or when translating the workflow to Lumerical FDTD (addmode, addmodeexpansion, DFT/power monitors, transmission/S-parameters).
---

# fdtdx ModeSource / ModeOverlapDetector 速用指南（可迁移到 Lumerical）

## 适用场景（触发词）
- 用户在学/用 `fdtdx.ModePlaneSource`、`fdtdx.ModeSource`、`fdtdx.ModeOverlapDetector`
- 想理解 “TE/TM vs TE-like/TM-like”、模式注入、模式重叠、S 参数/透射测量
- 想把 fdtdx 的模式源/模式探测流程迁移到 **Lumerical FDTD**（脚本：`addmode`、`addmodeexpansion`、`adddftmonitor`/`addpower`）

## 核心概念（只记这几条）
- **模式（mode）**：波导截面的本征解；在传播中横向形状保持不变（集成光子里多为 hybrid，“TE-like/TM-like”按主电场分量约定）。
- **ModeSource**：先求模态，再把模态场分布作为激励注入；用于“像真实波导那样”激励单模/多模器件。
- **ModeOverlapDetector**：把仿真得到的场与目标模态做**重叠积分**，得到该模的耦合系数/功率；用于更干净的透射、反射、模态转换效率、S 参数。

## fdtdx 工作流（照做的顺序）
1. **搭场景**：材料 + 几何 + 边界（PML）+ 仿真体积。
2. **放源**：`ModePlaneSource(mode_index=0, filter_pol="te"/"tm", direction='+')`
   - `partial_grid_shape` 在传播方向通常取 1，其它方向覆盖模式截面（可用 `None` 自动扩展）。
3. **（关键）apply_params**：运行前必须 `fdtdx.apply_params(...)`，它会完成模式计算等预处理。
4. **先确认模式对不对**：查看源内部计算出的 `_E/_H` 分量分布（是否 TE-like/TM-like）。
5. **跑时域**：`fdtdx.run_fdtd`（必要时用能量探测器看传播是否正常）。
6. **放重叠探测器**：`ModeOverlapDetector(...)` 放在离源有距离的位置。
   - **探测器开窗**：用 `OnOffSwitch` 在后段/稳态段采样（例如从 0.75*time 开始），避免过渡态污染频域结果。
7. **计算 overlap**：`compute_overlap(...)`，通常关注 \|overlap\| 或功率比例。

## 常见坑与排查
- **overlap < 1** 常见原因：网格/数值色散、PML 吸收、源非完美匹配、探测器太近（未稳态）、仿真时间不够、模式阶次/偏振选错。
- **偏振理解**：集成波导一般用 “TE-like/TM-like”，不要拿严格 TE/TM 的定义硬套。
- **时间不够**：确保光从源传播到探测器并进入稳态，再做频域/重叠。

## 迁移到 Lumerical FDTD（对照表）
### 源（ModeSource）
- fdtdx `ModePlaneSource` ≈ Lumerical **mode source**：脚本命令 `addmode; set("injection axis", ...); set("center wavelength", ...); set("wavelength span", ...);`
- 要点：
  - 在波导截面注入（与传播方向垂直的平面）。
  - 用合适的中心波长/带宽覆盖目标波段。

### 探测（ModeOverlapDetector）
两种层级：
1. **简单透射**：`adddftmonitor`/`addpower` 读总功率，再算 \(T=P_{out}/P_{in}\)。
2. **按模式分解（更像 overlap）**：优先用 **mode expansion monitor**（脚本命令 `addmodeexpansion;`）
   - 目标：提取“基模/指定模”的前向/后向功率，得到更干净的透射、反射与模态转换（类 S 参数）。

## 在现有脚本系统里怎么落地（最短路径建议）
- 已用 `addmode` + `adddftmonitor` 的项目：
  - 若只要总透射：保持 `adddftmonitor`，用足够的频率点（例如 150）。
  - 若要“只看基模透射/忽略辐射和高阶模”：在输出端新增 `addmodeexpansion`，并在后处理优先读 mode expansion 的结果（基模前向功率 / 输入功率）。

## 结果呈现建议
- 同时输出两条曲线最稳：
  - **Total T**（功率监视器法）
  - **Fundamental-mode T**（mode expansion / overlap 法）
这样可以直观看到辐射/高阶模导致的差异。

