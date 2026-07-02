# reference.md — 实现细节

## peak_bracketed 参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `half_win_um` | 0.006 | 半高搜索半窗 (µm)，约 ±6 nm |
| `shoulder_um` | 0.003 | 肩台区离峰中心偏移 (µm) |
| `search_half_window_um` | 0.05 | 自动找峰窗口 (µm) |

## transmission_peak 参数

| 参数 | 默认 (m) | 说明 |
|------|----------|------|
| `peak_search_half_window` | 5e-9 | 找峰窗口半宽 |
| `local_half_window` | 20e-9 | 基线 min(T) 窗口半宽 |

## 与 Khafagy 2025 FOM

灵敏度论文常用 **S_λ = Δλ/ΔC**，与 Q 不同。Q 用于谐振腔尖锐度；勿把 Q 与 S_λ/FOM 混在同一公式里。
