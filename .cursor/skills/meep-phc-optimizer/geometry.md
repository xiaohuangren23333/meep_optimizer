# Quadratic-taper cavity geometry

## Side view (x propagation, symmetric)

```
        Mirror          Taper (t²)         Defect         Taper (t²)          Mirror
   ←—— N_mirror ——→  ←—— N_taper ——→   (no hole)   ←—— N_taper ——→   ←—— N_mirror ——→
   (a_m, rx_m, ry_m)  a,rx,ry interpolate   gap    mirror params      uniform holes
```

## Parameter sets

| Symbol | Role | Optimized in |
|--------|------|--------------|
| `a_m, rx_m, ry_m` | Mirror unit cell | Phase 1 |
| `a_c, rx_c, ry_c` | Defect / taper endpoint | Phase 2 steps 1–4 |
| `N_taper` | Taper hole count each side | Phase 2 step 5 |
| `N_mirror` | Mirror holes each side | Phase 2 step 6 |

Total holes each side: `N_mirror + N_taper` (defect center has no hole).

## Taper table example

`N_taper=3`, mirror `(a_m,rx_m,ry_m)=(0.44,0.14,0.28)`, defect `(a_c,rx_c,ry_c)=(0.42,0.12,0.24)`:

| i | t | a_i | rx_i | ry_i |
|---|---|-----|------|------|
| 1 | 1/3 | 0.4278 | 0.1267 | 0.2533 |
| 2 | 2/3 | 0.4311 | 0.1333 | 0.2667 |
| 3 | 1   | 0.4400 | 0.1400 | 0.2800 |

Formula: `param_i = param_c + (param_m - param_c) * (i/N_taper)²`

## Cell size

3D (`run_3d_cavity.py`):

```
total_holes = 2 * (N_mirror + N_taper)
sx = 2*dpml + max(total_holes * a_m, 4.0) + 2*pad
sy = 2*dpml + w_wg + 2*pad
sz = 2*dpml + h_total + 1.5
```

Reference simulation uses **same sx, sy, sz** with geometry = sub + slab + ridge only.

## Copy-paste: mirror bandgap geometry (3D)

See `optimize_3d_ridge.build_geom(a, rx, ry, N)` — uniform elliptical holes along x.

## Copy-paste: full cavity

See `run_3d_cavity.build_cavity_3d(mirror, cavity)` — do not duplicate geometry logic in new files; import and call.
