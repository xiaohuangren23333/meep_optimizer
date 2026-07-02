# Meep PHC Reference

## Official tutorials (read order)

1. [Resonant Modes and Transmission in a Waveguide Cavity](https://meep.readthedocs.io/en/latest/Python_Tutorials/Resonant_Modes_and_Transmission_in_a_Waveguide_Cavity/) — periodic holes, cavity mode, transmission spectrum
2. [Basics — transmittance spectrum](https://meep.readthedocs.io/en/latest/Python_Tutorials/Basics/) — `GaussianSource`, `add_flux`, normalization
3. [Python User Interface](https://meep.readthedocs.io/en/latest/Python_User_Interface/) — `Simulation`, `Source`, `FluxRegion`, `Ellipsoid`, `PML`
4. [bend-flux.py](https://github.com/NanoComp/meep/blob/master/python/examples/bend-flux.py) — save/load flux for R/T normalization pattern

## Normalization (Meep docs pattern)

Official holey-waveguide cavity tutorial: run **structure with holes**, then run **same setup without holes** (or empty mirror region), divide fluxes:

```
T(ω) = flux_structure(ω) / flux_reference(ω)
```

Critical details from Meep bend-flux example:

- Use **identical** `cell_size`, sources, `FluxRegion`, `fcen`, `df`, `nfreq`
- Reference geometry = waveguide only (no air holes)
- For reflection setups use `get_flux_data` + `load_minus_flux_data`; **this project uses transmission-only** (single flux plane at output)

## Elliptical hole in Meep

`mp.Ellipsoid` `size=(2*rx, 2*ry, depth)` defines axis-aligned ellipse with semi-axes rx, ry.

Mirror hole at period index i (centered lattice):

```python
hole_start = -(N-1)*a / 2.0
x = hole_start + i * a
mp.Ellipsoid(material=mp.Medium(index=n_air),
             center=mp.Vector3(x, 0, z_center),
             size=mp.Vector3(2*rx, 2*ry, hole_depth))
```

## Quadratic taper (二次渐变)

For taper index `i ∈ {1,…,N_taper}` from defect toward mirror:

```
t = i / N_taper
a(t)  = a_c  + (a_m  - a_c)  * t²
rx(t) = rx_c + (rx_m - rx_c) * t²
ry(t) = ry_c + (ry_m - ry_c) * t²
```

Hole center spacing uses local period `a(t)`; positions in `run_3d_cavity.py`:

- Left taper: `x = -(N_taper - i + 0.5) * ai`
- Right taper: `x = +(N_taper - i + 0.5) * ai`

**Do not** replace with linear taper unless explicitly requested — this repo standard is **quadratic**.

## Simulation parameters (project default)

```python
n_wg, n_sub = 2.18, 1.44
w_wg = 1.5          # μm
h_slab, h_ridge = 0.2, 0.2
h_total = 0.4
lambda_min, lambda_max = 1.25, 1.75   # μm
fcen = 1.0 / 1.55
df = (1/lambda_min - 1/lambda_max) / 2
resolution_3d = 20
dpml = 1.0
pad = 2.0
```

Source (3D ridge):

```python
mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ey,
          center=mp.Vector3(src_x, 0, h_total/2),
          size=mp.Vector3(0, w_wg, h_total))
```

Flux monitor (output, before PML):

```python
mon_x = sx/2 - dpml - 0.5
mp.FluxRegion(center=mp.Vector3(mon_x, 0, h_total/2),
              size=mp.Vector3(0, 2*w_wg, 2*h_total))
```

Field decay (cavity): `stop_when_fields_decayed(50, mp.Ey, pt, 1e-5)`  
Bandgap scan: `1e-4`

## Q extraction (two methods in repo)

### A. Lorentzian fit (`spectrum_analysis.fit_peak`)

- Model: `offset + A*γ²/((λ-λ₀)²+γ²)`
- Q = λ₀ / (2γ), constrained γ > 0
- Peak must fall inside mirror bandgap JSON range

### B. FWHM bracket (`optical-resonance-q`)

- For raw FDTD spectra with wide stopband shoulders
- Script: `.cursor/skills/optical-resonance-q/scripts/q_peak.py --method peak_bracketed`

Use both when validating suspicious Q (2D vs 3D gap).

## Repository entry points

```bash
source scripts/cloud_env.sh   # Cloud / after bootstrap
python run.py mirror-v2       # Phase 1
python run.py cavity-3d-v2    # Phase 2
python run.py phase1-2
python run.py cavity-3d       # single design validation
```

## Output artifacts

| Path | Content |
|------|---------|
| `results_3d/best_mirror_phase1.json` | a, rx, ry, N, gap, T_min |
| `results/best_cavity_design.json` | a_c, rx_c, ry_c, Nt, Nm, Q, λ₀, T_peak |
| `results/cavity_3d/cavity_3d_flux.csv` | wavelength, flux_h, flux_ref, T |
| `results/cavity_3d/cavity_3d_analysis.json` | validation + fit metadata |

## Harminv (future Phase 3)

For mode Q independent of Lorentzian fit, Meep provides `Harminv` after long run with narrowband excitation near resonance. Not yet wired in this repo — prefer extending `run_3d_cavity.py` rather than new ad-hoc scripts.
