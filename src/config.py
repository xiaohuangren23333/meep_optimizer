"""Shared physical and simulation parameters for PHC optimization."""

# Materials
n_wg = 2.18
n_air = 1.0
n_sub = 1.44

# 3D ridge waveguide geometry (μm)
w_wg = 1.5
h_slab = 0.2
h_ridge = 0.2
h_total = h_slab + h_ridge

# Spectral range
lambda_min = 1.25
lambda_max = 1.75
lambda_target = 1.55
fcen = 1.0 / lambda_target
df = (1.0 / lambda_min - 1.0 / lambda_max) / 2.0

# FDTD settings
resolution_2d = 30
resolution_3d = 20
dpml = 1.0
pad = 2.0

# Frequency samples
nfreq_bandgap_3d = 800
nfreq_cavity_2d = 2000
nfreq_cavity_3d = 2000

# Optimization targets
Q_target = 1000
T_peak_target = 0.9
T_min_target = 0.1
