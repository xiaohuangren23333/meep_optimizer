"""Shared physical and simulation parameters for PHC optimization."""

# Materials
n_wg = 2.18
n_air = 1.0
n_sub = 1.44

# Fabrication minimum feature size
min_feature_nm = 200
min_feature_um = min_feature_nm / 1000.0  # 0.2 μm

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


def check_hole_geometry(a, rx, ry, w_wg=w_wg, min_um=min_feature_um):
    """Validate a single elliptical hole against the minimum feature size.

    The 200 nm minimum feature applies to the smallest fabricable dimension:
      - hole DIAMETER (2*rx, 2*ry) — a hole of diameter 200 nm is printable
      - solid wall between adjacent holes (a - 2*rx)
      - remaining ridge width beside the hole (w_wg - 2*ry)
    (Previously this incorrectly required the semi-axis >= 200 nm, i.e. a
    400 nm diameter, which forced the period into the leaky 2nd-order gap.)
    """
    tol = 1e-6
    violations = []
    if 2 * rx + tol < min_um:
        violations.append(f"hole_dx={2*rx*1000:.0f}nm<{min_feature_nm}nm")
    if 2 * ry + tol < min_um:
        violations.append(f"hole_dy={2*ry*1000:.0f}nm<{min_feature_nm}nm")
    wall_x = a - 2 * rx
    if wall_x + tol < min_um:
        violations.append(f"wall_x={wall_x*1000:.0f}nm<{min_feature_nm}nm")
    wall_y = w_wg - 2 * ry
    if wall_y + tol < min_um:
        violations.append(f"wall_y={wall_y*1000:.0f}nm<{min_feature_nm}nm")
    return len(violations) == 0, violations


def check_periodic_geometry(a, rx, ry, w_wg=w_wg):
    """Validate periodic hole array geometry."""
    return check_hole_geometry(a, rx, ry, w_wg=w_wg)


def check_cavity_geometry(a_m, rx_m, ry_m, a_c, rx_c, ry_c, N_taper, w_wg=w_wg):
    """Validate mirror, defect, and taper holes in a cavity design."""
    for label, a, rx, ry in (
        ("mirror", a_m, rx_m, ry_m),
        ("defect", a_c, rx_c, ry_c),
    ):
        ok, violations = check_hole_geometry(a, rx, ry, w_wg=w_wg)
        if not ok:
            return False, [f"{label}: {v}" for v in violations]

    for i in range(1, N_taper + 1):
        t = i / N_taper
        ai = a_c + (a_m - a_c) * t**2
        rxi = rx_c + (rx_m - rx_c) * t**2
        ryi = ry_c + (ry_m - ry_c) * t**2
        ok, violations = check_hole_geometry(ai, rxi, ryi, w_wg=w_wg)
        if not ok:
            return False, [f"taper_{i}: {v}" for v in violations]
    return True, []


def filter_valid_scan_values(values, validator):
    """Keep only scan values that satisfy a geometry validator."""
    valid = []
    for val in values:
        ok, _ = validator(val)
        if ok:
            valid.append(val)
    return valid


def max_hole_rx(a, min_um=min_feature_um):
    """Maximum elliptical x semi-axis for periodic lattice period a (wall >= min)."""
    return (a - min_um) / 2.0


def max_hole_ry(w_wg=w_wg, min_um=min_feature_um):
    """Maximum elliptical y semi-axis for ridge width w_wg (remaining wall >= min)."""
    return (w_wg - min_um) / 2.0


def min_hole_semi_axis(min_um=min_feature_um):
    """Minimum semi-axis so hole diameter >= minimum feature size."""
    return min_um / 2.0
