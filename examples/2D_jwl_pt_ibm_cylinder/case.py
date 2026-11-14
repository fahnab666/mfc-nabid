#!/usr/bin/env python3
"""2D JWL-products slab driving a shock over an immersed cylinder with JWL/ambient PT relaxation.

Environment: MFC_IBM_MODE = static | prescribed | free   (moving_ibm 0 | 1 | 2)
Synthetic coefficients as in examples/1D_jwl_pt_relaxation; a solver verification case, not a material calibration.
"""

import json
import math
import os

mode = os.environ.get("MFC_IBM_MODE", "static").strip().lower()
if mode not in ("static", "prescribed", "free"):
    raise ValueError("MFC_IBM_MODE must be 'static', 'prescribed' or 'free'")
moving = {"static": 0, "prescribed": 1, "free": 2}[mode]

A, B, R1, R2, omega, rho0 = 2e10, 1e8, 4.0, 1.1, 0.3, 1000.0
cv_j, cv_a, gamma_a, pi_a, t0 = 1000.0, 750.0, 0.4, 0.0, 300.0
eps = 1e-6


def jwl_pressure(rho_j, temp):
    pref = A * math.exp(-R1 * rho0 / rho_j) + B * math.exp(-R2 * rho0 / rho_j)
    tref = t0 * (rho_j / rho0) ** omega
    return pref + omega * rho_j * cv_j * (temp - tref)


def jwl_density_at(pres, temp):
    lo, hi = 1e-3, 1e4
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if jwl_pressure(mid, temp) < pres:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


# Products slab (rho_j = 800, T = 1800 K) against ambient air (p = 1e5 Pa, rho = 1.2).
rj_hi, t_hi = 800.0, 1800.0
p_hi = jwl_pressure(rj_hi, t_hi)
ra_hi = (p_hi + pi_a) / (gamma_a * cv_a * t_hi)
p_lo, ra_lo = 1.0e5, 1.2
t_lo = (p_lo + pi_a) / (gamma_a * ra_lo * cv_a)
rj_lo = jwl_density_at(p_lo, t_lo)

nx, ny, half, x_slab = 200, 100, 0.25, 0.30

case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "y_domain%beg": -half,
    "y_domain%end": half,
    "m": nx - 1,
    "n": ny - 1,
    "p": 0,
    "dt": 2.0e-7,
    "t_step_start": 0,
    "t_step_stop": 600,
    "t_step_save": 150,
    "num_patches": 2,
    "model_eqns": 2,
    "num_fluids": 2,
    "relax": "T",
    "relax_model": 5,
    "mpp_lim": "F",
    "mixture_err": "F",
    "alt_soundspeed": "T",
    "time_stepper": 3,
    "weno_order": 3,
    "weno_eps": 1e-16,
    "mp_weno": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "fd_order": 2,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "ib": "T",
    "num_ibs": 1,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "F",
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5 * x_slab,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": x_slab,
    "patch_icpp(1)%length_y": 2 * half,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p_hi,
    "patch_icpp(1)%alpha_rho(1)": (1 - eps) * rj_hi,
    "patch_icpp(1)%alpha_rho(2)": eps * ra_hi,
    "patch_icpp(1)%alpha(1)": 1 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(2)%geometry": 3,
    "patch_icpp(2)%x_centroid": x_slab + 0.5 * (1.0 - x_slab),
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%length_x": 1.0 - x_slab,
    "patch_icpp(2)%length_y": 2 * half,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p_lo,
    "patch_icpp(2)%alpha_rho(1)": eps * rj_lo,
    "patch_icpp(2)%alpha_rho(2)": (1 - eps) * ra_lo,
    "patch_icpp(2)%alpha(1)": eps,
    "patch_icpp(2)%alpha(2)": 1 - eps,
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": 0.55,
    "patch_ib(1)%y_centroid": 0.0,
    "patch_ib(1)%radius": 0.05,
    "patch_ib(1)%slip": "T",
    "patch_ib(1)%moving_ibm": moving,
    "fluid_pp(1)%eos": "jwl",
    "fluid_pp(1)%jwl_a": A,
    "fluid_pp(1)%jwl_b": B,
    "fluid_pp(1)%jwl_r1": R1,
    "fluid_pp(1)%jwl_r2": R2,
    "fluid_pp(1)%jwl_omega": omega,
    "fluid_pp(1)%jwl_rho0": rho0,
    "fluid_pp(1)%jwl_t0": t0,
    "fluid_pp(1)%cv": cv_j,
    "fluid_pp(2)%eos": "stiffened_gas",
    "fluid_pp(2)%gamma": 1 / gamma_a,
    "fluid_pp(2)%pi_inf": (1 + 1 / gamma_a) * pi_a,
    "fluid_pp(2)%cv": cv_a,
}
if moving == 1:
    case["patch_ib(1)%vel(1)"] = 100.0
if moving == 2:
    case["patch_ib(1)%vel(1)"] = 0.0
    case["patch_ib(1)%mass"] = 20.0

print(json.dumps(case))
