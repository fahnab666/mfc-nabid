#!/usr/bin/env python3
"""Uniform, initially unequal-temperature JWL/ambient mixture with a known PT equilibrium."""

import json
import math

# Synthetic coefficients; this is a solver verification case, not a material calibration.
A, B, R1, R2, omega, rho0 = 2e10, 1e8, 4.0, 1.1, 0.3, 1000.0
cv_j, cv_a, gamma_a, pi_a, t0 = 1000.0, 750.0, 0.4, 1e7, 300.0
rho_j, temperature, y_j = 800.0, 1800.0, 0.4
pref = A * math.exp(-R1 * rho0 / rho_j) + B * math.exp(-R2 * rho0 / rho_j)
ecold = (A / R1 * math.exp(-R1 * rho0 / rho_j) + B / R2 * math.exp(-R2 * rho0 / rho_j)) / rho0
tref = t0 * (rho_j / rho0) ** omega
pressure = pref + omega * rho_j * cv_j * (temperature - tref)
rho_a = (pressure + pi_a) / (gamma_a * cv_a * temperature)
rho = 1 / (y_j / rho_j + (1 - y_j) / rho_a)
mj, ma = rho * y_j, rho * (1 - y_j)
energy = mj * (ecold + cv_j * (temperature - tref)) + ma * (cv_a * temperature + pi_a / rho_a)

# Perturb the initial phase volumes while keeping the target masses and total energy.
alpha_j = 0.8 * mj / rho_j
initial_rho_j = mj / alpha_j
initial_pref = A * math.exp(-R1 * rho0 / initial_rho_j) + B * math.exp(-R2 * rho0 / initial_rho_j)
initial_ecold = (A / R1 * math.exp(-R1 * rho0 / initial_rho_j) + B / R2 * math.exp(-R2 * rho0 / initial_rho_j)) / rho0
gamma_mix = alpha_j / omega + (1 - alpha_j) / gamma_a
pi_mix = mj * initial_ecold - alpha_j * initial_pref / omega + (1 - alpha_j) * (1 + gamma_a) * pi_a / gamma_a
initial_pressure = (energy - pi_mix) / gamma_mix

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "m": 63,
            "n": 0,
            "p": 0,
            "dt": 1e-7,
            "t_step_start": 0,
            "t_step_stop": 10,
            "t_step_save": 5,
            "num_patches": 1,
            "model_eqns": "5eq",
            "num_fluids": 2,
            "relax": "T",
            "relax_model": 5,
            "mpp_lim": "F",
            "mixture_err": "F",
            "alt_soundspeed": "T",
            "time_stepper": "rk3",
            "weno_order": 3,
            "weno_eps": 1e-16,
            "mp_weno": "F",
            "riemann_solver": "hllc",
            "wave_speeds": "direct",
            "avg_state": "arithmetic",
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "format": "silo",
            "precision": "double",
            "parallel_io": "T",
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "T_wrt": "T",
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": initial_pressure,
            "patch_icpp(1)%alpha_rho(1)": mj,
            "patch_icpp(1)%alpha_rho(2)": ma,
            "patch_icpp(1)%alpha(1)": alpha_j,
            "patch_icpp(1)%alpha(2)": 1 - alpha_j,
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
    )
)
