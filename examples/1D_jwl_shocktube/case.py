#!/usr/bin/env python3
# Single-fluid JWL shock tube. Both states sit on or above the principal isentrope of a
# synthetic (non-calibrated) parameter set, so the run exercises the JWL reference
# curve, its inversion, and the frozen sound speed without any material catalog.
import json
import math

# Numerical setup
Nx = 400
dx = 1.0 / (1.0 * (Nx + 1))

Tend, Nt = 5.0e-05, 500
mydt = Tend / (1.0 * Nt)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0e00,
            "x_domain%end": 1.0e00,
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": mydt,
            "t_step_start": 0,
            "t_step_stop": int(Nt),
            "t_step_save": int(math.ceil(Nt / 10.0)),
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": "5eq",
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": "rk3",
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": "hllc",
            "wave_speeds": "direct",
            "avg_state": "arithmetic",
            "bc_x%beg": -3,
            "bc_x%end": -3,
            # Formatted Database Files Structure Parameters
            "format": "silo",
            "precision": "double",
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: compressed products at the reference density
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.25,
            "patch_icpp(1)%length_x": 0.5,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": 2.0e10,
            "patch_icpp(1)%alpha_rho(1)": 1600.0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Patch 2: expanded products
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": 3.0e09,
            "patch_icpp(2)%alpha_rho(1)": 1200.0,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Fluids Physical Parameters: user-supplied synthetic JWL coefficients
            "fluid_pp(1)%eos": "jwl",
            "fluid_pp(1)%jwl_A": 5.0e11,
            "fluid_pp(1)%jwl_B": 8.0e09,
            "fluid_pp(1)%jwl_R1": 4.5,
            "fluid_pp(1)%jwl_R2": 1.2,
            "fluid_pp(1)%jwl_omega": 0.3,
            "fluid_pp(1)%jwl_rho0": 1600.0,
        }
    )
)
