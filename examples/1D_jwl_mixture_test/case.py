#!/usr/bin/env python3
import json
import os

eps = 1.0e-8

# Mixture closure to benchmark: 0 isobaric, 1 Kuhl, 2 p-T equilibrium, 3 Rocflu blend.
jwl_mix_type = int(os.environ.get("JWL_MIX_TYPE", "0"))

rho_jwl = 1630.0
rho_air = 1.225

jwl = {
    "A": 3.712e11,
    "B": 3.231e9,
    "R1": 4.15,
    "R2": 0.95,
    "omega": 0.30,
    "rho0": rho_jwl,
    "E0": 1.0089e10,
    "air_e0": 2.5575e5,
    "air_rho0": rho_air,
    "air_gamma": 0.4,
}

print(
    json.dumps(
        {
            "run_time_info": "T",
            # Domain
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "m": 399,
            "n": 0,
            "p": 0,
            "dt": 5.0e-8,
            "t_step_start": 0,
            "t_step_stop": 600,
            "t_step_save": 150,
            # Numerics
            "num_patches": 2,
            "model_eqns": 2,
            "num_fluids": 2,
            "jwl_mix_type": jwl_mix_type,
            "mpp_lim": "T",
            "mixture_err": "T",
            "alt_soundspeed": "F",
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            # Output
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "F",
            "rho_wrt": "T",
            "pres_wrt": "T",
            "c_wrt": "T",
            "parallel_io": "F",
            # Patch 1: mostly air background
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": 101325.0,
            "patch_icpp(1)%alpha_rho(1)": eps * rho_jwl,
            "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_air,
            "patch_icpp(1)%alpha(1)": eps,
            "patch_icpp(1)%alpha(2)": 1.0 - eps,
            # Patch 2: hot JWL detonation-products driver. p must exceed the JWL cold pressure at rho0
            # (~7.1 GPa here) so the products carry positive internal energy / temperature.
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.15,
            "patch_icpp(2)%length_x": 0.30,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": 1.2e10,
            "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_jwl,
            "patch_icpp(2)%alpha_rho(2)": eps * rho_air,
            "patch_icpp(2)%alpha(1)": 1.0 - eps,
            "patch_icpp(2)%alpha(2)": eps,
            # Fluid 1: JWL products
            "fluid_pp(1)%eos": 2,
            "fluid_pp(1)%gamma": 1.0 / 0.4,
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%cv": 613.5,
            "fluid_pp(1)%jwl_A": jwl["A"],
            "fluid_pp(1)%jwl_B": jwl["B"],
            "fluid_pp(1)%jwl_R1": jwl["R1"],
            "fluid_pp(1)%jwl_R2": jwl["R2"],
            "fluid_pp(1)%jwl_omega": jwl["omega"],
            "fluid_pp(1)%jwl_rho0": jwl["rho0"],
            "fluid_pp(1)%jwl_E0": jwl["E0"],
            "fluid_pp(1)%jwl_air_e0": jwl["air_e0"],
            "fluid_pp(1)%jwl_air_rho0": jwl["air_rho0"],
            "fluid_pp(1)%jwl_air_gamma": jwl["air_gamma"],
            # Fluid 2: ideal-gas air (cv set so p-T equilibrium has R_air = (gamma-1)*cv = 0.4*717.5 = 287 J/kg/K)
            "fluid_pp(2)%eos": 1,
            "fluid_pp(2)%gamma": 1.0 / 0.4,
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%cv": 717.5,
        }
    )
)
