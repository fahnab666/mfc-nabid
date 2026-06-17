#!/usr/bin/env python3
import argparse
import json

parser = argparse.ArgumentParser(description="Compact JWL closure-cost benchmark")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--mix-type", type=int, choices=(0, 1, 2, 3), default=0)
args = parser.parse_args()

eps = 1.0e-8
rho_jwl = 1630.0
rho_air = 1.225

params = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "m": 399,
    "n": 0,
    "p": 0,
    "dt": 5.0e-8,
    "t_step_start": 0,
    "t_step_stop": 200,
    "t_step_save": 200,
    "num_patches": 2,
    "model_eqns": 2,
    "num_fluids": 2,
    "jwl_mix_type": args.mix_type,
    "mpp_lim": "T",
    "mixture_err": "T",
    "time_stepper": 3,
    "recon_type": 1,
    "weno_order": 3,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "rho_wrt": "T",
    "pres_wrt": "T",
    "c_wrt": "T",
    "parallel_io": "F",
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%length_x": 1.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%pres": 101325.0,
    "patch_icpp(1)%alpha_rho(1)": eps * rho_jwl,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_air,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1.0 - eps,
    "patch_icpp(2)%geometry": 1,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": 0.15,
    "patch_icpp(2)%length_x": 0.3,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%pres": 1.2e10,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_jwl,
    "patch_icpp(2)%alpha_rho(2)": eps * rho_air,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    "fluid_pp(1)%eos": 2,
    "fluid_pp(1)%gamma": 2.5,
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%cv": 613.5,
    "fluid_pp(1)%jwl_A": 3.712e11,
    "fluid_pp(1)%jwl_B": 3.231e9,
    "fluid_pp(1)%jwl_R1": 4.15,
    "fluid_pp(1)%jwl_R2": 0.95,
    "fluid_pp(1)%jwl_omega": 0.30,
    "fluid_pp(1)%jwl_rho0": rho_jwl,
    "fluid_pp(1)%jwl_E0": 1.0089e10,
    "fluid_pp(1)%jwl_air_e0": 2.5575e5,
    "fluid_pp(1)%jwl_air_rho0": rho_air,
    "fluid_pp(1)%jwl_air_gamma": 0.4,
    "fluid_pp(2)%eos": 1,
    "fluid_pp(2)%gamma": 2.5,
    "fluid_pp(2)%pi_inf": 0.0,
    "fluid_pp(2)%cv": 717.5,
}

print(json.dumps(params))
