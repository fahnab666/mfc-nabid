"""Nearly dry binary impact: restitution and momentum from a full MFC run."""

import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument("--dt", type=float, default=0.00025)
args, _ = parser.parse_known_args()
case = {
    "m": 95,
    "n": 47,
    "p": 0,
    "x_domain%beg": -0.6,
    "x_domain%end": 0.6,
    "y_domain%beg": -0.3,
    "y_domain%end": 0.3,
    "model_eqns": "5eq",
    "num_fluids": 1,
    "num_patches": 1,
    "time_stepper": "rk3",
    "dt": args.dt,
    "t_step_start": 0,
    "t_step_stop": round(0.025 / args.dt),
    "t_step_save": 1,
    "weno_order": 3,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "ib": "T",
    "num_ibs": 2,
    "ib_state_wrt": "T",
    "run_time_info": "T",
    "many_ib_patch_parallelism": "T",
    "ib_neighborhood_radius": 0,
    "collision_model": 1,
    "collision_time": 0.01,
    "coefficient_of_restitution": 0.9,
    "ib_coefficient_of_friction": 0.0,
    "fd_order": 2,
    "parallel_io": "T",
    "precision": "double",
    "format": "silo",
    "prim_vars_wrt": "T",
    "fluid_pp(1)%gamma": 2.5,
    "fluid_pp(1)%pi_inf": 0.0,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": 1.2,
    "patch_icpp(1)%length_y": 0.6,
    "patch_icpp(1)%alpha_rho(1)": 1e-9,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%pres": 1e-6,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
}
for i, sign in [(1, -1), (2, 1)]:
    for key, val in {
        "geometry": 2,
        "x_centroid": sign * 0.1025,
        "y_centroid": 0.0,
        "radius": 0.1,
        "mass": 1.0,
        "moving_ibm": 2,
        "slip": "T",
        "vel(1)": -sign * 0.5,
        "vel(2)": 0.0,
        "vel(3)": 0.0,
    }.items():
        case[f"patch_ib({i})%{key}"] = val
print(json.dumps(case))
