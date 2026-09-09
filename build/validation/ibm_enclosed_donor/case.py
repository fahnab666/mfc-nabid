"""Exercise the IBM image-point fallback inside overlapping solid bodies."""

import json

case = {
    "run_time_info": "T",
    "x_domain%beg": -0.5,
    "x_domain%end": 0.5,
    "y_domain%beg": -0.5,
    "y_domain%end": 0.5,
    "z_domain%beg": -0.5,
    "z_domain%end": 0.5,
    "m": 47,
    "n": 47,
    "p": 47,
    "dt": 1.0e-4,
    "t_step_start": 0,
    "t_step_stop": 2,
    "t_step_save": 2,
    "num_patches": 1,
    "model_eqns": "5eq",
    "num_fluids": 1,
    "time_stepper": "rk3",
    "weno_order": 3,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "bc_z%beg": -3,
    "bc_z%end": -3,
    "ib": "T",
    "num_ibs": 2,
    "many_ib_patch_parallelism": "T",
    "fd_order": 2,
    "parallel_io": "T",
    "precision": "double",
    "format": "silo",
    "prim_vars_wrt": "T",
    "fluid_pp(1)%gamma": 2.5,
    "fluid_pp(1)%pi_inf": 0.0,
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": 1.0,
    "patch_icpp(1)%length_y": 1.0,
    "patch_icpp(1)%length_z": 1.0,
    "patch_icpp(1)%alpha_rho(1)": 1.0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%pres": 1.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
}

for index, center in enumerate((-0.05, 0.05), start=1):
    case[f"patch_ib({index})%geometry"] = 8
    case[f"patch_ib({index})%x_centroid"] = center
    case[f"patch_ib({index})%y_centroid"] = 0.0
    case[f"patch_ib({index})%z_centroid"] = 0.0
    case[f"patch_ib({index})%radius"] = 0.15
    case[f"patch_ib({index})%mass"] = 1.0
    case[f"patch_ib({index})%moving_ibm"] = 1
    case[f"patch_ib({index})%slip"] = "F"

print(json.dumps(case))
