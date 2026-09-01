#!/usr/bin/env python3
"""Small deterministic 2D static-IBM + one-block-AMR smoke case.

This case is intentionally not a physical benchmark.  It isolates the AMR
fine-block and immersed-boundary setup on a laptop-friendly grid:

* 64 x 64 coarse cells with one fixed 24 x 24-cell AMR block (2:1);
* one stationary circular no-slip body with an off-grid-line centre;
* quiescent ideal gas, so the expected field remains uniform;
* 20 fixed-size steps and only initial/final output.

Run locally on CPU:
    ./mfc.sh run examples/2D_ibm_amr_local_smoke/case.py --no-gpu

Run on a CUDA/NVHPC host with one GPU:
    ./mfc.sh run examples/2D_ibm_amr_local_smoke/case.py --gpu=acc -n 1 -g 0

The intentionally small block is not a six-rank scaling case.  IBM AMR blocks
are owned whole, so use the full Daoud case (or a separately enlarged block)
when validating six GPUs.
"""

import json

case = {
    "run_time_info": "T",
    "format": "silo",
    "precision": "double",
    "parallel_io": "F",
    "prim_vars_wrt": "T",
    "rho_wrt": "T",
    "pres_wrt": "T",
    "vel_wrt(1)": "T",
    "vel_wrt(2)": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "y_domain%beg": 0.0,
    "y_domain%end": 1.0,
    "m": 63,
    "n": 63,
    "p": 0,
    "dt": 1.0e-4,
    "t_step_start": 0,
    "t_step_stop": 20,
    "t_step_save": 20,
    "num_patches": 1,
    "model_eqns": "5eq",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%y_centroid": 0.5,
    "patch_icpp(1)%length_x": 1.0,
    "patch_icpp(1)%length_y": 1.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 1.0,
    "patch_icpp(1)%alpha_rho(1)": 1.0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 2.5,
    "fluid_pp(1)%pi_inf": 0.0,
    "ib": "T",
    "num_ibs": 1,
    "fd_order": 2,
    "viscous": "F",
    # Deliberately asymmetric: no image point lies on a cell boundary.
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": 0.488,
    "patch_ib(1)%y_centroid": 0.513,
    "patch_ib(1)%radius": 0.10,
    "patch_ib(1)%slip": "F",
    "amr": "T",
    "amr_ref_ratio": 2,
    "amr_max_level": 1,
    "amr_max_blocks": 1,
    "amr_max_grid_size": 24,
    "amr_subcycle": "T",
    "amr_block_beg(1)": 20,
    "amr_block_end(1)": 43,
    "amr_block_beg(2)": 20,
    "amr_block_end(2)": 43,
    "amr_regrid_int": 0,
}

print(json.dumps(case, indent=4))
