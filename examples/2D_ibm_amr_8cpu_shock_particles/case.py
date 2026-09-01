#!/usr/bin/env python3
"""Moderate eight-rank CPU verification: shock, four IBM bodies, and AMR.

This is deliberately more substantial than the local smoke tests: a Mach-1.5
shock crosses four fixed no-slip cylinders over 250 RK3 steps on a 512 x 512
coarse grid.  One static 64 x 64-cell level-1 AMR block resolves the particle
cluster at twice the base resolution.  The block fits an 8-rank 2D split
(fine extent is 127 x 127 cells), so it is suitable for a local MacBook run.

Force-driven particle motion is intentionally excluded: moving_ibm=2 with AMR
is not yet a validated MFC configuration.  This case verifies the supported
multi-body static-IBM AMR path with non-trivial shock dynamics.

Run:
    ./mfc.sh run examples/2D_ibm_amr_8cpu_shock_particles/case.py \
        --no-gpu -n 8 -t pre_process simulation
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
    "schlieren_wrt": "T",
    "schlieren_alpha(1)": 10.0,
    "ib_state_wrt": "T",
    # 512² coarse mesh; the AMR block below is 64² coarse cells.
    "x_domain%beg": 0.0,
    "x_domain%end": 1.0,
    "y_domain%beg": 0.0,
    "y_domain%end": 1.0,
    "m": 511,
    "n": 511,
    "p": 0,
    "dt": 4.0e-4,
    "t_step_start": 0,
    "t_step_stop": 250,
    "t_step_save": 50,
    # Single-fluid Mach-1.5 shock tube.
    "num_patches": 2,
    "model_eqns": "5eq",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": "arithmetic",
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "viscous": "T",
    "fd_order": 2,
    "bc_x%beg": -17,
    "bc_x%end": -3,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    # Ambient gas.
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.5,
    "patch_icpp(1)%y_centroid": 0.5,
    "patch_icpp(1)%length_x": 1.0,
    "patch_icpp(1)%length_y": 1.0,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": 1.0,
    "patch_icpp(1)%alpha_rho(1)": 1.4,
    "patch_icpp(1)%alpha(1)": 1.0,
    # Post-shock state, initially left of x=0.38.
    "patch_icpp(2)%geometry": 3,
    "patch_icpp(2)%x_centroid": 0.19,
    "patch_icpp(2)%y_centroid": 0.5,
    "patch_icpp(2)%length_x": 0.38,
    "patch_icpp(2)%length_y": 1.0,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%vel(1)": 0.6944,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": 2.4583,
    "patch_icpp(2)%alpha_rho(1)": 2.6069,
    "patch_icpp(2)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 2.5,
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%Re(1)": 2.5e5,
    # Four fixed circles inside the fine block, offset from grid symmetry.
    "ib": "T",
    "num_ibs": 4,
    "many_ib_patch_parallelism": "T",
    "patch_ib(1)%geometry": 2,
    "patch_ib(1)%x_centroid": 0.462,
    "patch_ib(1)%y_centroid": 0.466,
    "patch_ib(1)%radius": 0.017,
    "patch_ib(1)%slip": "F",
    "patch_ib(2)%geometry": 2,
    "patch_ib(2)%x_centroid": 0.512,
    "patch_ib(2)%y_centroid": 0.482,
    "patch_ib(2)%radius": 0.017,
    "patch_ib(2)%slip": "F",
    "patch_ib(3)%geometry": 2,
    "patch_ib(3)%x_centroid": 0.477,
    "patch_ib(3)%y_centroid": 0.531,
    "patch_ib(3)%radius": 0.017,
    "patch_ib(3)%slip": "F",
    "patch_ib(4)%geometry": 2,
    "patch_ib(4)%x_centroid": 0.529,
    "patch_ib(4)%y_centroid": 0.539,
    "patch_ib(4)%radius": 0.017,
    "patch_ib(4)%slip": "F",
    # Static AMR, refined 2:1 around the particle cluster.
    "amr": "T",
    "amr_ref_ratio": 2,
    "amr_max_level": 1,
    "amr_max_blocks": 1,
    "amr_max_grid_size": 64,
    "amr_subcycle": "T",
    "amr_block_beg(1)": 224,
    "amr_block_end(1)": 287,
    "amr_block_beg(2)": 224,
    "amr_block_end(2)": 287,
    "amr_regrid_int": 0,
}

print(json.dumps(case, indent=4))
