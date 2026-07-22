#!/usr/bin/env python3
# Benchmark jwl_eos_weno_order_3_riemann_solver_2
# Additional Benchmarked Features
# - fluid_pp(1)%eos : jwl (generalized Mie-Gruneisen mixture accumulators)
# - model_equations : 2
# - time_stepper : 3
# - weno_order : 3
# - riemann_solver : 2
#
# Synthetic O(1) JWL coefficients (no material catalog): a ball of JWL products
# expands into an ideal-gas ambient, so every Riemann face and conversion call
# runs the mg_mixture accumulator path, with pure-phase, mixed, and
# eps-vanished-phase cells all present. Patch states match the admissible set
# used by the "JWL EOS" regression tests (e > 0, c^2 > 0 against the principal
# isentrope p_s(rho=1) ~ 2.97).

import argparse
import json
import math

parser = argparse.ArgumentParser(
    prog="5eq_jwl_weno3_hllc",
    description="This MFC case was created for the purposes of benchmarking MFC.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--gbpp", type=int, metavar="MEM", default=16, help="Adjusts the problem size per rank to fit into [MEM] GB of GPU memory per GPU.")
parser.add_argument("--steps", type=int, default=None, help="Override t_step_stop/t_step_save.")

ARGS = vars(parser.parse_args())
DICT = ARGS["mfc"]

size = 1 if DICT["gpu"] else 0

ppg = 8000000 / 16.0
procs = DICT["nodes"] * DICT["tasks_per_node"]
ncells = math.floor(ppg * procs * ARGS["gbpp"])
s = math.floor((ncells / 2.0) ** (1 / 3))
Nx, Ny, Nz = 2 * s, s, s

eps = 1e-06

# Products ball: radius and centroid sit on the y/z symmetry planes so the
# reflective boundaries complete the sphere.
r0 = 0.5

# Fastest signal is the products' frozen sound speed at the ball state
# (p = 10, rho ~ 1, omega = 0.3, isentrope terms O(10)), c ~ 6; the domain is
# 4 x 2 x 2, so a conservative acoustic CFL fixes dt independent of resolution
# via the x cell width.
dx = 4.0 / Nx
c_max = 8.0
cfl = 0.3
dt = cfl * dx / c_max

n_steps = ARGS["steps"] if ARGS["steps"] is not None else int(2 * (5 * size + 5))

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -2.0,
            "x_domain%end": 2.0,
            "y_domain%beg": 0.0,
            "y_domain%end": 2.0,
            "z_domain%beg": 0.0,
            "z_domain%end": 2.0,
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "cyl_coord": "F",
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": n_steps,
            "t_step_save": n_steps,
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 2,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 3,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -2,
            "bc_y%end": -3,
            "bc_z%beg": -2,
            "bc_z%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "F",
            # Patch 1: ambient (ideal air, fluid 2)
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 1.0,
            "patch_icpp(1)%z_centroid": 1.0,
            "patch_icpp(1)%length_x": 4.0,
            "patch_icpp(1)%length_y": 2.0,
            "patch_icpp(1)%length_z": 2.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0 * eps,
            "patch_icpp(1)%alpha_rho(2)": 0.5 * (1 - eps),
            "patch_icpp(1)%alpha(1)": eps,
            "patch_icpp(1)%alpha(2)": 1 - eps,
            # Patch 2: JWL products ball (fluid 1)
            "patch_icpp(2)%geometry": 8,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%z_centroid": 0.0,
            "patch_icpp(2)%radius": r0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(2)%pres": 10.0,
            "patch_icpp(2)%alpha_rho(1)": 1.0 * (1 - eps),
            "patch_icpp(2)%alpha_rho(2)": 0.5 * eps,
            "patch_icpp(2)%alpha(1)": 1 - eps,
            "patch_icpp(2)%alpha(2)": eps,
            # Fluids Physical Parameters: user-supplied synthetic JWL coefficients
            "fluid_pp(1)%eos": 4,
            "fluid_pp(1)%jwl_A": 50.0,
            "fluid_pp(1)%jwl_B": 8.0,
            "fluid_pp(1)%jwl_R1": 4.5,
            "fluid_pp(1)%jwl_R2": 1.2,
            "fluid_pp(1)%jwl_omega": 0.3,
            "fluid_pp(1)%jwl_rho0": 1.0,
            "fluid_pp(2)%gamma": 1.0e00 / (1.4 - 1.0e00),
            "fluid_pp(2)%pi_inf": 0.0,
        }
    )
)
