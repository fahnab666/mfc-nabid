#!/usr/bin/env python3
# Generic products-expansion: a slab of high-pressure JWL products expands into ideal
# air across a diffuse interface, exercising the mixed-cell pressure-equilibrium
# closure. The JWL coefficients are synthetic and load from the material file
# products.yaml in this directory through the external material-file mechanism.
import json
import math
import os
import sys

case_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(case_dir, "..", "..", "toolchain"))
from mfc.materials import material_fluid

eps = 1.0e-06

# Numerical setup
Nx = 400
dx = 1.0 / (1.0 * (Nx + 1))

Tend, Nt = 2.0e-05, 800
mydt = Tend / (1.0 * Nt)

case = {
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
    "num_fluids": 2,
    "mpp_lim": "F",
    "mixture_err": "T",
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
    # Patch 1: pressurized products slab
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": 0.15,
    "patch_icpp(1)%length_x": 0.3,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%pres": 2.0e10,
    "patch_icpp(1)%alpha_rho(1)": 1600.0 * (1.0 - eps),
    "patch_icpp(1)%alpha_rho(2)": 1.2 * eps,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    # Patch 2: ambient air
    "patch_icpp(2)%geometry": 1,
    "patch_icpp(2)%x_centroid": 0.65,
    "patch_icpp(2)%length_x": 0.7,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%pres": 1.0e05,
    "patch_icpp(2)%alpha_rho(1)": 1600.0 * eps,
    "patch_icpp(2)%alpha_rho(2)": 1.2 * (1.0 - eps),
    "patch_icpp(2)%alpha(1)": eps,
    "patch_icpp(2)%alpha(2)": 1.0 - eps,
    # Fluid 2: ideal air
    "fluid_pp(2)%gamma": 1.0e00 / (1.4 - 1.0e00),
    "fluid_pp(2)%pi_inf": 0.0,
}

# Fluid 1: JWL products from the provenance-tagged material file in this directory
case.update(material_fluid(1, "products.yaml", case_dir))

print(json.dumps(case))
