#!/usr/bin/env python3
"""Daoud IBM-AMR case at 10 L0 / 20 L1 cells per particle diameter.

The physical domain, shock, particle cloud, collision settings, fixed AMR
corridor, and 1 ms end time are inherited from case_amr.py.  Only the mesh and
time step change from the original 15/30 cells-per-diameter configuration.
"""

import contextlib
import io
import json
import math
import runpy
from pathlib import Path

source = Path("/Users/fahadnabid/Downloads/case_amr.py")
with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
    namespace = runpy.run_path(source)

case = namespace["case"].copy()
cells_per_diameter_l0 = 10.0
dx_target = namespace["PARTICLE_DIAMETER"] / cells_per_diameter_l0
nx = int(round((namespace["X_END"] - namespace["X_BEG"]) / dx_target))
ny = int(round((namespace["Y_END"] - namespace["Y_BEG"]) / dx_target))
dx = (namespace["X_END"] - namespace["X_BEG"]) / nx
dy = (namespace["Y_END"] - namespace["Y_BEG"]) / ny


def first_cell_at_or_after(coord, domain_beg, spacing):
    return max(0, int(math.ceil((coord - domain_beg) / spacing - 0.5)))


def last_cell_at_or_before(coord, domain_beg, spacing, count):
    return min(count - 1, int(math.floor((coord - domain_beg) / spacing - 0.5)))


bx0 = first_cell_at_or_after(namespace["L1_X0"], namespace["X_BEG"], dx)
bx1 = last_cell_at_or_before(namespace["L1_X1"], namespace["X_BEG"], dx, nx)
by0 = first_cell_at_or_after(namespace["L1_Y0"], namespace["Y_BEG"], dy)
by1 = last_cell_at_or_before(namespace["L1_Y1"], namespace["Y_BEG"], dy, ny)
ref_ratio = 2
max_signal_speed = namespace["MAX_SIGNAL_SPEED"]
dt_nominal = namespace["CFL"] * min(dx, dy) / (ref_ratio * max_signal_speed)
num_steps = int(math.ceil(namespace["T_END"] / dt_nominal))
dt = namespace["T_END"] / num_steps

case.update(
    {
        "m": nx - 1,
        "n": ny - 1,
        "dt": dt,
        "t_step_stop": num_steps,
        "t_step_save": max(1, int(round(1.0e-5 / dt))),
        "amr_block_beg(1)": bx0,
        "amr_block_end(1)": bx1,
        "amr_block_beg(2)": by0,
        "amr_block_end(2)": by1,
        "collision_time": 20.0 * dt / 4.0,
    }
)

print(json.dumps(case, indent=4))
