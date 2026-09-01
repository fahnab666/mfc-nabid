#!/usr/bin/env python3
"""Short 2D five-particle fixed-block AMR test for eight CPU ranks.

The Daoud shock/curtain geometry is retained, with 5 cells per diameter on
L0 and a static 2:1 AMR corridor (10 cells per diameter on L1).
"""

import contextlib
import io
import json
import math
import runpy
from pathlib import Path

source = Path("/Users/fahadnabid/Downloads/case_amr.py")
with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
    ns = runpy.run_path(source)

case = ns["case"].copy()
cells_per_diameter_l0 = 5.0
dx_target = ns["PARTICLE_DIAMETER"] / cells_per_diameter_l0
nx = int(round((ns["X_END"] - ns["X_BEG"]) / dx_target))
ny = int(round((ns["Y_END"] - ns["Y_BEG"]) / dx_target))
dx = (ns["X_END"] - ns["X_BEG"]) / nx
dy = (ns["Y_END"] - ns["Y_BEG"]) / ny


def first_cell_at_or_after(coord, domain_beg, spacing):
    return max(0, int(math.ceil((coord - domain_beg) / spacing - 0.5)))


def last_cell_at_or_before(coord, domain_beg, spacing, count):
    return min(count - 1, int(math.floor((coord - domain_beg) / spacing - 0.5)))


bx0 = first_cell_at_or_after(ns["L1_X0"], ns["X_BEG"], dx)
bx1 = last_cell_at_or_before(ns["L1_X1"], ns["X_BEG"], dx, nx)
by0 = first_cell_at_or_after(ns["L1_Y0"], ns["Y_BEG"], dy)
by1 = last_cell_at_or_before(ns["L1_Y1"], ns["Y_BEG"], dy, ny)
dt_nominal = ns["CFL"] * min(dx, dy) / (2 * ns["MAX_SIGNAL_SPEED"])
short_steps = 20
dt = dt_nominal

case.update(
    {
        "m": nx - 1,
        "n": ny - 1,
        "dt": dt,
        "t_stop": short_steps * dt,
        "t_save": short_steps * dt,
        "t_step_stop": short_steps,
        "t_step_save": short_steps,
        "particle_cloud(1)%num_particles": 5,
        "particle_cloud(1)%moving_ibm": 0,
        "collision_model": 0,
        "amr_block_beg(1)": bx0,
        "amr_block_end(1)": bx1,
        "amr_block_beg(2)": by0,
        "amr_block_end(2)": by1,
        "amr_ref_ratio": 2,
        "amr_max_level": 1,
        "amr_max_blocks": 1,
        "amr_regrid_int": 0,
    }
)

print(json.dumps(case, indent=4))
