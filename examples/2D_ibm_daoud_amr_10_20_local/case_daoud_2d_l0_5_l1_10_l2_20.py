#!/usr/bin/env python3
"""2D Daoud curtain analogue: L0=5, L1=10, L2=20 cells per diameter."""

import contextlib
import io
import json
import math
import runpy
from pathlib import Path

with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
    ns = runpy.run_path("/Users/fahadnabid/Downloads/case_amr.py")

case = ns["case"].copy()
d = ns["PARTICLE_DIAMETER"]
x_beg, x_end = -0.080, 0.120
y_beg, y_end = -0.006, 0.006
l0_cells_per_diameter = 5.0
dx_target = d / l0_cells_per_diameter
nx = int(round((x_end - x_beg) / dx_target))
ny = int(round((y_end - y_beg) / dx_target))
dx = (x_end - x_beg) / nx
dy = (y_end - y_beg) / ny


def first_cell(coord, beg, spacing):
    return max(0, int(math.ceil((coord - beg) / spacing - 0.5)))


def last_cell(coord, beg, spacing, count):
    return min(count - 1, int(math.floor((coord - beg) / spacing - 0.5)))


def nested_block(lo, hi, beg, spacing, count):
    i0 = first_cell(lo, beg, spacing)
    i1 = last_cell(hi, beg, spacing, count)
    # The static L2 builder uses a quarter-width inset on each side.
    extent = 4 * ((i1 - i0 + 1) // 4)
    return i0, i0 + extent - 1


# L2 is the centre half of this L1 block: x about [-1, 45] mm and y about
# [-1.25, 1.25] mm. L1 provides a 10-cells/D transition to the coarse ends.
bx0, bx1 = nested_block(-0.024, 0.068, x_beg, dx, nx)
by0, by1 = nested_block(-0.0025, 0.0025, y_beg, dy, ny)
if 2 * ((bx1 - bx0 + 1) // 2) > bx1 - bx0 + 1 or 2 * ((by1 - by0 + 1) // 2) > by1 - by0 + 1:
    raise RuntimeError("L1 dimensions must support the nested L2 block")

ref_ratio = 2
cfl = 0.30
dt_nominal = cfl * min(dx, dy) / ns["MAX_SIGNAL_SPEED"]
t_end = 1.0e-4
steps = int(math.ceil(t_end / dt_nominal))
dt = t_end / steps
save_every = max(1, int(round(1.0e-5 / dt)))

case.update(
    {
        "x_domain%beg": x_beg,
        "x_domain%end": x_end,
        "y_domain%beg": y_beg,
        "y_domain%end": y_end,
        "m": nx - 1,
        "n": ny - 1,
        "dt": dt,
        "cfl_adap_dt": "F",
        "t_stop": t_end,
        "t_save": save_every * dt,
        "t_step_stop": steps,
        "t_step_save": save_every,
        "patch_icpp(1)%x_centroid": 0.5 * (x_beg + x_end),
        "patch_icpp(1)%y_centroid": 0.5 * (y_beg + y_end),
        "patch_icpp(1)%length_x": x_end - x_beg,
        "patch_icpp(1)%length_y": y_end - y_beg,
        "patch_icpp(2)%x_centroid": 0.5 * (x_beg + ns["SHOCK_X"]),
        "patch_icpp(2)%y_centroid": 0.5 * (y_beg + y_end),
        "patch_icpp(2)%length_x": ns["SHOCK_X"] - x_beg,
        "patch_icpp(2)%length_y": y_end - y_beg,
        "amr_block_beg(1)": bx0,
        "amr_block_end(1)": bx1,
        "amr_block_beg(2)": by0,
        "amr_block_end(2)": by1,
        "amr_ref_ratio": ref_ratio,
        "amr_max_level": 2,
        "amr_max_blocks": 2,
        "amr_regrid_int": 0,
        "amr_subcycle": "T",
        "collision_time": 5.0 * dt,
    }
)

print(json.dumps(case, indent=4))
