#!/usr/bin/env python3
"""Eight-rank 10-minute Daoud curtain interaction test at 10/20 cells per D."""

import contextlib
import io
import json
import math
import runpy
from pathlib import Path

with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(Path(__file__).with_name("case.py"))["case"]

x_beg, x_end = -0.025, 0.060
y_beg, y_end = -0.002, 0.002
d = 115.0e-6
dx_target = d / 10.0
nx = int(round((x_end - x_beg) / dx_target))
ny = int(round((y_end - y_beg) / dx_target))
dx = (x_end - x_beg) / nx
dy = (y_end - y_beg) / ny


def first_cell(coord, beg, spacing):
    return max(0, int(math.ceil((coord - beg) / spacing - 0.5)))


def last_cell(coord, beg, spacing, count):
    return min(count - 1, int(math.floor((coord - beg) / spacing - 0.5)))


shock_x = -1.0e-4
bx0 = first_cell(-0.0015, x_beg, dx)
bx1 = last_cell(0.0035, x_beg, dx, nx)
by0 = first_cell(-0.0010, y_beg, dy)
by1 = last_cell(0.0010, y_beg, dy, ny)
steps = 160

case.update(
    {
        "x_domain%beg": x_beg,
        "x_domain%end": x_end,
        "y_domain%beg": y_beg,
        "y_domain%end": y_end,
        "m": nx - 1,
        "n": ny - 1,
        "cfl_adap_dt": "F",
        "patch_icpp(1)%x_centroid": 0.5 * (x_beg + x_end),
        "patch_icpp(1)%y_centroid": 0.5 * (y_beg + y_end),
        "patch_icpp(1)%length_x": x_end - x_beg,
        "patch_icpp(1)%length_y": y_end - y_beg,
        "patch_icpp(2)%x_centroid": 0.5 * (x_beg + shock_x),
        "patch_icpp(2)%y_centroid": 0.5 * (y_beg + y_end),
        "patch_icpp(2)%length_x": shock_x - x_beg,
        "patch_icpp(2)%length_y": y_end - y_beg,
        "amr_block_beg(1)": bx0,
        "amr_block_end(1)": bx1,
        "amr_block_beg(2)": by0,
        "amr_block_end(2)": by1,
        "t_step_stop": steps,
        "t_step_save": 80,
        "t_stop": steps * case["dt"],
        "t_save": 80 * case["dt"],
    }
)

print(json.dumps(case, indent=4))
