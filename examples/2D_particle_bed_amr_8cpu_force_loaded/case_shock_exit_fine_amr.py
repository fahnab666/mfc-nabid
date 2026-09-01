#!/usr/bin/env python3
"""Long moving-particle shock run with fine tiled L1 AMR."""

import contextlib
import io
import json
import runpy
from pathlib import Path

source = Path(__file__).with_name("case_4x_area_amr_domain.py")
with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(source)["case"]

# Start the shock upstream of the bed.
case["amr_subcycle"] = "F"
case["patch_icpp(2)%x_centroid"] = 0.15
case["patch_icpp(2)%length_x"] = 0.30

# Cover the bed and the downstream trajectory to the outflow boundary.
case.update(
    {
        "amr_block_beg(1)": 80,
        "amr_block_end(1)": 207,
        "amr_block_beg(2)": 80,
        "amr_block_end(2)": 175,
        "amr_max_blocks": 8,
        "amr_max_grid_size": 96,
        "amr_regrid_int": 8,
        "amr_tag_eps": 0.05,
        "amr_buf": 16,
        "t_step_stop": 1600,
        "t_step_save": 200,
    }
)

print(json.dumps(case, indent=4))
