#!/usr/bin/env python3
"""Force-driven IBM particle bed with a generic tiled L1 AMR domain."""

import contextlib
import io
import json
import runpy
from pathlib import Path

source = Path(__file__).with_name("case.py")
with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(source)["case"]

# Double each side of the original 32x32 coarse-cell block: 64x64 cells,
# centered at the bed. This is four times the refined area.
case["amr_subcycle"] = "T"
case["particle_cloud(1)%moving_ibm"] = 2
case["particle_cloud(1)%num_particles"] = 4
case["particle_cloud(1)%length_x"] = 0.09
case["particle_cloud(1)%length_y"] = 0.09
case["particle_cloud(1)%radius"] = 0.015
case.update(
    {
        "amr_block_beg(1)": 80,
        "amr_block_end(1)": 175,
        "amr_block_beg(2)": 80,
        "amr_block_end(2)": 175,
        "amr_max_blocks": 64,
        "amr_max_grid_size": 32,
        "t_step_stop": 80,
        "t_step_save": 40,
    }
)

print(json.dumps(case, indent=4))
