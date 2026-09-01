#!/usr/bin/env python3
"""Force-driven IBM particle bed with a four-times wider L1 AMR domain."""

import contextlib
import io
import json
import runpy
from pathlib import Path

source = Path(__file__).with_name("case.py")
with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(source)["case"]

# Original block is 32x32 Level-0 cells: [112,143] in both directions.
# Expand fourfold about the bed center (128,128): [64,191].
case.update(
    {
        "amr_block_beg(1)": 80,
        "amr_block_end(1)": 191,
        "amr_block_beg(2)": 80,
        "amr_block_end(2)": 191,
        "amr_max_blocks": 16,
        "amr_max_grid_size": 32,
        "t_step_stop": 80,
        "t_step_save": 20,
    }
)

print(json.dumps(case, indent=4))
