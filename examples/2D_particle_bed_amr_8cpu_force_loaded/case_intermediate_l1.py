#!/usr/bin/env python3
"""Intermediate 48x48 Level-0-cell L1 block diagnostic for eight ranks."""

import contextlib
import io
import json
import runpy
from pathlib import Path

with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(Path(__file__).with_name("case.py"))["case"]
case.update(
    {
        "amr_block_beg(1)": 104,
        "amr_block_end(1)": 151,
        "amr_block_beg(2)": 104,
        "amr_block_end(2)": 151,
        "t_step_stop": 40,
        "t_step_save": 20,
    }
)
print(json.dumps(case, indent=4))
