#!/usr/bin/env python3
"""Short serial initialization check for the nested Daoud AMR mesh."""

import contextlib
import io
import json
import runpy
from pathlib import Path

with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(Path(__file__).with_name("case_daoud_2d_l0_5_l1_10_l2_20.py"))["case"]

case["t_step_stop"] = 4
case["t_step_save"] = 4
case["t_stop"] = 4 * case["dt"]
case["t_save"] = case["t_stop"]
print(json.dumps(case, indent=4))
