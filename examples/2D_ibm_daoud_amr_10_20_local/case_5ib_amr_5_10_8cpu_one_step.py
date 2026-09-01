#!/usr/bin/env python3
"""One-step completion test for static particle-cloud plus fixed AMR."""

import contextlib
import io
import runpy
from pathlib import Path

source = Path(__file__).with_name("case_5ib_amr_5_10_8cpu.py")
with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(source)["case"]

case.update({"t_stop": case["dt"], "t_save": case["dt"], "t_step_stop": 1, "t_step_save": 1})
print(__import__("json").dumps(case, indent=4))
