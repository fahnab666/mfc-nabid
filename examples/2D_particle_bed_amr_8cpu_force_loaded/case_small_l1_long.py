#!/usr/bin/env python3
"""Longer force-driven IBM run using the validated small L1 AMR block."""

import contextlib
import io
import json
import runpy
from pathlib import Path

with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(Path(__file__).with_name("case.py"))["case"]
case["t_step_stop"] = 80
case["t_step_save"] = 20
print(json.dumps(case, indent=4))
