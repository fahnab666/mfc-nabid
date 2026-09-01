#!/usr/bin/env python3
"""Runnable eight-step local performance benchmark for the 10/20 Daoud case."""

import contextlib
import io
import json
import runpy
from pathlib import Path

source = Path(__file__).with_name("case.py")
with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(source)["case"]

benchmark_steps = 8
case["t_step_stop"] = benchmark_steps
case["t_step_save"] = benchmark_steps
case["t_stop"] = benchmark_steps * case["dt"]
case["t_save"] = case["t_stop"]

print(json.dumps(case, indent=4))
