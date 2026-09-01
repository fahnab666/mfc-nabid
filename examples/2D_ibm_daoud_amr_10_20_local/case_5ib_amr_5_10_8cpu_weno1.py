#!/usr/bin/env python3
"""Diagnostic variant: disable grid-dependent WENO coefficients during AMR setup."""

import contextlib
import io
import runpy
from pathlib import Path

source = Path(__file__).with_name("case_5ib_amr_5_10_8cpu.py")
with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(source)["case"]

case.update({"weno_order": 1, "mapped_weno": "F", "mp_weno": "F"})
print(__import__("json").dumps(case, indent=4))
