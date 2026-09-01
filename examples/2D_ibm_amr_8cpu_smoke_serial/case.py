#!/usr/bin/env python3
"""Eight-rank IBM + AMR smoke case using serial per-rank output.

This reuses the distributed case's deterministic physics and turns off the
parallel-I/O path so AMR post-processing reads each rank's own restart data.
"""

import contextlib
import io
import json
import runpy
from pathlib import Path

source = Path(__file__).parents[1] / "2D_ibm_amr_8cpu_smoke" / "case.py"
with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(source)["case"]

case["parallel_io"] = "F"
print(json.dumps(case, indent=4))
