#!/usr/bin/env python3
"""Eight-rank Daoud curtain interaction test sized for a short local run."""

import contextlib
import io
import json
import runpy
from pathlib import Path

with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(Path(__file__).with_name("case.py"))["case"]

# Keep the 10-cells/D coarse mesh and fixed 20-cells/D L1 corridor. Move the
# shock close enough to reach the 2 mm curtain during this local test.
shock_x = -5.0e-4
steps = 200
case.update(
    {
        "patch_icpp(2)%x_centroid": 0.5 * (case["x_domain%beg"] + shock_x),
        "patch_icpp(2)%length_x": shock_x - case["x_domain%beg"],
        "t_step_stop": steps,
        "t_step_save": 100,
        "t_stop": steps * case["dt"],
        "t_save": 100 * case["dt"],
    }
)

print(json.dumps(case, indent=4))
