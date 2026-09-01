#!/usr/bin/env python3
"""Short 2D fixed-block AMR check for the Daoud curtain geometry.

Uses the 10/20 cells-per-diameter mesh and exact static AMR corridor from
case_10_20_local, with bodies held fixed so this isolates block AMR + IBM
from unvalidated force-driven AMR coupling.
"""

import contextlib
import io
import json
import runpy
from pathlib import Path

source = Path(__file__).with_name("case.py")
with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(source)["case"]

short_stop = 8 * case["dt"]
case.update(
    {
        "particle_cloud(1)%moving_ibm": 0,
        "collision_model": 0,
        "t_stop": short_stop,
        "t_save": short_stop,
        "t_step_stop": 8,
        "t_step_save": 8,
    }
)

print(json.dumps(case, indent=4))
