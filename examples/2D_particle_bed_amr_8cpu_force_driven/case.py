#!/usr/bin/env python3
"""Force-driven counterpart of the minimal particle-bed + AMR baseline.

This intentionally changes only particle_cloud(1)%moving_ibm to 2 and shortens
the run to 40 steps.  It is a diagnostic case for the currently unvalidated
force-driven particle-cloud + AMR path.
"""

import contextlib
import io
import json
import runpy
from pathlib import Path

source = Path(__file__).parents[1] / "2D_particle_bed_amr_8cpu_static" / "case.py"
with contextlib.redirect_stdout(io.StringIO()):
    case = runpy.run_path(source)["case"]

case["particle_cloud(1)%moving_ibm"] = 2
case["t_step_stop"] = 40
case["t_step_save"] = 20
case["collision_model"] = 1
case["coefficient_of_restitution"] = 0.9
case["collision_time"] = 0.01
case["ib_coefficient_of_friction"] = 0.1

print(json.dumps(case, indent=4))
