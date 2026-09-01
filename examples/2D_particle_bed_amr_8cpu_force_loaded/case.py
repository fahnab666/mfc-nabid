#!/usr/bin/env python3
"""Short, shock-loaded force-driven particle-bed + AMR diagnostic."""

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
# Initial shock ends at x=0.46, so it meets the x=0.5 particle bed here.
case["patch_icpp(2)%x_centroid"] = 0.23
case["patch_icpp(2)%length_x"] = 0.46
case["collision_model"] = 1
case["coefficient_of_restitution"] = 0.9
case["collision_time"] = 0.01
case["ib_coefficient_of_friction"] = 0.1

print(json.dumps(case, indent=4))
