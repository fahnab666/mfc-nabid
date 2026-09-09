"""Run the repository's six-sphere Mach-3 reproduction in an isolated directory."""

import json
from pathlib import Path

repo = Path(__file__).resolve().parents[3]
case = json.loads((repo / "benchmarks/ibm_moving_validation/daoud.json").read_text())
# Match the effective contact resolution reported for the failed production run.
# Restitution accuracy is checked separately by the 40-step binary-impact case.
case["collision_steps_per_contact"] = 4
print(json.dumps(case))
