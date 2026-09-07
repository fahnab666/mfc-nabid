"""Check saved single-ideal-gas states using the actual MPI checkpoint grid."""

import argparse
import json
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("directory", type=Path)
parser.add_argument("--impact", action="store_true")
parser.add_argument("--stop", type=float, required=True)
args = parser.parse_args()
data = args.directory / "restart_data"
edges = [np.fromfile(data / f"lustre_{axis}_cb.dat", dtype=np.float64) for axis in "xy"]
zpath = data / "lustre_z_cb.dat"
is3d = zpath.exists() and zpath.stat().st_size > 8
if is3d:
    edges.append(np.fromfile(zpath, dtype=np.float64))
shape = tuple(len(e) - 1 for e in edges)
coords = np.meshgrid(*[(e[1:] + e[:-1]) / 2 for e in edges], indexing="ij", sparse=True)
rows = []
for path in sorted(data.glob("ib_state_*.dat"), key=lambda p: int(p.stem.split("_")[-1])):
    label = path.stem.split("_")[-1]
    ib = np.fromfile(path, dtype=np.float64).reshape(-1, 20)
    state = np.fromfile(data / f"lustre_{label}.dat", dtype=np.float64)
    q = state.reshape((*shape, len(shape) + 3), order="F")
    rho = q[..., 0]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        pressure = 0.4 * (q[..., len(shape) + 1] - np.sum(q[..., 1 : len(shape) + 1] ** 2, axis=-1) / (2 * rho))
    fluid = np.ones(shape, dtype=bool)
    for particle in ib:
        distance2 = 0.0
        for axis, grid in enumerate(coords):
            delta = grid - particle[16 + axis]
            if axis > 0:
                width = edges[axis][-1] - edges[axis][0]
                delta -= np.rint(delta / width) * width
            distance2 = distance2 + delta**2
        fluid &= distance2 > particle[19] ** 2
    gaps = []
    for i, particle in enumerate(ib):
        for other in ib[:i]:
            delta = particle[16:19] - other[16:19]
            for axis in range(1, len(shape)):
                width = edges[axis][-1] - edges[axis][0]
                delta[axis] -= np.rint(delta[axis] / width) * width
            gaps.append((np.linalg.norm(delta) - particle[19] - other[19]) / (particle[19] + other[19]))
    invalid = ~np.isfinite(pressure) | (pressure <= 0) | (rho <= 0) | ~np.all(np.isfinite(q), axis=-1)
    rows.append(
        {
            "save": int(label),
            "time_s": float(ib[0, 0]),
            "nonfinite_values": int((~np.isfinite(q)).sum()),
            "invalid_fluid_cells": int((invalid & fluid).sum()),
            "invalid_solid_cells": int((invalid & ~fluid).sum()),
            "min_fluid_pressure": float(pressure[fluid].min()),
            "min_fluid_density": float(rho[fluid].min()),
            "minimum_pair_gap_D": float(min(gaps)),
            "max_particle_speed": float(np.linalg.norm(ib[:, 7:10], axis=1).max()),
            "particles_finite": bool(np.isfinite(ib).all()),
        }
    )
summary = {
    "checkpoints": len(rows),
    "last_time_s": rows[-1]["time_s"],
    "max_overlap_D": max(0.0, -min(r["minimum_pair_gap_D"] for r in rows)),
    "all_finite": all(r["nonfinite_values"] == 0 and r["particles_finite"] for r in rows),
    "all_fluid_states_positive": all(r["invalid_fluid_cells"] == 0 for r in rows),
    "max_invalid_solid_cells": max(r["invalid_solid_cells"] for r in rows),
}
summary["reached_stop"] = summary["last_time_s"] >= args.stop * (1 - 1e-8)
if args.impact:
    summary["restitution"] = float(ib[1, 7] - ib[0, 7])
    summary["momentum_x"] = float(ib[:, 7].sum())
print(json.dumps({"summary": summary, "checkpoints": rows}, indent=2))
assert summary["all_finite"] and summary["all_fluid_states_positive"] and summary["reached_stop"], summary
if args.impact:
    assert abs(summary["restitution"] - 0.9) < 0.01, summary
    assert abs(summary["momentum_x"]) < 1e-6, summary
