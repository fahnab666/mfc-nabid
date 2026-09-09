"""Check preservation of the exact uniform solution across a periodic seam."""

from pathlib import Path

import numpy as np

root = Path(__file__).parent / "restart_data"
edges = [np.fromfile(root / f"lustre_{axis}_cb.dat", dtype=np.float64) for axis in "xyz"]
shape = tuple(len(edge) - 1 for edge in edges)
last = max(root.glob("lustre_[0-9]*.dat"), key=lambda path: int(path.stem.split("_")[-1]))
label = last.stem.split("_")[-1]
q = np.fromfile(last, dtype=np.float64).reshape((*shape, 6), order="F")
ib = np.fromfile(root / f"ib_state_{label}.dat", dtype=np.float64).reshape(-1, 20)
rho = q[..., 0]
mom = q[..., 1:4]
pressure = 0.4 * (q[..., 4] - np.sum(mom * mom, axis=-1) / (2.0 * rho))
coords = np.meshgrid(*[(edge[:-1] + edge[1:]) / 2.0 for edge in edges], indexing="ij", sparse=True)
distance2 = 0.0
for axis, grid in enumerate(coords):
    delta = grid - ib[0, 16 + axis]
    width = edges[axis][-1] - edges[axis][0]
    delta -= np.rint(delta / width) * width
    distance2 = distance2 + delta**2
fluid = distance2 > ib[0, 19] ** 2
velocity = mom / rho[..., None]
metrics = {
    "checkpoint": last.name,
    "time": float(ib[0, 0]),
    "particle_centroid_x": float(ib[0, 16]),
    "particle_velocity_x": float(ib[0, 7]),
    "force_norm": float(np.linalg.norm(ib[0, 1:4])),
    "all_finite": bool(np.isfinite(q[fluid]).all() and np.isfinite(ib).all()),
    "max_density_error": float(np.max(abs(rho[fluid] - 1.0))),
    "max_pressure_error": float(np.max(abs(pressure[fluid] - 1.0))),
    "max_velocity_error": float(np.max(abs(velocity[fluid] - np.array([1.0, 0.0, 0.0])))),
}
print(metrics)
assert metrics["all_finite"]
assert metrics["max_density_error"] < 1.0e-10
assert metrics["max_pressure_error"] < 1.0e-10
assert metrics["max_velocity_error"] < 1.0e-10
assert metrics["force_norm"] < 1.0e-10
assert abs(metrics["particle_velocity_x"] - 1.0) < 1.0e-12
