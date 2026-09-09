"""Quantitative checks for the reduced MFC-paper static-sphere case."""

from pathlib import Path

import numpy as np

root = Path(__file__).parent / "restart_data"
edges = [np.fromfile(root / f"lustre_{axis}_cb.dat", dtype=np.float64) for axis in "xyz"]
shape = tuple(len(edge) - 1 for edge in edges)
last = max(root.glob("lustre_[0-9]*.dat"), key=lambda path: int(path.stem.split("_")[-1]))
q = np.fromfile(last, dtype=np.float64).reshape((*shape, 6), order="F")
rho = q[..., 0]
mom = q[..., 1:4]
energy = q[..., 4]
pressure = 0.4 * (energy - np.sum(mom * mom, axis=-1) / (2.0 * rho))
velocity = mom / rho[..., None]
centers = [(edge[:-1] + edge[1:]) / 2.0 for edge in edges]
x, y, z = np.meshgrid(*centers, indexing="ij", sparse=True)
solid = (x + 3.0e-3) ** 2 + y**2 + z**2 <= 0.05**2
fluid = ~solid
rho_inf = 0.2199
p_inf = 10918.2549
u_inf = 527.2
mach_inf = u_inf / np.sqrt(1.4 * p_inf / rho_inf)
iy = int(np.argmin(abs(centers[1])))
iz = int(np.argmin(abs(centers[2])))
line = rho[:, iy, iz]
upstream = centers[0] < -0.053
compressed = np.flatnonzero(upstream & (line > 1.25 * rho_inf))
shock_x = float(centers[0][compressed[0]]) if compressed.size else float("nan")
standoff_D = (-0.053 - shock_x) / 0.1
rho_sym_y = np.linalg.norm(rho - rho[:, ::-1, :]) / np.linalg.norm(rho)
rho_sym_z = np.linalg.norm(rho - rho[:, :, ::-1]) / np.linalg.norm(rho)

metrics = {
    "checkpoint": last.name,
    "inflow_mach": float(mach_inf),
    "fluid_finite": bool(np.isfinite(q[fluid]).all()),
    "min_fluid_density": float(rho[fluid].min()),
    "min_fluid_pressure": float(pressure[fluid].min()),
    "max_density_ratio": float(rho[fluid].max() / rho_inf),
    "max_pressure_ratio": float(pressure[fluid].max() / p_inf),
    "max_speed": float(np.linalg.norm(velocity[fluid], axis=-1).max()),
    "bow_shock_standoff_D_threshold": standoff_D,
    "density_symmetry_y_rel_l2": float(rho_sym_y),
    "density_symmetry_z_rel_l2": float(rho_sym_z),
}
print(metrics)
assert metrics["fluid_finite"]
assert metrics["min_fluid_density"] > 0.0
assert metrics["min_fluid_pressure"] > 0.0
assert abs(metrics["inflow_mach"] - 2.0) < 0.01
assert metrics["max_density_ratio"] > 1.5
assert np.isfinite(standoff_D) and standoff_D > 0.0
assert max(rho_sym_y, rho_sym_z) < 1.0e-10
