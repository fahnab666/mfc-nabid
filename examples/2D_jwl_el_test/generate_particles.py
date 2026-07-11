#!/usr/bin/env python3
# Generate the Euler-Lagrange particle input file for the 2D JWL blast case.
# Places a cloud of particles in an annulus around the high-pressure core so the
# expanding blast sweeps a broad band of them. Writes input/lag_particles.dat,
# an 8-column formatted file: x y z vx vy vz R0 T (one particle per line).
import os

import numpy as np

N = 10000  # must be <= lag_params%nParticles_glb in case.py
R0 = 1.0e-5  # particle radius [m] (mass = 4/3 pi R0^3 * particle_pp%rho0ref_particle)
r_inner, r_outer = 0.014, 0.042  # annulus bounds [m]; stay off the [-0.05,0.05] boundaries

rng = np.random.default_rng(42)
r = np.sqrt(rng.uniform(r_inner**2, r_outer**2, N))  # area-uniform sampling
th = rng.uniform(0.0, 2.0 * np.pi, N)
x, y = r * np.cos(th), r * np.sin(th)
z = np.zeros(N)
vx = np.zeros(N)
vy = np.zeros(N)
vz = np.zeros(N)
rad = np.full(N, R0)
T = np.full(N, 300.0)

os.makedirs("input", exist_ok=True)
np.savetxt(
    "input/lag_particles.dat",
    np.column_stack([x, y, z, vx, vy, vz, rad, T]),
    fmt="%.10e",
)
print(f"Wrote {N} particles to input/lag_particles.dat")
