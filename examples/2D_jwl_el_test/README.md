# 2D JWL products blast with Euler-Lagrange particles

A single-material JWL (PETN products) blast coupled to the Euler-Lagrange particle
solver. A high-pressure circular core (30 GPa) expands into lower-pressure products
(10 GPa) and sweeps a cloud of 10,000 Lagrangian particles (two-way coupled: drag,
added mass, pressure force, collisions).

This case exercises the JWL Mie-Grüneisen products EOS (`fluid_pp%eos = 2`) together
with `particles_lagrange`. It is sized to run on 8 CPU ranks on a laptop:

- `t_stop = 2e-5` finishes in ~2 min (1165 steps)
- `t_stop = 1e-4` finishes in ~11 min (5343 steps)

both well under 15 minutes on an 8-core machine (400x400 grid, adaptive CFL).

## Running

```bash
# 1. Generate the particle input file (writes input/lag_particles.dat)
cd examples/2D_jwl_el_test && python3 generate_particles.py && cd -

# 2. Run pre_process + simulation + post_process on 8 ranks
./mfc.sh run examples/2D_jwl_el_test/case.py -n 8
```

## Notes

- Boundaries use extrapolation (`-3`). Characteristic/CBC boundaries (`-6`) are
  prohibited with `eos_jwl` because they bypass the JWL closure.
- `particles_lagrange` requires at least 2D (`n > 0`) and reads an 8-column
  formatted file `input/lag_particles.dat`: `x y z vx vy vz R0 T`.
- The JWL parameters are PETN products (A, B, R1, R2, omega, rho0, Cv, Q). `jwl_Q`
  is the detonation energy per unit mass; MFC derives `jwl_E0 = rho0 * Q` internally.
