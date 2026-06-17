# JWL Closure Modes Benchmark

This benchmark compares the cost of the selectable JWL mixture closures on the same compact 1D JWL/air blast setup used by the registered closure goldens.

Run the baseline isobaric closure and the p-T equilibrium closure with the same build options:

```bash
./mfc.sh run benchmarks/jwl_closure_modes/case.py -- --mix-type 0
./mfc.sh run benchmarks/jwl_closure_modes/case.py -- --mix-type 2
```

Mode `2` performs a bounded scalar (bisection) solve for the products volume fraction in pressure and energy recovery, so its per-cell recovery work is strictly larger than the closed-form mode `0`.

## Reference timing

Measured on an Apple M4 Pro, single MPI rank, CPU RelDebug build (`gfortran 15.2.0`), `m = 399`, 200 steps, `riemann_solver = 2`, `weno_order = 3`. Per-step wall time from the run-time-information output:

| Mode | Closure | s/step (3 runs) |
| ---- | ------- | --------------- |
| 0    | isobaric            | 4.23e-3, 4.63e-3, 4.58e-3 |
| 2    | p-T equilibrium     | 4.41e-3, 4.13e-3, 4.11e-3 |

On this compact 1D setup the two modes are within run-to-run noise: the bounded root-find converges in a few iterations per cell and its cost is amortized below the per-step variance dominated by reconstruction and the Riemann solve. This does **not** mean the bisection is free — it adds a per-cell iterative solve in the pressure/energy recovery path that scales with iteration count, so larger grids, stiffer states, or tighter tolerances can surface a measurable difference. The benchmark exists to make that cost observable on the target machine, not to claim it is negligible.
