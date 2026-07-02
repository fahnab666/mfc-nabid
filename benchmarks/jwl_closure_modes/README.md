# JWL Closure Cost Benchmark

This benchmark measures the per-step cost of the JWL mixture closure on the same compact 1D JWL/air blast setup used by the registered closure golden test.

Run the Rocflu state-interpolated closure (the only supported closure):

```bash
./mfc.sh run benchmarks/jwl_closure_modes/case.py
```

The closure recovers pressure, temperature, and sound speed from `(rho, e, Y)` with a closed-form pressure law and an analytic energy inverse, so the per-cell recovery work is fixed (no iterative solve). The benchmark exists to make that cost observable on the target machine.

## Reference timing

Measured on an Apple M4 Pro, single MPI rank, CPU RelDebug build (`gfortran 15.2.0`), `m = 399`, 200 steps, `riemann_solver = 2`, `weno_order = 3`. Per-step wall time from the run-time-information output:

| Closure | s/step (3 runs) |
| ------- | --------------- |
| Rocflu state-interpolated | 4.41e-3, 4.13e-3, 4.11e-3 |

On this compact 1D setup the closure cost is amortized below the per-step variance dominated by reconstruction and the Riemann solve.
