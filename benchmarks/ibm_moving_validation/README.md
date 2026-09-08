# Moving IBM and contact validation

This benchmark checks fixes for the six-particle Mach-3 Daoud failure and the
particle time integration error in a binary impact.
Run commands from the repository root. Use separate output directories for
comparisons; MFC writes checkpoints beside the input file.

See [the HPC test plan](HPC_TEST_PLAN.md) for staged multi-GPU correctness,
convergence, physics-validation, and scaling tests with acceptance criteria.

## Solver changes

- Particle velocity, position, and orientation use the same current-stage/step-start
  RK weights as the fluid. Position and orientation use the incoming stage velocity.
- Replacing moving-solid momentum preserves internal energy, including the JWL
  reactant/product reference offset. No new EOS approximation is introduced.
- Fresh-cell interpolation uses physical distances on stretched meshes. When
  the immediate neighborhood has no established-fluid donor, the search expands
  within the allocated halo and reports an error if that search also fails.
- IBM markers cover the full allocated halo. Image interpolation excludes solid
  donors, including the exact-cell-center case. An all-solid local stencil uses
  the nearest fluid donor in the bounded halo search; if none exists, it aborts.
- Contacts use geometry and minimum-image distances, with ownership by the
  smaller global particle ID's rank. The scan is over local-neighborhood pairs,
  O(N_local^2), and uses no marker-pair list or per-step allocation. This favors
  correctness for the current 698-sphere setup; a spatial-bin acceleration can
  be benchmarked separately for much larger local particle counts.
- Adaptive dt is capped by collision_time/collision_steps_per_contact (default
  20), plus a quarter-cell bound on translation and sphere/circle surface speed.
  Fixed-dt runs must choose adequate contact and motion resolution explicitly.
- Automatic MPI neighborhood sizing uses two bounding radii and the smallest
  directional rank width. It remains a startup coverage estimate.
- Packing supports per-axis periodic overrides and uses equal-width periodic
  hash bins so a short remainder bin cannot hide a near-seam neighbor.

The linear spring/dashpot law, restitution, friction coefficient, and physical
contact duration are retained. No contact subcycling or new lubrication law is
introduced. The nearest-donor fallback is first order; spatial/time convergence
and GPU/distributed production performance remain separate validation tasks.

## Full binary impact

Copy `impact.py` to a separate run directory before executing:

```sh
OMP_NUM_THREADS=1 ./mfc.sh run PATH/impact.py -t pre_process simulation --no-gpu --no-build -n 2
python3 benchmarks/ibm_moving_validation/check_outputs.py PATH --impact --stop 0.025
```

The low-density gas makes hydrodynamic impulse negligible relative to particle
inertia. The initial relative speed is 1, restitution target 0.9, contact time
0.01, and dt 0.00025 (40 steps/contact). Refine dt to check restitution and peak
overlap convergence; this does not substitute for a fluid-loaded impact study.

## Daoud reproduction

`daoud.json` freezes the original six sphere centers from save 0, and the
original gas, grid, boundary, and collision inputs. Copy it to a separate run
directory and run on eight ranks:

```sh
OMP_NUM_THREADS=1 ./mfc.sh run PATH/daoud.json -t pre_process simulation --no-gpu --no-build -n 8
python3 benchmarks/ibm_moving_validation/check_outputs.py PATH --stop 8e-7
```

This preserves the failed run's particle realization while testing dynamics.
The checkpoint checker reads the saved grid directly and reports fluid and
solid validity separately. It is specific to these single-ideal-gas Cartesian
cases, with periodic y/z. Saved overlap extrema are sampled, not stage maxima.

The production case remains a separate 1.40-billion-cell GPU validation problem.

## Case optimization

Build and run each comparison in a separate directory using the same case:

```sh
./mfc.sh build -t simulation --no-gpu --case-optimization --input PATH/daoud.json -j 8
OMP_NUM_THREADS=1 ./mfc.sh run PATH/daoud.json -t simulation --case-optimization --no-gpu --no-build -n 8
```

Generate initial data with `pre_process` first, or copy the same restart into
both directories. Compare checkpoint times as well as values: the final step
can be shortened to reach `t_stop`, and a final output can overwrite an earlier
output with the same integer `mytime/t_save` index.

## CPU validation (2026-09-07)

With gfortran and MPI in double precision:

- The Daoud reproduction reached 0.80 microseconds on eight ranks, resuming
  from the 0.44-microsecond checkpoint after fixing the donor search. All 41
  saved states were finite with positive fluid pressure/density; maximum sampled
  overlap was 0.001318 diameters. This was a restarted run, not a fresh continuous
  run of the final source.
- Binary-impact restitution was 0.901522 at 40 steps/contact and 0.900782 at
  80 steps/contact, approaching the specified 0.9; net x momentum was zero.
  Before the RK correction, the 40-step result was 0.850371.
- A 60-particle packing check produced identical initial states on one and two
  ranks. Exhaustive y/z-periodic pair distances were at least 0.115113 for a
  required separation of 0.115.
- At the same restarted Daoud checkpoint, case optimization differed from the
  ordinary build by at most 1.10e-11 in the fluid state and 1.01e-12 in particle
  state, using `abs(optimized - ordinary)/max(1, abs(ordinary))`.
- The four JWL golden regressions passed through all three executables.
- All 56 IBM golden regressions passed without updating reference files.

GPU execution, the other supported compilers, and the production mesh have
not been validated locally.
