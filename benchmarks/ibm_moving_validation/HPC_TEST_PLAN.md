# Multi-GPU particle-shock validation and improvement plan

This plan qualifies moving immersed boundaries, particle contacts, and
particle-shock calculations before larger HPC production runs. Start from
commit `3713317e` or a documented descendant. The first priority is a fresh,
continuous run of the final solver: the successful local Daoud result used a
restart at approximately 0.44 microseconds.

The thresholds below are proposed engineering targets, not established solver
accuracy or universal DNS criteria. Preserve existing golden-test tolerances.
Advance only after the preceding stage passes. Do not refresh golden files to
hide an unexplained difference.

## Available cases and current evidence

- [impact.py](impact.py): nearly dry, two-circle impact; contact time 0.01,
  restitution target 0.9, default timestep 0.00025 (40 steps/contact).
- [daoud.json](daoud.json): six fixed initial sphere locations, Mach-3 flow,
  viscous no-slip IBM, adaptive timestepping, final time 0.8 microseconds.
- [check_outputs.py](check_outputs.py): checkpoint analysis for these
  single-ideal-gas Cartesian cases. It is not a JWL pressure checker.
- [README.md](README.md): CPU validation results and implementation limitations.

Additional stress cases and diagnostics listed below still need implementation.
The older `benchmarks/1D_jwl_mixture_closure_validation/README.md` describes the
previous weighted closure; its reported results must not be treated as validation
of the current pressure-temperature equilibrium implementation.

## EOS integration gate

The upstream integration includes the state-dependent EOS framework from
[PR #1811](https://github.com/MFlowCode/MFC/pull/1811), alongside this fork's
existing equilibrium JWL closure. These are different mixture models:

| Selector | Behavior |
|---|---|
| `jwl_pt` (2) | Existing pressure-temperature equilibrium closure, program burn, afterburn and JWL reaction progress |
| `jwl` (4) | Per-phase JWL coefficients, frozen mixture acoustics, five- and six-equation models |
| `mie_gruneisen` (3), `vinet` (5) | PR #1811 reference-curve families |
| `ideal_gas` (6) | Zero-stiffness gas; use the name when importing upstream cases |

Numeric `2` remains the existing JWL selector in this fork. Do not copy upstream
numeric ideal-gas selectors without translating them. The `jwl_pt` and per-phase
families cannot be mixed in one case. Use `jwl_wrt` for equilibrium temperature
and `T_wrt` for per-phase temperature. Reaction parameters also belong to their
respective models; switching the EOS name alone does not convert a burn case.

The existing capabilities retain their established locations:

| Location | Connection |
|---|---|
| `src/common/m_jwl.fpp` | Equilibrium pressure, temperature, energy and sound speed; initialized/finalized by `m_variables_conversion` for all three executables |
| `src/common/m_variables_conversion.fpp` | Primitive/conservative conversion and dispatch to the equilibrium closure or PR reference-curve helpers |
| `src/simulation/m_jwl_sources.fpp` | Existing program-burn, afterburn and reaction sources, called from `m_rhs` |
| `src/simulation/m_riemann_solver_{hll,hllc,lf}.fpp` | Equilibrium face energy/sound-speed calls; reactive progress uses HLLC |
| `src/simulation/m_ibm.fpp` | Equilibrium ghost/fresh-cell energy reconstruction, including reaction progress; per-phase cells use the shared EOS helpers |
| `src/post_process/m_start_up.fpp` | Existing `jwl_wrt` equilibrium diagnostics alongside per-phase `T_wrt` |
| `toolchain/mfc/params/definitions.py`, `case_validator.py`, `case.py` | Parameter bindings, model compatibility and case-optimization selection |

The obsolete inline Riemann include is removed; its equilibrium calls are wired
directly into the refactored solvers. Existing uppercase JWL coefficient names
remain accepted alongside the PR's lowercase names.

On each compiler/backend, run both JWL paths before particle-shock production:

```sh
./mfc.sh test --only Kernel JWL --no-examples --gpu acc -a -j 1
./mfc.sh test --only eos=jwl --no-examples --gpu acc -a -j 1
./mfc.sh test --only eos=mie_gruneisen --no-examples --gpu acc -a -j 1
./mfc.sh test --only eos=vinet --no-examples --gpu acc -a -j 1
./mfc.sh test --only 'Reactive Burn' --no-examples --gpu acc -a -j 1
```

Repeat with `--case-optimization`, and use `--gpu mp` for the OpenMP target lane.
Also run the registered JWL isentropic-release and Mie-Gruneisen acoustic/impact
convergence tests (`./mfc.sh test -l` lists their IDs). Compare each closure with
its own reference results; equality between the two mixture models is not an
acceptance criterion. GPU results and DNS qualification remain pending until
these gates and the particle/mesh convergence stages below pass.

## 1. Freeze the configuration

Begin in double precision, with fast math and case optimization disabled.
Record the following for every run:

| Category | Required record |
|---|---|
| Source | Commit, dirty/clean status, complete input file |
| Software | Compiler, MPI, GPU backend, module versions |
| Hardware | GPU model, nodes, ranks per node, GPU binding |
| Geometry | Particle realization, seed, dimensions, concentration |
| Numerics | Mesh, decomposition, actual timestep history, output interval |
| Results | Physical stop time, exit status, diagnostic summary |

Use one MPI rank per scheduler-visible GPU device initially. Follow the site's
device convention when GPUs expose independently scheduled partitions.
`-n` is ranks per node; `-N` is nodes. Choose the module slug and scheduler
template from `toolchain/modules` and `toolchain/templates`. Confirm that local
subdomains satisfy the solver's decomposition and halo requirements.

## 2. CPU/GPU and MPI correctness

Fresh-cell reconstruction follows the established practice of supplying a fluid
state when a moving solid uncovers a cell. The EOS consistency fix rebuilds total
and six-equation phase energies from the reconstructed pressure and densities;
equilibrium JWL also retains its reaction-dependent energy reference. Test
`D116E4F1` moves a circle through a two-fluid MG/ideal-gas field on two ranks.
Restoring the old fresh-cell block changes final phase-energy values by about
0.0029, well above the existing comparison tolerance.

[Mittal et al. (2008)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2834215/) discuss
fresh-cell treatment, and
[Brahmachary et al. (2018)](https://doi.org/10.1002/fld.4479) demonstrate
inverse-distance-based reconstruction for high-speed compressible flow.
These support the approach, not an accuracy claim for this exact stencil.
The current reconstruction does not enforce a swept-volume conservation balance.
[Seo and Mittal (2011)](https://doi.org/10.1016/j.jcp.2011.06.003) connect geometric
conservation errors to pressure oscillations. Therefore measure fluid/particle
mass, momentum and energy budgets, pressure/force oscillations, and mesh/timestep
convergence before DNS qualification. Positive interpolation weights alone do
not guarantee thermodynamic admissibility for a nonlinear mixture EOS.

Run regressions inside an allocation with the site's supported compiler/backend.
Use the ordinary build first, then repeat selected cases with case optimization.

| Test | Matrix | Measurements | Initial pass target |
|---|---|---|---|
| Existing regressions | IBM and JWL on supported GPU backend | Golden comparisons | All relevant tests pass |
| Binary impact | CPU reference; 1, 2, 4 GPUs where valid | Restitution, momentum, overlap, trajectories | Restitution within 0.005 of 0.9 at 40 steps/contact |
| Fresh Daoud | 1, 2, 4, 8 GPUs where valid | Fluid and particle histories | Continuous completion to 0.8 microseconds; no invalid states |
| Case optimization | Off/on, same case and decomposition | Matching-time fields and particles | Short-run normalized differences around 1e-8 or smaller |
| Repeatability | Three repeats of one multi-GPU case | Forces and trajectories | No systematic drift or rank-boundary jumps |

Compare identical physical times, not filenames alone. A final save can overwrite
an earlier checkpoint sharing the integer `mytime/t_save` index.
Normalize each quantity using a fixed physical reference scale, and record that
scale. Do not require bitwise equality after GPU reductions. For long chaotic
runs, compare statistics and uncertainty instead of pointwise trajectories.

Stop on changed initial packing, lost or duplicate particle IDs, duplicate/missing
contacts, rank-dependent force jumps, nonfinite states, or unexplained golden
differences.

## 3. Moving-boundary and contact stress cases

These additional cases isolate failure modes before shock-curtain runs.

| Case | Setup | Failure mode tested |
|---|---|---|
| Prescribed translation | Sphere moving through uniform flow across ranks | Fresh-cell reconstruction and ownership |
| Periodic translation | Cross faces, then edges and corners | Periodic markers and ownership |
| Partitioned pair contact | Shift the contact relative to rank boundaries | Missing or double-counted forces |
| Simultaneous contact | Symmetric line of three spheres | Both contacts detected; symmetry retained |
| Oblique impact | Tangential velocity and friction | Torque sign and rotation |
| Unequal spheres | Different radii and masses using explicit patches | Effective mass and torque assumptions |
| Opening a tight cluster | Nearly touching particles move apart | Donors in newly uncovered fluid |
| Empty particle ranks | Particles occupy only part of the domain | Zero-sized arrays and communication |

Compare one rank with multiple decompositions. For isolated contact, check
equal-and-opposite forces and total momentum. For inelastic/frictional contact,
account for spring energy and dissipation; kinetic energy should not be constant.
The current collision law assumes circles/spheres. General IBM geometry support
does not establish correct collision handling for arbitrary STL bodies.

## 4. Time and mesh convergence

Keep physical geometry, material properties, and contact duration fixed.
First refine time on one mesh. For the Daoud contact duration:

| collision_steps_per_contact | Contact timestep cap (seconds) |
|---:|---:|
| 20 | 2.47878e-10 |
| 40 | 1.23939e-10 |
| 80 | 6.19695e-11 |

These are caps, not guaranteed timesteps. Record whether contact, fluid CFL, or
body motion controls each step. For fixed-timestep impact, change `dt` directly;
the adaptive contact parameter does not refine that case.

Compare individual particle impulse, final velocity, peak overlap, contact
duration, shock arrival, pressure impulse, curtain centroid, and curtain thickness.
Target medium-to-fine differences below 1% for impulse and mean particle motion,
with decreasing differences. Inspect pressure peaks separately for sampling error.

Then use approximately 18, 27, and 40 cells per particle diameter, measured in
the particle region. Keep the same physical domain and exact particle centers.
Choose a timestep small enough that temporal error does not mask spatial error.
Target medium-to-fine changes below 2% for integral quantities and below 5% for
resolved pressure peaks, with consistent convergence. Investigate non-monotonic
results rather than forcing a Richardson extrapolation or grid convergence index.

At fixed 3D domain size, 27 to 40 cells/diameter requires about 3.25 times the
cells. Runtime can approach 4.8 times if the timestep also shrinks with cell
width; contact-limited runs may scale differently. Start with separate time and
mesh sweeps rather than an expensive full Cartesian product.

## 5. Conservation and reconstruction diagnostics

The current checkpoint checker does not provide all these diagnostics. Add
compact reductions and histories before the large campaigns.

| Diagnostic | Sampling |
|---|---|
| Minimum density, pressure, sound-speed squared | Every step |
| Maximum overlap and contact count | Every particle RK stage |
| Fresh cells requiring expanded donor search | Every stage, aggregated per rank |
| Image points using nearest-fluid fallback | Every stage, aggregated per rank |
| Particle count and unique global IDs | Startup, migration, restart |
| Fluid/particle momentum and boundary impulse | Every step or short interval |
| Energy, boundary flux, source contributions | Short interval |
| Maximum owned/neighborhood particle counts | Periodically |

Use scalar reductions rather than full fields every stage. Compute conservation
over the physical fluid region, accounting for boundary fluxes and moving-boundary
effects. Including artificial solid-interior states in a whole-array sum does not
give a physical conservation budget.

For JWL/reactive cases, use the active EOS and include energy-reference and
reaction-source accounting. Do not apply the ideal-gas checkpoint pressure formula.
A provisional conservation target is residual below 0.1% of the relevant
transported quantity, decreasing under refinement. Persistent residuals correlated
with donor fallback would prioritize reconstruction/conservation improvements.

## 6. Validate shock and particle physics

| Configuration | Primary quantities |
|---|---|
| Shock without particles | Shock speed, jump states, pressure history |
| Shock over one fixed sphere | Surface pressure, drag and impulse |
| Shock accelerating one free sphere | Acceleration, velocity, force impulse |
| Two spheres at varying separation | Shielding, gap flow, interaction forces |
| Small fixed curtain | Reflected/transmitted shocks and pressure impulse |
| Small moving curtain | Expansion, particle velocity distribution, momentum transfer |
| Full curtain | Experimental observables and uncertainty |

Match Mach number, particle diameter/density, curtain thickness, volume fraction,
lateral coverage, and boundary conditions to the chosen reference experiment.
Matching diameter alone is insufficient. Sandia's 115-micrometer glass-particle
curtain measurements are a candidate reference, not an already matched case.

Validate viscous near-contact behavior separately. Published interface-resolved
collision work includes lubrication corrections for gaps smaller than the mesh.
Test whether such a correction is needed in the intended regime before adding it.

## 7. Particle count and GPU scaling

Progress through 6, 60, 200, 698, then several thousand particles, subject to
memory and local-count limits. Separate increasing domain size at fixed
concentration from increasing concentration in a fixed domain.

The current `src/common/m_constants.fpp` defines `num_local_ibs_max = 2000` and
`num_ib_patches_max_namelist = 54000`. These are implementation limits to audit,
not demonstrated operating capacities. Log both owned and neighborhood counts;
global particle count alone does not establish that local storage is adequate.

| Study | Hold fixed | Vary |
|---|---|---|
| Strong scaling | Global mesh, particles, physical duration | GPU count |
| Weak scaling | Cells and particle loading per GPU | Domain size and GPU count |
| Contact scaling | Mesh and decomposition in a controlled setup | Local particle count |
| I/O scaling | Simulation workload | Output frequency |

Start with 1, 2, and 4 nodes, then expand. The six-particle case is a correctness
test, not a representative scaling benchmark. Measure a representative steady
timing window after startup and report checkpoint cost separately.

Record median timestep time, fluid update time, IBM reconstruction time, contact
time, MPI communication/wait time, per-rank imbalance, peak device memory, and
checkpoint time. Repeat performance measurements three times.

For fixed work, strong-scaling efficiency relative to P0 GPUs is
`E(P) = T(P0) * P0 / (T(P) * P)`. Investigate efficiency below 70%, rank-time
imbalance above 20%, or contact work exceeding 10-15% of timestep time. These are
profiling triggers, not correctness failures. If contact cost grows quadratically
with local particle count, spatial bins become a justified optimization.

## 8. Restarts, uncertainty, and DNS qualification

Compare continuous runs with restarts before shock arrival, during acceleration,
during contact, and near an ownership transfer. First retain the decomposition;
then test changed GPU counts as a separate restart-portability experiment.
Check IDs, positions, velocities, force continuity, and physical time.

Adaptive restart time currently starts from `n_start*t_save`, although the saved
state may occur slightly after that nominal time. Measure this timing difference
explicitly before attributing all restart differences to the physical solver.

Before calling production results DNS, demonstrate resolution of post-shock
boundary layers, wakes, and relevant turbulent scales as well as convergence of
integral quantities. Cells per diameter alone are insufficient. Report local
post-shock Reynolds numbers and the smallest resolved flow scales.

For random curtains, start with three realizations to estimate variability, then
increase the ensemble until confidence intervals for the selected observables
are sufficiently narrow.

## Launch procedure

Run all project commands from the repository root. Replace the uppercase
placeholders below. `acc` is the NVIDIA/OpenACC example; select `mp` where the
site/compiler supports OpenMP offload. Do not assume the same module slug,
template, or device count across systems.

```sh
source ./mfc.sh load -c SYSTEM_SLUG -m g
./mfc.sh precheck -j 8

./mfc.sh run RUN/daoud.json \
  -t pre_process simulation \
  --gpu acc \
  -e batch -c SYSTEM_TEMPLATE \
  -N 2 -n 4 \
  -a ACCOUNT -w 01:00:00 \
  --case-optimization --dry-run
```

Inspect the generated job and GPU binding, then remove `--dry-run` to submit.
Omit `--case-optimization` for the initial baseline; enable it for the paired
comparison. Use a separate run directory for every variant. The normal run
command builds the required configuration; use `--no-build` only after building
that exact configuration. A dry run is a submission check, not simulation evidence.

Inside an appropriately sized allocation, the regression commands are:

```sh
./mfc.sh test --only IBM --no-examples --gpu acc -j 1
./mfc.sh test --only JWL --no-examples --gpu acc -a -j 1
```

Confirm launcher and device visibility first. `-j` controls concurrent jobs, not
GPU allocation. Keep it at one initially so simultaneous tests do not compete for
the same device. Inspect `./mfc.sh test --help` and the site's launch setup before
increasing concurrency.

## First allocation and improvement decisions

Prioritize stages 1-4 on small cases. The most valuable early results are fresh-run
stability, decomposition independence, timestep convergence, and donor statistics.

| Observation | Next focused improvement |
|---|---|
| Rank-dependent contacts or trajectory jumps | Ownership, periodicity, MPI exchange |
| Conservation residual follows donor fallback | Reconstruction and moving-boundary budgets |
| Contact outcomes change with timestep | Contact resolution; assess subcycling only after convergence |
| Mesh-converged near-contact physics disagrees with reference | Contact/lubrication model assessment |
| Contact runtime grows quadratically | Spatial binning with identical pair-force checks |
| Large rank imbalance | Particle-aware decomposition/communication profiling |
| Optimized and ordinary builds disagree | Case-specialization and device-state audit |

## References

- [NASA: examining spatial grid convergence](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html).
- [Costa et al. (2015): collision model for fully resolved particle-laden flows](https://arxiv.org/abs/1506.01880).
- [Sandia: flash X-ray measurements of shock-induced particle-curtain dispersal](https://www.sandia.gov/research/publications/details/flash-x-ray-measurements-on-the-shock-induced-dispersal-of-a-dense-particle-2015-12-01/).
- [MFC getting started](https://mflowcode.github.io/documentation/getting-started.html).
  For this fork, prefer local `./mfc.sh --help`, module definitions, and templates
  when they differ from upstream documentation.
