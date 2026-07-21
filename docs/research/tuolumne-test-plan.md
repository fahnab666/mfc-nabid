# Tuolumne Test and Improvement Plan for the EOS PR Ladder (PR1 to PR5)

Target machine: LLNL Tuolumne (Cray EX, AMD MI300A APUs, gfx942, Flux scheduler).
MFC already knows the machine: slug `tuo` in `toolchain/modules` (cpe/25.03 + rocm/6.3.1,
`craype-accel-amd-gfx942`, `HSA_XNACK=0`) and batch template `toolchain/templates/tuo.mako`.
Branch under test: `feature/eos-combined` (PR1 interface + PR2 selector + PR3 MG module +
PR4a JWL + PR5 material files, plus the review-fix and performance commits).

Everything below has already passed locally on macOS gfortran CPU: 618/618 golden tests
bit-identical, both JWL examples stable, 7/7 bench cases at grind parity with master.
Tuolumne closes the four lanes a laptop cannot: Cray ftn, AMD flang, GPU offload on MI300A,
and multi-node MPI scaling. These are exactly the maintainer's open questions.

## Phase 0: environment and baselines

1. Load modules per mode before every build (the source is required):
   `source ./mfc.sh load -c tuo -m c` for CPU work, `-m g` for GPU work.
2. Clone twice or use two build dirs: once at the merge base
   (master `0857ace`) for baselines, once at `feature/eos-combined`. Never run two
   `./mfc.sh` commands concurrently in one tree; it corrupts the build staging.
3. Record `build/lock.yaml` for every configuration before trusting any timing.
   A stale debug flag in the lock invalidated one local benchmark this week.

## Phase 1: build matrix (correct compilation is itself a test)

Build all three targets in each configuration, on both trees:

| # | Compiler | Flags | What it gates |
|---|----------|-------|---------------|
| 1 | Cray ftn (cpe) | (default CPU) | Cray CPU correctness |
| 2 | Cray ftn | `--debug` | bounds and IEEE checks on new code |
| 3 | Cray ftn | `--gpu mp` | OpenMP offload, the primary GPU lane |
| 4 | amdflang | `--gpu mp` | AMD offload; validates the cross-TU `mixture_closure` GPU_UPDATE fix |
| 5 | Cray ftn | `--single` | single precision discipline |
| 6 | Cray ftn | `--mixed` | validates the `f_bulk_modulus` stp/wp fixes; this exact configuration failed to compile before commit 362793d2 |

Configuration 6 is the highest-value quick win: it is a compile-only check that the two
mixed-precision review findings are truly closed on compilers other than gfortran. The
`cray_inline` directives on `f_bulk_modulus`, `f_validate_state`, `s_mg_jwl_reference`, and
`s_accumulate_mixture_properties` are Cray-only and have never been exercised by a Cray
compiler; watch build 1 and 3 logs for INLINEALWAYS diagnostics.

## Phase 2: correctness

1. Full golden suite, CPU, Cray ftn: `./mfc.sh test -j 32`. Expect 618/618. Any diff on a
   pre-existing test is a bug in the stack, not noise; do not regenerate goldens.
2. Full golden suite, GPU (inside a Flux allocation):
   `./mfc.sh test -j 4 --gpu mp -g 0 1 2 3`. The three JWL tests (976FEEC6 and the two
   example tests) run the mg_mixture path on device for the first time; to smoke them
   first: `./mfc.sh test --gpu mp -g 0 -o 976FEEC6 8F3DAFFE 5B8C364F`.
3. CPU versus GPU agreement on the two JWL examples: run
   `examples/1D_jwl_shocktube` and `examples/1D_jwl_products_expansion` in both modes and
   diff the final D data to tolerance. This is the direct check that the JWL accumulator,
   the device `exp()` calls, and the `GPU_UPDATE` set behave identically on gfx942.
4. AMD flang staleness regression check (finding from the review workflow): on build 4,
   run the JWL shocktube and confirm finite, nonzero pressure output. Before commit
   09fa0a17 a stale device copy of `mixture_closure` could silently zero the mixture
   state on exactly this compiler.
5. Vanished-phase robustness: run the products-expansion example with `mpp_lim = T`
   added to the case. Before 09fa0a17 this could NaN at the material interface. Expect a
   clean run with finite output.
6. Toolchain unit tests on the login node:
   `python -m pytest mfc/test_mg_eos.py mfc/test_materials.py -q` from `toolchain/`.

## Phase 3: performance (the maintainer's question, answered on his hardware class)

1. CPU bench, both trees, identical locks:
   `./mfc.sh bench --mem 1 -j 32 -o base.yaml` then the same on the combined tree;
   compare with `./mfc.sh bench_diff base.yaml combined.yaml`. Expect parity; the one
   real CPU regression (+8.3 percent HLLC from the accumulator dispatch) was found and
   fixed locally (commit 07fc90ac), and this run confirms the fix holds under Cray ftn,
   whose inliner differs from gfortran.
2. GPU bench, both trees, single MI300A: same bench pair under `-m g` with `--gpu mp`
   builds. This is the first GPU timing of the stack anywhere. The stiffened hot path
   should be identical to master; the only new device code behind `if (mg_mixture)` is
   never taken in the bench cases.
3. Weak scaling: HLLC bench case at 1, 2, 4, and 8 ranks (up to 2 nodes x 4 APUs) with
   fixed memory per rank, submitted through Flux:
   `./mfc.sh run benchmarks/5eq_rk3_weno3_hllc/case.py -e batch -c tuo -N <nodes> -n <ranks_per_node> -w 00:15:00 -a <bank> -t pre_process simulation --output-summary <artifact>.yaml`.
   Compare parallel efficiency base versus combined. The stack changes no MPI code, so
   any efficiency delta indicates a problem.
4. JWL-specific grind: time the JWL shocktube at production resolution on one APU and
   report the grind-time ratio of a JWL cell versus a stiffened cell (expect a modest
   constant factor from the two `exp()` evaluations per cell per stage; record it so the
   PR body can state the cost of enabling the feature honestly).
5. Record everything in yaml artifacts and keep the lock files; the PR1 reopen argument
   is strongest as a table of same-lock grind numbers on the reviewers' own hardware class.

## Phase 4: what to report back into the PRs

1. PR1 (#1663): the CPU bench table already posted, extended with Cray CPU, Cray GPU, and
   AMD GPU columns from Phase 3, plus the `--mixed` compile fix. This directly answers
   "is it fast, does it scale, on hardware I trust".
2. PR2/PR3: build matrix results only; they have no hot-path exposure.
3. PR4a/PR5: the CPU versus GPU agreement runs and the JWL cost factor.

## Phase 5: maintainability, mergeability, and scalability audits

These need no GPU time; run them anywhere.

Maintainability
1. `./mfc.sh precheck` must stay green on every commit of the stack (already enforced
   locally by the hook; re-run once on Tuolumne to catch platform-dependent lint).
2. Extension-point drill: add a dummy EOS family end to end (enum value, one `case` arm
   in `s_mg_mixture_variables`, validator rejection) and count the touched files. The
   target is three; more means a dispatch site has leaked outside the sanctioned points
   (`select case (eos_fl(i))` and the `@:accumulate_mixture` macro). Revert the drill
   commit afterwards.
3. Grep audit that no solver file computes EOS algebra outside `m_thermodynamics` and
   `m_eos_mie_gruneisen`, and that `mg_mixture` dispatch appears only via the macro and
   the two conversion sites.
4. Doc freshness: `lint_docs` runs inside precheck; confirm `case.md`, `contributing.md`,
   and this plan still reference existing symbols after any rebase.

Mergeability
1. Rebase the ladder onto current upstream master (`MFlowCode/MFC`) in a scratch clone
   and record the conflict surface per PR. The ladder is only submittable if PR1 rebases
   clean; later rungs may conflict with each other by design (they stack).
2. Per-rung independence check: build and run the golden suite at each rung boundary
   (PR1 tip, PR1+2, PR1+2+3, full stack) so every PR is individually green and the
   maintainer can merge them one at a time. The rung tips exist as local branches;
   re-verify after any rebase.
3. Diff budget per rung: report net LOC and files touched per PR against the request in
   each PR body. Anything that grew during the fix cycle gets rechecked against scope.
4. Golden-file discipline: the stack must add goldens only for new JWL tests and change
   zero pre-existing goldens; `git diff master --stat tests/` is the one-line check.

Scalability (beyond the Phase 3 node scaling)
1. Feature scaling: the extension-point drill above is the code-scalability test; the
   Phase 3 JWL cell-cost factor is the physics-cost test.
2. Fluid-count scaling: run the JWL mixture test at num_fluids = 2 and 3 (JWL + two
   ambients) and confirm grind scales linearly with the accumulator loop, no worse.
3. Problem-size scaling on one APU: the JWL shocktube at 4x and 16x resolution;
   grind per cell should be flat, showing the JWL path is compute-bound like the rest
   of the solver rather than memory-pathological.

## Improvement plan by PR

PR1 (thermodynamics interface)
1. Done this cycle: `cray_inline` on `f_bulk_modulus` and `f_validate_state`; mixed
   precision argument kinds; NaN branch parity in `f_validate_state`.
2. Remaining: fold the IBM and Lagrangian-bubble coupling exceptions behind the interface
   once a bit-identical routing exists (documented exceptions today). Low priority until
   a reviewer asks.

PR2 (EOS selector)
1. Known gap: extend the `s_check_eos` loop bound to the bubbles_euler extra fluid slot.
2. Improvement: readable case-file names for every enum value, with validator tests per
   rejected (eos, feature) pair.

PR3 (Mie-Gruneisen module)
1. The generic MG family is currently exercised only through JWL. Add
   `eos_mie_gruneisen` as a selectable family: one new `case` arm in
   `s_mg_mixture_variables` calling `s_mg_stiffened_reference` style providers, plus a
   golden test showing MG-with-linear-curve equals stiffened gas to machine precision at
   runtime, not just in the unit test.

PR4a (JWL)
1. The `select case (eos_fl(i))` and the `@:accumulate_mixture` macro are the extension
   points; keep them the only dispatch sites.
2. Implement the reserved closures (`composition_weighted`, `modified_composition_weighted`,
   `kuhl`) as new `mixture_closure` values; each is a separate small PR on the ladder.
3. Case-optimization: under `--case-optimization`, `mg_mixture` is a known constant at
   build time; baking it removes even the scalar branch from the flux kernels. Worth a
   line in the PR body, not code yet.
4. Water-scale ambients need the N-constituent pressure-equilibrium closure; documented
   gap, out of scope for the ladder.

PR5 (material files)
1. Add `$MFC_MATERIALS_DIR` documentation and a curated open-literature catalog repo
   (external, per the export-control ledger; the repo itself stays synthetic-only).
2. Extend `EOS_FAMILIES` with `mie_gruneisen` when PR3's family goes live, and later a
   `table` family whose data loads through the same provenance-gated mechanism.

## Execution plan for 1-hour sessions, 2 nodes, 4 APUs per node

All runs go through the APU nodes: CPU-mode runs use the Zen 4 cores of the same MI300A
nodes, GPU-mode runs use the accelerators. Only compilation and post-processing happen on
the login nodes. Prebuild everything before each session so the hour starts computing.

On the login node before any session (no allocation needed):
1. Phase 0 setup, then all six builds of the matrix on both trees. The `--mixed` and
   `--debug` configurations are compile-only gates and are finished the moment they build.
2. Toolchain pytest.
3. Prebuild the GPU binaries (`-m g`, `--gpu mp`, Cray ftn and amdflang) so no session
   burns time compiling.

Session 0 (1 node, CPU mode): the CPU golden suite (618 tests) under Cray ftn with high
`-j` across the node's cores; it takes about ten minutes locally at `-j 8`, so one hour
on a Tuolumne node covers it with room for the Phase 5 rung-boundary suite repeats.

Session 1 (1 node is enough, correctness): GPU golden subset including the three JWL
tests; both JWL examples on one APU; dump final D data. If time remains, the amdflang
binary repeat of the JWL shocktube (stale-selector check) and the `mpp_lim` run. All are
small 1D cases; this fits in well under an hour.

Session 2 (1 node, single-APU performance): GPU bench pair, base then combined, one APU,
identical locks. The seven cases at `--mem 1` are sized for a single device; interleave
base and combined per case if the hour looks tight so partial results are still paired.
Compare offline with `bench_diff` after the session.

Session 3 (2 nodes, scaling): HLLC bench case at 1, 2, 4, and 8 ranks (up to the full
2 nodes x 4 APUs) with fixed memory per rank, base and combined. Eight short runs; submit
them as one Flux script so scheduling gaps do not eat the hour. This replaces the
multi-day weak-scaling ladder; 8 APUs across 2 nodes already answers "does it scale" at
the granularity the PR discussion needs, and the JWL cell-cost timing (Phase 3 item 4)
fits at the end of this session on one rank.

CPU versus GPU agreement diffs, bench_diff tables, and the Phase 4 report are all offline
post-processing of session artifacts. Any red result stops the line for that lane; report
it with the log rather than working around it.

## Appendix A: artifact convention

One directory per session under `$HOME/eos_ladder_runs/`, named
`s<session>_<tree>_<mode>/`, e.g. `s2_combined_gpu/`. Every run writes an
`--output-summary` yaml into it, and every session directory gets a copy of the tree's
`build/lock.yaml` and the `git rev-parse HEAD` it ran. A result without its lock and SHA
is not evidence; the lock file is what caught an invalid benchmark locally.

## Appendix B: session command blocks

Login node, once per tree (base at 0857ace, combined at feature/eos-combined):

    source ./mfc.sh load -c tuo -m c
    ./mfc.sh build -j 16                          # CPU Release
    ./mfc.sh build -j 16 --debug                  # compile gate
    ./mfc.sh build -j 16 --mixed                  # compile gate (fix 362793d2)
    ./mfc.sh build -j 16 --single                 # compile gate
    source ./mfc.sh load -c tuo -m g
    ./mfc.sh build -j 16 --gpu mp                 # Cray ftn offload
    (amdflang build per toolchain/modules AMD entries, if configured for tuo)

Session 0, CPU suite (1 node, interactive allocation):

    flux alloc -N 1 -t 60m
    source ./mfc.sh load -c tuo -m c
    ./mfc.sh test -j 32 2>&1 | tee $ART/s0_suite.log     # expect 618 passed, 0 failed

Session 1, GPU correctness (1 node):

    flux alloc -N 1 -t 60m
    source ./mfc.sh load -c tuo -m g
    ./mfc.sh test --gpu mp -g 0 -o 976FEEC6 8F3DAFFE 5B8C364F   # JWL on device first
    ./mfc.sh test -j 4 --gpu mp -g 0 1 2 3                       # then the full GPU set
    ./mfc.sh run examples/1D_jwl_shocktube/case.py -t pre_process simulation
    ./mfc.sh run examples/1D_jwl_products_expansion/case.py -t pre_process simulation
    (repeat the shocktube with the amdflang binary, then with mpp_lim = T in the case)

Session 2, single-APU bench (1 node, run per tree, interleave if tight):

    ./mfc.sh bench -m 1 -j 8 -o $ART/s2_<tree>_gpu.yaml
    # offline: ./mfc.sh bench_diff s2_base_gpu.yaml s2_combined_gpu.yaml

Session 3, scaling (2 nodes): submit all eight runs (1, 2, 4, 8 ranks x both trees) as
batch jobs at session start so they pack the hour:

    ./mfc.sh run benchmarks/5eq_rk3_weno3_hllc/case.py -e batch -c tuo \
        -N <nodes> -n <ranks_per_node> -w 00:12:00 -a <bank> \
        -t pre_process simulation --output-summary $ART/s3_<tree>_r<ranks>.yaml

## Appendix C: pass and fail criteria

| Check | Pass | Action on fail |
|---|---|---|
| Compile matrix (6 configs x 2 trees) | all build | capture the full build log; the config name is the finding |
| CPU suite | 618 passed, 0 failed | do not regenerate goldens; diff is the bug report |
| GPU suite | same pass set as CPU | isolate the first failing UUID with `-o <uuid>` |
| CPU vs GPU JWL examples | final D data agree to golden tolerance | dump both, bisect field by field |
| amdflang JWL shocktube | finite nonzero pressure everywhere | suspect the mixture_closure device copy first (commit 09fa0a17) |
| mpp_lim products run | completes, finite output | suspect the vanished-phase floor (commit 09fa0a17) |
| GPU bench pair | every case within 5 percent, no one-sided trend | rerun the outlier alternating, as done for hypo_hll and HLLC locally |
| Scaling | efficiency delta base vs combined under 2 percent at 8 ranks | check GPU_UPDATE placement around MPI in the diff |

## Appendix D: failure triage kit

Capture before leaving any failed session: the run directory, `build/lock.yaml`, the git
SHA, `module list` output, and for GPU failures `rocm-smi` and the value of `HSA_XNACK`.
On Cray, rebuild the failing target with `--debug` and rerun once; the bounds checker
turns most silent wrong answers into a named line number. If a GPU result is wrong but
CPU is right on the same binary tree, the first suspects in this stack are, in order:
device copies of `mixture_closure`, `eos_fl`, and the `jwl_*` arrays (all updated at
m_variables_conversion line ~350), then the `@:accumulate_mixture` dispatch sites.

## Appendix E: final audit findings and platform watchpoints

A four-dimension adversarially verified audit of the final branch state (macro rewrite,
select case dispatch, GPU data coverage, checker gates) produced the following. The two
code defects and the checker gap are fixed on the branch; the watchpoints are Tuolumne
work items and are folded into the session checklists below.

Fixed on the branch before HPC testing
1. `acoustic_source` was not prohibited for JWL fluids; the source builds Tait
   coefficients from the stiffened `gammas`/`pi_infs` arrays, which a JWL fluid leaves
   at the -1e6 default, injecting NaN sources past every checker. Now prohibited in
   `src/simulation/m_checker.fpp` and `case_validator.py`.
2. The probe (`probe_wrt`) sound speed passed the shortcut enthalpy
   `((gamma+1)p + pi_inf)/rho` while `s_compute_speed_of_sound` subtracts `qv/rho`;
   exact for qv = 0, wrong for JWL where the qv slot carries the reference-energy
   accumulator. The shortcut now includes `+ qv` at all three probe sites in
   `m_data_output.fpp`, which is bit-identical for every qv = 0 case.
3. `relativity` with JWL was rejected only by the Python validator; now also a Fortran
   `@:PROHIBIT`, so hand-edited input files cannot bypass it.

Watchpoints for the Tuolumne sessions (verified as well founded, not defects)
1. Cray ftn GPU builds get no INLINEALWAYS: `cray_inline` emits its directive only on
   Cray CPU builds (`parallel_macros.fpp`); under OpenACC/OpenMP offload,
   `s_mg_mixture_variables` (now a nested select case plus a call to
   `s_mg_jwl_reference`) is a plain `routine seq`, and inlining into the flux kernels is
   at the compiler's discretion. Session 2 must therefore compare the JWL burn-case grind
   on Cray GPU explicitly, not assume the gfortran CPU parity carries over. Add to
   session 2: one JWL shocktube timing per tree.
2. amdflang cross-TU declare-target: `mg_mixture` is declared in
   `m_global_parameters_common`, updated to device in `m_variables_conversion` (~line
   362), and read inside flux kernels compiled in the Riemann solver TUs via the
   `@:accumulate_mixture` macro. If the device copy is stale under amdflang OpenMP
   offload, JWL cells silently run the stiffened accumulator (wrong answers, no crash),
   which the CPU golden suite cannot catch. Session 1 must include the CPU versus GPU
   diff of a JWL example ON THE AMDFLANG BINARY specifically; agreement there closes
   this watchpoint. If it fails, the fix pattern is the local re-update used for
   `Re_size` in `m_riemann_state`.
3. Cray codegen of nested select case inside a device `routine seq` has historic
   fragility; the single-fluid JWL shocktube versus the gfortran golden (session 1) is
   the sufficient probe, and a `--debug` Cray rebuild is the first triage step.

Audit dimensions that came back clean: all 14 `@:accumulate_mixture` expansions match the
original call sites including the `num_fluids - 1` and `alpha_lim` variants; the
`select case (eos_fl(i))` refactor is semantically identical to the prior if/else for
every enum value on all three targets; no raw acc/omp pragmas entered with the refactor.
