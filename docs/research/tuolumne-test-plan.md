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
2. Full golden suite, GPU: `./mfc.sh test -j 8 --gpu` (check `test --help` for the current
   GPU flag form on this MFC revision). The three JWL tests (976FEEC6 and the two example
   tests) run the mg_mixture path on device for the first time.
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
3. Weak scaling: HLLC bench case at 1, 2, 4, 8 nodes (4 APUs per node) with fixed memory
   per rank, submitted through Flux:
   `./mfc.sh run benchmarks/5eq_rk3_weno3_hllc/case.py -e batch -N <nodes> -n 4 ...`.
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

## Execution order on Tuolumne

Day 1: Phase 0 and builds 1 to 6 (compile matrix), Phase 2 items 1 and 6.
Day 2: Phase 2 items 2 to 5 (GPU correctness), Phase 3 items 1 and 2.
Day 3: Phase 3 items 3 to 5 (scaling), assemble the Phase 4 report.
Any red result stops the line for that lane; report it with the log rather than working
around it.
