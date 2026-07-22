# EOS PR Ladder: Refined Map and Implementation Plan

This is the implementation-state companion to the maintainer's ten-PR sequence for
incremental EOS support. It records, per rung, what the canonical scope is, what is
implemented on `feature/eos-combined`, what was learned during implementation, and what
remains before submission. The hardware validation campaign for the implemented rungs is
a separate document, `tuolumne-test-plan.md` in this directory; this file does not
restate it.

Governing principles, from the maintainer and preserved throughout:

1. Thermodynamics, reaction kinetics, and material data stay separate.
2. Pyrometheus remains the gas-chemistry backend and is not forced to represent
   phenomenological explosive burn laws.
3. The repository ships equation families, never a catalog of military-relevant
   materials. Every in-repo coefficient set is synthetic.
4. Each rung changes no answers except where its scope says so, and proves it with the
   golden suite.

## State summary

| Rung | Canonical scope | State | Where |
|---|---|---|---|
| PR0 | Probe enthalpy qv correction (split out, see PR1 notes) | Identified, unsubmitted | commit `d04ade6a` (partial) |
| PR1 | Central thermodynamics interface, no changed answers | Submitted | upstream #1663; tip `feature/thermodynamics-interface` |
| PR2 | Explicit `fluid_pp(:)%eos` selector, existing adapters only | Submitted | upstream #1664 |
| PR3 | Generic Mie-Gruneisen backend | Implemented, unsubmitted | scope on `feature/eos-combined` |
| PR4 | Nonreactive JWL as a reference-curve specialization | Implemented, unsubmitted | scope on `feature/eos-combined` |
| PR5 | External material-file mechanism | Implemented, unsubmitted | scope on `feature/eos-combined` |
| PR6 to PR10 | Reaction separation and families | Not started | design constraints below |

A warning that shapes the whole submission plan: the historical rung-tip branches
(`feature/eos-mie-gruneisen`, `feature/eos-jwl`) are stale. Twenty-one commits of
hardening, review fixes, and performance work landed on `feature/eos-combined` after
those tips were cut, including the vanished-phase gate, the amdflang firstprivate
idiom, the mg_mixture gate restructure, and the HLL/LF golden coverage. Submitting the
old tips would ship unhardened code. The rung branches are therefore rebuilt from the
combined diff at submission time, per the reconstruction procedure below; the old tips
are historical markers only.

`feature/eos-combined` carries all five implemented rungs plus the review-fix and
performance commits (40 commits, roughly +2900 net LOC over master `0857ace6`, about
half of which is tests, benchmarks, and documentation). Verification state: 620/620
golden tests bit-identical for all pre-JWL paths, drift-canceled grind parity with
master on every stock benchmark case (median pairwise deltas +0.1 to +0.4 percent), and
a measured JWL cell cost of +23 percent grind over stiffened HLLC (89.2 versus 72.4 ns
per grid point per equation per RHS on the reference laptop).

## Rung reconstruction: fix assignment and procedure

Every hardening commit on the combined branch has exactly one home rung. The assignment
below is the authority when the rung branches are rebuilt; a fix that ships in the
wrong rung either breaks the no-answer-change contract of an early rung or arrives too
late to protect the rung that needs it.

| Combined-branch work | Home rung | Reason |
|---|---|---|
| Six-op interface module, adapter wiring, `cray_inline` hints (`cf8b0819`) | PR1 | interface scope |
| Mixed-precision `f_bulk_modulus` call-site kinds, `f_validate_state` NaN parity (`362793d2`) | PR1 | correctness of the relocated helpers |
| Probe sound-speed enthalpy `+ qv` term (part of `d04ade6a`) | PR0, standalone | behavior change for any qv nonzero case; see below |
| Selector enum, name mapping, rejection checkers | PR2 | selector scope |
| MG leaf module, accumulator identity, `s_mg_mixture_variables` | PR3 | MG scope |
| Dispatch hoist and `@:accumulate_mixture` macro (`07fc90ac`, `9899810f`) | PR3 | the accumulator's performance contract |
| mg_mixture gate keyed on non-legacy families (`a460b17c`) | PR3 | the extension contract |
| JWL reference curve, parameters, checkers, examples | PR4 | JWL scope |
| Vanished-phase gate (`fcf286bd`, refining `09fa0a17`) | PR4 | JWL isentrope robustness |
| amdflang firstprivate host copy `mg_mixture_loc` (`702c88f4`) | PR4 | first rung whose device path reads the gate in kernels |
| Device refresh of `mixture_closure` from the defining TU (`09fa0a17`) | PR3 | the gate's own device copy |
| JWL feature gates incl. acoustic_source, relativity (`d04ade6a` remainder), sim_data (`b9f33ce6`) | PR4 | JWL scope |
| HLL/LF mixture goldens (`721d26ed`), example goldens, suite case | PR4 | JWL test coverage |
| `5eq_jwl_weno3_hllc` bench case and `bench.yaml` entry (`c0eeeff0`) | PR4 | CI would fail benching a JWL case before PR4 merges |
| `materials.py`, `test_materials.py` | PR5 | material-file scope |
| Comment and documentation corrections (`de3d52d0`, `96c02284`) | with their subject's rung | |
| Tuolumne plan, this document | fork only | not upstream material |

PR0, the probe qv correction, is deliberately split out and submitted first as a tiny
standalone fix. The simulation probe output path builds the sound-speed enthalpy
without the qv term while the main path (`s_compute_enthalpy`) and the post-process
energy loop both include it; the correction changes probe output for any case with
nonzero qv, stiffened cases included, and no golden covers probes with qv. Hiding a
behavior change inside PR1, whose entire claim is "no changed answers", would hand a
reviewer a reason to distrust the ladder. As its own two-line PR with the two
authoritative call sites quoted, it is instead an easy merge.

Reconstruction procedure, applied per rung as its predecessor merges:

1. Branch from current upstream master (predecessor merged).
2. Apply the rung's scope as a fresh, small commit series extracted from the combined
   diff using the file and assignment tables, not by cherry-picking the 40-commit
   history; the history interleaves rungs and its intermediate states are not
   individually verified.
3. Regenerate the rung's new goldens on the rebuilt tip and confirm they match the
   combined-branch values; test UUIDs are stable (CRC32 of the trace string) so the
   identifiers carry over.
4. Build all three targets, run the full suite, run precheck, and for PR3 onward run
   the drift-canceled paired benchmark against the rebuilt tip's own base.
5. Diff the rebuilt tip against the corresponding file states on `feature/eos-combined`
   apart from intentional rebase drift; any unexplained difference is a porting error.

The combined branch remains the integration proof (all rungs coexist and pass
together) and the Tuolumne test article. It is never itself submitted.

## PR1: thermodynamics interface (submitted, #1663)

Canonical scope: centralize `pressure`, `internal_energy`, `temperature`, `sound_speed`,
`thermodynamic_derivatives`, `validate_state` behind one interface with
`stiffened_gas` and `pyrometheus` adapters; bit-identical results.

Implemented as `src/common/m_thermodynamics.fpp` (251 lines). Learnings folded back
into the branch:

1. Fypp macros, not procedure pointers, are the dispatch mechanism that survives all
   four CI compilers and both GPU backends. Polymorphic or pointer-based dispatch
   inside device kernels is exactly the construct that diverges across nvfortran, Cray,
   and amdflang offload; textual expansion keeps the hot path inlinable everywhere.
   The `@:accumulate_mixture` macro in `src/simulation/include/inline_riemann.fpp` is
   the pattern.
2. Dispatch placement is a measured performance matter, not style. Putting a runtime
   branch inside the leaf accumulator blocked gfortran inlining and cost 8 percent on
   HLLC grind; hoisting the branch to the call sites restored parity. Any reviewer
   question about interface overhead is answered by the drift-canceled bench table.
3. Small pure helpers carry `cray_inline` hints. These have never been exercised by a
   Cray compiler locally; the Tuolumne build matrix gates them.

Remaining: fold the Tuolumne Cray CPU, Cray GPU, and AMD GPU bench columns into the PR
body when available, and land PR0 first so this PR's no-answer-change claim is exact.

## PR2: explicit EOS selector (submitted, #1664)

Canonical scope: stable enumeration (`stiffened_gas`, `ideal_gas_mixture`,
`mie_gruneisen`, `jwl`, `table`), readable names in case files, integer internal,
stiffened gas default, unsupported combinations rejected explicitly, no promise of
arbitrary per-cell family mixing.

Implemented exactly so: constants in `src/common/m_constants.fpp`, per-fluid `%eos` in
the derived types, name-to-enum mapping in the toolchain, rejection in
`m_checker_common.fpp` and `case_validator.py`. `mie_gruneisen` and `table` are
enumerated but rejected by the checkers until their rungs land, which keeps the
enumeration stable from day one.

## PR3: generic Mie-Gruneisen (implemented, unsubmitted)

Canonical scope: reference pressure and energy curves, Gruneisen coefficient,
pressure/energy inversion, analytic sound speed and derivatives, synthetic manufactured
tests. JWL then becomes a reference curve, not a special path.

Implemented as `src/common/m_eos_mie_gruneisen.fpp` (a 126-line leaf) with
manufactured tests in `toolchain/mfc/test_mg_eos.py` (CI-run pytest, 203 lines). The
central identity every family plugs into is the accumulator triple: given a family's
(Gamma, p_ref, dp_ref, e_ref, de_ref) at the phasic density,

    gamma  += alpha / Gamma
    pi_inf += alpha * ((rho_i * dp_ref - p_ref) / Gamma - p_ref)
    qv     += alpha_rho * (e_ref + rho_i * de_ref - dp_ref / Gamma)

which reduces exactly to the legacy stiffened slot arithmetic for a linear reference
curve. That reduction is the manufactured test and the reason all-stiffened cells stay
bit-identical. The identity reuses the existing gamma, pi_inf, and qv mixture slots in
their existing algebraic roles, so no downstream consumer of the mixture state changes.

Implementation refinements beyond the canonical text:

1. The five-equation model always has mixture cells, so the mixture rule cannot wait
   for PR4. `s_mg_mixture_variables` in `m_variables_conversion.fpp` owns it, with one
   `select case (eos_fl(i))` arm per family.
2. The extension contract is deliberately loud: the `mg_mixture` gate keys on "any
   fluid outside the two legacy-slot families", so a new family activates the
   generalized path with no gate edit, and a family missing its case arm contributes
   nothing and fails immediately with a zero mixture gamma rather than silently running
   stiffened arithmetic. Adding a family touches three files: the constants enum, the
   case arm, and the checker allowlist.
3. This rung exposes no new user-facing choice. `fluid_pp(:)%eos = 'mie_gruneisen'`
   stays checker-rejected: the backend exists for families to build on, and exposing a
   bare generic family without a validated closure, parameter checks, and tests of its
   own invites misuse. The first consumer is JWL in PR4; a user-facing generic MG
   (for example a linear-Hugoniot solid) becomes its own later rung with its own
   checkers when something needs it. Because nothing user-visible changes, this PR's
   acceptance is the same as PR1's: full suite bit-identical plus the manufactured
   tests.

Remaining: rebuild per the reconstruction procedure after PR2 merges; the module is a
leaf so rebase conflicts should be near zero.

## PR4: nonreactive JWL (implemented, unsubmitted)

Canonical scope: the JWL reference curve, derivatives and inversion, user-supplied
parameters, generic products-expansion and shock-tube cases, no named-explosive catalog,
no reaction source.

Implemented on the PR3 base: `s_mg_jwl_reference` evaluates

    E1 = A exp(-R1 V),  E2 = B exp(-R2 V),  V = rho0 / rho
    p_ref = E1 + E2,    dp_ref = (V / rho) (R1 E1 + R2 E2)
    e_ref = E1 / (rho0 R1) + E2 / (rho0 R2),  de_ref = p_ref / rho^2,  Gamma = omega

and the accumulator identity above does the rest; there is no JWL-specific solver path.
User parameters are `fluid_pp(i)%jwl_A/B/R1/R2/omega/rho0`, positivity- and
ordering-checked (`R1 > R2 > 0`) in `s_check_jwl_inputs`. In-repo cases use synthetic
coefficients only (examples: A = 5e11, B = 8e9, R1 = 4.5, R2 = 1.2, omega = 0.3,
rho0 = 1600; suite and benchmark: the O(1) nondimensional set).

Hardening shipped with this rung, all part of its review story:

1. Vanished-phase gate: a phase driven below `sgm_eps` by reconstruction overshoot has
   an ill-defined phasic density; its reference-curve contribution is skipped while its
   mass still enters the mixture density. This trades a NaN for a bias of order
   `sgm_eps` and is bit-identical on every current test cell.
2. amdflang cross-TU staleness: the Riemann kernels read `mg_mixture` through a
   firstprivate host-local copy, following the proven `Re_size_loc` idiom already in
   those kernels.
3. Feature gates: JWL rejects bubbles, phase change, elasticity, MHD, surface tension,
   IGR, immersed boundaries, acoustic sources, relativity, and post-process `sim_data`,
   each with a named checker message. Gates are lifted rung by rung, never implicitly.
4. Coverage: golden tests for HLLC, HLL, and Lax-Friedrichs mixture accumulation, two
   example-based goldens, and the `5eq_jwl_weno3_hllc` benchmark case that prices the
   feature (+23 percent grind on JWL cells, from the two exponentials per fluid per
   accumulation). The `bench.yaml` entry must ride in this PR and no earlier, since
   upstream benchmark CI would otherwise run a case the solver rejects.

Remaining: rebuild after PR3 merges; attach the CPU versus GPU agreement runs and the
measured cost factor from Tuolumne.

## PR5: external material files (implemented, unsubmitted)

Canonical scope: Cantera-style material file with `family`, `parameters`, and
`provenance` (`citation`, `release_status`); search by explicit path, case directory,
then configured public directory; analytic parameters stay runtime values.

Implemented as `toolchain/mfc/materials.py` (a 126-line loader) with unit tests in
`test_materials.py` (148 lines): schema validation, provenance required,
`release_status` gated, loader resolves into ordinary `fluid_pp` runtime values so the
Fortran side needs no change. The repository keeps the loader and synthetic
demonstration files; every calibrated real-explosive set stays external by
construction. The canonical schema sketch includes an energy-release parameter `Q`;
the nonreactive ladder has no consumer for it, so the schema accepts only parameters
its family consumes and `Q` enters the schema with PR8, keeping file validation strict
instead of silently carrying dead keys.

Remaining: rebuild after PR4 merges; consider the CI grep gate that no real-explosive
coefficient block appears under `src/` or committed `examples/`.

## Anticipated review objections, per rung

Preparing the answer before the question is most of why the combined branch exists.

| Rung | Likely objection | Prepared answer |
|---|---|---|
| PR0 | is the qv term correct? | the two authoritative paths (`s_compute_enthalpy`, post energy loop) both include qv; the probe path is the outlier |
| PR1 | interface overhead in the hot path | drift-canceled bench table: +0.1 to +0.4 percent, within noise, after the dispatch-placement fix |
| PR1 | why macros, not objects | four compilers and three GPU configurations cannot all inline dynamic dispatch in kernels; measured 8 percent penalty from one misplaced branch |
| PR1 | pre-existing defects visible in relocated code | upstream issue filed beforehand (see below); relocation is byte-identical to master |
| PR2 | dead enum values | stability of the enumeration is the point; checkers reject them loudly until implemented |
| PR3 | why is the mixture rule in `m_variables_conversion` and not the leaf | inlining contract; the leaf stays pure so the accumulator inlines, the dispatch lives at call sites |
| PR3 | what stops a half-added family | zero-gamma loud failure by construction, three-file extension contract, extension drill in the test plan |
| PR4 | why no wave-speed or solver changes | the accumulator identity feeds the existing solvers; JWL is state evaluation only, which is the point of PR3 |
| PR4 | cost of the feature | measured: +23 percent grind on JWL cells, zero on stiffened cases; quoted in the PR body |
| PR4 | robustness at interfaces | vanished-phase gate with stated bias bound; eps-vanished phases exercised in suite, examples, and the bench case |
| PR5 | is this an explosives library | no: loader plus synthetic demonstrations; provenance and release gating are how real data stays external |

## Submission sequence and dependencies

The chain is strict: PR0, then PR1 through PR5, each submitted only after its
predecessor merges, each rebuilt onto current master at submission time per the
reconstruction procedure. Tuolumne results slot into PR bodies as follows: PR1 gets the
three-platform bench table, PR2 and PR3 get build-matrix results, PR4 and PR5 get the
CPU/GPU agreement runs and the JWL cost factor.

Before PR1 review advances, the two pre-existing upstream defects found during review
(conserved-variable slot misuse and a scalar-stress energy loop in the relocated
`s_compute_pressure`) are filed as an upstream issue with the byte-identical evidence
against master. Reviewers reading PR1 will see the relocated routine and may attribute
the defects to the move; a pre-filed issue converts that moment from a liability into
demonstrated diligence, and keeps their fixes out of the ladder where they would break
the bit-identity claims.

## PR6 to PR10: reactions (not started, design constraints recorded)

These rungs are future work; what follows pins the architecture so nothing in PR1 to
PR5 has to be undone.

PR6 separates thermodynamic backend, reaction mechanism, and reaction integrator, with
`reaction_rates`, `advance_reactions`, and `reaction_timestep_constraint` as the
integrator surface. Nothing in the implemented rungs couples EOS selection to reactions,
so this split starts clean.

PR7 hardens the Pyrometheus path (mechanism identity in build hashes, provenance,
reactor tests, CPU/GPU agreement) before any new EOS family couples to it.

PR8 adds the generic progress variable for condensed materials: lambda, a supplied rate
law, reactant EOS, product EOS, energy release. Pyrometheus is deliberately not the
vehicle for this; program burn, ignition-and-growth, and pressure-triggered condensed
reactions use the same integrator with a different mechanism backend. The material-file
schema gains `Q` here, with its provenance fields mandatory.

PR9 is the hard thermodynamics: reactant and product EOS coupled through reaction
progress with consistent mixture energy, pressure equilibrium, the reaction energy
offset, and intermediate-state sound speeds. It is explicitly its own rung so it cannot
hide inside a reaction source term. The generalized accumulator structure from PR3 is
the intended home: a reacting cell is a mixture whose composition weights move with
lambda.

PR10 adds specific families one PR each (programmed burn, pressure-dependent progress,
ignition-and-growth, Pyrometheus-coupled afterburn where appropriate), synthetic
parameters only, calibration external through the PR5 mechanism.

Reference implementations for the reaction rungs exist in the frozen research tree
(`m_jwl_sources.fpp`: exact exponential progress update, front-CFL guard, solid-cell
gating) and are ported and decomposed when their rung starts, never merged wholesale.
