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
| PR1 | Central thermodynamics interface, no changed answers | Submitted | upstream #1663; tip `feature/thermodynamics-interface` |
| PR2 | Explicit `fluid_pp(:)%eos` selector, existing adapters only | Submitted | upstream #1664 |
| PR3 | Generic Mie-Gruneisen backend | Implemented, unsubmitted | tip `feature/eos-mie-gruneisen` |
| PR4 | Nonreactive JWL as a reference-curve specialization | Implemented, unsubmitted | tip `feature/eos-jwl` |
| PR5 | External material-file mechanism | Implemented, unsubmitted | on `feature/eos-combined` |
| PR6 to PR10 | Reaction separation and families | Not started | design notes below |

`feature/eos-combined` carries all five implemented rungs plus the review-fix and
performance commits (40 commits, roughly +2900 net LOC over master `0857ace6`, half of
which is tests, benchmarks, and documentation). Verification state: 620/620 golden tests
bit-identical for all pre-JWL paths, drift-canceled grind parity with master on every
stock benchmark case (median pairwise deltas +0.1 to +0.4 percent), and a measured JWL
cell cost of +23 percent grind over stiffened HLLC (89.2 versus 72.4 ns per grid point
per equation per RHS on the reference laptop).

## PR1: thermodynamics interface (submitted, #1663)

Canonical scope: centralize `pressure`, `internal_energy`, `temperature`, `sound_speed`,
`thermodynamic_derivatives`, `validate_state` behind one interface with
`stiffened_gas` and `pyrometheus` adapters; bit-identical results.

Implemented as `src/common/m_thermodynamics.fpp`. Learnings folded back into the branch:

1. Fypp macros, not procedure pointers, are the dispatch mechanism that survives all
   four CI compilers and both GPU backends. Textual expansion keeps the hot path
   inlinable; the `@:accumulate_mixture` macro in
   `src/simulation/include/inline_riemann.fpp` is the pattern.
2. Dispatch placement is a measured performance matter, not style. Putting a runtime
   branch inside the leaf accumulator blocked gfortran inlining and cost 8 percent on
   HLLC grind; hoisting the branch to the call sites restored parity. Any reviewer
   question about interface overhead is answered by the drift-canceled bench table.
3. Small pure helpers carry `cray_inline` hints. These have never been exercised by a
   Cray compiler locally; the Tuolumne build matrix gates them.

Remaining: fold the Tuolumne Cray CPU, Cray GPU, and AMD GPU bench columns into the PR
body when available.

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

Implemented as `src/common/m_eos_mie_gruneisen.fpp` with manufactured tests in
`toolchain/mfc/test_mg_eos.py` (CI-run pytest, 203 lines). The central identity every
family plugs into is the accumulator triple: given a family's
(Gamma, p_ref, dp_ref, e_ref, de_ref) at the phasic density,

    gamma  += alpha / Gamma
    pi_inf += alpha * ((rho_i * dp_ref - p_ref) / Gamma - p_ref)
    qv     += alpha_rho * (e_ref + rho_i * de_ref - dp_ref / Gamma)

which reduces exactly to the legacy stiffened slot arithmetic for a linear reference
curve. That reduction is the manufactured test and the reason all-stiffened cells stay
bit-identical.

Implementation refinement beyond the canonical text: the five-equation model always has
mixture cells, so the mixture rule cannot wait for PR4. `s_mg_mixture_variables` in
`m_variables_conversion.fpp` owns it, with one `select case (eos_fl(i))` arm per family.
The extension contract is deliberately loud: the `mg_mixture` gate keys on "any fluid
outside the two legacy-slot families", so a new family activates the generalized path
with no gate edit, and a family missing its case arm contributes nothing and fails
immediately with a zero mixture gamma rather than silently running stiffened arithmetic.
Adding a family touches three files: the constants enum, the case arm, and the checker
allowlist.

Remaining: submit after PR2 merges; rebase is expected clean since the module is a leaf.

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

Hardening added during review, all part of this rung's story:

1. Vanished-phase gate: a phase driven below `sgm_eps` by reconstruction overshoot has
   an ill-defined phasic density; its reference-curve contribution is skipped while its
   mass still enters the mixture density. This trades a NaN for a bias of order
   `sgm_eps` and is bit-identical on every current test cell.
2. amdflang cross-TU staleness: the Riemann kernels read `mg_mixture` through a
   firstprivate host-local copy, following the proven `Re_size_loc` idiom.
3. Feature gates: JWL rejects bubbles, phase change, elasticity, MHD, surface tension,
   IGR, immersed boundaries, and post-process `sim_data`, each with a named checker
   message. Gates are lifted rung by rung, never implicitly.
4. Coverage: golden tests for HLLC, HLL, and Lax-Friedrichs mixture accumulation, two
   example-based goldens, and the `5eq_jwl_weno3_hllc` benchmark case that prices the
   feature (+23 percent grind on JWL cells, from the two exponentials per fluid per
   accumulation).

Remaining: submit after PR3; attach the CPU versus GPU agreement runs and the measured
cost factor from Tuolumne.

## PR5: external material files (implemented, unsubmitted)

Canonical scope: Cantera-style material file with `family`, `parameters`, and
`provenance` (`citation`, `release_status`); search by explicit path, case directory,
then configured public directory; analytic parameters stay runtime values.

Implemented as `toolchain/mfc/materials.py` with unit tests in `test_materials.py`
(148 lines): schema validation, provenance required, `release_status` gated, loader
resolves into ordinary `fluid_pp` runtime values so the Fortran side needs no change.
The repository keeps the loader and synthetic demonstration files; every calibrated
real-explosive set stays external by construction.

Remaining: submit after PR4; consider the CI grep gate that no real-explosive
coefficient block appears under `src/` or committed `examples/`.

## Submission sequence and dependencies

The chain is strict: PR1 then PR2 then PR3 then PR4 then PR5, each submitted only after
its predecessor merges, each rebased onto current master at submission time. The rung
tips exist as local branches and are re-verified (build plus golden suite) after any
rebase. Tuolumne results slot into PR bodies as follows: PR1 gets the three-platform
bench table, PR2 and PR3 get build-matrix results, PR4 and PR5 get the CPU/GPU
agreement runs and the JWL cost factor. Two pre-existing upstream defects found during
review (conserved-variable slot misuse and a scalar-stress energy loop in the relocated
`s_compute_pressure`) are byte-identical to master and will be reported upstream
separately, not fixed in this ladder.

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
reactions use the same integrator with a different mechanism backend.

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
