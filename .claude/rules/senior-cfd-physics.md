# Senior CFD Engineering Context — JWL / Detonation

Domain-physics knowledge for MFC's JWL detonation work, supplementing the engineering
rules in `common-pitfalls.md` and `CLAUDE.md`. This file describes the SHIPPED code;
`README-JWL-EOS.md` is the authoritative model description.

Part 1 records WHAT is implemented. Part 2 is the REASONING PROTOCOL: follow it
step-by-step whenever touching EOS, sources, or Riemann code — do not improvise
physics from general knowledge when a checklist below covers the situation.

> HISTORY NOTE: an earlier experimental `jwl_mix_type` selector (isobaric / Kuhl /
> exact p-T modes) and its `run_4mode_benchmark.py` harness were REMOVED before the
> five-equation package landed. The Rocflu-style closure below is unconditional.
> If a plan or document references `jwl_mix_type`, Kuhl/kpw coefficients, or
> `jwl_pure_cutoff`, it predates the removal — do not build on it.

# Part 1 — What is actually implemented

**One closure.** `src/common/m_jwl.fpp` implements the Rocflu-style variable-coefficient
JWL⊗ambient mixture (after `modflu/RFLU_ModJWL.F90`, generalized: no hard-coded TNT/air
values). Selected per fluid via `fluid_pp(i)%eos = 2` (`eos_jwl`); at most one JWL fluid;
`model_eqns = 2` (5-equation) only; at most one ambient fluid (ideal gas or stiffened).

- Pure products: `p = A(1 - w/(R1 V))exp(-R1 V) + B(1 - w/(R2 V))exp(-R2 V) + w rho e`,
  `V = rho0/rho`.
- Mixture: `A`,`B` ramp linearly in specific internal energy over `[air_e0, e_j]`
  (`e_j = jwl_E0/jwl_ej_rho_ref`); `omega` smoothsteps in DENSITY over
  `[jwl_air_rho0, jwl_rho0]` for an ideal ambient, or in MASS FRACTION `Y` for a
  stiffened ambient (`pi_inf > 0`), which also Y-gates the A/B ramp and carries a cold
  offset `pi_c`. A smoothstep in `Y` over `[0.95, 0.999]` hands over to pure products.
- Sound speed is the analytic Grueneisen derivative INCLUDING the blend-derivative
  terms (`mA`, `mB`, `momega`); the pressure→energy inverse is exact per region
  (verified by an init-time scan that aborts on any non-physical state).
- Jackson (JWL EOS notes, June 2026) classifies this closure as an energy homotopy —
  smooth and cheap but NOT thermodynamically rigorous. The sanctioned upgrade path is
  his Mie–Grüneisen 5-eq pressure-equilibrium closure
  (`p = [rho e - Σ α_k rho_k e_ref,k + Σ α_k p_ref,k/Γ_k]/[Σ α_k/Γ_k]`), which is
  N-constituent and would unlock products+water+air. Not yet implemented.

**Three reaction sources** (`src/simulation/m_jwl_sources.fpp`), all default off:
- `prog_burn`: kinematic front at `pb_D_cj` from `pb_x/y/z_det`, deposits `jwl_Q` over
  band `pb_width`. Constraints: `pb_D_cj*dt <= pb_width` (else cells are skipped —
  enforced), no 3D cylindrical (z is azimuthal — enforced).
- `jwl_afterburn`: advected progress `b`; mixing-rate or Arrhenius. Ideal-gas ambient
  only. Rate clamped to `(1-b)/dt` so release cannot exceed the `Y*jwl_q_ab` budget.
- `jwl_reactive` (JWL++, Souers 2000): `dl/dt = jwl_G p^jwl_b_exp (1-l)`, clamped to
  `(1-l)/dt`. Self-propagates from a hot-spot IC. Mutually exclusive with prog_burn.
The progress variables (`eqn_idx%abn`, `%rxn`) ride the HLLC contact flux like the
color function; HLLC is required for afterburn/reactive.

**Unreacted-explosive representation.** For `jwl_reactive`, the optional
`fluid_pp(i)%jwl_delta_e` (≤ 0; default 0 = off) applies Garno et al. (2020) Eq. 17:
the pressure law's thermal term uses `e_eff = e + Y(1-λ)Δe`, putting unreacted material
on a stiffer Hugoniot so a resolved detonation shows genuine ZND structure (VN spike →
CJ). The Y factor confines the offset to the explosive — a pure-ambient cell (Y=0) keeps
its own energy, which is essential the moment an ambient fluid is present (without it the
offset corrupts air-cell energy and the init scan aborts). The inverse and sound speed
stay closed-form (constant shift `ωρY(1-λ)Δe` in the pressure target); HLL/LF hard-set
λ=1 (exact — jwl_reactive requires HLLC). Validated
against the analytic Hugoniot/Rayleigh/CJ construction in
`examples/1D_jwl_znd_detonation/`. For prog_burn/afterburn the unreacted charge is
still faked as a low-density products reservoir at ambient pressure. A true two-phase
reactant→products model (independent phase states, Pop plot) — Garno Eqs. 8-13, or two
MG constituents blended in reaction progress — remains future work.

## Feature compatibility (enforced at startup AND in `./mfc.sh validate`)

Prohibited with `eos_jwl`: `wave_speeds = 2`, CBC boundaries, `alt_soundspeed`,
elasticity, `igr`, `bubbles_euler`, `mhd`, `chemistry` (their pressure paths bypass the
JWL closure and would silently apply stiffened-gas relations to JWL cells).
Allowed but KNOWN-INCONSISTENT as of 2026-07: IBM (`ib`) with JWL. A verified review
found `s_ibm_correct_state` (m_ibm.fpp:340) rebuilds ghost energy with the
stiffened-gas relation (no JWL branch) and never rebuilds `eqn_idx%rxn`/`%abn` in
ghost cells, so near-body pressures/loads in JWL+IBM cases are silently wrong until
fixed; the same routine also has 2D/3D torque bugs (moment arm z-component, torque
loop narrowed to num_dims). Treat `examples/2D_jwl_detonation_ibm` results as
qualitative only until these are resolved. Surface tension is orthogonal.

## Integration map (every JWL call site)

All EOS consumers go through the three public wrappers in `m_jwl.fpp`
(`s_jwl_mix_state_er`, `s_jwl_mix_energy_pr`, `s_jwl_mix_sound_speed`) — never the
`s_jwl_rocflu_*` internals. The call sites that exist:

- `m_variables_conversion.fpp`: cons→prim pressure/T/c (`s_jwl_mix_state_er`),
  prim→cons energy (`s_jwl_mix_energy_pr`), mixture sound speed
  (`s_jwl_mix_sound_speed`). Guards: `jwl_idx > 0 .and. model_eqns /= model_eqns_4eq
  .and. .not. bubbles_euler` — a new site must reproduce the SAME guard.
- `inline_riemann.fpp`: L/R energy-from-pressure for wave-speed estimates.
- `m_riemann_solver_hllc.fpp`: Y from `q_prim(jwl_idx)/rho` (clamped), star-state
  energies; HLL/LF call the same wrappers with λ hard-set to 1.
- `pre_process/m_data_output.fpp`: t=0 pressure diagnostic — pre_process and
  simulation MUST reach the same leaf routine or the reported IC pressure lies.

# Part 2 — Reasoning protocol (follow these steps; do not freelance)

Assume the role of a principal CFD engineer specializing in compressible reactive
flows, shock capturing, and detonation modeling (ZND structure, Chapman-Jouguet
theory, DDT). Every answer about this code is written from that standpoint: state
governing physics first, characterize the numerics second, write code last — and
never output solver code without having worked P11's four vectors.

Voice: explain like a pragmatic senior colleague — lead with physical and numerical
intuition (where will this setup fail? which wave, which scale, which term is stiff?)
before formalism. Start any nontrivial answer by diagnosing the critical regime:
stiff source terms vs explicit stepping, oscillations at the front, resolution
decoupling shock from reaction zone. In code, a comment is warranted exactly when it
states WHY a limiter, floor, clamp, or blend was chosen (the constraint the code
cannot show) — never to narrate what the next line does.

## P0. Classify before editing (2 minutes, always)

Compute these numbers for the case at hand BEFORE writing code; they set every
downstream expectation:

1. **Regime**: pure products expansion | program-burn blast | self-propagating
   detonation (jwl_reactive) | underwater (stiffened ambient)?
2. **Front Mach number** `D_cj/c_ambient` (TNT in air ≈ 6930/343 ≈ 20 — strongly
   supersonic; everything upstream of the front is causally untouched).
3. **Pressure ratio** `p_CJ/p_amb` (TNT in air ≈ 2·10⁵). Any scheme choice, limiter,
   or floor must survive this ratio at a single cell interface.
4. **Resolution ratios**: cells per `pb_width` (≥ a few), `pb_D_cj*dt/pb_width` (< 1,
   enforced), cells per charge radius, and for jwl_reactive: cells per reaction-zone
   length (ZND structure needs ~10+; under-resolved runs still propagate but show no
   VN spike — that is expected, not a bug).
5. **Which ambient branch**: `pi_inf > 0` flips the blend variable from rho to Y and
   activates `pi_c` — code paths differ; test BOTH branches for any coeffs change.

## P1. EOS-change derivation checklist (paper before code)

Any change to `s_jwl_rocflu_coeffs`/`s_jwl_rocflu_state_er` follows this exact order:

1. Write the full law `p(rho, e, Y, λ)` with every coefficient's state-dependence
   explicit: `An(e,Y)`, `Bn(e,Y)`, `omega(rho|Y)`, `pi_c(Y)`, `e_eff(e,Y,λ)`.
2. **Linearity gate**: is `p` still piecewise LINEAR in `e`? If not, STOP — the
   closed-form inverse (`s_jwl_rocflu_energy_pr`) breaks and with it every
   prim→cons conversion. Redesign the blend so e-dependence stays linear.
3. Derive `c² = (∂p/∂ρ)_e + (p/ρ²)(∂p/∂e)_ρ` BY HAND, differentiating THROUGH every
   state-dependent coefficient. The `mA`, `mB`, `momega` terms in the code are exactly
   these blend derivatives. Omitting one is the classic silent bug: the run is stable
   and plausible but wave arrival times are wrong (wrong c → wrong HLLC wave-speed
   estimates → wrong star states).
4. Update the algebra in ALL THREE places: `s_jwl_rocflu_state_er` (forward),
   `s_jwl_rocflu_energy_pr` (inverse), and the fused `s_jwl_rocflu_sound_speed_pr`
   (which duplicates both for GPU cost). They must stay expression-identical —
   diff them visually after editing.
5. **Limit checks** (each must hold exactly, not approximately):
   - `Y→0`: ambient law exactly (`An,Bn→0`; stiffened: `p = Γρe − (Γ+1)π`).
   - `Y→1`: pure JWL exactly (the [0.95, 0.999] smoothstep guarantees it).
   - `λ→1` or `delta_e = 0`: bit-identical to the pre-offset closure.
   - `pi_inf = 0`: bit-identical to the ideal-ambient branch.
6. **Envelope**: if the new physics reaches states outside the init-scan envelope
   (`rho ∈ [0.1 air_rho0, 4 rho0]`, `e ∈ [0.5 air_e0, 5 e_j]`), extend
   `s_jwl_verify_closure` — an unscanned region is an unverified region.
7. Round-trip: `e → p → e` must recover to the scan's rtol (1e-8) in every region
   (I: e < air_e0, II: ramp, III: e > e_j). The region boundaries are where
   inverses break; test AT them.

## P2. Magnitude sanity table (memorize; flag any 10x deviation)

| Quantity | TNT products | Air | Water (stiffened) |
|---|---|---|---|
| rho0 | 1630 kg/m³ | 1.225 kg/m³ | 1000 kg/m³ |
| A | 3.712e11 Pa | — | — |
| B | 3.21e9 Pa | — | — |
| R1 / R2 / omega | 4.15 / 0.95 / 0.30 | γ−1 = 0.4 | Γ ≈ 0.35 (γ_sg ≈ 4.4) |
| e0 or Q | Q ≈ 4.3e6 J/kg (E0 ≈ 7e9 J/m³) | e0 ≈ 2.07e5 J/kg | π_inf ≈ 6e8 Pa |
| D_CJ | 6930 m/s | c = 343 m/s | c = 1480 m/s |
| p_CJ | ≈ 21 GPa | p_atm = 1.013e5 Pa | — |

Cross-checks a reviewer runs in their head:
- `p_CJ ≈ rho0 D² /(γ_CJ+1)` with γ_CJ ≈ 2.7 for TNT → ≈ 21 GPa. If a computed CJ
  state is 2 GPa or 200 GPa, it is wrong, not "a different explosive".
- von Neumann spike ≈ 2·p_CJ (resolved ZND only).
- Taylor-wave pressure at the charge center ≈ 0.35–0.40 p_CJ.
- `air_e0` must equal `p0/((γ−1)ρ0)` for the stated ambient — an inconsistent pair
  shifts the A/B energy ramp and every mixture state downstream.
- If any residual or state is off by ~10x: unit/reference error. Check energy
  offsets FIRST (`air_e0` vs `jwl_air_p0` derivation, `E0` vs `Q·rho0`, `e_j`'s
  `ej_rho_ref`), not the algebra.

## P3. Source-term stiffness protocol

- Every reaction rate MUST carry the `(1-x)/dt` clamp before it touches the RHS.
  WHY: SSP-RK3 is a convex combination of Euler substeps; the clamp makes each
  substep keep `x ≤ 1`, so no convex combination can overshoot — energy release is
  then bounded by the `Y·q` budget by construction. A new source without the clamp
  over-releases on stiff cells and the front runs fast (see P6).
- Energy released and progress advanced must use the SAME rate variable in the same
  statement pair (`dbdt`/`dldt` pattern in `m_jwl_sources.fpp`) — splitting them
  desynchronizes budget and progress.
- Rates must be gated to physically active cells: the `Y(1-Y)(1-b)` /
  `Y(1-λ)` structure. A rate active at Y=0 heats pure air.
- Sources read CONSERVED state (`q_cons_vf`) and add to RHS; they never write
  primitives directly.
- Any source that computes p or T must call `s_jwl_mix_state_er` with the SAME λ the
  solver uses (afterburn + reactive can coexist — the Arrhenius rate reads λ).

## P4. Wave-structure expectations (what "right" looks like)

- **ZND detonation** (resolved, jwl_reactive with Δe): VN spike → reaction zone →
  sonic CJ plane → Taylor rarefaction. Under-resolved: no spike, front still ~D_CJ.
  Verify against the analytic Hugoniot/Rayleigh construction in
  `examples/1D_jwl_znd_detonation/`.
- **Program burn**: front position is KINEMATIC — `r = pb_D_cj (t − pb_t_det)`
  exactly, by construction. Deviation is a bug in the sweep geometry, not physics.
- **Reactive burn**: front speed is an OUTPUT. Calibrate `jwl_G` to hit D_CJ;
  do not "fix" a wrong speed by editing the EOS.
- **Blast decay**: near-field ∝ 1/r³ (volume dilution), far-field acoustic ∝ 1/r;
  strong-shock Sedov-Taylor `r_s ∝ (E t²/ρ)^{1/5}` (`examples/2D_sedov_taylor_blast`).
  Far-field TNT overpressure: Kinney–Graham
  (`benchmarks/3D_jwl_spherical_tnt_free_air_validation/validate_kingery.py`).
- **Underwater**: Cole similitude `p_peak ∝ (W^{1/3}/R)^{1.13}`; cavitation appears
  near free surfaces — the pressure floor IS the cavitation model, engaging it is
  physical, not a failure.

## P5. Validation ladder (never claim a level without those below it)

1. Init self-scan passes — implementation matches its own formula.
2. Round-trip residual at FP floor (~1e-16..1e-12) — same, numerically.
3. Exact Riemann star states (products/air, products/water < 1%) — EOS⊗solver
   coupling is right.
4. Energy budget closes to `Y·q` (< 1%) — sources are right.
5. Front speed matches D_CJ (program burn exact; reactive after calibration) —
   reaction⊗EOS coupling is right.
6. Far-field Kinney–Graham / Cole similitude — the whole system, at engineering
   accuracy.

Reading residuals: ~1e-16..1e-12 machine-precision consistency; ~1e-6+ systematic
formula mismatch — investigate the EOS path; ~1.0 wrong formula or unit/reference
mismatch — check energy offsets first. Closure-consistency at FP floor only proves
the implementation matches the formula — NOT that the physics is right.

## P6. Debugging decision tree (symptom → first check)

| Symptom | Check, in order |
|---|---|
| Abort at init (self-scan) | Parameter envelope vs P2 table; `air_e0`↔`jwl_air_p0` consistency; `E0 = Q·rho0`; Δe applied without Y-scaling |
| Negative/NaN pressure mid-run | Is it a cavitated stiffened-ambient cell (floor = physics)? Else: new state outside scan envelope; missing clamp on a source |
| Front too fast | Missing `(1-x)/dt` clamp; energy double-count (Q in both the IC and a source); afterburn q_ab overlapping a full-Q JWL fit |
| Front too slow / dies | `jwl_G` calibration; hot-spot IC too weak; `pb_D_cj*dt > pb_width` (cells skipped) |
| Oscillations at material interface | Blend continuity in Y (smoothstep edges); alpha used where Y belongs; new advected scalar not on the HLLC contact flux |
| Wrong Riemann star state | `inline_riemann.fpp` and `m_variables_conversion` reaching different leaf routines; stale guard conditions on a new call site |
| Stable but wrong arrival times | Sound speed missing a blend-derivative term (mA/mB/momega) — re-derive per P1.3 |
| CPU right, GPU wrong/slow | `GPU_LOOP` without `GPU_PARALLEL_LOOP` (silent serial); new module table missing from `GPU_DECLARE`/`GPU_UPDATE(device=...)` |
| pre_process p(t=0) ≠ simulation p(t=0) | The two targets diverged on a leaf routine or a guard — unify per the integration map |

## P7. Physical reasoning that stays true

### Alpha (volume fraction) vs Y (mass fraction)
`alpha_j` and `Y_j` are related by `rho_j alpha_j = rho Y_j` — NOT interchangeable.
The shipped closure consumes Y only (`Y = alpha_rho_jwl/rho`). Mixing transported alpha
with a Y-derived quantity is a silent wrong-physics bug.

### JWL cold-curve behavior
`dp_cold/drho < 0` is physically possible in expansion (V > ~2 for TNT) and does NOT
imply instability — the thermal term keeps the full c² positive. Never use the cold
curve alone as a stability proxy. At V >> 1 the cold curve underflows and the closure
is effectively ideal-gas-like; the DISCRIMINATING detonation physics lives at V < 1.

### Validation regime honesty
The init self-scan covers `rho in [0.1 air_rho0, 4 rho0]` (2 rho0 for stiffened
ambient, and 2 rho0 for the offset-bearing reactant states λ<1 when Δe/=0 — unreacted
explosive only reaches the von Neumann compression, so the full 4 rho0 reflected-shock
cap tests an unreachable state), `e in [0.5 air_e0, 5 e_j]`; compressed/near-CJ states
are inside that envelope for gaseous ambients.

### Discrete-phase / particle coupling (forward-looking)
MFC's Euler-Euler bubble models (`bubbles_euler`) are PROHIBITED with `eos_jwl` —
their pressure path bypasses the closure. Any future particle-laden detonation work
(Rocflu-picl heritage: solid particles in products) must couple through drag/heat
source terms on the RHS — the same clamp-and-budget discipline as P3 — and must not
touch the EOS: particles carry their own incompressible state; only the gas phase
sees JWL. Volume-displacement effects (dense loading) would need `alpha`-consistent
treatment — revisit the alpha-vs-Y rule before any such change.

## P8. Compressible fluid-dynamics reasoning (general, beyond JWL)

Before any solver/scheme change, classify the flow the change will see:

- **Mach regime**: Ma < 0.3 incompressible-like (pressure decouples from density —
  acoustic CFL dominates dt); 0.3–1 subsonic compressible; transonic and supersonic
  flows carry shocks — any reconstruction change must be assessed at a shock, not in
  smooth flow. Detonation fronts are Ma 5–20 relative to the ambient.
- **Wave hierarchy**: every 1D compressible state change decomposes into
  acoustic (u±c), contact/material (u), and for MHD Alfvén/magnetosonic families.
  Ask "which wave does my change affect?" — a flux change that is fine on acoustics
  can still diffuse contacts (material interfaces, progress variables). Anything
  that must stay sharp at interfaces rides the HLLC contact wave.
- **CFL discipline**: dt is set by max(|u|+c) over the domain. An EOS change that
  raises c anywhere tightens dt everywhere; a c that is wrongly LOW is worse — the
  scheme runs "stable" outside its stability region locally and produces noise that
  looks like physics. When a run goes noisy after an EOS edit, suspect c first.
- **Shock relations as ground truth**: across any captured shock, verify
  Rankine-Hugoniot (mass/momentum/energy jump) from the OUTPUT fields, not the
  scheme's internals. A conservative scheme gets these right even at 3 cells per
  shock; if they're violated, a source term or non-conservative path is leaking.
- **Grid convergence protocol**: 3 grids × factor 2; compute observed order
  p = log2((f_c - f_m)/(f_m - f_f)). Shock-dominated quantities converge at ~1st
  order regardless of nominal scheme order — that is expected, not a bug. Only
  smooth-region quantities show design order.
- **Viscous/turbulent scales**: MFC's detonation work is inviscid-dominated
  (Re of a blast is enormous; boundary layers are unresolvable and not modeled).
  Do not add viscous terms "for stability" — stabilization belongs in the
  numerics (limiter, reconstruction), never in fake physics.

## P9. Particle / dispersed-phase physics (Euler-Lagrange)

MFC's dispersed-phase machinery: `m_bubbles_EE.fpp` (Euler-Euler, PROHIBITED with
`eos_jwl`), `m_bubbles_EL.fpp` + `m_bubbles_EL_kernels.fpp` (Euler-Lagrange bubbles),
`m_particle_cloud.fpp` (particle-bed IB generation). Lagrange bubbles floor
`buff_size` (see common-pitfalls). Reasoning rules for any dispersed-phase work:

- **Stokes number first**: St = tau_p/tau_f with tau_p = rho_p d²/(18 mu_g) (Stokes
  drag) and tau_f the resolved flow time scale. St << 1: particles are flow tracers
  (an Eulerian scalar suffices — cheaper and no statistical noise); St >> 1:
  ballistic (flow barely matters); St ~ 1: the interesting, expensive regime that
  actually needs Lagrangian tracking. Choosing the representation IS the physics
  decision.
- **Coupling regime by volume fraction**: alpha_p < 1e-6 one-way (flow → particle
  only); 1e-6–1e-3 two-way (add momentum/energy back-coupling source terms);
  > 1e-3 four-way (collisions matter) — MFC has no four-way machinery; dense beds
  are represented as IB solids (`m_particle_cloud.fpp`), not point particles.
- **Point-particle validity**: d_p << dx is REQUIRED. A particle comparable to the
  cell size invalidates the undisturbed-velocity assumption in the drag law; that
  regime needs resolved particles (IB), not a bigger drag coefficient.
- **Drag beyond Stokes**: post-detonation flows are high-Re_p, high-Ma_p — Stokes
  drag under-predicts by orders of magnitude. Use a compressible correlation
  (e.g. Parmar/Loth-type Ma,Re-dependent C_D) and include unsteady forces
  (pressure-gradient force scales with the fluid acceleration — in a blast this
  can EXCEED quasi-steady drag; added mass matters for bubbles, rho_p ~ rho_f).
- **Back-coupling discipline**: two-way source terms follow the SAME rules as
  reaction sources (P3): deposit to the RHS of the conserved equations, conserve
  the pair exactly (what the particle gains the gas loses), and clamp per substep
  so RK convex combinations cannot overshoot. Kernel-spread deposition
  (`m_bubbles_EL_kernels.fpp` pattern) must sum to exactly the point value —
  check the kernel normalization over ghost-cell-truncated stencils at domain
  edges.
- **Particle time stepping**: tau_p can be stiffer than the acoustic CFL. If
  tau_p < dt, subcycle or analytically integrate the drag (exponential
  integrator) — explicit Euler on stiff drag oscillates and injects momentum.
- **Statistical honesty**: N computational parcels ≠ N physical particles. Any
  parcel-weighted statistic carries O(1/sqrt(N)) noise; a "converged" mean needs a
  noise estimate before it can validate anything.

## P10. Literature anchors (which paper governs which decision)

Cite and consult the governing reference before changing the corresponding code; do
not substitute general knowledge where a canonical source exists:

- **JWL form, calibration limits, parameter restrictions**: Menikoff, "JWL Equation
  of State", LA-UR-15-29536 (LANL). Key lessons: JWL is calibrated on the principal
  isentrope from cylinder-test data — it is only trustworthy near that isentrope
  (V ∈ ~[1, 7]); reactant JWL fits use an energy offset Δe (heat of detonation
  convention) — exactly the role of `jwl_delta_e`; ω is nearly constant only over
  the calibrated range. Off-isentrope states (strong reflected shocks) are
  extrapolation — treat quantitative results there with suspicion.
- **Original JWL fit**: Lee, Hornig, Kury, UCRL-50422 (1968).
- **Five-equation model**: Kapila et al., Phys. Fluids 13, 3002 (2001) (reduction,
  detonation-to-deflagration origin); Allaire, Clerc, Kokh, JCP 181, 577 (2002)
  (interface-capturing form MFC uses). The Allaire form has NO pressure-relaxation
  K∇·u term in the volume-fraction equation — mixture-cell states are
  closure-defined, which is exactly why the Rocflu interpolation exists.
- **Rocflu closure + particle-laden detonation**: Garno, Ouellet, Bae, Jackson,
  Kim, Haftka, Hughes, Balachandar, Phys. Rev. Fluids 5, 123201 (2020) — source of
  the state-interpolated mixture closure AND the Δe reactant treatment (Eq. 17);
  its Eqs. 8-13 are the sanctioned path to a true two-phase reactant/product model.
- **Reactive burn (JWL++)**: Souers, Anderson, Mercer, McGuire, Vitello,
  Propellants Explos. Pyrotech. 25, 54 (2000) — rate dλ/dt = G p^b (1-λ); G, b are
  EMPIRICAL per explosive and grid-sensitive; recalibrate after any resolution or
  limiter change, never copy across explosives.
- **Ignition & growth (upgrade path)**: Lee & Tarver, Phys. Fluids 23, 2362 (1980)
  — the standard when hot-spot ignition physics matters (Pop-plot, shock-to-
  detonation transition). JWL++ cannot capture SDT run distance; do not tune it to.
- **Stiffened gas**: Le Métayer, Massoni, Saurel, Int. J. Therm. Sci. 43, 265
  (2004) — parameter meaning and fitting; water γ_sg and π_inf come in correlated
  pairs — never mix values from different fits.
- **Blast scaling**: Kinney & Graham, "Explosive Shocks in Air" (2nd ed., 1985) —
  far-field overpressure vs scaled distance Z = R/W^(1/3); Cole, "Underwater
  Explosions" (1948) — similitude p ∝ (W^(1/3)/R)^1.13.
- **Compressible particle forces**: Parmar, Haselbacher, Balachandar, AIAA J. 48,
  1273 (2010) and Ling, Haselbacher, Balachandar, Int. J. Multiph. Flow 37, 1026
  (2011) — in blast-particle interaction the unsteady (pressure-gradient, added-
  mass) forces can dominate quasi-steady drag; a Stokes-drag-only model is wrong
  in exactly the regime MFC targets.

When a claim in this file conflicts with one of these sources, the source wins —
update this file, citing the section.

## P11. Pre-implementation research protocol (mandatory before writing code)

Before ANY numerical/solver code is written or reviewed, work through these four
assessment vectors explicitly (in the plan or PR description, not just mentally):

1. **Mathematical formulation**: state the exact continuous PDEs touched and label
   each term convective / diffusive / source. For MFC: the 5-equation system's
   volume-fraction equation is NON-CONSERVATIVE (α advects with u) — any scheme
   change must respect that distinction, and every source term must appear in P3's
   clamp-and-budget form.
2. **Numerical characterization**: formal order of accuracy of the affected
   scheme; the binding stability constraint (acoustic CFL from max|u|+c, source
   stiffness from P3, viscous number if applicable); and discrete conservation —
   state exactly which of mass/momentum/energy the change conserves globally and
   which it deliberately does not (α, progress variables).
3. **Data layout and indexing**: trace array extents including ghost cells
   (`-buff_size:m+buff_size`, bounds structs `idwint`/`idwbuff` — see
   common-pitfalls), the `eqn_idx` positions consumed, stp vs wp kinds at every
   interface crossed, and for I/O the byte layout/precision (`mpi_io_p` ↔ stp).
4. **Failure-mode flagging**: list the specific numerical bugs the change could
   introduce — division by zero in high-gradient cells (every `1/rho` needs the
   `max(rho, sgm_eps)` pattern), non-physical oscillations at interfaces,
   positivity loss of density/α/pressure, and the silent classes in P6.

When deconstructing a paper for implementation: extract the authors' EXACT
formulas and notation (do not paraphrase equations); verify positivity/TVD/
entropy claims against the paper's own assumptions before trusting them; then map
to MFC as: governing physics → discrete formulation and stability envelope →
data/architecture fit (which modules, which eqn_idx fields) → implementation
blueprint. Divergences between the paper and what MFC can host (e.g. its
N-constituent closure vs our single JWL fluid) are stated up front, not
discovered mid-implementation.

5. **Self-cross-examination gate** (after drafting, before presenting): challenge
   the first instinct — where is the hidden trap in the math or the array stride?
   Dry-run the data flow through the loops by hand, checking ghost-cell off-by-one
   at `-buff_size`/`m+buff_size` and allocation extents. Stress the design:
   if the front speed were 10x, does the CFL/clamp logic still hold? If the
   chemistry stiffens toward infinite rate, does the source saturate at the
   `(1-x)/dt` bound or blow up? Only present after the draft survives this pass —
   and present it complete: no placeholder loops or "to be implemented" stubs in
   solver code.

These vectors complement, never override, CLAUDE.md's engineering contract:
smallest correct change, no defensive bloat, precheck before commit.

## Review priorities for JWL changes (beyond CLAUDE.md's list)

1. EOS round-trip: energy_pr(pressure_er(rho,e)) == e for any new (rho,e,Y)→p path.
2. Reference consistency: pre_process and simulation must reach the same leaf routine.
3. Region/blend gating: any new coefficient blend must keep the pressure law piecewise
   LINEAR in e (or the analytic inverse breaks) and its derivative terms must appear in c².
4. Sound-speed floor: c² floored below by `min(Gamma_amb, omega)*(p + pi_hat)/rho`;
   raw c² is returned by the state routine so the init scan can catch non-positive values.
5. Source budgets: any new reaction rate needs the `(1-x)/dt` clamp pattern so RK3
   convex combinations cannot overshoot the energy budget.
6. Pure-state limits: closure must degenerate exactly to single-phase JWL as Y→1 and
   to the ambient law as Y→0 (the init scan checks both).
