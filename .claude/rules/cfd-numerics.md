# Numerical Methods Reference — MFC Schemes and Worked Examples

How to reason about MFC's discrete machinery. Physics equations:
`cfd-equations.md`; reasoning protocol: `senior-cfd-physics.md`. Parameter names
below are case-file parameters (`./mfc.sh params <name>` for details).

## 1. What MFC actually implements (choose from these, don't invent)

- **Reconstruction** (`recon_type`): WENO (`weno_order` = 1/3/5/7, variants
  `mapped_weno`, `wenoz` (+`wenoz_q`), `teno` (+`teno_CT`)) or MUSCL
  (`muscl_order`, `muscl_lim`). WENO-Z: sharper shocks, less dissipation than
  classic WENO-JS; TENO: even sharper, cutoff-sensitive (`teno_CT` too small →
  oscillations, too large → dissipative).
- **Riemann solver** (`riemann_solver`): HLL, HLLC, exact, and LF variants.
  HLLC restores the contact wave — REQUIRED for material interfaces, color
  function, and the JWL progress variables; HLL/LF smear contacts (and hard-set
  λ=1 on the JWL path). `wave_speeds=2` is prohibited with JWL.
- **Time stepping** (`time_stepper`): SSP-RK1/2/3. SSP property = each stage is a
  convex combination of Euler steps → anything bounded per Euler step (positivity,
  the (1-x)/dt clamps) stays bounded for the full step. Adaptive `cfl_dt` targets
  `cfl_target`; watch source terms that assume constant dt (P6 prog_burn trap).
- Only `src/simulation` is GPU-accelerated; the reconstruction and Riemann kernels
  are the hot loops — anything added inside them is O(cells × faces × stages).

## 2. Scheme-selection reasoning

- Formal order is a SMOOTH-flow property. At captured shocks every scheme is
  ~1st order; grid convergence of shock-borne quantities at ~O(Δx) is correct
  behavior. Choose higher order for the smooth regions (Taylor wave, acoustics),
  not to "sharpen the shock".
- Dissipation ladder (most → least): LF > HLL > HLLC > exact. More dissipation =
  more robust and more smeared. Debug strategy: if a case survives with LF but
  blows up with HLLC, the problem is a genuinely unphysical state, not the solver.
- Interface sharpness is set by the CONTACT treatment, not the reconstruction
  order: WENO7+HLL still smears Y across ~many cells; WENO3+HLLC keeps it tight.
- Anything advected that must correlate with Y (progress variables λ, b; color
  function) must ride the SAME contact flux — mixed treatments de-correlate the
  fields and create phantom mixture states (P7 alpha-vs-Y class).

## 3. Stability reasoning beyond the CFL number

- The acoustic CFL uses max(|u|+c) — for detonations that maximum lives AT the
  front (VN spike: both u and c peak). An EOS error that under-predicts c there
  looks like a "stability improvement" (larger allowed dt) — treat any dt gain
  after an EOS edit as a red flag, not a win (P6).
- Source stiffness is a separate limit: explicit sources need τ_source ≳ dt or a
  clamp/analytic integration (P3). The (1-x)/dt clamp converts an unconditionally
  unstable explicit rate into a saturating one — that is WHY it exists.
- Positivity: reconstruction can produce negative ρ/p at strong gradients even
  when cell averages are fine. MFC floors via sgm_eps patterns — a floor firing
  every step is not "working", it is masking an unphysical state; find the state.

## 4. Worked example A — applying the protocol to an EOS change

Task: "make the JWL A-coefficient blend quadratic in e for smoothness."
- P0: touches the closure → affects every mixture cell in every JWL case.
- P1.2 linearity gate: p would become quadratic in e → the closed-form
  (rho,p,Y)→e inverse breaks → every prim→cons conversion breaks.
- Verdict: REJECT as posed. Offer the compliant alternative: keep An linear in e,
  smooth the BLEND WEIGHT in Y or rho instead (those enter coefficients, not e),
  which preserves the analytic inverse. This rejection-with-alternative IS the
  correct output; writing the requested code would be the failure.

## 5. Worked example B — debugging with the ladder

Symptom: reactive-burn front runs 8% above D_CJ on refinement.
- P6 row "front too fast" → check clamp (present), energy double-count (none).
- P5 ladder bottom-up: init scan passes; round-trip at FP floor; star states OK
  → EOS and solver are sound; the discrepancy is in the SOURCE calibration.
- P10 Souers: G, b are empirical AND grid-sensitive — 8% fast on a finer grid is
  the documented behavior, fixed by recalibrating jwl_G at the production
  resolution, not by touching the EOS or the solver.
- [MEMORY LOG]: JWL++ rate constants are resolution-bound; recalibrate on grid
  change; never absorb rate error into EOS parameters.

## 6. Worked example C — new advected scalar checklist

Adding any new per-cell transported quantity (following eqn_idx%abn's pattern):
1. Register in `eqn_idx` (Fortran builder AND `QPVF_IDX_VARS` in toolchain
   `case.py` — mismatch = silent wrong index, see common-pitfalls).
2. Advect on the HLLC contact flux like the color function (Section 2).
3. Rebuild it in IBM ghost cells in `s_ibm_correct_state` (the shipped %rxn/%abn
   MISS this — P6 known issue; do not copy that omission).
4. Include in MPI halo exchange sizing and GPU_UPDATE around it.
5. Clamp its source with (1-x)/dt (P3); gate to physically active cells.
6. Extend the checkers (`m_checker*.fpp` + `case_validator.py`) for its
   feature-compatibility constraints.
