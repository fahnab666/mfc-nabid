# PR4 design notes: JWL via generalized mixture accumulators

Working notes for the nonreactive-JWL PR. Not part of the PR itself unless wanted.

## The central identity

MFC's five-equation pressure recovery is p = (rho e_int - pi_inf - qv)/gamma with

- gamma  = sum_i alpha_i/(gamma_i - 1)
- pi_inf = sum_i alpha_i gamma_i pi_inf_i/(gamma_i - 1)
- qv     = sum_i alpha_i rho_i qv_i

These are the constant-reference-curve special case of the Mie-Gruneisen mixture rule.
For MG fluids with a single cell pressure, e_i = e_ref_i + (p - p_ref_i)/(Gamma_i rho_i)
summed over phases gives

  p = [rho e_int - sum_i alpha_i (rho_i e_ref_i - p_ref_i/Gamma_i)] / (sum_i alpha_i/Gamma_i)

so the SAME formula holds with generalized accumulators:

- gamma  -> sum_i alpha_i/Gamma_i                     (identical for stiffened gas)
- pi_inf + qv -> sum_i alpha_i (rho_i e_ref_i - p_ref_i/Gamma_i)

Checked: for stiffened gas (Gamma = gamma_s - 1, p_ref = -pi_inf_i, e_ref = pi_inf_i/rho_i
+ qv_i) this reproduces MFC's pi_inf and qv accumulators exactly, term by term.

Energy (prim to cons) is the same identity inverted: energy = gamma*pres + pi_inf + qv,
also unchanged.

## Sound speed

MFC's non-alt five-equation path is c^2 = (H - 0.5|u|^2 - qv/rho)/gamma. Deriving the
frozen mixture sound speed for general MG (uniform compression at frozen alpha,
d(rho e) = rho h drho/rho along the isentrope):

  c^2 = [ h - sum_i y_i (e_ref_i + rho_i e_ref_i' - p_ref_i'/Gamma_i) ] / gamma

with y_i = alpha_i rho_i/rho and primes d/drho_i. The bracket term generalizes the qv
slot: for stiffened gas e_ref + rho e_ref' - p_ref'/Gamma = (pi/rho + qv) - pi/rho - 0
= qv_i, recovering MFC's -qv/rho term exactly. Checked against c^2 = gamma_s(p+pi)/rho
for a single stiffened fluid: exact.

For JWL the bracket is e_s + p_s/rho_p - p_s'/omega (nonzero, density dependent), so the
sound-speed path needs one generalized quantity in the qv slot:

  qv_acc = sum_i alpha_i rho_i (e_ref_i + rho_i e_ref_i' - p_ref_i'/Gamma_i)

which equals the existing qv for all-stiffened cells (bit-identical regression guard) and
carries the JWL correction otherwise. Threading: qv is already passed everywhere the
sound speed needs it (Riemann L/R states, conversion), so no new plumbing.

CAUTION: qv_acc coincides with qv in the pressure/energy formulas ONLY for stiffened
fluids. For MG fluids pressure needs sum alpha_i(rho_i e_ref_i - p_ref_i/Gamma_i) in the
pi_inf+qv slot while sound speed needs qv_acc in the qv slot; these differ. Decide in
implementation whether to (a) fold the pressure form into pi_inf and keep qv pure for
sound speed, or (b) add one explicit argument. Option (a) preferred: pi_inf absorbs
sum alpha_i(rho_i e_hat_ref_i - p_ref_i/Gamma_i) (e_hat = e_ref without qv), qv absorbs
the sound-speed bracket. Verify pi_inf is not used by the (H-...)/gamma sound path
(it is not today) and qv is not used by the pressure path other than subtraction
(it is only subtracted). Audit every consumer of gamma/pi_inf/qv mixture variables
before committing to (a).

## Where the code changes

1. m_derived_types: fluid_pp gains jwl_A, jwl_B, jwl_R1, jwl_R2, jwl_omega, jwl_rho0.
2. definitions.py: register the six per-fluid params; eos choices already include 4.
3. Defaults: dflt_real for all six.
4. m_mpi_proxy: extend the hand-listed fluid_pp broadcast loops.
5. m_eos_mie_gruneisen: add s_mg_jwl_reference (V = rho0/rho; E1 = A exp(-R1 V);
   E2 = B exp(-R2 V); p_ref = E1+E2; e_ref = E1/(rho0 R1)+E2/(rho0 R2);
   dp_ref = (V/rho)(R1 E1 + R2 E2); de_ref = p_ref/rho^2; Gamma = omega).
6. m_variables_conversion mixture-variable routines: MG branch computing the
   generalized accumulators from alpha_i, alpha_i rho_i via the providers. This is the
   whole solver-side change; Riemann/conversion/pre_process inherit it through the
   existing arguments.
7. Checker/validator: accept eos = jwl (model_eqns 2, one JWL fluid + one
   stiffened/ideal ambient, num_fluids = 2, A,B,R1,R2,omega,rho0 positive, R1 > R2;
   same incompatibility parity as the eos selector). mixture_closure enum: value 1
   pressure_equilibrium implemented; composition_weighted and kuhl reserved+rejected
   (the harness for future mixing models).
8. Cases: generic products-expansion + shock tube, synthetic coefficients; one golden.
9. misc/mg_eos_verification.py: Test C for the JWL provider (already has the oracle)
   plus mixture-rule reduction checks.

## Regression strategy

All-stiffened cells hit the constant branch and must be bit-identical (existing goldens
untouched). The MG branch activates only when a fluid has eos /= stiffened_gas, which is
impossible in every existing test (validator rejected it until this PR).

## PR4a first-draft approaches (decision record)

PR4a = the first shippable JWL PR: reference curve, the one default closure
(pressure_equilibrium), and the mixture_closure harness so later closures add without
call-site changes. Three approaches were weighed for where the closure lives.

Approach 1 (chosen): generalize the accumulators in place. The closure is implemented
inside the two existing mixture-variable routines in m_variables_conversion
(s_convert_species_to_mixture_variables at line 113 and _acc at line 181, the only two
accumulator sites on the merged base). If any fluid selects MG/JWL, the loop accumulates
the generalized slots from the per-fluid providers instead of the constant
gammas/pi_infs/qvs arrays. Pros: zero call-site changes; Riemann, conversion,
pre_process, and post inherit through the existing gamma_K/pi_inf_K/qv_K arguments; the
all-stiffened path is the untouched original loop, so bit-identity is structural, not
tested-for. Cons: needs rho_i = alpha_rho_i/max(alpha_i, sgm_eps) inside the loop, and
every consumer of the accumulator outputs must be audited for the density-dependence the
constant curves never had (increment 3).

Approach 2 (rejected): a standalone m_jwl closure module called from the pressure and
sound-speed sites, jwl_idx-gated (the old simple-tree style). Rejected because scattered
per-site dispatch is exactly what the maintainer ladder removes, and it duplicates the
centralization PR1 just landed.

Approach 3 (rejected for now): dispatch per-op inside m_thermodynamics on %eos. Rejected
because MFC's mixture closure lives in the accumulators, not the ops: the ops already
consume gamma/pi_inf/qv and stay correct once those are generalized. Revisit only if
review asks for the adapter-registry shape.

Harness sketch (the "keep other models addable" contract):

- m_constants.fpp: mixture_closure_pressure_equilibrium = 1,
  mixture_closure_composition_weighted = 2, mixture_closure_kuhl = 3.
- One global integer `mixture_closure`, default pressure_equilibrium, registered in
  definitions.py (mixing is a property of the mixture, so it is global, not per-fluid).
- The only select case sits at the top of the accumulator branch:
    select case (mixture_closure)
    case (mixture_closure_pressure_equilibrium)
        call s_mixture_pressure_equilibrium_coeffs(alpha_K, alpha_rho_K, gamma_K, pi_inf_K, qv_K)
    end select
  A future closure adds one case line and one routine with the same signature
  (alpha_i, alpha_i rho_i in; gamma/pi_inf/qv slots out); Y_i and cv_i are derivable
  from the arguments, so composition_weighted and kuhl need no signature change.
- Validator and checker reject values 2 and 3 with a named message until their PRs land
  (same reserved-value pattern PR2 uses for eos = mie_gruneisen/jwl/table).

Guard: a module logical set once at init (any fluid with eos = mie_gruneisen or jwl)
keeps stiffened-only runs on the original loop with a single well-predicted branch; the
two exp() per mixed cell only ever run when a JWL fluid is selected.

## Implementation sequence (increments, one commit each)

Order chosen so every increment builds, prechecks, and keeps the full golden suite
bit-identical until the last two, which add new cases without touching old ones.

1. Provider. `s_mg_jwl_reference(A, B, R1, R2, omega, rho0, rho, gamma_mg, p_ref,
   dp_ref, e_ref, de_ref)` in `m_eos_mie_gruneisen.fpp`, sharing the two exponentials
   across all four outputs. Extend `toolchain/mfc/test_mg_eos.py` (the Python oracle
   `jwl_reference` already exists there): provider-vs-oracle equality, isentrope
   consistency de_ref = p_ref/rho^2, FD checks. No solver change; suite untouched.
2. Parameters. `fluid_pp` members jwl_A/B/R1/R2/omega/rho0 in `m_derived_types.fpp`;
   `_r()` + `_nv()` registrations and TYPED_DECLS rows in `definitions.py`; dflt_real
   defaults in `s_assign_default_values_to_user_inputs`; the hand-listed fluid_pp
   broadcast loops in each target's `m_mpi_proxy`. Dead data until increment 4; suite
   untouched.
3. Accumulator audit + slot decision. Before code: grep every consumer of the
   gamma/pi_inf/qv mixture variables in all three targets and confirm option (a) holds
   (pi_inf never enters the (H-...)/gamma sound path; qv only ever subtracted in the
   pressure/energy path). Record the audit table here. Then generalize the
   mixture-variable routines in `m_variables_conversion`: an MG branch computes the
   generalized accumulators via the providers, guarded by `any(fluid_pp(:)%eos ==
   mie_gruneisen .or. jwl)`; the stiffened path is the untouched existing code, not a
   rewrite. Full suite must be bit-identical (the branch is unreachable).
4. Selection + harness. Checker/validator accept eos = jwl (model_eqns 2, exactly one
   JWL fluid plus one stiffened/ideal ambient, A,B,rho0 > 0, R1 > R2 > 0, omega > 0);
   global `mixture_closure` parameter, default `pressure_equilibrium`, the reserved
   values rejected at validation with a clear message; the one dispatch site in the
   mixture-variable routine. case_validator PHYSICS_DOCS entries. Suite still
   bit-identical (no case selects jwl).
5. Cases + goldens. `1D_jwl_products_expansion` (products at synthetic JWL state
   expanding into ambient air; exercises the mixed-cell closure) and
   `1D_jwl_shocktube` (single-material JWL shock tube; exercises pure-phase
   pressure/energy/sound speed). Synthetic coefficients only (the test_mg_eos.py set:
   A=5e11, B=8e9, R1=4.5, R2=1.2, omega=0.3, rho0=1600). One targeted regression test
   in `toolchain/mfc/test/cases.py` per case; generate goldens.
6. Verification ladder, recorded in the PR body: pure-phase limits exact (alpha -> 0/1
   reduces to the single-fluid formulas); full existing suite bit-identical; init-scan
   c^2 > 0 over the case envelope; analytic sound speed vs finite difference on the
   compiled path; the two new goldens stable across a rebuild.

GPU note: the accumulator loops are existing GPU_PARALLEL_LOOP regions in simulation;
the provider is seq-routine only, same as PR3, so no new data movement. The JWL branch
adds two exp() calls per mixed cell per evaluation; acceptable because the stiffened
branch is unchanged for stiffened-only runs.

Out of scope (own PRs later): composition_weighted and kuhl closures (harness values
reserved), reaction sources, named-explosive data, m_thermodynamics adapter unification
(follows PR1 review).
