# JWL EOS Notes

This page documents the Jones-Wilkins-Lee (JWL) equation-of-state support used by the five-equation model in MFC.

## Pure JWL Products

For products density `rho`, specific internal energy `e`, and relative volume `V = rho0/rho`, MFC evaluates

```text
p = A (1 - omega/(R1 V)) exp(-R1 V)
  + B (1 - omega/(R2 V)) exp(-R2 V)
  + omega rho e.
```

The cold-curve pressure is evaluated through a shared helper so pressure recovery, energy recovery, and sound-speed code use the same expression.

## Mixture Closures

The JWL mixture closure applies to five-equation mixtures with one JWL products fluid and one non-JWL ambient fluid, which may be an ideal gas (air) or a stiffened gas (e.g. water; see "Stiffened-Gas Ambient" below). Whenever a JWL fluid is present, MFC applies the Rocflu single-fluid closure from `modflu/RFLU_ModJWL.F90`. It ramps the JWL `A` and `B` coefficients linearly in specific internal energy between the ambient-gas energy `jwl_air_e0` and the products reference energy `e_j = jwl_E0/jwl_ej_rho_ref`, ramps `omega` with mixture density via a smoothstep between `jwl_air_rho0` and `jwl_rho0`, and blends the heat capacity as the mass-weighted `cv = Y*cv_products + (1 - Y)*cv_air`. The transition from these mixture coefficients to pure-products JWL coefficients is a smoothstep in products mass fraction `Y` over `[0.95, 0.999]`, replacing the earlier hard cutoff at `Y = 0.99`; because the smoothstep weight depends only on `Y`, the analytic pressure-to-energy inverse remains exact under the blend.

This closure is not TNT-specific. MFC reads `jwl_A`, `jwl_B`, `jwl_R1`, `jwl_R2`, `jwl_omega`, `jwl_rho0`, `cv`, and either `jwl_Q` or `jwl_E0` from the JWL fluid, plus `jwl_air_rho0` and either `jwl_air_e0` or `jwl_air_p0` for the co-existing gas, so any explosive products model with a valid JWL parameter set can be used. `jwl_Q` is the specific detonation energy in J/kg; internally MFC derives `jwl_E0 = jwl_rho0*jwl_Q`. If both are provided, they must be consistent. The ambient-gas Grüneisen coefficient is taken from the ideal-gas fluid's own `gamma` (`Gamma_air = 1/gamma`; with a single JWL fluid the JWL fluid's own `gamma` is used), and the optional `jwl_ej_rho_ref` sets the products-energy reference density (default `jwl_rho0`). Unlike Rocflupicl's case-specific implementation, MFC does not hard-code TNT density, TNT energy scaling, ambient density, ambient energy, or air heat capacity inside the EOS. The reference parameters must satisfy `jwl_rho0 > jwl_air_rho0` and `jwl_E0/jwl_ej_rho_ref > jwl_air_e0`.

## Stiffened-Gas Ambient (Underwater/Condensed)

When the non-JWL fluid has `pi_inf > 0`, MFC treats it as a stiffened gas and switches the blend variable: `omega`, the `A`/`B` energy ramp, and the cold stiffness `pi_c = (Gamma + 1)*pi_hat` all blend in products mass fraction `Y` (smoothstep over [0, 1]) instead of density, because a density ramp of a stiff liquid's large Grueneisen coefficient makes the mixture sound speed non-convex. The pure-ambient limit then recovers `p = Gamma*rho*e - (Gamma + 1)*pi_inf` exactly at any density, so pure-water shocks follow the exact stiffened-gas law. The true stiffness is derived from the ambient fluid's own inputs (`pi_inf_true = fluid_pp%pi_inf/(fluid_pp%gamma + 1)`); no new case parameters are needed, and `jwl_air_p0` resolves the ambient energy through the stiffened form `e0 = (p0*gamma + pi_inf)/rho0` (MFC-convention `gamma`/`pi_inf`). Because `pi_c` is independent of both `e` and `rho`, the analytic pressure-to-energy inverse remains exact and the Grueneisen sound-speed identity is unchanged. The pressure floor doubles as the standard cavitation (pressure-cutoff) model, and the sound-speed safety floor becomes `min(Gamma, omega)*(p + pi_hat)/rho`. The ideal-gas path (`pi_inf = 0`) is bit-identical to the closure described above. Note the closure's temperature remains `T = e/cv` in the pure-ambient limit, which for a stiffened liquid counts the cold (stiffness) energy as thermal; `T` is a diagnostic on this path (pressure and sound speed are `cv`-free) and is not a physical liquid temperature. A 1D products/water shock tube (`examples/1D_jwl_underwater_shocktube`) reproduces the exact two-material Riemann star state (JWL isentrope vs. water shock Hugoniot) to under 1%.

Finite pressure, temperature, and energy floors are applied only after explicit finite checks. The state routine returns the raw squared sound speed; the public wrappers bound it below by `Gamma_air*p/rho`, a safety floor that sits below any physical mixture sound speed (pure air gives `(Gamma_air + 1)*p/rho`, pure products `(omega + 1)*p/rho` plus cold-curve terms) so it engages only on unphysical states. NaNs are otherwise intentionally preserved so bad states are visible during debugging instead of being converted into plausible-looking floor values.

## Reaction Sources

Three optional detonation-energy sources sit on top of the closure (all default off, so a plain products-expansion run is unaffected):

- **Program burn** (`prog_burn`): a Rocflu-style kinematic front expands from the detonation point (`pb_x_det`/`pb_y_det`/`pb_z_det`, time `pb_t_det`) at speed `pb_D_cj`, depositing the detonation energy `jwl_Q` over a reaction zone of width `pb_width`. No extra field is transported; each swept cell receives exactly `jwl_Q` per unit explosive mass.
- **Afterburn** (`jwl_afterburn`): products-air combustion energy `jwl_q_ab` released through an advected progress variable, either mixing-rate (`jwl_ab_model = 1`, time scale `jwl_ab_tau`) or Arrhenius (`jwl_ab_model = 2`, default; `jwl_ab_A`, `jwl_ab_theta`, `jwl_ab_n`). Requires an ideal-gas ambient. Since `jwl_q_ab` adds to `jwl_Q`, use a detonation-only JWL fit to avoid double counting.
- **Reactive burn** (`jwl_reactive`, JWL++ after Souers 2000): a self-propagating detonation driven by the local pressure, `dλ/dt = jwl_G·p^jwl_b_exp·(1−λ)`, releasing `jwl_Q` as the explosive reacts. Mutually exclusive with `prog_burn`; `jwl_G`/`jwl_b_exp` are calibrated per explosive to reproduce its CJ velocity.

## Validation Scope

The production closure is covered by registered golden tests: the Rocflu products-air closure (with a homogeneous 50/50 slab exercising the density/energy interpolation), a stiffened-gas (water) ambient, the program-burn plus afterburn sources, and the JWL++ reactive burn. Beyond the golden regression, verified physics includes: the underwater products/water shock tube against the exact two-material Riemann star state (<1%); the afterburn and JWL++ energy budgets closing to their `Y·q` bounds (<1%); the program-burn front speed matching `pb_D_cj`; and the 3D free-air TNT blast tracking the Kinney-Graham overpressure curve in the far field (`benchmarks/3D_jwl_spherical_tnt_free_air_validation/validate_kingery.py`).

At initialization MFC scans the closure over the configured material's `(rho, e, Y)` envelope and aborts if any state yields a non-positive or non-finite squared sound speed, or if the pressure-to-energy round-trip fails to recover the input energy. An inconsistent parameter set therefore fails fast at startup rather than producing silently wrong states during the run.

The closure follows Rocflu's pressure, temperature, inverse-energy, and sound-speed formulas, but replaces its case-specific hard-coded air values and explosive energy divisor with the corresponding MFC material inputs. Its inverse energy selects the exact low-, blended-, or high-energy branch of that pressure law; this removes the legacy fallback's pressure/energy round-trip mismatch.
