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

The JWL mixture closure applies to five-equation JWL/ideal-gas mixtures with one JWL products fluid and one non-JWL ideal-gas fluid. Whenever a JWL fluid is present, MFC applies the Rocflu single-fluid closure from `modflu/RFLU_ModJWL.F90`. It ramps the JWL `A` and `B` coefficients linearly in specific internal energy between the ambient-gas energy `jwl_air_e0` and the products reference energy `e_j = jwl_E0/jwl_ej_rho_ref`, ramps `omega` with mixture density via a smoothstep between `jwl_air_rho0` and `jwl_rho0`, and blends the heat capacity as the mass-weighted `cv = Y*cv_products + (1 - Y)*cv_air`. The transition from these mixture coefficients to pure-products JWL coefficients is a smoothstep in products mass fraction `Y` over `[0.95, 0.999]`, replacing the earlier hard cutoff at `Y = 0.99`; because the smoothstep weight depends only on `Y`, the analytic pressure-to-energy inverse remains exact under the blend.

This closure is not TNT-specific. MFC reads `jwl_A`, `jwl_B`, `jwl_R1`, `jwl_R2`, `jwl_omega`, `jwl_rho0`, `cv`, and either `jwl_Q` or `jwl_E0` from the JWL fluid, plus `jwl_air_rho0` and either `jwl_air_e0` or `jwl_air_p0` for the co-existing gas, so any explosive products model with a valid JWL parameter set can be used. `jwl_Q` is the specific detonation energy in J/kg; internally MFC derives `jwl_E0 = jwl_rho0*jwl_Q`. If both are provided, they must be consistent. The ambient-gas Grüneisen coefficient is taken from the ideal-gas fluid's own `gamma` (`Gamma_air = 1/gamma`; with a single JWL fluid the JWL fluid's own `gamma` is used), and the optional `jwl_ej_rho_ref` sets the products-energy reference density (default `jwl_rho0`). Unlike Rocflupicl's case-specific implementation, MFC does not hard-code TNT density, TNT energy scaling, ambient density, ambient energy, or air heat capacity inside the EOS. The reference parameters must satisfy `jwl_rho0 > jwl_air_rho0` and `jwl_E0/jwl_ej_rho_ref > jwl_air_e0`.

Finite pressure, temperature, and energy floors are applied only after explicit finite checks. The state routine returns the raw squared sound speed; the public wrappers bound it below by the ideal-gas value `Gamma_air*p/rho`. NaNs are otherwise intentionally preserved so bad states are visible during debugging instead of being converted into plausible-looking floor values.

## Validation Scope

The production closure is covered by a registered golden test. The regression includes a homogeneous 50/50 products-air slab so the Rocflu density and energy interpolation is exercised in addition to its endpoint branches.

At initialization MFC scans the closure over the configured material's `(rho, e, Y)` envelope and aborts if any state yields a non-positive or non-finite squared sound speed, or if the pressure-to-energy round-trip fails to recover the input energy. An inconsistent parameter set therefore fails fast at startup rather than producing silently wrong states during the run.

The closure follows Rocflu's pressure, temperature, inverse-energy, and sound-speed formulas, but replaces its case-specific hard-coded air values and explosive energy divisor with the corresponding MFC material inputs. Its inverse energy selects the exact low-, blended-, or high-energy branch of that pressure law; this removes the legacy fallback's pressure/energy round-trip mismatch.
