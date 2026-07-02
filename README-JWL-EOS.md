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

The JWL mixture closure applies to five-equation JWL/ideal-gas mixtures with one JWL products fluid and one non-JWL ideal-gas fluid. Whenever a JWL fluid is present, MFC applies the Rocflu single-fluid closure from `modflu/RFLU_ModJWL.F90`. It interpolates the JWL `A` and `B` coefficients with specific internal energy and interpolates `omega` and heat capacity with mixture density. The closure uses adapted mixture coefficients through products mass fraction `Y <= 0.99` and pure-products JWL coefficients for `Y > 0.99`.

This closure is not TNT-specific. MFC reads `jwl_A`, `jwl_B`, `jwl_R1`, `jwl_R2`, `jwl_omega`, `jwl_rho0`, `cv`, and either `jwl_Q` or `jwl_E0` from the JWL fluid, so any explosive products model with a valid JWL parameter set can be used. `jwl_Q` is the specific detonation energy in J/kg; internally MFC derives `jwl_E0 = jwl_rho0*jwl_Q`. If both are provided, they must be consistent. Unlike Rocflupicl's case-specific implementation, MFC does not hard-code TNT density, TNT energy scaling, ambient density, ambient energy, or air heat capacity inside the EOS. The reference parameters must satisfy `jwl_rho0 > jwl_air_rho0` and `jwl_E0/jwl_rho0 > jwl_air_e0`.

Finite pressure, temperature, energy, and sound-speed floors are applied only after explicit finite checks. NaNs are intentionally preserved so bad states are visible during debugging instead of being converted into plausible-looking floor values.

## Validation Scope

The production closure is covered by a registered golden test. The regression includes a homogeneous 50/50 products-air slab so the Rocflu density and energy interpolation is exercised in addition to its endpoint branches.

The closure follows Rocflu's pressure, temperature, inverse-energy, and sound-speed formulas, but replaces its case-specific hard-coded air values and explosive energy divisor with the corresponding MFC material inputs. Its inverse energy selects the exact low-, blended-, or high-energy branch of that pressure law; this removes the legacy fallback's pressure/energy round-trip mismatch.
