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

The `jwl_mix_type` selector is available only for five-equation JWL/ideal-gas mixtures with one JWL products fluid and one non-JWL ideal-gas fluid.

- `0`, `isobaric`: closed-form mechanical-equilibrium closure. This is the default validation path.
- `1`, `kuhl`: the supplied Kuhl/Khasainov piece-wise caloric model and additive pressure `p_air + p_products`.
- `2`, `ptequil`: the supplied `model_exact.f90` pressure-temperature equilibrium model. MFC eliminates the four unknowns in the reference matrix to a bracketed products-volume-fraction solve while retaining the same pressure, temperature, volume, and energy equations. Air uses the Kuhl caloric table shifted by `jwl_air_e0`; products require positive `cv`.
- `3`, `rocflu`: the Rocflu single-fluid closure from `modflu/RFLU_ModJWL.F90`. It interpolates the JWL coefficients with specific internal energy and the Gruneisen coefficient and heat capacity with mixture density, with ideal-air and pure-products endpoint branches.

Finite pressure, temperature, energy, and sound-speed floors are applied only after explicit finite checks. NaNs are intentionally preserved so bad states are visible during debugging instead of being converted into plausible-looking floor values.

## Validation Scope

Closure selectability is covered by registered golden tests for `jwl_mix_type = 0, 1, 2, 3`. The mode-3 regression includes a homogeneous 50/50 products-air slab so the Rocflu density and energy interpolation is exercised in addition to its endpoint branches.

Mode `3` follows Rocflu's pressure, temperature, and sound-speed formulas, but replaces its case-specific hard-coded air values and explosive energy divisor with the corresponding MFC material inputs. Its inverse energy selects the exact low-, blended-, or high-energy branch of that pressure law; this removes the legacy fallback's pressure/energy round-trip mismatch.
