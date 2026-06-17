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
- `1`, `kuhl`: Kuhl/Khasainov temperature-form additive closure. It requires positive `cv` on both products and air.
- `2`, `ptequil`: pressure-temperature equilibrium closure. It solves a bounded scalar root for the products volume fraction and is substantially more expensive than mode `0`.
- `3`, `rocflu`: Garno/Rocflu-style single-fluid blend. Its sound speed is evaluated by the Rocflu-specific Gruneisen form rather than by phasic volume fractions.

Finite pressure, temperature, energy, and sound-speed floors are applied only after explicit finite checks. NaNs are intentionally preserved so bad states are visible during debugging instead of being converted into plausible-looking floor values.

## Validation Scope

The exact-reference validation in this branch is scoped to five-equation JWL cases: a 1D pure-JWL shock tube and a compact 3D quasi-1D repeat of the same published Shyue-style Riemann reference. Closure selectability is covered by registered golden tests for `jwl_mix_type = 0, 1, 2, 3`.

The p-T equilibrium mode remains selectable and physically distinct. Its root-find cost is measured by `benchmarks/jwl_closure_modes`; it is not optimized away by the case-optimization gates.
