# JWL pressure-temperature relaxation

This uniform synthetic mixture starts with unequal constituent temperatures.
Its conserved masses and energy correspond to a known equilibrium at 1800 K.
Pre-processing and each completed simulation step relax the phase volumes while
preserving both masses, momentum and total energy.

```bash
./mfc.sh run examples/1D_jwl_pt_relaxation/case.py -n 2
```

The example uses `relax = T`, `relax_model = 5`, the five-equation model, and one
JWL fluid with an ideal or stiffened ambient. Other cases keep master's existing
EOS and relaxation paths. The hyperbolic step retains master's frozen mixture
sound speed; the scalar equilibrium sound speed is used to reject unstable PT
roots. This is a relaxation model, rather than an always-on equilibrium EOS in
the Riemann solver.

The scalar solve uses independent phase volumes, precision-dependent convergence
bounds and the same JWL caloric reference temperature as the phase EOS. A failed
solve aborts with a diagnostic instead of silently selecting a different closure.
Immersed bodies, EL coupling, elasticity and heat conduction are not supported
with this PT option.
