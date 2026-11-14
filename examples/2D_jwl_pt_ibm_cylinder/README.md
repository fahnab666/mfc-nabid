# 2D JWL PT relaxation with an immersed cylinder

A high-pressure JWL-products slab drives a shock over an immersed cylinder with
JWL/ambient pressure-temperature relaxation (`relax = T`, `relax_model = 5`).
The PT projection runs after each full RK step on fluid cells only; cells
inside the body and on its ghost layer (`ib_markers /= 0`) are skipped.

```bash
MFC_IBM_MODE=static     ./mfc.sh run examples/2D_jwl_pt_ibm_cylinder/case.py -n 4
MFC_IBM_MODE=prescribed ./mfc.sh run examples/2D_jwl_pt_ibm_cylinder/case.py -n 4
MFC_IBM_MODE=free       ./mfc.sh run examples/2D_jwl_pt_ibm_cylinder/case.py -n 4
```

`MFC_IBM_MODE` selects `moving_ibm = 0` (fixed), `1` (prescribed 100 m/s in x)
or `2` (free flight, mass 20). Moving JWL bodies need positive `cv` and
`jwl_t0`, which this case sets.

Coefficients are synthetic (as in `examples/1D_jwl_pt_relaxation`). This is a
solver-coupling verification case, not a calibrated material or a physical
validation. Use `relax = F` (and drop `relax_model`) for the standard
mixture-closure baseline on the same initial state.
