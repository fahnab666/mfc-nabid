# 3D JWL Quasi-1D Reference

This case repeats the pure-JWL shock tube from
`examples/1D_jwl_single_material_shocktube` in a thin 3D domain whose initial
data are uniform in the transverse (`y`, `z`) directions. Its purpose is a
**consistency check**: the 3D solver path (3D reconstruction and 3D `y`/`z`
Riemann fluxes) must reproduce the 1D JWL solution on the centerline and must
not introduce any transverse variation that is not present in the data.

This case verifies 3D solver-path consistency under transverse extrusion. It is
not the primary 3D validation case; for a genuinely multidimensional JWL blast
with gauge comparison see `examples/3D_jwl_spherical_blast_validation`.

It is *not* a standalone reproduction of the published exact-Riemann star state.
That is the job of the 1D case, which compares against the analytic Riemann
solution at high resolution (the self-similar star state for this problem is
`p* = 11.53 GPa`, `u* = 2026 m/s`). At the compact resolution and short time
used here the self-similar plateau is not resolved, so the centerline is
compared against the **1D MFC solution on the same grid and at the same time**,
not against the analytic star state.

## Setup

- Domain: `x in [0, 1] m`, `y, z in [0, 0.02] m`; grid `m = 399`, `n = p = 9`
  (400 x-cells, 10x10 transverse cells). The x-discretization matches the 1D
  sibling case exactly.
- Reconstruction: `recon_type = 2` (MUSCL), `muscl_order = 2`, `muscl_lim = 2`;
  `riemann_solver = 2` (HLLC); `time_stepper = 3` (RK3).
- Boundary conditions: `bc_{x,y,z}%{beg,end} = -3` (ghost-cell extrapolation) on
  all six faces.
- Time: `dt = 1e-8 s`, `t_step_stop = 600`, i.e. final time `t = 6e-6 s`,
  matching the 1D sibling.
- Initial states (two `geometry = 9` box patches spanning the full transverse
  extent), single JWL fluid `rho0 = 1770 kg/m^3`, `u = 0`:
  - Left  (`x < 0.5`): `p = 30 GPa`.
  - Right (`x > 0.5`): `p = 10 GPa`.
- JWL parameters: `A = 371.2 GPa`, `B = 3.231 GPa`, `R1 = 4.15`, `R2 = 0.95`,
  `omega = 0.30`, `rho0 = 1770 kg/m^3`, `E0 = 1.0089e10 J/m^3`.

## Run

```bash
./mfc.sh run examples/3D_jwl_quasi1d_reference/case.py -n 4
./mfc.sh run examples/3D_jwl_quasi1d_reference/case.py -n 4 -t post_process
```

Single-rank also works; with 4 ranks the domain decomposes in `x`
(`(m+1)/4 = 100 >= num_stcls_min * muscl_order`).

## Validation metrics

Comparing the 3D centerline (transverse cell nearest `y = z = 0.01`) against the
1D sibling solution on the identical `x`-grid at `t = 6e-6 s`, and measuring the
transverse spread `max_yz(.) - min_yz(.)` at each `x`:

| Quantity                         | Result        |
| -------------------------------- | ------------- |
| Centerline pressure rel. L1 vs 1D | `1.6e-16`     |
| Centerline velocity rel. L1 vs 1D | `2.0e-15`     |
| Max transverse spread, `p`        | `0` Pa        |
| Max transverse spread, `u`        | `0` m/s       |
| Max transverse spread, `v`        | `0` m/s       |
| Max transverse spread, `w`        | `0` m/s       |

**Pass expectation:** the centerline matches the 1D solution to round-off
(rel. L1 `< 1e-12`), and the transverse spreads are at or near floating-point
zero. The transverse spreads here are exactly zero because the configuration is
symmetric in `y` and `z`, so the 3D update reproduces identical values in every
transverse cell.

## Limitation

This verifies the 3D JWL solver path by preserving a published 1D reference
solution under transverse extrusion. It is not claimed as full multidimensional
blast validation.
