@page lso-filter-testing LSO Filter Verification Plan

# LSO Filter Verification Plan

This page defines the acceptance tests for the LSO filtered statistical products and
Euler-Lagrange closure output. A successful build or a visually smooth field is not a
correctness test. Production use requires the numerical, MPI, GPU, and physics checks below.

## Supported scope

The simulation forms each nonlinear statistical product from the unfiltered state and then
filters the product. Post-process reads those simulation products. An optional post-process
pass can widen the filter by filtering the products and the gas mask again.

The 11 product blocks use `lso_R_gas` for ideal and stiffened gases and the common EOS
temperature interface for state-dependent families. Closure reconstruction is currently
limited to one calorically perfect ideal or stiffened gas. Chemistry, state-dependent
EOS closures, and stretched grids remain outside the closure acceptance scope. Serial LSO
output and `lso_down_sample_factor > 1` are supported; downsampling uses coarse-grid MPI
views and requires divisibility in every active grid direction.

## Numerical invariants

Every implementation change must verify:

1. Constant preservation: filtering a constant field returns the same field to roundoff.
2. Kernel normalization: each pass satisfies
   `a0 + 2*(a1 + a2 + a3 + a4) = 1`.
3. Design convergence: filter generation fails when the requested transfer-function error
   is not reached within `lso_max_passes`.
4. Phase partition: `phi_p + filter(g) = 1` to the output precision.
5. Product ordering: `rho_uu` equals `filter(g*rho*u*u)`, not a product formed from filtered
   density and momentum.
6. Mask normalization: filtered conserved variables equal
   `filter(g*q)/max(filter(g), lso_w_floor)`.
7. Widening equivalence: the post-process pass gives
   `filter2(filter1(g*X))` for every product and uses
   `filter2(filter1(g))` for normalization.
8. Uniform-flow limit: `R_sg`, `Q_T`, `E_ku`, `W_tau_u`, `R_mu_sg`, and `R_lam_sg` are zero
   to roundoff away from physical boundaries.
9. Manufactured nonuniform field: each closure agrees with an independent NumPy evaluation
   using the saved filter coefficients.

The negative version of the product-ordering test is required: constructing `rho_uu` after
filtering must fail the manufactured-field comparison.

## Functional cases

Run these cases with bounds checking:

```shell
./mfc.sh build -j 8 --debug
./mfc.sh run examples/1D_lso_mach3_shocktube/case.py -n 2
./mfc.sh run examples/1D_lso_mach3_shocktube/case.py -n 2 -t post_process
```

The one-dimensional case verifies nonperiodic shock-tube boundaries, nonzero sub-filter
stress near the shock, finite closure output, and post-process reading of `lso_stat` data.
The initial restart has no filtered companion file, so the initial database must contain
the ordinary flow fields and omit LSO statistics and closures. Every later saved step must
contain exactly 11 statistical blocks and the dimension-appropriate closure count.

The reduced three-dimensional IBM case must contain at least one resolved sphere and use
periodic transverse boundaries. It verifies `phi_p`, stationary `phi_p*u_p = 0`, mask
normalization, and all 3D tensor components. A moving-body variant must verify that the
saved in-situ particle velocity, rather than the case-file initial velocity, reaches
post-process output.

## MPI tests

Run the manufactured and reduced IBM cases with one, two, and four ranks. Compare every
filtered conserved field, mask, statistical product, and closure after assembling the
global arrays. The expected result is bitwise equality when the decomposition is unchanged
by I/O packing, otherwise agreement within the configured storage-precision tolerance.

Use a debug build for the two-rank shock-tube test. It specifically protects the global
QOI-layout initialization, compact MPI reads into buffered fields, four-cell filter halos,
and closure output split into communication-sized chunks.

Include decompositions that cut through the shock and through an immersed body. Verify
periodic wrapping in y and z and nonperiodic extrapolation in x independently.

## GPU and compiler matrix

The production gate consists of:

| Lane | Required result |
| --- | --- |
| gfortran CPU debug and release | Build and numerical checks pass |
| nvfortran OpenACC | Build and numerical checks pass |
| nvfortran OpenMP target | Build and numerical checks pass |
| Cray OpenACC and OpenMP target | Build and numerical checks pass |
| AMD flang OpenMP target | Build and numerical checks pass |
| Intel ifx CPU | Build and numerical checks pass |
| Case optimization, CPU and available GPU lanes | Same supported-case results |

CPU success is not evidence of GPU race freedom. Compare device output with the CPU
reference for every statistical field, including the viscous products.

## Validation tests

Unit tests must reject:

- statistics without filtered output;
- statistical output without parallel I/O;
- nonpositive `filter_sigma`;
- stretched grids;
- non-divisible LSO downsampling;
- closures without statistics;
- closures for chemistry, multiple fluids, state-dependent EOS families, or nonpositive
  `lso_R_gas`.

## Release gate

Run, in order:

```shell
./mfc.sh format -j 8
./mfc.sh precheck -j 8
./mfc.sh build -j 8 --debug
./mfc.sh test -j 8 --only <LSO test UUIDs>
```

No golden may be generated until the independent numerical checks pass. Generate only the
LSO goldens from a clean checkout of the reviewed commit. Record compiler, backend, rank
count, precision, case-optimization state, and the maximum error for each invariant.

<div style='text-align:center; font-size:0.75rem; color:#888; padding:16px 0 0;'>Page last updated: 2026-09-17</div>
