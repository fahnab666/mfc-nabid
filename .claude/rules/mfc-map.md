# MFC Codebase Map — Where Physics Lives and How a Step Flows

Orientation for working anywhere in MFC beyond the JWL corner. This maps modules
to responsibilities; it does not restate CLAUDE.md (workflow/rules) or
`docs/documentation/` (authoritative API docs — freshness-checked by precheck).
When this map and the code disagree, the code wins — update this file.

## 1. The three executables and their contract

- **pre_process**: grid (`m_grid`), IC patches (`m_icpp_patches`,
  `m_patch_geometries`, `m_assign_variables`, analytic ICs compiled via Fypp),
  perturbation seeding, writes the initial conserved fields.
- **simulation**: the solver (only GPU target). Reads pre_process output,
  advances in time, writes restart/output data.
- **post_process**: reads simulation output, derives variables
  (`m_derived_variables`), writes silo/HDF5 for VisIt/ParaView.
Contract: all three share `src/common/` — an EOS or index change there must be
tested on all three; pre and sim MUST reach the same leaf routines for the same
physics (P6 t=0 pressure-mismatch class).

## 2. Simulation timestep data flow (who calls whom)

`p_main` → `m_time_steppers` (SSP-RK stages) → each stage calls `m_rhs`:
1. Ghost/boundary fill: `m_boundary_common` (+ CBC via `m_cbc`/`m_compute_cbc`)
2. cons→prim: `m_variables_conversion` (EOS dispatch lives HERE — stiffened gas,
   JWL, MHD; any new EOS hooks in here)
3. Reconstruction to faces: `m_weno` / `m_muscl` / `m_thinc` (interface
   sharpening)
4. Riemann fluxes: `m_riemann_solvers` dispatching to `_hll`/`_hllc`/`_hlld`
   (MHD)/`_lf`, shared state setup in `m_riemann_state`,
   `include/inline_riemann.fpp`
5. Source terms, each its own module: `m_jwl_sources`, `m_bubbles_EE`,
   `m_bubbles_EL`(+`_kernels`), `m_acoustic_src`, `m_body_forces`,
   `m_hypoelastic`/`m_hyperelastic`, `m_surface_tension`, `m_viscous`,
   `m_chemistry` (common), `m_qbmm`
6. Corrections after the update: `m_pressure_relaxation` (multiphase),
   `m_ibm` (`s_ibm_correct_state` — ghost-cell rebuild; see P6 known JWL issues),
   `m_phase_change`
7. `m_igr` is an alternative information-geometric regularization path (bypasses
   WENO+Riemann; incompatible with JWL).
MPI halos: `m_mpi_common` + per-target `m_mpi_proxy`; GPU_UPDATE(host) before
send / (device) after receive.

## 3. Physics feature matrix (what combines with what)

Feature flags are case parameters (~180 registered in
`toolchain/mfc/params/definitions.py`; search with `./mfc.sh params`):
- `model_eqns`: 1 (Γ/π 4-eq), 2 (5-eq Allaire — JWL lives here only), 3 (6-eq
  with p-relaxation), 4 (4-eq). Index layout (`eqn_idx`) depends on it.
- Bubbles: `bubbles_euler` (EE, sub-grid) vs `bubbles_lagrange` (EL parcels) —
  mutually distinct machinery; EE prohibited with JWL; EL floors buff_size.
- `ib` (IBM): static/moving solids, STL/model support (`m_model`,
  `m_compute_levelset`, `m_ib_patches`, `m_particle_cloud` beds).
- `mhd` (`m_riemann_solver_hlld`, powell terms), `chemistry` (multispecies,
  `m_chemistry`), elasticity (`hypo/hyperelasticity`), `surface_tension`
  (color function c), `relax` phase change, `acoustic_source`, `igr`.
- Compatibility is ENFORCED in `m_checker_common` (shared), per-target
  `m_checker`, and mirrored in `toolchain/mfc/case_validator.py` — a new
  constraint goes in BOTH (Fortran aborts at runtime, Python at `validate`).
  The JWL prohibition list is the worked example (senior-cfd-physics Part 1).

## 4. Boundary conditions

`bc_x/y/z%beg/%end` integer codes: periodic, reflective, extrapolation,
slip/no-slip walls, CBC characteristic in/outflow (`m_cbc` — subsonic
non-reflecting; prohibited with JWL), plus per-patch BC overrides
(`num_bc_patches`, `m_boundary_common`). Dirichlet-style forcing goes through
patches/sources, not hand-edited ghost cells — ghost logic lives in ONE place
(`m_boundary_common`) so all three executables agree.

## 5. Toolchain map (Python, `./mfc.sh`)

- `toolchain/mfc/params/`: definitions.py (source of truth), generators emit
  Fortran declarations/namelists/broadcasts at build time (see common-pitfalls
  for what stays manual).
- `case.py` files ARE the input format (Python dict → namelist);
  `case_validator.py` = pre-flight physics constraints; `toolchain/mfc/case.py`
  holds `QPVF_IDX_VARS` (analytic-IC variable → eqn_idx map).
- `toolchain/mfc/test/cases.py` generates the golden-file suite;
  `toolchain/mfc/lint_source.py` is the forbidden-pattern list precheck runs.
- Batch/HPC: `toolchain/templates/*.mako` (job scripts), `toolchain/modules`
  (per-cluster module slugs for `./mfc.sh load`).

## 6. Orientation heuristics (how to find anything fast)

- A physics term → grep its parameter name in `definitions.py` first (gives the
  canonical spelling), then grep that in `src/` to find the consuming module.
- An index question → read `s_initialize_eqn_idx`
  (`m_global_parameters_common.fpp`) — never assume positions.
- A "who calls this" question → the module naming convention is the call graph:
  `s_initialize_X_module` is called from the target's `m_start_up`, compute
  entry points from `m_rhs` (sim) or `m_derived_variables` (post).
- A tolerance/golden-file question → `docs/documentation/testing.md`, then
  `toolchain/mfc/test/`.
- Before assuming a capability exists or not: `./mfc.sh params <keyword>` and
  `ls examples/ | grep -i <keyword>` — examples are the living feature list.
