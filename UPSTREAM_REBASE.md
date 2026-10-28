# Rebasing EL, JWL and LSO onto upstream master

The feature branch uses `MFlowCode/MFC` as upstream and `fahnab666/mfc-nabid`
as origin. Keep upstream changes in the branch's ancestry and feature changes
in commits above that base. Rebase onto `upstream/master` rather than merging
master repeatedly. Upstream master itself does not receive these fork features.

## October 5 baseline

The October 11, 2026 rebase starts at upstream commit `ce6194eb`.
The previous feature tip was `7db05b7b`; its three commits were replayed without
squashing. The original history is retained on
`backup/el-jwl-iso-before-rebase-20261011`.

Two files conflicted. In `src/simulation/m_data_output.fpp`, retain both
upstream's `f_probe_owned` helper and the feature's `s_write_lso_stat_file`.
In `src/simulation/m_ibm.fpp`, retain upstream's reacting-surface arrays,
locals and GPU reductions, adding the JWL locals `rho_IP_q` and `rho_GP_q`
to both the declaration and the GPU private list. The burn-test and
`viz_wrt` commits replayed unchanged. Particle-cloud IB initialization retains
upstream's inert thermal and reacting-surface defaults.

LSO statistics for a state-dependent EOS now require positive heat capacity
and reference temperature. JWL statistics use the phase-temperature evaluator,
which needs `jwl_t0`; the validator previously accepted inputs without it.
Filtering without statistics retains its existing requirements. A validator
regression rejects missing or zero `jwl_t0` and accepts a positive value.

## Next rebase

Run from a clean checkout of the feature branch. Save any local work first.
Use a unique backup name for each rebase. The expected remote SHA protects
against overwriting another contributor's changes.

```bash
git switch feature/el-jwl-iso
git fetch origin
git fetch upstream
git status --short
git merge --ff-only origin/feature/el-jwl-iso
expected_remote=$(git rev-parse origin/feature/el-jwl-iso)
git branch backup/el-jwl-iso-before-next-rebase
git rebase upstream/master
```

If there are conflicts, inspect both changes, stage the resolved files and run
`git rebase --continue`. Use `git rebase --abort` to return to the saved tip.
Never select one side for an entire shared solver file without reviewing the
feature and upstream changes. Feature hooks share EOS, MPI, IBM, startup and
parameter code with upstream, so future conflict-free rebases cannot be guaranteed.

Validate before publishing:

Finish building before running tests, and do not change build configurations
while a suite is running. Nested test commands read the shared `build/lock.yaml`
and can otherwise switch to missing binaries mid-run.

On macOS with Homebrew coreutils installed, prepend
`/opt/homebrew/opt/coreutils/libexec/gnubin` to `PATH` for precheck. Upstream's
scheduler-monitor tests use GNU `stat`; BSD `stat` leaves their file-stability
polling loop waiting indefinitely.

```bash
./mfc.sh precheck -j 8
./mfc.sh build -t pre_process simulation post_process --no-gpu -j 8
./mfc.sh test --only eos=jwl --no-examples --no-build --no-gpu -j 2
./mfc.sh test --only 'Reactive Burn' --no-examples --no-build --no-gpu -j 2
./mfc.sh test --only E60E4AF8 70846B2E --no-examples --no-build --no-gpu -j 2
./mfc.sh test --only 'LSO Filter' --no-examples --no-build --no-gpu -j 2
./mfc.sh test --only IBM --no-chemistry --no-examples --no-build --no-gpu -j 2
git merge-base --is-ancestor upstream/master HEAD
git diff --check upstream/master...HEAD -- . ':!tests/*/golden-metadata.txt'
git push origin HEAD:feature/el-jwl-iso \
  --force-with-lease="refs/heads/feature/el-jwl-iso:$expected_remote"
```

Include EL particle and reactive-burn regressions when the rebase touches those
paths. Run GPU and chemistry coverage on the appropriate machines when affected.
Do not regenerate goldens to make an unexplained regression pass.
The whitespace check excludes generated golden metadata, whose existing
environment reports contain trailing spaces.

After a published history rewrite, an existing checkout may still have the old
history. Preserve local work and use a fresh checkout, or deliberately move the
clean local branch to the new remote tip. Do not merge the old feature tip back
into the rebased branch.

## Validation of the October 5 rebase

On macOS with GNU Fortran 15, double precision and CPU MPI:

- All seven precheck gates passed, using GNU coreutils on `PATH`.
- All three solver targets built from the rebased sources.
- All 61 non-chemistry IBM regressions passed across all three solver targets.
- The selected EOS, burn, LSO and periodic IBM regressions passed 21 of 23 cases.
  `5179D69D` and `3D70B98C`, the stiff bounded-burn cases, failed. Their goldens
  predate upstream's burn update: they came from a fork burn that always applied a
  bounded exponential update, while upstream integrates `rburn%substeps = 0` as an
  explicit RHS source. The two cases now set `rburn%substeps = 1`, which selects
  upstream's operator-split burn, and their goldens were regenerated from it. The
  reactant is spent in one step, the product holds the full initial mass and the
  energy is unchanged. No solver source changed.
- With that change the full CPU suite without chemistry passes 605 of 607 cases
  before the golden regeneration, the two failures being those stiff cases, and all
  12 reactive-burn cases pass after it. Every golden owned by upstream matches.
- A reduced two-rank EL/LSO smoke case completed with exit code zero and produced
  nonempty restart fields at steps 0, 1 and 2 and LSO fields at steps 1 and 2.
  Post-processing those fields on two ranks also completed with exit code zero.
- Upstream master is an ancestor. All 36 upstream-only changed files match it
  byte for byte, and range-diff shows only the IBM declaration adaptation.

GPU compilers and the separate chemistry configuration were not tested locally.

## Validation of the October 11 rebase

On macOS with GNU Fortran 15, double precision and CPU MPI:

- All seven precheck gates passed after the validator improvement, using GNU
  coreutils on `PATH`. All 115 case-validator unit tests passed; the new LSO
  regression failed before the fix and passed afterward.
- All three targets built in release mode and with GNU debug bounds checking.
- All 609 non-chemistry regression cases passed across all three targets against
  unchanged goldens. The full release run passed 608 cases in 14m 23.4s;
  `FBF808BE` failed before solver launch because the PMIx listener thread could
  not start. Its isolated retry passed in 2.7s with unchanged inputs.
- Upstream's reacting-surface regression `F52F0D4C` passed across all three
  targets with GNU debug bounds checking and unchanged goldens in 8m 22.9s.
- A fresh two-rank, two-step EL/LSO case completed pre-processing, simulation
  and post-processing. Restart fields exist at steps 0, 1 and 2, filtered fields
  at steps 1 and 2, and Silo output from post-processing.
- A fresh three-rank visualization case wrote four frames over 30 selected
  cells. Density, velocity and pressure match the restart-derived primitive
  fields bit for bit after conversion to float32. Frame times equal fixed-step
  times exactly.
- Upstream master is an ancestor. All 59 upstream-only changed files match it
  byte for byte. Before the validator and documentation follow-ups, all 61
  feature-only changed files matched the old feature tip. Range-diff records
  only the adaptations required by the two conflicts; the burn-test and
  visualization commits are unchanged.

No golden files were regenerated. GPU execution was not tested locally.
Completed regression commands:

```bash
OMP_NUM_THREADS=1 ./mfc.sh test --no-examples --no-chemistry --no-build \
  --no-gpu --no-debug --no-reldebug -a -j 8
OMP_NUM_THREADS=1 ./mfc.sh test --only FBF808BE --no-examples --no-build \
  --no-gpu --no-debug --no-reldebug -a -j 1
OMP_NUM_THREADS=1 ./mfc.sh test --only F52F0D4C --no-build \
  --no-gpu --debug --no-reldebug -a -j 1
```

Completed suite logs are `build/validation/rebase-20261011-release-regressions.log`,
`rebase-20261011-release-retry.log` and `rebase-20261011-reacting-surface-final.log`.
Build, precheck and fresh smoke artifacts share the `rebase-20261011-` prefix.

## October 27 rebase

The branch is rebased onto upstream master `8c2bec88`, including the merged
isothermal immersed-boundary implementation from PR #2038 (`d01ec1fb`). The
earlier copy of that PR, `e876419f`, is dropped because its contents are upstream.
The previous feature tip is retained locally as
`backup/el-jwl-iso-before-upbase-20261027` at `e876419f`.

The JWL same-temperature density helper remains: it supplies the moving-body
pressure correction, rather than a second prescribed-wall-temperature path.
The MPI conflict retains upstream's full-array-section arguments together with
LSO's initial-field selection. The two additional LSO MPI writes use the same
full-array-section convention as upstream's corresponding writes.

Before these MPI argument and documentation follow-ups, all 70 feature-only
files match the previous feature tip, all 49 upstream-only files match master,
and the IBM and EOS sources and isothermal golden match the previous tip.
Range-diff records only the required EOS, IBM and MPI adaptations; the burn,
visualization, LSO validation and previous documentation commits are unchanged.


Validation against `8c2bec88`, using GNU Fortran 15 and CPU MPI:

- All seven precheck gates passed. All three targets built in release and debug.
- The complete release suite without examples or chemistry passed 613 of 614
  cases across all three targets in 22m 8.6s. `B43894C4`, the moving JWL circle,
  has a tolerance mismatch in both release and debug: the largest absolute
  discrepancy among failing values is `7.97e-4` in total energy. Its EOS and IBM
  source files match the previous feature tip exactly. The cause of the mismatch
  has not been established; no golden files were regenerated.
- Upstream's isothermal regression `5BBDFB29` passed in release and debug.
- Fresh two-rank EL/LSO pre-processing, simulation and post-processing completed.
  Restart fields at steps 0, 1 and 2, filtered fields at steps 1 and 2, and Silo
  output were verified.
- Fresh three-rank visualization produced four frames over 30 selected cells.
  Primitive fields match restart-derived float32 values bit for bit, and frame
  times match the fixed-step times exactly.

Logs and smoke artifacts use the `build/validation/upbase-20261027-` prefix.
GPU execution and chemistry were not tested.


## October 28 upstream follow-up

The final rebase base is upstream master `732209d0`. This adds upstream's
second-order IBM implementation from PR #1958 (`4e96e0da`) and its coverage-map
refresh. The rebase applied without conflicts. The feature delta relative to
this base retains the JWL moving-body density correction, while the isothermal
implementation continues to come from upstream PR #2038. The older copy
`e876419f` is absent from the final history.


Validation of the final base:

- All three release targets rebuilt, and all seven precheck gates passed after
  formatting the added MPI call.
- The affected non-chemistry IBM suite passed 63 of 64 cases across all three
  targets in 4m 12.2s, including isothermal case `5BBDFB29` and both new
  second-order cases, `731ACB1A` and `284ABF4C`. The sole failure remains
  `B43894C4`, with the same numerical discrepancy reported above.
- All 55 upstream-only files match master. Of 65 feature-only files, only this
  rebase report differs from the previous feature tip. The older isothermal
  copy is excluded, and the merged upstream implementation is an ancestor.
- No golden files were regenerated locally.

Final build, precheck, IBM and integrity logs use the
`build/validation/upbase-20261028-` prefix.
