# Rebasing EL, JWL and LSO onto upstream master

The feature branch uses `MFlowCode/MFC` as upstream and `fahnab666/mfc-nabid`
as origin. Keep upstream changes in the branch's ancestry and feature changes
in commits above that base. Rebase onto `upstream/master` rather than merging
master repeatedly. Upstream master itself does not receive these fork features.

## Current baseline

The October 5, 2026 rebase starts at upstream commit `3dc5b2f2`.
The previous feature tip was `79508187`. Its complete feature delta was
consolidated into one commit without changing its tree before the rebase.
The original history is retained on `backup/el-jwl-iso-before-rebase-20261005`.

The one conflict was in `src/simulation/m_ibm.fpp`: retain upstream's
`Ys_IP` declaration using `NUM_SPECIES`, along with the JWL locals
`rho_IP_q` and `rho_GP_q`. Both density locals must remain private in the
GPU loop. Upstream's periodic IBM, centroid and patch-marking fixes are retained.

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

## Validation of this rebase

On macOS with GNU Fortran 15, double precision and CPU MPI:

- All seven precheck gates passed, using GNU coreutils on `PATH`.
- All three solver targets built from the rebased sources.
- All 61 non-chemistry IBM regressions passed across all three solver targets.
- The selected EOS, burn, LSO and periodic IBM regressions passed 21 of 23 cases.
  `5179D69D` and `3D70B98C`, the stiff bounded-burn cases, failed. The original
  checkout's October 1 log, `build/validation/upstream-alignment-regressions.log`,
  records the same failures and matching values before this rebase. Reactive-burn,
  RHS and time-stepper sources are unchanged from the original feature tip.
  These failures remain unresolved; no goldens were regenerated.
- A reduced two-rank EL/LSO smoke case completed with exit code zero and produced
  nonempty restart fields at steps 0, 1 and 2 and LSO fields at steps 1 and 2.
  Post-processing those fields on two ranks also completed with exit code zero.
- Upstream master is an ancestor. All 36 upstream-only changed files match it
  byte for byte, and range-diff shows only the IBM declaration adaptation.

GPU compilers and the separate chemistry configuration were not tested locally.
