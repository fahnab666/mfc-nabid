---
name: run-mfc
description: Build, run, test, and smoke-check MFC (the Fortran CFD solver) — use when asked to run MFC, build the three targets, run a case.py, verify a change end-to-end, or check that the solver still works.
---

# Run MFC

MFC is a CLI scientific solver: three Fortran executables (pre_process,
simulation, post_process) driven entirely by `./mfc.sh` from the repo root.
No GUI — "running the app" means running a case end-to-end and checking the
field output. All paths below are relative to the repo root.

## Agent path: the smoke driver

```bash
.claude/skills/run-mfc/smoke.sh                    # Sod shock tube, 2 ranks
.claude/skills/run-mfc/smoke.sh examples/2D_jwl_detonation 4   # any case, N ranks
```

It builds pre_process+simulation (no-op when up to date), statically validates
the case, runs it with MPI, and exits 0 only if field-data output for multiple
time steps exists. Verified here: Sod (1000 steps) completes in ~2 s of
simulation time; a from-scratch build takes ~5–15 min.

## Build

```bash
./mfc.sh build -j $(getconf _NPROCESSORS_ONLN)     # all three targets
./mfc.sh build -t simulation -j 12                  # one target
```

First build bootstraps a Python venv and compiles dependencies — slow once,
cached in `build/` after.

## Run a case directly

```bash
./mfc.sh validate examples/1D_sodshocktube/case.py                 # seconds, no run
./mfc.sh run examples/1D_sodshocktube/case.py -t pre_process simulation -n 2
```

Progress prints per time step; a final banner shows `Exit Code: 0` and total
time. Omit `-t …` to also run post_process.

## Where output lands

Inside the case directory:
- `restart_data/lustre_<step>.dat` — conserved fields (cases with
  `parallel_io: 'T'`, which is most of them). This is the file set the smoke
  driver asserts on.
- `D/*.dat` — serial-IO cases only. Do NOT assume `D/` is populated; the Sod
  case leaves it empty.
- post_process adds silo/HDF5 under the case dir for VisIt/ParaView.

## Test suite

```bash
./mfc.sh test -j 8 --only <filter>   # targeted golden-file tests
./mfc.sh test -l                     # list all tests
```

Full suite is hours; always filter. Golden files are tolerance-compared —
an unexplained diff is a bug report, not noise (see CLAUDE.md).

## Gotchas

- `ls case/restart_data/lustre_*.dat case/D/*.dat` exits nonzero when either
  glob is empty even if the other matched — guard with `|| true` in `set -e`
  scripts (the smoke driver does).
- `./mfc.sh run` reuses stale output silently; delete `restart_data/`, `D/`,
  `p_all/` before a run you intend to assert on.
- Case parameter errors surface fastest via `./mfc.sh validate` (static) —
  don't debug them through a full run.
- Run everything via `./mfc.sh`, never cmake/compilers directly; on HPC
  clusters `source ./mfc.sh load -c <slug> -m <g|c>` first (see CLAUDE.md).

## Troubleshooting

- Smoke fails with `expected >=2 output files, found 0` after a green
  simulation banner → the case writes somewhere unexpected; check
  `parallel_io` in the case.py and look for the newest files under the case
  dir (`find <case_dir> -newer <marker>`).
- Build errors after editing `toolchain/mfc/params/*.py` → re-run
  `./mfc.sh build` (regenerates Fortran declarations); a NEW file under
  `params/` needs one cmake reconfigure (see .claude/rules/common-pitfalls.md).
