# Moving IBM correction validation

These results validate commit `38f22faa` (`Fix moving IBM ghost-state coupling`).
The generated restart and Silo output files are intentionally excluded from Git.

## Results

| Case | Configuration | Acceptance result |
| --- | --- | --- |
| `paper_ibm_static` | Reduced version of the documented Mach-2 helium flow over a fixed sphere, 128 x 64 x 64, 250 steps | Passed: all fluid states finite and positive; maximum density ratio 2.6816; maximum pressure ratio 5.0056; transverse density symmetry error below 4.8e-11 |
| `ibm_comoving_periodic` | Moving no-slip sphere and uniform fluid translating together across a periodic seam, 96 x 48 x 48, 200 steps | Passed: density, pressure, and velocity errors below 1.3e-15; force norm 9.4e-19 |
| `impact.py` | Two-body collision with restitution coefficient 0.9 | Passed: zero net momentum drift; measured restitution coefficient 0.90152 |
| `daoud_six` | Six-sphere Mach-3 shock/contact reproduction, 0.8 microseconds, 680 adaptive steps | Passed: all 41 saved fluid and particle states finite; positive fluid density and pressure; no invalid solid cells or saved particle overlap |

## Commands

From the repository root:

```console
OMP_NUM_THREADS=1 ./mfc.sh run build/validation/paper_ibm_static/case.py -t pre_process simulation --no-gpu --no-build -n 4
./build/venv/bin/python3 build/validation/paper_ibm_static/check.py

OMP_NUM_THREADS=1 ./mfc.sh run build/validation/ibm_comoving_periodic/case.py -t pre_process simulation --no-gpu --no-build -n 2
./build/venv/bin/python3 build/validation/ibm_comoving_periodic/check.py

OMP_NUM_THREADS=1 ./mfc.sh run build/validation/daoud_six/case.py -t pre_process simulation --case-optimization --no-gpu --no-build -n 8
./build/venv/bin/python3 benchmarks/ibm_moving_validation/check_outputs.py build/validation/daoud_six --stop 8e-7
```

The paper-derived case is a reduced local physics/stability validation, not a grid-convergence reproduction of the full published case. The full production checkpoint still requires validation on its target HPC/GPU configuration.
