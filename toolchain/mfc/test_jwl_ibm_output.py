"""Opt-in moving-wall/MPI regression: MFC_TEST_JWL_IBM_OUTPUT=1."""

import json
import math
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest


@pytest.mark.skipif(os.environ.get("MFC_TEST_JWL_IBM_OUTPUT") != "1", reason="requires MFC builds and MPI execution")
@pytest.mark.parametrize("ranks,velocity,solver", [(1, 100, 5), (2, 100, 5), (4, 100, 5), (2, 0, 5), (2, 100, 7)])
def test_equilibrium_sphere_has_no_force_across_mpi_partition(tmp_path, ranks, velocity, solver):
    root = Path(__file__).resolve().parents[2]
    # Synthetic material, common p/T and velocity. The sphere crosses the x=.5 MPI face.
    pressure = 2e8 * math.exp(-3.5) + 4e6 * math.exp(-0.8) + 0.4 * 1000 * 750 * (1500 - 300)
    ambient_density = pressure / (0.4 * 750 * 1500)
    case = {
        "m": 31,
        "n": 31,
        "p": 31,
        "dt": 2e-6,
        "t_step_start": 0,
        "t_step_stop": 20,
        "t_step_save": 20,
        "num_patches": 1,
        "model_eqns": 2,
        "num_fluids": 2,
        "relax": "T",
        "relax_model": solver,
        "mpp_lim": "F",
        "mixture_err": "F",
        "alt_soundspeed": "T",
        "time_stepper": 3,
        "weno_order": 3,
        "weno_eps": 1e-16,
        "riemann_solver": 2,
        "wave_speeds": 1,
        "avg_state": 2,
        "fd_order": 2,
        "ib": "T",
        "num_ibs": 1,
        "parallel_io": "T",
        "ib_state_wrt": "T",
        "ib_force_wrt": "T",
        "patch_icpp(1)%geometry": 9,
        "patch_icpp(1)%pres": pressure,
        "patch_icpp(1)%alpha_rho(1)": 300.0,
        "patch_icpp(1)%alpha_rho(2)": 0.7 * ambient_density,
        "patch_icpp(1)%alpha(1)": 0.3,
        "patch_icpp(1)%alpha(2)": 0.7,
        "patch_ib(1)%geometry": 8,
        "patch_ib(1)%x_centroid": 0.45,
        "patch_ib(1)%y_centroid": 0.5,
        "patch_ib(1)%z_centroid": 0.5,
        "patch_ib(1)%radius": 0.125,
        "patch_ib(1)%slip": "T",
        "patch_ib(1)%moving_ibm": int(velocity != 0),
        "patch_ib(1)%vel(1)": velocity,
        "fluid_pp(1)%eos": "jwl",
        "fluid_pp(1)%jwl_a": 2e8,
        "fluid_pp(1)%jwl_b": 4e6,
        "fluid_pp(1)%jwl_r1": 3.5,
        "fluid_pp(1)%jwl_r2": 0.8,
        "fluid_pp(1)%jwl_omega": 0.4,
        "fluid_pp(1)%jwl_rho0": 1000.0,
        "fluid_pp(1)%jwl_t0": 300.0,
        "fluid_pp(1)%cv": 750.0,
        "fluid_pp(2)%eos": "stiffened_gas",
        "fluid_pp(2)%gamma": 2.5,
        "fluid_pp(2)%pi_inf": 0.0,
        "fluid_pp(2)%cv": 750.0,
    }
    for i, axis in enumerate("xyz", 1):
        case.update({f"{axis}_domain%beg": 0.0, f"{axis}_domain%end": 1.0, f"bc_{axis}%beg": -3, f"bc_{axis}%end": -3})
        case.update({f"patch_icpp(1)%{axis}_centroid": 0.5, f"patch_icpp(1)%length_{axis}": 1.0, f"patch_icpp(1)%vel({i})": velocity if i == 1 else 0.0})
    path = tmp_path / "case.json"
    path.write_text(json.dumps(case))
    result = subprocess.run(
        [str(root / "mfc.sh"), "run", str(path), "-t", "pre_process", "simulation", "--mpi", "--no-gpu", "--no-single", "-n", str(ranks), "-j", "8"],
        cwd=root,
        capture_output=True,
        text=True,
        timeout=900,
        check=False,
        env={**os.environ, "OMP_NUM_THREADS": "1"},
    )
    assert result.returncode == 0, result.stdout + result.stderr
    force = np.loadtxt(tmp_path / "D/ib_forces.dat", ndmin=2)
    assert np.isfinite(force).all()
    assert np.max(np.linalg.norm(force[:, 2:5], axis=1)) / (pressure * np.pi * 0.125**2) < 1e-10
    assert np.max(np.linalg.norm(force[:, 5:8], axis=1)) / (pressure * np.pi * 0.125**3) < 1e-10
    initial = np.fromfile(tmp_path / "restart_data/lustre_0.dat").reshape(8, 32**3)
    final = np.fromfile(tmp_path / "restart_data/lustre_20.dat").reshape(8, 32**3)
    assert np.isfinite(final).all()
    # Zero transverse momentum needs an absolute acoustic-momentum scale.
    scale = np.maximum(np.abs(initial), 1.0)
    scale[2:5] = math.sqrt(pressure * (300 + 0.7 * ambient_density))
    assert np.max(np.abs(final - initial) / scale) < 1e-10
