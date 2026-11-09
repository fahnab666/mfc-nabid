"""Opt-in, full-stage MPI checks for JWL PT relaxation: MFC_TEST_JWL_PT_OUTPUT=1."""

import contextlib
import io
import json
import os
import re
import runpy
import subprocess
from pathlib import Path

import numpy as np
import pytest


@pytest.mark.skipif(os.environ.get("MFC_TEST_JWL_PT_OUTPUT") != "1", reason="requires MFC builds and MPI execution")
@pytest.mark.parametrize("single,ranks,swapped", [(False, 1, False), (False, 2, True), (True, 2, False)])
def test_nonuniform_pt_preserves_conserved_fields_and_equalizes_temperature(tmp_path, single, ranks, swapped):
    from .viz.silo_reader import read_silo_file

    root = Path(__file__).resolve().parents[2]
    with contextlib.redirect_stdout(io.StringIO()) as output:
        runpy.run_path(str(root / "examples/1D_jwl_pt_relaxation/case.py"))
    case = json.loads(output.getvalue())
    case["precision"] = "single" if single else "double"
    case.update({"num_patches": 2, "bc_x%beg": -2, "bc_x%end": -2})
    patch = {k: v for k, v in case.items() if k.startswith("patch_icpp(1)")}
    case.update({k.replace("patch_icpp(1)", "patch_icpp(2)"): v for k, v in patch.items()})
    case.update({"patch_icpp(1)%x_centroid": 0.25, "patch_icpp(2)%x_centroid": 0.75})
    case.update({"patch_icpp(1)%length_x": 0.5, "patch_icpp(2)%length_x": 0.5})
    case["patch_icpp(2)%pres"] *= 1.1
    if swapped:

        def swap(match):
            return match[1] + str(3 - int(match[2])) + match[3]

        case = {re.sub(r"(fluid_pp\(|%alpha(?:_rho)?\()(1|2)(\))", swap, k): v for k, v in case.items()}
    path = tmp_path / "case.json"
    path.write_text(json.dumps(case))
    result = subprocess.run(
        [str(root / "mfc.sh"), "run", str(path), "--mpi", "--no-gpu", "--no-debug", "--no-reldebug", "--single" if single else "--no-single", "-n", str(ranks), "-j", "8"],
        cwd=root,
        capture_output=True,
        text=True,
        timeout=900,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    initial = np.fromfile(tmp_path / "restart_data/lustre_0.dat", dtype=np.float32 if single else np.float64).reshape(6, 64).astype(float)
    jidx, aidx = (1, 0) if swapped else (0, 1)
    rtol = 5e-4 if single else 2e-8
    for step in (0, 5, 10):
        q = np.fromfile(tmp_path / f"restart_data/lustre_{step}.dat", dtype=np.float32 if single else np.float64).reshape(6, 64).astype(float)
        assert np.isfinite(q).all()
        rho = q[0] + q[1]
        aj, aa = q[4 + jidx], q[4 + aidx]
        rj, ra = q[jidx] / aj, q[aidx] / aa
        ej = case[f"fluid_pp({jidx + 1})%jwl_a"] * np.exp(-4 * 1000 / rj)
        eb = 1e8 * np.exp(-1.1 * 1000 / rj)
        ecold = (ej / 4 + eb / 1.1) / 1000
        gamma = aj / 0.3 + aa / 0.4
        pi = q[jidx] * ecold - aj * (ej + eb) / 0.3 + aa * 3.5e7
        pressure = (q[3] - q[2] ** 2 / (2 * rho) - pi) / gamma
        tj = 300 * (rj / 1000) ** 0.3 + (pressure - ej - eb) / (0.3 * rj * 1000)
        ta = (pressure + 1e7) / (0.4 * ra * 750)
        np.testing.assert_allclose(tj, ta, rtol=rtol)
        np.testing.assert_allclose(aj + aa, 1, rtol=0, atol=1e-7)
        # Reflecting walls: the hyperbolic update and PT projection conserve both masses and total energy.
        np.testing.assert_allclose(q[[0, 1, 3]].sum(axis=1), initial[[0, 1, 3]].sum(axis=1), rtol=5e-6 if single else 1e-12)
        for rank in range(ranks):
            fields = read_silo_file(str(tmp_path / f"silo_hdf5/p{rank}/{step}.silo")).variables
            np.testing.assert_allclose(fields["T1"], fields["T2"], rtol=rtol)
    assert np.ptp(q[2]) > 1e-3, "test must evolve a nonuniform state"
