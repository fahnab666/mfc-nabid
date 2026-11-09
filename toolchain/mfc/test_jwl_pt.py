"""Manufactured equilibrium and conservation checks against the actual scalar Fortran solver."""

import ctypes
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from .case_validator import CaseValidator
from .params.definitions import CONSTRAINTS
from .test.case import TestCase as RegressionCase

ROOT = Path(__file__).resolve().parents[2]
# Synthetic coefficients, with no calibrated material data.
MAT = np.array([2e10, 1e8, 4.0, 1.1, 0.3, 1000.0, 0.0, 1000.0, 750.0, 0.4, 0.0, 0.0, 0.0])


@pytest.fixture(scope="module", params=[32, 64])
def solve(request, tmp_path_factory):
    compiler = shutil.which("gfortran-15") or shutil.which("gfortran")
    if not compiler:
        pytest.skip("requires GNU Fortran for the scalar EOS checks")
    bits = request.param
    folder = tmp_path_factory.mktemp(f"jwl-pt-{bits}")
    # CPU build of the production source: strip only Fypp includes and GPU annotations.
    source = (ROOT / "src/common/m_jwl_pt.fpp").read_text()
    source = "\n".join(line for line in source.splitlines() if not line.startswith("#:") and not line.lstrip().startswith("$:"))
    (folder / "core.f90").write_text(source)
    (folder / "precision.f90").write_text(f"module m_precision_select\nuse iso_fortran_env\ninteger,parameter::wp=real{bits}\nend module\n")
    (folder / "wrapper.f90").write_text(
        f"""subroutine pt_eval(x,y) bind(C)
use iso_c_binding
use m_jwl_pt
implicit none
real(c_{'float' if bits == 32 else 'double'}),intent(in)::x(17)
real(c_{'float' if bits == 32 else 'double'}),intent(out)::y(6)
type(jwl_pt_parameters)::mat
logical::ok
mat=jwl_pt_parameters(x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17))
call s_jwl_pt_equilibrium(x(1),x(2),x(3),x(4),mat,y(1),y(2),y(3),y(4),y(5),ok)
y(6)=merge(1.,0.,ok)
end subroutine
"""
    )
    library = folder / "pt.so"
    subprocess.run(
        [compiler, "-shared", "-fPIC", "-O2", "-ffree-line-length-none", "-fcheck=all", "precision.f90", "core.f90", "wrapper.f90", "-o", str(library)],
        cwd=folder,
        check=True,
        capture_output=True,
        text=True,
    )
    lib = ctypes.CDLL(str(library))
    func = lib.pt_eval
    pointer = ctypes.POINTER(ctypes.c_float if bits == 32 else ctypes.c_double)
    func.argtypes = [pointer, pointer]
    dtype = np.float32 if bits == 32 else np.float64

    def call(mj, ma, energy, mat=MAT, guess=0.5):
        args = np.array([mj, ma, energy, guess, *mat], dtype=dtype)
        result = np.empty(6, dtype=dtype)
        func(args.ctypes.data_as(pointer), result.ctypes.data_as(pointer))
        return result.astype(float)

    call.bits = bits
    call.rtol = 5e-4 if bits == 32 else 2e-8
    return call


def manufactured(rj, temp, yj, mat):
    a, b, r1, r2, omega, rho0, t0, cj, ca, ga, pi, qj, qa = mat
    ea, eb = a * np.exp(-r1 * rho0 / rj), b * np.exp(-r2 * rho0 / rj)
    tref = t0 * (rj / rho0) ** omega
    pressure = ea + eb + omega * rj * cj * (temp - tref)
    ra = (pressure + pi) / (ga * ca * temp)
    rho = 1 / (yj / rj + (1 - yj) / ra)
    mj, ma = rho * yj, rho * (1 - yj)
    ec = (ea / r1 + eb / r2) / rho0
    energy = mj * (ec + cj * (temp - tref) + qj) + ma * (ca * temp + pi / ra + qa)
    return mj, ma, energy, pressure, mj / rj


@pytest.mark.parametrize("stiffness,t0", [(0.0, 0.0), (1e7, 0.0), (0.0, 300.0), (1e7, 300.0)])
def test_manufactured_equilibria(solve, stiffness, t0):
    mat = MAT.copy()
    mat[6], mat[10] = t0, stiffness
    rng = np.random.default_rng(7231)
    for _ in range(80):
        rj = 10 ** rng.uniform(1, 3.5)
        temp = t0 * (rj / mat[5]) ** mat[4] + 10 ** rng.uniform(2.5, 4.5)
        yj = rng.uniform(0.001, 0.999)
        mj, ma, energy, pressure, alpha = manufactured(rj, temp, yj, mat)
        actual = solve(mj, ma, energy, mat)
        assert actual[5] == 1, (rj, temp, yj, actual)
        np.testing.assert_allclose(actual[:3], [pressure, temp, alpha], rtol=solve.rtol)
        assert actual[4] > 0
        assert actual[2] + actual[3] == pytest.approx(1, abs=1e-7)


def test_precision_and_pressure_unit_regressions(solve):
    mj, ma, energy, pressure, _ = manufactured(855.8975858530579, 718.6692532510846, 0.008478677689621544, MAT)
    for scale in (1.0, 1e-10, 1e5):
        mat = MAT.copy()
        mat[[0, 1, 7, 8, 10, 11, 12]] *= scale
        actual = solve(mj, ma, energy * scale, mat)
        assert actual[5] == 1
        assert actual[0] / scale == pytest.approx(pressure, rel=solve.rtol)


def test_endpoints_and_cutoff_temperature(solve):
    mat = MAT.copy()
    mat[10] = 1e7
    for yj in (0.0, 0.999e-9, 1.001e-9, 0.1, 1 - 1e-12, 1.0):
        mj, ma, energy, pressure, _ = manufactured(1e5 / (0.3 * 1000 * 300), 300.0, yj, mat)
        actual = solve(mj, ma, energy, mat)
        assert actual[5] == 1, (yj, actual)
        np.testing.assert_allclose(actual[:2], [pressure, 300], rtol=solve.rtol)


def test_equilibrium_sound_speed_matches_isentropic_pressure_derivative(solve):
    mj, ma, energy, _, _ = manufactured(800, 1800, 0.4, MAT)
    rho = mj + ma
    base = solve(mj, ma, energy)
    assert base[5] == 1
    h = 0.005 if solve.bits == 32 else 1e-5
    states = []
    for sign in (-1, 1):
        # Fixed mass fractions; de = p/rho^2 drho is the independent isentropic perturbation.
        e = energy / rho + sign * h * base[0] / rho
        result = solve(mj * (1 + sign * h), ma * (1 + sign * h), rho * (1 + sign * h) * e)
        assert result[5] == 1
        states.append(result[0])
    derivative = (states[1] - states[0]) / (2 * rho * h)
    assert base[4] == pytest.approx(derivative, rel=0.01 if solve.bits == 32 else 1e-6)


def test_no_physical_equilibrium_is_an_explicit_failure(solve):
    assert solve(1, 1, -1e6)[5] == 0


def test_bad_volume_guess_recovers_the_physical_root(solve):
    mj, ma, energy, pressure, _ = manufactured(800, 1800, 0.4, MAT)
    for guess in (1e-5, 1 - 1e-5):
        actual = solve(mj, ma, energy, guess=guess)
        assert actual[5] == 1
        assert actual[0] == pytest.approx(pressure, rel=solve.rtol)


def test_unresolved_stiff_pressure_is_not_accepted(solve):
    if solve.bits == 64:
        pytest.skip("targets the single-precision stiffness cancellation limit")
    mat = MAT.copy()
    mat[10] = 1e7
    mj, ma, energy, _, _ = manufactured(1e-6, 300, 0.5, mat)
    assert solve(mj, ma, energy, mat)[5] == 0


def test_formation_energy_is_conserved(solve):
    mat = MAT.copy()
    mat[11:] = [2e5, -1e5]
    mj, ma, energy, pressure, _ = manufactured(800, 1800, 0.4, mat)
    actual = solve(mj, ma, energy, mat)
    assert actual[5] == 1
    assert actual[0] == pytest.approx(pressure, rel=solve.rtol)


def pt_case():
    eos = CONSTRAINTS["fluid_pp(1)%eos"]["names"]
    return {
        "num_fluids": 2,
        "model_eqns": 2,
        "riemann_solver": 2,
        "wave_speeds": 1,
        "relax": "T",
        "relax_model": 5,
        "fluid_pp(1)%eos": eos["jwl"],
        "fluid_pp(1)%jwl_a": 2e10,
        "fluid_pp(1)%jwl_b": 1e8,
        "fluid_pp(1)%jwl_r1": 4.0,
        "fluid_pp(1)%jwl_r2": 1.1,
        "fluid_pp(1)%jwl_omega": 0.3,
        "fluid_pp(1)%jwl_rho0": 1000,
        "fluid_pp(1)%cv": 1000,
        "fluid_pp(2)%eos": eos["ideal_gas"],
        "fluid_pp(2)%gamma": 2.5,
        "fluid_pp(2)%cv": 750,
    }


@pytest.mark.parametrize("change", [{"relax_model": 6}, {"model_eqns": 3}, {"num_fluids": 3}, {"fluid_pp(1)%cv": 0}, {"fluid_pp(1)%jwl_b": -1}, {"ib": "T"}])
def test_validator_rejects_unsupported_pt_inputs(change):
    validator = CaseValidator(pt_case() | change)
    validator.check_eos_selector()
    assert validator.errors


def test_validator_accepts_jwl_pt():
    validator = CaseValidator(pt_case())
    validator.check_eos_selector()
    assert validator.errors == []


def test_regression_case_does_not_inherit_an_inactive_gamma():
    params = RegressionCase("JWL PT fixture", pt_case()).params
    validator = CaseValidator(params)
    validator.check_eos_selector()
    assert validator.errors == []
    assert "fluid_pp(1)%gamma" not in params


def test_regression_case_keeps_explicit_invalid_eos_inputs():
    params = RegressionCase("JWL PT fixture", pt_case() | {"fluid_pp(1)%gamma": 2.5}).params
    validator = CaseValidator(params)
    validator.check_eos_selector()
    assert any("gamma is not read" in error for error in validator.errors)
