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
real(c_{'float' if bits == 32 else 'double'}),intent(in)::x(18)
real(c_{'float' if bits == 32 else 'double'}),intent(out)::y(6)
type(jwl_pt_parameters)::mat
logical::ok
mat=jwl_pt_parameters(x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17))
if (nint(x(18))==7) then
 call s_jwl_pt_equilibrium_matrix(x(1),x(2),x(3),x(4),mat,y(1),y(2),y(3),y(4),y(5),ok)
else
 call s_jwl_pt_equilibrium(x(1),x(2),x(3),x(4),mat,y(1),y(2),y(3),y(4),y(5),ok)
end if
y(6)=merge(1.,0.,ok)
end subroutine
subroutine cw_eval(x,y) bind(C)
use iso_c_binding
use m_jwl_pt
implicit none
real(c_{'float' if bits == 32 else 'double'}),intent(in)::x(18)
real(c_{'float' if bits == 32 else 'double'}),intent(out)::y(2)
type(jwl_pt_parameters)::mat
logical::ok
mat=jwl_pt_parameters(x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17))
call s_pt_cw_guess(x(1),x(2),x(3),mat,y(1),ok)
y(2)=merge(1.,0.,ok)
end subroutine

subroutine wall_eval(x,y) bind(C)
use iso_c_binding
use m_jwl_pt
implicit none
real(c_{'float' if bits == 32 else 'double'}),intent(in)::x(16)
real(c_{'float' if bits == 32 else 'double'}),intent(out)::y(2)
type(jwl_pt_parameters)::mat
logical::ok
mat=jwl_pt_parameters(x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16))
call s_jwl_wall_density(x(1),x(2),x(3),mat,y(1),ok)
y(2)=merge(1.,0.,ok)
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
    cw_func = lib.cw_eval
    cw_func.argtypes = [pointer, pointer]
    dtype = np.float32 if bits == 32 else np.float64

    def call(mj, ma, energy, mat=MAT, guess=0.5, solver="fast"):
        args = np.array([mj, ma, energy, guess, *mat, 7 if solver == "matrix" else 5], dtype=dtype)
        result = np.empty(6, dtype=dtype)
        func(args.ctypes.data_as(pointer), result.ctypes.data_as(pointer))
        return result.astype(float)

    call.bits = bits
    call.rtol = 5e-4 if bits == 32 else 2e-8

    def cw_guess(mj, ma, energy, mat=MAT):
        args = np.array([mj, ma, energy, 0.5, *mat, 5], dtype=dtype)
        result = np.empty(2, dtype=dtype)
        cw_func(args.ctypes.data_as(pointer), result.ctypes.data_as(pointer))
        return result.astype(float)

    call.cw_guess = cw_guess

    wall_func = lib.wall_eval
    wall_func.argtypes = [pointer, pointer]

    def wall(rho, p_from, p_to, mat=MAT):
        args = np.array([rho, p_from, p_to, *mat], dtype=dtype)
        result = np.empty(2, dtype=dtype)
        wall_func(args.ctypes.data_as(pointer), result.ctypes.data_as(pointer))
        return result.astype(float)

    call.wall = wall
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


def test_wall_density_manufactured_isotherms(solve):
    """The production inverse recovers known densities without a pressure-ratio surrogate."""
    rng = np.random.default_rng(20261116)
    checked = 0
    for t0 in (0.0, 300.0):
        mat = MAT.copy()
        mat[6] = t0
        a, b, r1, r2, omega, rho0, _, cv = mat[:8]

        def pressure(rho, temp):
            return a * np.exp(-r1 * rho0 / rho) + b * np.exp(-r2 * rho0 / rho) + omega * rho * cv * (temp - t0 * (rho / rho0) ** omega)

        for _ in range(200):
            rho = rho0 * 10 ** rng.uniform(-1.3, 1)
            target = rho * rng.uniform(0.5, 2)
            temp = t0 * (max(rho, target) / rho0) ** omega + a / (omega * min(rho, target) * cv)
            p_from, p_to = pressure(rho, temp), pressure(target, temp)
            actual = solve.wall(rho, p_from, p_to, mat)
            assert actual[1] == 1, (rho, target, temp, actual)
            assert actual[0] == pytest.approx(target, rel=solve.rtol)
            checked += 1
    assert checked == 400


@pytest.mark.parametrize("t0,rho,target,temp", [(666.6666666666666, 2000.0, 1000.0, 1213.005273848596), (0.0, 10000.0, 5000.0, 2 / 3)])
def test_wall_density_stage_one_counterexamples(solve, t0, rho, target, temp):
    mat = MAT.copy()
    mat[:8] = [2e8, 4e6, 3.5, 0.8, 0.4, 1000, t0, 750]

    def pressure(r):
        return mat[0] * np.exp(-mat[2] * mat[5] / r) + mat[1] * np.exp(-mat[3] * mat[5] / r) + mat[4] * r * mat[7] * (temp - t0 * (r / mat[5]) ** mat[4])

    actual = solve.wall(rho, pressure(rho), pressure(target), mat)
    assert actual[1] == 1
    assert actual[0] == pytest.approx(target, rel=solve.rtol)


@pytest.mark.parametrize("rho,p_from,p_to", [(0.0, 1.0, 1.0), (1.0, -1.0, 1.0), (1.0, 1.0, 0.0), (1.0, float("nan"), 1.0)])
def test_wall_density_rejects_invalid_states(solve, rho, p_from, p_to):
    assert solve.wall(rho, p_from, p_to)[1] == 0


def test_wall_density_rejects_unreachable_pressure(solve):
    mat = MAT.copy()
    mat[:8] = [2e8, 4e6, 3.5, 0.8, 0.4, 1000.0, 666.6666666666666, 750.0]
    assert solve.wall(2000.0, 237436068.87423152, 2 * 237436068.87423152, mat)[1] == 0


@pytest.mark.parametrize("stiffness,solver", [(0.0, "fast"), (1e7, "fast"), (0.0, "matrix"), (1e7, "matrix")])
def test_reconstructed_wall_mixture_is_a_pt_fixed_point(solve, stiffness, solver):
    """The two-phase ghost state must retain common pressure/temperature under the branch's PT rule."""
    mat = MAT.copy()
    mat[:8] = [2e8, 4e6, 3.5, 0.8, 0.4, 1000.0, 300.0, 750.0]
    mat[10] = stiffness
    a, b, r1, r2, omega, rho0, t0, cv = mat[:8]
    for rho in (500.0, 1000.0, 2000.0):
        for ratio in (0.5, 1.0, 1.5):
            target = rho * ratio
            temp = t0 * (max(rho, target) / rho0) ** omega + 2000.0

            def pressure(r):
                return a * np.exp(-r1 * rho0 / r) + b * np.exp(-r2 * rho0 / r) + omega * r * cv * (temp - t0 * (r / rho0) ** omega)

            pres = pressure(target)
            wall = solve.wall(rho, pressure(rho), pres, mat)
            assert wall[1] == 1
            for alpha in (0.05, 0.5, 0.95):
                ra = (pres + stiffness) / (mat[9] * mat[8] * temp)
                mj, ma = alpha * wall[0], (1 - alpha) * ra
                ec = (a * np.exp(-r1 * rho0 / wall[0]) / r1 + b * np.exp(-r2 * rho0 / wall[0]) / r2) / rho0
                energy = mj * (ec + cv * (temp - t0 * (wall[0] / rho0) ** omega)) + ma * (mat[8] * temp + stiffness / ra)
                actual = solve(mj, ma, energy, mat, alpha, solver)
                assert actual[5] == 1, (rho, ratio, alpha, stiffness, solver, actual)
                np.testing.assert_allclose(actual[:3], [pres, temp, alpha], rtol=max(solve.rtol, 1e-4 if solver == "fast" else 0.0))
                rj_eq, ra_eq = mj / actual[2], ma / actual[3]
                ec_eq = (a * np.exp(-r1 * rho0 / rj_eq) / r1 + b * np.exp(-r2 * rho0 / rj_eq) / r2) / rho0
                energy_eq = mj * (ec_eq + cv * (actual[1] - t0 * (rj_eq / rho0) ** omega)) + ma * (mat[8] * actual[1] + stiffness / ra_eq)
                assert energy_eq == pytest.approx(energy, rel=solve.rtol)


def test_fast_solver_bounds_the_constituent_pressure_difference(solve):
    """Every Fast PT return, the CW gate included, leaves |p_j - p_a| within the gate's tau**2 of the pressure."""
    a, b, r1, r2, omega, rho0, _, cj, ca, ga = MAT[:10]
    rng = np.random.default_rng(20261116)
    bound = max(1e-8, 64 * np.finfo(np.float32 if solve.bits == 32 else np.float64).eps)
    for _ in range(2000):
        rj, temp, yj = 10 ** rng.uniform(-0.5, 3), 10 ** rng.uniform(2.4, 3.7), 10 ** rng.uniform(-4, 0) * rng.uniform(0.1, 1)
        mj, ma, energy, _, _ = manufactured(rj, temp, yj, MAT)
        pres, temp, aj, aa, _, ok = solve(mj, ma, energy, MAT, 0.0)
        assert ok == 1
        pj = a * np.exp(-r1 * rho0 * aj / mj) + b * np.exp(-r2 * rho0 * aj / mj) + omega * mj / aj * cj * temp
        assert abs(pj - ga * ma / aa * ca * temp) <= 4 * bound * pres


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


def test_matrix_solver_manufactured_equilibria(solve):
    rng = np.random.default_rng(815)
    for _ in range(30):
        rj = 10 ** rng.uniform(1, 3.5)
        temp = 300 * (rj / MAT[5]) ** MAT[4] + 10 ** rng.uniform(2.5, 4.2)
        yj = rng.uniform(0.05, 0.95)
        mj, ma, energy, pressure, alpha = manufactured(rj, temp, yj, MAT)
        actual = solve(mj, ma, energy, MAT, alpha * 0.8, solver="matrix")
        assert actual[5] == 1, (rj, temp, yj, actual)
        np.testing.assert_allclose(actual[:3], [pressure, temp, alpha], rtol=solve.rtol)


def test_fast_and_matrix_solvers_agree(solve):
    rng = np.random.default_rng(191)
    errors = np.zeros(4)
    for _ in range(40):
        mat = MAT.copy()
        mat[6] = rng.choice([0.0, 300.0])
        mat[10] = rng.choice([0.0, 1e7])
        mat[11:] = rng.uniform(-2e5, 2e5, 2)
        rj = 10 ** rng.uniform(1.5, 3.2)
        temp = mat[6] * (rj / mat[5]) ** mat[4] + 10 ** rng.uniform(2.7, 3.8)
        yj = rng.uniform(0.03, 0.97)
        mj, ma, energy, _, _ = manufactured(rj, temp, yj, mat)
        fast = solve(mj, ma, energy, mat, solver="fast")
        matrix = solve(mj, ma, energy, mat, solver="matrix")
        assert fast[5] == matrix[5] == 1, (rj, temp, yj, fast, matrix)
        errors = np.maximum(errors, np.abs(fast[:4] - matrix[:4]) / np.maximum(np.abs(fast[:4]), 1.0))
    assert np.max(errors) <= solve.rtol, errors


def test_matrix_solver_resolves_pressure_in_water_like_ambient(solve):
    if solve.bits == 32:
        pytest.skip("p << pi_inf is below single-precision resolution")
    # pi_inf >> p makes rho e mostly ambient cold energy; the tolerance must still be relative to p.
    mat = MAT.copy()
    mat[8], mat[9], mat[10] = 588.0, 3.4, 6e8
    rng = np.random.default_rng(1116)
    for _ in range(40):
        rj, temp, yj = 10 ** rng.uniform(1, 3.2), 10 ** rng.uniform(2.5, 3.6), rng.uniform(0.05, 0.95)
        mj, ma, energy, pressure, alpha = manufactured(rj, temp, yj, mat)
        actual = solve(mj, ma, energy, mat, alpha * 0.8, solver="matrix")
        assert actual[5] == 1, (rj, temp, yj, actual)
        assert actual[0] == pytest.approx(pressure, rel=1e-9)
        assert actual[2] + actual[3] == pytest.approx(1, abs=4 * np.finfo(float).eps)
        # A common shift of the formation energies is the same physical state.
        shifted = mat.copy()
        shifted[11:] = 5e6
        moved = solve(mj, ma, energy + (mj + ma) * 5e6, shifted, alpha * 0.8, solver="matrix")
        assert moved[5] == 1 and moved[0] == pytest.approx(pressure, rel=1e-9)


def test_composition_weighted_guess_matches_ideal_limit(solve):
    mat = MAT.copy()
    mat[[0, 1, 6, 10]] = 0.0
    for yj in (0.01, 0.2, 0.8, 0.99):
        mj, ma, energy, _, alpha = manufactured(800, 1800, yj, mat)
        guess, valid = solve.cw_guess(mj, ma, energy, mat)
        assert valid == 1
        assert guess == pytest.approx(alpha, rel=solve.rtol)


def test_composition_weighted_guess_rejects_unphysical_state(solve):
    assert solve.cw_guess(1, 1, -1e6)[1] == 0


def test_composition_weighted_guess_with_offsets_and_stiffened_gas(solve):
    mat = MAT.copy()
    mat[6], mat[10], mat[11], mat[12] = 300.0, 1e7, 2e5, -1e5
    mj, ma, energy, _, _ = manufactured(800, 1800, 0.4, mat)
    rho, e = mj + ma, energy / (mj + ma)
    # Coefficients blended by the heat-capacity share w; the JWL reference temperature enters through the cold curve.
    weight = mj * mat[7] / (mj * mat[7] + ma * mat[8])
    omega, pihat = mat[9] + weight * (mat[4] - mat[9]), (1 - weight) * mat[10]
    e -= (mj * mat[11] + ma * mat[12]) / rho
    v, tref = mat[5] / rho, mat[6] * (rho / mat[5]) ** mat[4]
    ea, eb = mat[0] * np.exp(-mat[2] * v), mat[1] * np.exp(-mat[3] * v)
    pcw = weight * (ea * (1 - omega / (mat[2] * v)) + eb * (1 - omega / (mat[3] * v)) + (omega - mat[4]) * rho * mat[7] * tref) + omega * rho * e - (mat[9] + 1) * pihat
    tcw = (e - weight * ((ea / mat[2] + eb / mat[3]) / mat[5] - mat[7] * tref) - pihat / rho) * rho / (mj * mat[7] + ma * mat[8])
    expected = 1 - ma * mat[9] * mat[8] * tcw / (pcw + mat[10])
    guess, valid = solve.cw_guess(mj, ma, energy, mat)
    assert valid == 1
    assert guess == pytest.approx(expected, rel=solve.rtol)


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


@pytest.mark.parametrize(
    "change",
    [
        {"relax_model": 6},
        {"model_eqns": 3},
        {"num_fluids": 3},
        {"fluid_pp(1)%cv": 0},
        {"fluid_pp(1)%jwl_b": -1},
        {"particles_lagrange": "T"},
        {"ib": "T", "num_ibs": 1, "patch_ib(1)%moving_ibm": 2},
    ],
)
def test_validator_rejects_unsupported_pt_inputs(change):
    validator = CaseValidator(pt_case() | change)
    validator.check_eos_selector()
    assert validator.errors


def test_validator_accepts_jwl_pt():
    for mode in (5, 7):
        validator = CaseValidator(pt_case() | {"relax_model": mode})
        validator.check_eos_selector()
        assert validator.errors == []


@pytest.mark.parametrize(
    "change",
    [{"ib": "T"}, {"ib": "T", "num_ibs": 1, "patch_ib(1)%moving_ibm": 1}, {"ib": "T", "num_ibs": 1, "patch_ib(1)%moving_ibm": 2, "fluid_pp(1)%jwl_t0": 300.0}],
)
def test_validator_accepts_jwl_pt_with_immersed_bodies(change):
    """Stationary, prescribed and free-flight bodies are allowed; a free-flight body also needs the JWL reference temperature."""
    validator = CaseValidator(pt_case() | change)
    validator.check_eos_selector()
    assert validator.errors == []


@pytest.mark.parametrize("change", [{"relax": "F"}, {"fluid_pp(1)%eos": CONSTRAINTS["fluid_pp(1)%eos"]["names"]["ideal_gas"]}])
def test_validator_reserves_matrix_mode_for_jwl_relaxation(change):
    validator = CaseValidator(pt_case() | {"relax_model": 7} | change)
    validator.check_eos_selector()
    assert validator.errors


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
