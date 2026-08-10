"""Manufactured verification for the Mie-Grueneisen functions in src/common/m_helper_basic.fpp.

The backend has no golden case yet (its mixture sound-speed term and the Riemann interface
are gated on the EOS de-duplication in MFlowCode/MFC#1708), so these tests mirror the Fortran
formulas line for line and run with the toolchain unit tests. They cover the stiffened-gas
reduction of the pressure relation, the analytic sound speed against a finite-difference of the
thermodynamic identity, and the reference-curve derivative against a finite difference. Physics
follows a constant-Grueneisen Mie-Grueneisen EOS with a linear-Us Hugoniot (Menikoff and Plohr,
Rev. Mod. Phys. 61, 75, 1989). All parameters here are synthetic, not calibrated material data.
"""


def mg_eta(rho, rho0):
    return (rho - rho0) / rho


def mg_p_hugoniot(p0, rho0, c0, s, eta):
    return p0 + rho0 * c0**2 * eta / (1.0 - s * eta) ** 2


def mg_reference(rho, rho0, c0, s, p0, gamma):
    eta = mg_eta(rho, rho0)
    if eta >= 0.0:
        p_h = mg_p_hugoniot(p0, rho0, c0, s, eta)
        return p_h * (1.0 - gamma * eta / (2.0 * (1.0 - eta))) - p0 * gamma * eta / (2.0 * (1.0 - eta))
    return p0 + c0**2 * (rho - rho0)


def mg_dref_drho(rho, rho0, c0, s, p0, gamma):
    eta = mg_eta(rho, rho0)
    if eta >= 0.0:
        p_h = mg_p_hugoniot(p0, rho0, c0, s, eta)
        return c0**2 * (1.0 - eta) * (1.0 - (1.0 + gamma / 2.0) * eta) * (1.0 + s * eta) / (1.0 - s * eta) ** 3 - gamma / (2.0 * rho0) * (p_h + p0)
    return c0**2


def mg_pressure(rho, e, e0, gamma, f_ref):
    return gamma * rho * (e - e0) + f_ref


def mg_internal_energy(rho, pres, e0, gamma, f_ref):
    return e0 + (pres - f_ref) / (gamma * rho)


def mg_temperature(e, e0, cv, t0):
    return t0 + (e - e0) / cv


def mg_sound_speed_sq(rho, e, e0, pres, gamma, dref):
    return gamma * (e - e0) + dref + gamma * pres / rho


# Synthetic (gamma, rho0, c0, s, p0) reference curves and evaluation states.
MG_CASES = [
    (0.4, 1.0, 1.0, 1.0, 0.1),
    (1.0, 2.0, 1.5, 1.2, 0.2),
    (2.0, 0.5, 2.0, 0.0, 0.0),
]
STATES = [(rho, e) for rho in (0.7, 1.0, 1.6, 2.4) for e in (0.5, 2.0, 5.0)]

# (gamma_s, pi_inf, rho, e) for ideal air, a stiffened liquid, a stiff solid.
SG_CASES = [
    (1.4, 0.0, 1.2, 2.5),
    (4.4, 6.0, 3.0, 0.8),
    (3.0, 1.0, 2.0, 0.4),
]


def rel(a, b):
    return abs(a - b) / max(1.0, abs(b))


def test_stiffened_gas_reduction():
    """A constant reference curve f = -gamma_s*pi_inf recovers the stiffened-gas pressure and
    round-trips through the energy inverse. This is the off-value the mixture path reduces to."""
    for gamma_s, pi_inf, rho, e in SG_CASES:
        gamma = gamma_s - 1.0
        f_ref = -gamma_s * pi_inf
        pres = mg_pressure(rho, e, 0.0, gamma, f_ref)
        assert rel(pres, (gamma_s - 1.0) * rho * e - gamma_s * pi_inf) < 1.0e-13
        assert rel(mg_internal_energy(rho, pres, 0.0, gamma, f_ref), e) < 1.0e-13
        c2 = mg_sound_speed_sq(rho, e, 0.0, pres, gamma, 0.0)
        assert rel(c2, gamma_s * (pres + pi_inf) / rho) < 1.0e-13


def test_pressure_energy_roundtrip():
    """The pressure relation and its energy inverse agree on the reference curve."""
    for gamma, rho0, c0, s, p0 in MG_CASES:
        for rho, e in STATES:
            f_ref = mg_reference(rho, rho0, c0, s, p0, gamma)
            pres = mg_pressure(rho, e, 0.0, gamma, f_ref)
            assert rel(mg_internal_energy(rho, pres, 0.0, gamma, f_ref), e) < 1.0e-13


def test_sound_speed_identity():
    """The assembled A.38 sound speed equals the thermodynamic form dp/drho|e + p/rho^2 dp/de|rho,
    the two computed from independent finite differences of the pressure relation."""
    h = 1.0e-6
    for gamma, rho0, c0, s, p0 in MG_CASES:
        for rho, e in STATES:
            if abs(mg_eta(rho, rho0)) < 5.0e-6:  # skip the df/drho discontinuity at rho == rho0
                continue
            f_ref = mg_reference(rho, rho0, c0, s, p0, gamma)
            pres = mg_pressure(rho, e, 0.0, gamma, f_ref)
            dref = mg_dref_drho(rho, rho0, c0, s, p0, gamma)
            c2 = mg_sound_speed_sq(rho, e, 0.0, pres, gamma, dref)

            dp = rho * h
            p_hi = mg_pressure(rho + dp, e, 0.0, gamma, mg_reference(rho + dp, rho0, c0, s, p0, gamma))
            p_lo = mg_pressure(rho - dp, e, 0.0, gamma, mg_reference(rho - dp, rho0, c0, s, p0, gamma))
            dp_drho = (p_hi - p_lo) / (2.0 * dp)
            de = max(abs(e), 1.0) * h
            dp_de = (mg_pressure(rho, e + de, 0.0, gamma, f_ref) - mg_pressure(rho, e - de, 0.0, gamma, f_ref)) / (2.0 * de)
            c2_thermo = dp_drho + pres / rho**2 * dp_de
            assert rel(c2, c2_thermo) < 1.0e-6


def test_dref_vs_finite_difference():
    """df/drho matches a finite difference of the reference curve away from rho == rho0."""
    h = 1.0e-6
    for gamma, rho0, c0, s, p0 in MG_CASES:
        for rho, _ in STATES:
            if abs(mg_eta(rho, rho0)) < 5.0e-6:
                continue
            dp = rho * h
            fd = (mg_reference(rho + dp, rho0, c0, s, p0, gamma) - mg_reference(rho - dp, rho0, c0, s, p0, gamma)) / (2.0 * dp)
            assert rel(mg_dref_drho(rho, rho0, c0, s, p0, gamma), fd) < 1.0e-6


def test_admissibility():
    """The single-material sound speed stays real for admissible (positive-pressure) states;
    strong tension is legitimately outside the envelope."""
    for gamma, rho0, c0, s, p0 in MG_CASES:
        for rho, e in STATES:
            f_ref = mg_reference(rho, rho0, c0, s, p0, gamma)
            pres = mg_pressure(rho, e, 0.0, gamma, f_ref)
            if pres <= 0.0:
                continue
            dref = mg_dref_drho(rho, rho0, c0, s, p0, gamma)
            assert mg_sound_speed_sq(rho, e, 0.0, pres, gamma, dref) > 0.0


def test_temperature_reduces_to_reference():
    """Temperature is linear in internal energy at constant specific heat."""
    assert rel(mg_temperature(2.0, 0.0, 0.5, 300.0), 300.0 + 4.0) < 1.0e-13
