!>
!! @file m_eos_mie_gruneisen.fpp
!! @brief Contains module m_eos_mie_gruneisen

#:include 'macros.fpp'

!> @brief Generic Mie-Gruneisen equation of state, p = p_ref(rho) + Gamma*rho*(e - e_ref(rho)) with constant Gamma, split into a
!! reference-curve provider and a curve-agnostic assembler (Arienti et al., GALCIT FM99-8, 2004). Materials differ only in the
!! provider: stiffened gas ships here and the JWL principal isentrope follows with its consumer. Leaf module; no solver path calls
!! it yet.
module m_eos_mie_gruneisen

    use m_precision_select

    implicit none

    private
    public :: s_mg_stiffened_reference, s_mg_jwl_reference, f_mg_pressure, f_mg_internal_energy, f_mg_dp_drho, f_mg_dp_de, &
        & f_mg_sound_speed_sq

contains

    !> Stiffened gas as a Mie-Gruneisen reference curve: Gamma = gamma_s - 1, p_ref = -pi_inf, e_ref = pi_inf/rho, plus the density
    !! derivatives. Template for the JWL isentrope; ideal gas is the pi_inf = 0 sub-case.
    pure subroutine s_mg_stiffened_reference(gamma_s, pi_inf, rho, gamma_mg, p_ref, dp_ref, e_ref, de_ref)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: gamma_s, pi_inf, rho
        real(wp), intent(out) :: gamma_mg, p_ref, dp_ref, e_ref, de_ref

        gamma_mg = gamma_s - 1._wp
        p_ref = -pi_inf
        dp_ref = 0._wp
        e_ref = pi_inf/rho
        de_ref = -pi_inf/rho**2

    end subroutine s_mg_stiffened_reference

    !> JWL principal isentrope as a Mie-Gruneisen reference curve (Arienti et al. 2004, Eq. 39): with V = rho0/rho, p_ref =
    !! A*exp(-R1*V) + B*exp(-R2*V), e_ref its isentrope energy, Gamma = omega, and the density derivatives, sharing the two
    !! exponentials. Coefficients are user-supplied case parameters; no material data lives here.
    pure subroutine s_mg_jwl_reference(jwl_A, jwl_B, jwl_R1, jwl_R2, jwl_omega, jwl_rho0, rho, gamma_mg, p_ref, dp_ref, e_ref, &
                                       & de_ref)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: jwl_A, jwl_B, jwl_R1, jwl_R2, jwl_omega, jwl_rho0, rho
        real(wp), intent(out) :: gamma_mg, p_ref, dp_ref, e_ref, de_ref
        real(wp)              :: V, E1, E2

        V = jwl_rho0/rho
        E1 = jwl_A*exp(-jwl_R1*V)
        E2 = jwl_B*exp(-jwl_R2*V)
        gamma_mg = jwl_omega
        p_ref = E1 + E2
        dp_ref = (V/rho)*(jwl_R1*E1 + jwl_R2*E2)
        e_ref = E1/(jwl_rho0*jwl_R1) + E2/(jwl_rho0*jwl_R2)
        de_ref = p_ref/rho**2

    end subroutine s_mg_jwl_reference

    !> p = p_ref + Gamma*rho*(e - e_ref).
    pure function f_mg_pressure(gamma_mg, rho, e_int, p_ref, e_ref) result(pres)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: gamma_mg, rho, e_int, p_ref, e_ref
        real(wp)             :: pres

        pres = p_ref + gamma_mg*rho*(e_int - e_ref)

    end function f_mg_pressure

    !> Inverse of f_mg_pressure: e = e_ref + (p - p_ref)/(Gamma*rho).
    pure function f_mg_internal_energy(gamma_mg, rho, pres, p_ref, e_ref) result(e_int)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: gamma_mg, rho, pres, p_ref, e_ref
        real(wp)             :: e_int

        e_int = e_ref + (pres - p_ref)/(gamma_mg*rho)

    end function f_mg_internal_energy

    !> Thermodynamic derivative dp/drho at constant specific internal energy, in terms of pressure: dp_ref + (p - p_ref)/rho -
    !! Gamma*rho*de_ref, where dp_ref and de_ref are the reference-curve density derivatives.
    pure function f_mg_dp_drho(gamma_mg, rho, pres, p_ref, dp_ref, de_ref) result(dp_drho)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: gamma_mg, rho, pres, p_ref, dp_ref, de_ref
        real(wp)             :: dp_drho

        dp_drho = dp_ref + (pres - p_ref)/rho - gamma_mg*rho*de_ref

    end function f_mg_dp_drho

    !> Thermodynamic derivative dp/de at constant density: Gamma*rho, the defining property of the Gruneisen coefficient.
    pure function f_mg_dp_de(gamma_mg, rho) result(dp_de)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: gamma_mg, rho
        real(wp)             :: dp_de

        dp_de = gamma_mg*rho

    end function f_mg_dp_de

    !> Squared sound speed, assembled from the exposed derivatives per the frozen identity c^2 = dp/drho|_e + (p/rho^2)*dp/de|_rho
    !! (Arienti et al. 2004; Menikoff and Plohr 1989). Can be negative outside the admissible domain, so the caller tests it before
    !! taking a square root.
    pure function f_mg_sound_speed_sq(gamma_mg, rho, pres, p_ref, dp_ref, de_ref) result(c_sq)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: gamma_mg, rho, pres, p_ref, dp_ref, de_ref
        real(wp)             :: c_sq

        c_sq = f_mg_dp_drho(gamma_mg, rho, pres, p_ref, dp_ref, de_ref) + pres/rho**2*f_mg_dp_de(gamma_mg, rho)

    end function f_mg_sound_speed_sq

end module m_eos_mie_gruneisen
