!>
!! @file
!! @brief JWL EOS and Rocflu state-interpolated mixture closure.

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief JWL equation of state and Rocflu state-interpolated mixture closure.
!!
!! Closure path: (rho, e, Y) -> (p, T, c), with Y the JWL products mass fraction.
!! Mixture rule after Garno et al., Phys. Rev. Fluids 5, 123201 (2020) and the Rocflu
!! implementation in modflu/RFLU_ModJWL.F90: the JWL A and B coefficients ramp linearly
!! from zero at e = air_e0 to full strength at e = e_j = E0/rho0, while omega and cv ramp
!! linearly with mixture density between air and products. The blended coefficients are
!! used up to 99% products mass fraction; above that the pure-products JWL law is used.
module m_jwl

    use m_global_parameters

    implicit none

    private
    public :: s_initialize_jwl_module, s_finalize_jwl_module, s_jwl_mix_state_er, s_jwl_mix_energy_pr, s_jwl_mix_sound_speed, &
        & jwl_idx

    ! Simulation builds use m_global_parameters tables.
#ifndef MFC_SIMULATION
    real(wp), allocatable, public, dimension(:) :: jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s
    real(wp), allocatable, public, dimension(:) :: jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas
    $:GPU_DECLARE(create='[jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s]')
    $:GPU_DECLARE(create='[jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas]')
#endif

    integer  :: jwl_idx                  !< JWL fluid index.
    real(wp) :: jwl_cv_prod, jwl_cv_air  !< Products/air specific heats.
    $:GPU_DECLARE(create='[jwl_idx]')
    $:GPU_DECLARE(create='[jwl_cv_prod, jwl_cv_air]')

contains

    !> Floor x to `floor`; NaNs pass through unchanged.
    subroutine s_jwl_floor(x, floor)

        $:GPU_ROUTINE(function_name='s_jwl_floor',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(inout) :: x
        real(wp), intent(in)    :: floor

        if (x == x) then
            if (x < floor) x = floor
        end if

    end subroutine s_jwl_floor

    ! Rocflu state-interpolated closure

    !> Rocflu effective coefficients between ambient air and JWL products. Reference: thierrydaoud/Rocflupicl,
    !! modflu/RFLU_ModJWL.F90 (RFLU_JWL_P_ER, RFLU_JWL_E_PR, RFLU_JWL_T_PR, and RFLU_JWL_C_ER).
    !!
    !! An and Bn blend linearly from 0 to A/B as e crosses [air_e0, e_j]; omega and cv blend
    !! linearly from air to products values as rho crosses [air_rho0, rho0]. Blended
    !! coefficients are used for Y <= 0.99, pure JWL for Y > 0.99.
    subroutine s_jwl_rocflu_coeffs(rho, e, Y, A, B, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, cv_a, An, Bn, omega, cv, &
                                   & mA, mB, momega)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_coeffs',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, cv_a
        real(wp), intent(out) :: An, Bn, omega, cv, mA, mB, momega
        real(wp)              :: e_j, xi
        real(wp), parameter   :: y_products_cutoff = 0.99_wp

        e_j = E0/rho0

        if (Y > y_products_cutoff) then
            An = A
            Bn = B
            omega = omega0
            cv = cv_j
            mA = 0._wp
            mB = 0._wp
            momega = 0._wp
            return
        end if

        ! Blended energy coefficients (Region I/II/III).
        mA = 0._wp
        mB = 0._wp
        if (e <= air_e0) then
            An = 0._wp
            Bn = 0._wp
        else if (e >= e_j) then
            An = A
            Bn = B
        else
            mA = A/(e_j - air_e0)
            mB = B/(e_j - air_e0)
            An = mA*(e - air_e0)
            Bn = mB*(e - air_e0)
        end if

        ! Density-blended omega, cv, and their rho-derivative.
        xi = min(max((rho - air_rho0)/(rho0 - air_rho0), 0._wp), 1._wp)
        omega = air_gamma + xi*(omega0 - air_gamma)
        cv = cv_a + xi*(cv_j - cv_a)
        if (rho > air_rho0 .and. rho < rho0) then
            momega = (omega0 - air_gamma)/(rho0 - air_rho0)
        else
            momega = 0._wp
        end if

    end subroutine s_jwl_rocflu_coeffs

    !> Rocflu single-fluid state-interpolated closure: (rho, e, Y) -> (p, T, c²).
    subroutine s_jwl_rocflu_state_er(rho, e, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, cv_a, pres, T, &
                                     & c2)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_state_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, cv_a
        real(wp), intent(out) :: pres, T, c2
        real(wp)              :: rho_s, Y_s, An, Bn, omega, cv, mA, mB, momega, V, exp1, exp2, coef1, coef2

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)
        ! Single Garno cutoff (YP = 99%): the blended branch covers all Y <= 0.99 and
        ! decays exactly to the ideal-gas air law (p = omega*rho*e, T = e/cv, c2 = (omega+1)p/rho)
        ! as An, Bn -> 0 (e -> air_e0) and the cold-curve exponentials -> 0 (rho -> air_rho0).
        call s_jwl_rocflu_coeffs(rho_s, e, Y_s, A, B, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, cv_a, An, Bn, omega, &
                                 & cv, mA, mB, momega)
        V = rho0/rho_s
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        coef1 = (1._wp - omega/(R1*V))*exp1
        coef2 = (1._wp - omega/(R2*V))*exp2
        pres = An*coef1 + Bn*coef2 + omega*rho_s*e
        T = (pres - An*exp1 - Bn*exp2)/(omega*cv*rho_s)
        c2 = exp1*An*(R1*rho0/rho_s**2 - omega/rho_s - omega/(R1*rho0) - rho_s*momega/(R1*rho0)) + mA*pres*coef1/rho_s**2 &
                      & + exp2*Bn*(R2*rho0/rho_s**2 - omega/rho_s - omega/(R2*rho0) - rho_s*momega/(R2*rho0)) &
                      & + mB*pres*coef2/rho_s**2 + omega*(e + pres/rho_s) + momega*rho_s*e
        call s_jwl_floor(pres, sgm_eps)
        call s_jwl_floor(T, 1._wp)
        call s_jwl_floor(c2, sgm_eps)

    end subroutine s_jwl_rocflu_state_er

    !> Analytic inverse of the Rocflu pressure closure: (rho, p, Y) -> e.
    !!
    !! The coefficient model makes An(e) and Bn(e) piecewise linear in e. Evaluating
    !! coefficients at the Region-II midpoint gives the exact linear slope for the
    !! active branch, then the low- and high-energy branches correct the saturated
    !! offsets used by the forward pressure law.
    subroutine s_jwl_rocflu_energy_pr(rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, cv_a, e)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, cv_a
        real(wp), intent(out) :: e
        real(wp)              :: rho_s, Y_s, e_j, e_eval, An, Bn, omega, cv, mA, mB, momega, V, C1, C2

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        e_j = E0/rho0
        e_eval = 0.5_wp*(e_j + air_e0)
        call s_jwl_rocflu_coeffs(rho_s, e_eval, Y_s, A, B, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, cv_a, An, Bn, &
                                 & omega, cv, mA, mB, momega)
        V = rho0/rho_s
        C1 = (1._wp - omega/(R1*V))*exp(-R1*V)
        C2 = (1._wp - omega/(R2*V))*exp(-R2*V)

        ! Region II: exact linear inverse for the active coefficient branch.
        e = (pres + (mA*C1 + mB*C2)*e_eval - An*C1 - Bn*C2)/max(mA*C1 + mB*C2 + omega*rho_s, sgm_eps)
        if (e < air_e0) then
            ! Region I: subtract the low-energy coefficient offsets, which are nonzero in the pure-JWL branch.
            e = (pres - (An - mA*(e_eval - air_e0))*C1 - (Bn - mB*(e_eval - air_e0))*C2)/max(omega*rho_s, sgm_eps)
        else if (e > e_j) then
            ! Region III: An=A, Bn=B -> p = A*C1 + B*C2 + omega*rho*e.
            e = (pres - A*C1 - B*C2)/max(omega*rho_s, sgm_eps)
        end if
        call s_jwl_floor(e, 0._wp)

    end subroutine s_jwl_rocflu_energy_pr

    ! Public entry points: look up fluid jidx's parameters, then evaluate the closure.

    !> Full state from energy for fluid jidx: (rho, e, Y) -> (p, T, c).
    subroutine s_jwl_mix_state_er(rho, e, Y, jidx, pres, T, c)

        $:GPU_ROUTINE(function_name='s_jwl_mix_state_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: pres, T, c
        real(wp)              :: c2

        call s_jwl_rocflu_state_er(rho, e, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                   & jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                   & jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, pres, T, c2)
        c = sqrt(max(c2, sgm_eps))

    end subroutine s_jwl_mix_state_er

    !> Energy from pressure for fluid jidx: (rho, p, Y) -> e.
    subroutine s_jwl_mix_energy_pr(rho, pres, Y, jidx, e)

        $:GPU_ROUTINE(function_name='s_jwl_mix_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: e

        call s_jwl_rocflu_energy_pr(max(rho, sgm_eps), pres, min(max(Y, 0._wp), 1._wp), jwl_As(jidx), jwl_Bs(jidx), &
                                    & jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), &
                                    & jwl_air_e0s(jidx), jwl_air_rho0s(jidx), jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, e)

    end subroutine s_jwl_mix_energy_pr

    !> Sound speed for fluid jidx: invert to energy, then evaluate c.
    subroutine s_jwl_mix_sound_speed(rho, pres, Y, jidx, c)

        $:GPU_ROUTINE(function_name='s_jwl_mix_sound_speed',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: c
        real(wp)              :: e, p_m, T_m, c2

        call s_jwl_mix_energy_pr(rho, pres, Y, jidx, e)
        call s_jwl_rocflu_state_er(rho, e, min(max(Y, 0._wp), 1._wp), jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                   & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                   & jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, p_m, T_m, c2)
        c = sqrt(max(c2, sgm_eps))

    end subroutine s_jwl_mix_sound_speed

    !> Initialize JWL parameter tables.
    impure subroutine s_initialize_jwl_module

        use m_mpi_common, only: s_mpi_abort
        use m_helper_basic, only: f_approx_equal, f_is_default

        integer  :: i, n_jwl, air_idx
        real(wp) :: jwl_E0_from_Q

        @:ALLOCATE(jwl_As(1:num_fluids), jwl_Bs(1:num_fluids), jwl_R1s(1:num_fluids), jwl_R2s(1:num_fluids), &
                   & jwl_omegas(1:num_fluids), jwl_rho0s(1:num_fluids), jwl_E0s(1:num_fluids))
        @:ALLOCATE(jwl_air_e0s(1:num_fluids), jwl_air_rho0s(1:num_fluids), jwl_air_gammas(1:num_fluids))

        jwl_idx = 0
        n_jwl = 0
        air_idx = 0
        do i = 1, num_fluids
            if (fluid_pp(i)%eos == eos_jwl .and. .not. f_is_default(fluid_pp(i)%jwl_rho0)) then
                if (f_is_default(fluid_pp(i)%jwl_E0) .and. .not. f_is_default(fluid_pp(i)%jwl_Q)) then
                    fluid_pp(i)%jwl_E0 = fluid_pp(i)%jwl_rho0*fluid_pp(i)%jwl_Q
                else if (.not. f_is_default(fluid_pp(i)%jwl_E0) .and. f_is_default(fluid_pp(i)%jwl_Q)) then
                    fluid_pp(i)%jwl_Q = fluid_pp(i)%jwl_E0/fluid_pp(i)%jwl_rho0
                end if
            end if
            jwl_As(i) = fluid_pp(i)%jwl_A
            jwl_Bs(i) = fluid_pp(i)%jwl_B
            jwl_R1s(i) = fluid_pp(i)%jwl_R1
            jwl_R2s(i) = fluid_pp(i)%jwl_R2
            jwl_omegas(i) = fluid_pp(i)%jwl_omega
            jwl_rho0s(i) = fluid_pp(i)%jwl_rho0
            jwl_E0s(i) = fluid_pp(i)%jwl_E0
            jwl_air_e0s(i) = fluid_pp(i)%jwl_air_e0
            jwl_air_rho0s(i) = fluid_pp(i)%jwl_air_rho0
            jwl_air_gammas(i) = fluid_pp(i)%jwl_air_gamma
            if (fluid_pp(i)%eos == eos_jwl) then
                jwl_idx = i
                n_jwl = n_jwl + 1
                if (f_is_default(fluid_pp(i)%jwl_A) .or. f_is_default(fluid_pp(i)%jwl_B) .or. f_is_default(fluid_pp(i)%jwl_R1) &
                    & .or. f_is_default(fluid_pp(i)%jwl_R2) .or. f_is_default(fluid_pp(i)%jwl_omega) &
                    & .or. f_is_default(fluid_pp(i)%jwl_rho0) .or. f_is_default(fluid_pp(i)%jwl_E0)) then
                    call s_mpi_abort('fluid_pp%eos = eos_jwl requires jwl_A, jwl_B, jwl_R1, jwl_R2, ' &
                                     & // 'jwl_omega, jwl_rho0, and either jwl_Q or jwl_E0 to be set.')
                end if
                if (.not. f_is_default(fluid_pp(i)%jwl_Q)) then
                    jwl_E0_from_Q = fluid_pp(i)%jwl_rho0*fluid_pp(i)%jwl_Q
                    if (.not. f_approx_equal(fluid_pp(i)%jwl_E0, jwl_E0_from_Q, 1.e-8_wp)) then
                        call s_mpi_abort('fluid_pp%eos = eos_jwl requires jwl_E0 = jwl_rho0*jwl_Q when both jwl_E0 and jwl_Q are set.')
                    end if
                end if
                if (f_is_default(fluid_pp(i)%jwl_air_e0) .or. f_is_default(fluid_pp(i)%jwl_air_rho0) &
                    & .or. f_is_default(fluid_pp(i)%jwl_air_gamma)) then
                    call s_mpi_abort('fluid_pp%eos = eos_jwl requires jwl_air_e0, jwl_air_rho0, and jwl_air_gamma to all be set.')
                end if
                if (fluid_pp(i)%jwl_R1 <= 0._wp .or. fluid_pp(i)%jwl_R2 <= 0._wp .or. fluid_pp(i)%jwl_omega <= 0._wp &
                    & .or. fluid_pp(i)%jwl_rho0 <= 0._wp .or. fluid_pp(i)%jwl_E0 <= 0._wp .or. fluid_pp(i)%jwl_air_rho0 <= 0._wp &
                    & .or. fluid_pp(i)%jwl_air_gamma <= 0._wp) then
                    call s_mpi_abort('JWL parameters jwl_R1, jwl_R2, jwl_omega, jwl_rho0, jwl_Q/jwl_E0, jwl_air_rho0, ' &
                                     & // 'and jwl_air_gamma must be positive.')
                end if
            else if (air_idx == 0) then
                air_idx = i
            end if
        end do

        if (n_jwl > 1) then
            call s_mpi_abort('At most one fluid may use eos_jwl; found more than one.')
        end if

        if (jwl_idx > 0 .and. model_eqns /= model_eqns_5eq) then
            call s_mpi_abort('eos_jwl is only supported with model_eqns_5eq.')
        end if

        jwl_cv_prod = 0._wp
        jwl_cv_air = 0._wp
        if (jwl_idx > 0) jwl_cv_prod = fluid_pp(jwl_idx)%cv
        if (air_idx > 0) jwl_cv_air = fluid_pp(air_idx)%cv

        if (jwl_idx > 0) then
            if (f_is_default(jwl_cv_prod) .or. jwl_cv_prod <= 0._wp) then
                call s_mpi_abort('The Rocflu closure requires positive fluid_pp%cv for the JWL fluid.')
            end if
            if (num_fluids > 1 .and. (air_idx == 0 .or. f_is_default(jwl_cv_air) .or. jwl_cv_air <= 0._wp)) then
                call s_mpi_abort('The Rocflu closure requires positive fluid_pp%cv for the non-JWL air fluid.')
            end if
            if (jwl_rho0s(jwl_idx) <= jwl_air_rho0s(jwl_idx) .or. jwl_E0s(jwl_idx)/jwl_rho0s(jwl_idx) <= jwl_air_e0s(jwl_idx)) then
                call s_mpi_abort('The Rocflu closure requires increasing air-to-products reference density and energy.')
            end if
        end if

        $:GPU_UPDATE(device='[jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s]')
        $:GPU_UPDATE(device='[jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas, jwl_idx]')
        $:GPU_UPDATE(device='[jwl_cv_prod, jwl_cv_air]')

    end subroutine s_initialize_jwl_module

    !> Deallocate the per-fluid JWL parameter tables.
    impure subroutine s_finalize_jwl_module

        @:DEALLOCATE(jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s)
        @:DEALLOCATE(jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas)

    end subroutine s_finalize_jwl_module

end module m_jwl
