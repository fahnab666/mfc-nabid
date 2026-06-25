!>
!! @file
!! @brief JWL EOS and Rocflu state-interpolated mixture closure.

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief JWL equation of state and Rocflu state-interpolated mixture closure.
!!
!! Closure path: (rho_m, e_m, Y(:), component_ids(:)) -> (p_m, T_m, c_m).
!! Mixture rule: Rocflu state-interpolated closure (Garno/Stanley, after RFLU_ModJWL.F90).
!! The blended An/Bn/omega coefficients ramp linearly from ambient-air values at e = air_e0
!! to full JWL values at e = e_j = E0/rho0, with a C1-smooth cubic Hermite transition in
!! mass fraction Y over [0.95, 1.0] replacing the original hard cut at Y = 0.99.
module m_jwl

    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use m_global_parameters
    use m_constants, only: jwl_mix_type_rocflu

    implicit none

    private
    public :: s_initialize_jwl_module, s_finalize_jwl_module, s_jwl_mix_pressure_er, s_jwl_mix_energy_pr, s_jwl_mix_sound_speed, &
        & jwl_idx, s_jwl_mixture_closure

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

    !> Component IDs used to identify explosive products in the multi-fluid Y array. These integer tags match the species-table
    !! convention from the upstream case setup.
    integer, parameter :: kpw_products_id = 2
    integer, parameter :: kpw_products_alt_id = 4

contains

    ! Shared utilities

    !> Floor x to a positive value; NaNs pass through unchanged.
    subroutine s_jwl_floor_positive(x, floor)

        $:GPU_ROUTINE(function_name='s_jwl_floor_positive',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(inout) :: x
        real(wp), intent(in)    :: floor

        if (x == x) then
            if (x < floor) x = floor
        end if

    end subroutine s_jwl_floor_positive

    !> Floor x to zero; NaNs pass through unchanged.
    subroutine s_jwl_floor_nonnegative(x)

        $:GPU_ROUTINE(function_name='s_jwl_floor_nonnegative',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(inout) :: x

        if (x == x) then
            if (x < 0._wp) x = 0._wp
        end if

    end subroutine s_jwl_floor_nonnegative

    !> Total clamped mass fraction of the explosive products components.
    subroutine s_jwl_products_mass_fraction(Y, component_ids, ncomp, Y_exp)

        $:GPU_ROUTINE(function_name='s_jwl_products_mass_fraction',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: Y(ncomp)
        integer, intent(in)   :: component_ids(ncomp), ncomp
        real(wp), intent(out) :: Y_exp
        real(wp)              :: Y_k, Y_sum
        integer               :: k

        Y_sum = 0._wp
        Y_exp = 0._wp
        do k = 1, ncomp
            Y_k = min(max(Y(k), 0._wp), 1._wp)
            Y_sum = Y_sum + Y_k
            if (component_ids(k) == kpw_products_id .or. component_ids(k) == kpw_products_alt_id) Y_exp = Y_exp + Y_k
        end do
        if (Y_sum > sgm_eps) then
            Y_exp = min(max(Y_exp/Y_sum, 0._wp), 1._wp)
        else
            Y_exp = 0._wp
        end if

    end subroutine s_jwl_products_mass_fraction

    ! Rocflu state-interpolated closure

    !> Rocflu effective coefficients between ambient air and JWL products. Reference: thierrydaoud/Rocflupicl,
    !! modflu/RFLU_ModJWL.F90 (RFLU_JWL_P_ER, RFLU_JWL_E_PR, RFLU_JWL_T_PR, and RFLU_JWL_C_ER).
    !!
    !! An and Bn blend linearly from 0 to A/B as e crosses [air_e0, e_j]; omega and cv blend
    !! linearly from air values to JWL reference values as rho crosses [air_rho0, rho0].
    !! A C1-smooth cubic Hermite S-curve in Y over [0.95, 1.0] replaces the original hard
    !! cut at Y = 0.99, eliminating the pressure jump at the pure-products interface.
    subroutine s_jwl_rocflu_coeffs(rho, e, Y, A, B, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, cv_a, An, Bn, omega, cv, &
                                   & mA, mB, momega)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_coeffs',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, cv_a
        real(wp), intent(out) :: An, Bn, omega, cv, mA, mB, momega
        real(wp)              :: e_j, xi, phi
        real(wp), parameter   :: y_blend_lo = 0.95_wp

        e_j = E0/rho0

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

        ! C1 cubic Hermite Y-blend: ramps smoothly from mixture coefficients at Y = y_blend_lo
        ! to pure-JWL coefficients at Y = 1. Replaces the original hard cut at Y = 0.99.
        if (Y <= y_blend_lo) return
        phi = (Y - y_blend_lo)/(1._wp - y_blend_lo)
        phi = phi*phi*(3._wp - 2._wp*phi)
        An = (1._wp - phi)*An + phi*A
        Bn = (1._wp - phi)*Bn + phi*B
        omega = (1._wp - phi)*omega + phi*omega0
        cv = (1._wp - phi)*cv + phi*cv_j
        mA = (1._wp - phi)*mA
        mB = (1._wp - phi)*mB
        momega = (1._wp - phi)*momega

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
        if (Y_s <= 0.01_wp) then
            pres = air_gamma*rho_s*e
            T = e/cv_a
            c2 = (air_gamma + 1._wp)*pres/rho_s
        else
            call s_jwl_rocflu_coeffs(rho_s, e, Y_s, A, B, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, cv_a, An, Bn, &
                                     & omega, cv, mA, mB, momega)
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
        end if
        call s_jwl_floor_positive(pres, sgm_eps)
        call s_jwl_floor_positive(T, 1._wp)
        call s_jwl_floor_positive(c2, sgm_eps)

    end subroutine s_jwl_rocflu_state_er

    !> Analytic inverse of the Rocflu pressure closure: (rho, p, Y) -> e.
    !!
    !! Derived by substituting An = mA*e + (A - mA*e_j) into the pressure formula and
    !! solving for e. The denominator (mA*C1 + mB*C2 + omega*rho) is strictly positive.
    subroutine s_jwl_rocflu_energy_pr(rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, cv_a, e)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, cv_a
        real(wp), intent(out) :: e
        real(wp)              :: rho_s, Y_s, e_j, An, Bn, omega, cv, mA, mB, momega, V, C1, C2

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)
        if (Y_s <= 0.01_wp) then
            e = pres/(air_gamma*rho_s)
            call s_jwl_floor_nonnegative(e)
            return
        end if

        e_j = E0/rho0
        ! Evaluate coefficients at the midpoint energy; omega, mA, mB depend only on rho and Y,
        ! so the evaluation point for e does not affect the inversion.
        call s_jwl_rocflu_coeffs(rho_s, 0.5_wp*(e_j + air_e0), Y_s, A, B, omega0, rho0, E0, air_e0, air_rho0, air_gamma, cv_j, &
                                 & cv_a, An, Bn, omega, cv, mA, mB, momega)
        V = rho0/rho_s
        C1 = (1._wp - omega/(R1*V))*exp(-R1*V)
        C2 = (1._wp - omega/(R2*V))*exp(-R2*V)
        ! General Region-II linear inverse; correct for all Y including the blend zone.
        e = (pres + (mA*e_j - A)*C1 + (mB*e_j - B)*C2)/(mA*C1 + mB*C2 + omega*rho_s)
        if (e < air_e0) then
            ! Region I: An=0, Bn=0 -> p = omega*rho*e.
            e = pres/(omega*rho_s)
        else if (e > e_j) then
            ! Region III: An=A, Bn=B -> p = A*C1 + B*C2 + omega*rho*e.
            e = (pres - A*C1 - B*C2)/max(omega*rho_s, sgm_eps)
        end if
        call s_jwl_floor_nonnegative(e)

    end subroutine s_jwl_rocflu_energy_pr

    ! Black-box closure and public dispatch wrappers

    !> Black-box Rocflu closure: (rho, e, Y) -> (p, T, c).
    subroutine s_jwl_rocflu_closure(rho_m, e_m, Y, component_ids, ncomp, jidx, p_m, T_m, c_m)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_closure',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho_m, e_m, Y(ncomp)
        integer, intent(in)   :: component_ids(ncomp), ncomp, jidx
        real(wp), intent(out) :: p_m, T_m, c_m
        real(wp)              :: Y_s, c2

        call s_jwl_products_mass_fraction(Y, component_ids, ncomp, Y_s)
        call s_jwl_rocflu_state_er(rho_m, e_m, Y_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                   & jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                   & jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, p_m, T_m, c2)
        c_m = sqrt(max(c2, sgm_eps))

    end subroutine s_jwl_rocflu_closure

    !> Black-box mixture closure dispatcher: (rho, e, Y, alpha_j, component_ids) -> (p, T, c).
    subroutine s_jwl_mixture_closure(rho_m, e_m, Y, alpha_j, component_ids, ncomp, jidx, p_m, T_m, c_m)

        $:GPU_ROUTINE(function_name='s_jwl_mixture_closure',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho_m, e_m, Y(ncomp), alpha_j
        integer, intent(in)   :: component_ids(ncomp), ncomp, jidx
        real(wp), intent(out) :: p_m, T_m, c_m

        select case (jwl_mix_type)
        case (jwl_mix_type_rocflu)
            call s_jwl_rocflu_closure(rho_m, e_m, Y, component_ids, ncomp, jidx, p_m, T_m, c_m)
        case default
            ! Unreachable: s_check_jwl_inputs prohibits all modes except jwl_mix_type_rocflu.
            p_m = ieee_value(0._wp, ieee_quiet_nan)
            T_m = p_m
            c_m = p_m
        end select

    end subroutine s_jwl_mixture_closure

    !> Dispatch pressure-from-energy to the active JWL closure.
    subroutine s_jwl_mix_pressure_er(rho, e, Y, alpha_j, jidx, pres)

        $:GPU_ROUTINE(function_name='s_jwl_mix_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, alpha_j
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: pres
        real(wp)              :: T_dummy, c2_dummy

        select case (jwl_mix_type)
        case (jwl_mix_type_rocflu)
            call s_jwl_rocflu_state_er(rho, e, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                       & jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                       & jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, pres, T_dummy, c2_dummy)
        case default
            pres = ieee_value(0._wp, ieee_quiet_nan)
        end select

    end subroutine s_jwl_mix_pressure_er

    !> Dispatch energy-from-pressure to the active JWL closure.
    subroutine s_jwl_mix_energy_pr(rho, pres, Y, alpha_j, jidx, e)

        $:GPU_ROUTINE(function_name='s_jwl_mix_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: e

        select case (jwl_mix_type)
        case (jwl_mix_type_rocflu)
            call s_jwl_rocflu_energy_pr(max(rho, sgm_eps), pres, min(max(Y, 0._wp), 1._wp), jwl_As(jidx), jwl_Bs(jidx), &
                                        & jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), &
                                        & jwl_air_e0s(jidx), jwl_air_rho0s(jidx), jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, e)
        case default
            e = ieee_value(0._wp, ieee_quiet_nan)
        end select

    end subroutine s_jwl_mix_energy_pr

    !> Dispatch sound speed to the active closure.
    subroutine s_jwl_mix_sound_speed(rho, pres, Y, alpha_j, jidx, c)

        $:GPU_ROUTINE(function_name='s_jwl_mix_sound_speed',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: c
        real(wp)              :: e, p_m, T_m, c2

        select case (jwl_mix_type)
        case (jwl_mix_type_rocflu)
            call s_jwl_mix_energy_pr(rho, pres, Y, alpha_j, jidx, e)
            call s_jwl_rocflu_state_er(rho, e, min(max(Y, 0._wp), 1._wp), jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), &
                                       & jwl_R2s(jidx), jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), &
                                       & jwl_air_rho0s(jidx), jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, p_m, T_m, c2)
            c = sqrt(max(c2, sgm_eps))
        case default
            c = ieee_value(0._wp, ieee_quiet_nan)
        end select

    end subroutine s_jwl_mix_sound_speed

    !> Initialize JWL parameter tables.
    impure subroutine s_initialize_jwl_module

        use m_mpi_common, only: s_mpi_abort
        use m_helper_basic, only: f_is_default

        integer :: i, n_jwl, air_idx

        @:ALLOCATE(jwl_As    (1:num_fluids))
        @:ALLOCATE(jwl_Bs    (1:num_fluids))
        @:ALLOCATE(jwl_R1s   (1:num_fluids))
        @:ALLOCATE(jwl_R2s   (1:num_fluids))
        @:ALLOCATE(jwl_omegas(1:num_fluids))
        @:ALLOCATE(jwl_rho0s (1:num_fluids))
        @:ALLOCATE(jwl_E0s   (1:num_fluids))
        @:ALLOCATE(jwl_air_e0s    (1:num_fluids))
        @:ALLOCATE(jwl_air_rho0s  (1:num_fluids))
        @:ALLOCATE(jwl_air_gammas (1:num_fluids))

        jwl_idx = 0
        n_jwl = 0
        air_idx = 0
        do i = 1, num_fluids
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
                                     & // 'jwl_omega, jwl_rho0, and jwl_E0 to all be set.')
                end if
                if (f_is_default(fluid_pp(i)%jwl_air_e0) .or. f_is_default(fluid_pp(i)%jwl_air_rho0) &
                    & .or. f_is_default(fluid_pp(i)%jwl_air_gamma)) then
                    call s_mpi_abort('fluid_pp%eos = eos_jwl requires jwl_air_e0, jwl_air_rho0, and jwl_air_gamma to all be set.')
                end if
                if (fluid_pp(i)%jwl_R1 <= 0._wp .or. fluid_pp(i)%jwl_R2 <= 0._wp .or. fluid_pp(i)%jwl_omega <= 0._wp &
                    & .or. fluid_pp(i)%jwl_rho0 <= 0._wp .or. fluid_pp(i)%jwl_air_rho0 <= 0._wp &
                    & .or. fluid_pp(i)%jwl_air_gamma <= 0._wp) then
                    call s_mpi_abort('JWL parameters jwl_R1, jwl_R2, jwl_omega, jwl_rho0, jwl_air_rho0, ' &
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

        if (jwl_idx > 0 .and. jwl_mix_type == jwl_mix_type_rocflu) then
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
        $:GPU_UPDATE(device='[jwl_mix_type]')

    end subroutine s_initialize_jwl_module

    !> Deallocate the per-fluid JWL parameter tables.
    impure subroutine s_finalize_jwl_module

        @:DEALLOCATE(jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s)
        @:DEALLOCATE(jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas)

    end subroutine s_finalize_jwl_module

end module m_jwl
