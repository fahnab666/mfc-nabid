!>
!! @file
!! @brief Contains module m_jwl

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Jones-Wilkins-Lee (JWL) equation of state and JWL/ideal-gas mixture closures. Holds the per-fluid JWL parameter tables and
!! the pressure<->energy / sound-speed routines for every supported mixture closure (jwl_mix_type). Split out of
!! m_variables_conversion so JWL physics lives in one place.
module m_jwl

    use m_global_parameters
    use m_constants, only: jwl_mix_type_isobaric, jwl_mix_type_kuhl, jwl_mix_type_ptequil, jwl_mix_type_rocflu

    implicit none

    private
    public :: s_initialize_jwl_module, s_finalize_jwl_module, s_jwl_mixture_sound_speed_squared, &
        & s_jwl_rocflu_sound_speed_squared, s_jwl_mix_pressure_er, s_jwl_mix_energy_pr, jwl_idx

    ! Per-fluid JWL parameter tables. In simulation these live in m_global_parameters
    ! (shared with the solver hot path); for pre/post_process they are declared here.
#ifndef MFC_SIMULATION
    real(wp), allocatable, public, dimension(:) :: jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s
    real(wp), allocatable, public, dimension(:) :: jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas
    $:GPU_DECLARE(create='[jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s]')
    $:GPU_DECLARE(create='[jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas]')
#endif

    !> Mass-fraction cutoff below/above which a cell is treated as pure air/products: the Rocflu blend switches to the
    !! single-material EOS, and the p-T-equilibrium solve (whose root sits a distance ~(1-Y) from the alpha bracket end and is
    !! ill-conditioned as Y -> 0 or 1) is skipped.
    real(wp), parameter :: jwl_pure_cutoff = 1.e-4_wp
    !> Relative stopping tolerance for scalar closure solves. Scales with the working precision while retaining the historical
    !! double-precision threshold used by existing goldens.
    real(wp), parameter :: jwl_root_rel_tol = max(1.e-12_wp, 100._wp*epsilon(1._wp))
    !> Endpoint tolerance for treating states as pure material in dispatcher shortcuts.
    real(wp), parameter :: jwl_endpoint_tol = max(1.e-12_wp, 100._wp*epsilon(1._wp))
    integer             :: jwl_idx                  !< Index of the JWL fluid (0 if none)
    real(wp)            :: jwl_cv_prod, jwl_cv_air  !< Products/air specific heats for the p-T-equilibrium closure.
    $:GPU_DECLARE(create='[jwl_idx]')
    $:GPU_DECLARE(create='[jwl_cv_prod, jwl_cv_air]')

contains

    !> Floor a positive scalar only when it is finite and below the supplied floor value.
    !! @param[inout] x Scalar to guard; NaNs are intentionally left unchanged.
    !! @param[in] floor Small positive lower bound for finite values.
    subroutine s_jwl_floor_positive(x, floor)

        $:GPU_ROUTINE(function_name='s_jwl_floor_positive',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(inout) :: x
        real(wp), intent(in)    :: floor

        if (x == x) then
            if (x < floor) x = floor
        end if

    end subroutine s_jwl_floor_positive

    !> Floor a scalar to zero only when it is finite and negative.
    !! @param[inout] x Scalar to guard; NaNs are intentionally left unchanged.
    subroutine s_jwl_floor_nonnegative(x)

        $:GPU_ROUTINE(function_name='s_jwl_floor_nonnegative',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(inout) :: x

        if (x == x) then
            if (x < 0._wp) x = 0._wp
        end if

    end subroutine s_jwl_floor_nonnegative

    !> Compute the JWL cold pressure term with relative volume V = rho0/rho.
    !! @param[in] rho Material density.
    !! @param[in] A First JWL pressure coefficient.
    !! @param[in] B Second JWL pressure coefficient.
    !! @param[in] R1 First JWL exponential coefficient.
    !! @param[in] R2 Second JWL exponential coefficient.
    !! @param[in] omega JWL Gruneisen coefficient.
    !! @param[in] rho0 Reference products density.
    !! @param[out] pcold Cold-curve pressure contribution.
    subroutine s_jwl_pcold(rho, A, B, R1, R2, omega, rho0, pcold)

        $:GPU_ROUTINE(function_name='s_jwl_pcold',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, A, B, R1, R2, omega, rho0
        real(wp), intent(out) :: pcold
        real(wp)              :: V

        V = rho0/max(rho, sgm_eps)
        pcold = A*(1._wp - omega/(R1*V))*exp(-R1*V) + B*(1._wp - omega/(R2*V))*exp(-R2*V)

    end subroutine s_jwl_pcold

    !> Compute d(pcold)/d(rho) for the JWL EOS.
    !! @param[in] rho Material density.
    !! @param[in] A First JWL pressure coefficient.
    !! @param[in] B Second JWL pressure coefficient.
    !! @param[in] R1 First JWL exponential coefficient.
    !! @param[in] R2 Second JWL exponential coefficient.
    !! @param[in] omega JWL Gruneisen coefficient.
    !! @param[in] rho0 Reference products density.
    !! @param[out] dpcold_drho Density derivative of the cold-curve pressure.
    subroutine s_jwl_dpcold_drho(rho, A, B, R1, R2, omega, rho0, dpcold_drho)

        $:GPU_ROUTINE(function_name='s_jwl_dpcold_drho',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, A, B, R1, R2, omega, rho0
        real(wp), intent(out) :: dpcold_drho
        real(wp)              :: V, rho_safe

        rho_safe = max(rho, sgm_eps)
        V = rho0/rho_safe
        dpcold_drho = A*exp(-R1*V)*(V/rho_safe)*(R1 - omega/V - omega/(R1*V**2)) + B*exp(-R2*V)*(V/rho_safe)*(R2 - omega/V &
                            & - omega/(R2*V**2))

    end subroutine s_jwl_dpcold_drho

    !> Compute the JWL isentropic sound-speed squared using the Rocflu/Stanley form.
    !! @param[in] rho Material density.
    !! @param[in] pres Pressure.
    !! @param[in] A First JWL pressure coefficient.
    !! @param[in] B Second JWL pressure coefficient.
    !! @param[in] R1 First JWL exponential coefficient.
    !! @param[in] R2 Second JWL exponential coefficient.
    !! @param[in] omega JWL Gruneisen coefficient.
    !! @param[in] rho0 Reference products density.
    !! @param[out] c2 Isentropic sound speed squared.
    subroutine s_jwl_sound_speed_squared(rho, pres, A, B, R1, R2, omega, rho0, c2)

        $:GPU_ROUTINE(function_name='s_jwl_sound_speed_squared',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, A, B, R1, R2, omega, rho0
        real(wp), intent(out) :: c2
        real(wp)              :: rho_safe, pcold, dpcold_drho, e

        ! c^2 = dp/drho|_e + (p/rho^2)*dp/de|_rho with p = pcold(rho) + omega*rho*e, so
        ! c^2 = dpcold/drho + omega*e + omega*p/rho, with e recovered from p and pcold.
        rho_safe = max(rho, sgm_eps)
        call s_jwl_pcold(rho_safe, A, B, R1, R2, omega, rho0, pcold)
        call s_jwl_dpcold_drho(rho_safe, A, B, R1, R2, omega, rho0, dpcold_drho)

        e = (pres - pcold)/(omega*rho_safe)
        c2 = dpcold_drho + omega*e + omega*pres/rho_safe
        call s_jwl_floor_positive(c2, sgm_eps)

    end subroutine s_jwl_sound_speed_squared

    !> Sound-speed squared for mixed JWL/ideal-gas states using the frozen (mass-weighted) mixture rule c^2 = sum_k Y_k*c_k^2, each
    !! phase evaluated at its own density rho_k = (Y_k*rho)/alpha_k. The frozen estimate is smooth and monotone between the phase
    !! values (avoiding Wood's interface dip) and only feeds the Riemann wave-speed estimate, not the pressure closure.
    subroutine s_jwl_mixture_sound_speed_squared(rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, &
        & air_gamma, c2)

        $:GPU_ROUTINE(function_name='s_jwl_mixture_sound_speed_squared',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: c2
        real(wp)              :: rho_safe, Y_safe, a_j, a_a, rho1, rho2, c2_air, c2_jwl

        rho_safe = max(rho, sgm_eps)
        Y_safe = min(max(Y, 0._wp), 1._wp)
        a_j = min(max(alpha_j, 0._wp), 1._wp)
        a_a = 1._wp - a_j

        if (a_j <= sgm_eps) then
            c2 = (air_gamma + 1._wp)*pres/rho_safe
            call s_jwl_floor_positive(c2, sgm_eps)
            return
        end if

        if (a_a <= sgm_eps) then
            call s_jwl_sound_speed_squared(rho_safe, pres, A, B, R1, R2, omega0, rho0, c2)
            call s_jwl_floor_positive(c2, sgm_eps)
            return
        end if

        rho1 = max(Y_safe*rho_safe/a_j, sgm_eps)
        rho2 = max((1._wp - Y_safe)*rho_safe/a_a, sgm_eps)

        call s_jwl_sound_speed_squared(rho1, pres, A, B, R1, R2, omega0, rho0, c2_jwl)
        c2_air = (air_gamma + 1._wp)*pres/rho2
        call s_jwl_floor_positive(c2_air, sgm_eps)
        call s_jwl_floor_positive(c2_jwl, sgm_eps)

        ! Frozen (mass-weighted) mixture sound speed: smooth and monotone between the phase values, avoiding Wood's interface dip
        c2 = Y_safe*c2_jwl + (1._wp - Y_safe)*c2_air
        call s_jwl_floor_positive(c2, sgm_eps)

    end subroutine s_jwl_mixture_sound_speed_squared

    !> JWL/ideal-gas mixture pressure from energy and density via the closed-form isobaric (mechanical-equilibrium) closure. Each
    !! phase is evaluated at its own density rho_k = (Y_k*rho)/alpha_k; both EOS are linear in pressure, so the common pressure that
    !! partitions the cell energy is exact (no iteration) and recovers the pure-air and pure-JWL limits.
    subroutine s_jwl_pressure_er(rho, e, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, pres)

        $:GPU_ROUTINE(function_name='s_jwl_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: pres
        real(wp)              :: rho_safe, Y_safe, a_j, a_a, rho1, rhoe, pref1, Kj

        rho_safe = max(rho, sgm_eps)
        Y_safe = min(max(Y, 0._wp), 1._wp)
        a_j = min(max(alpha_j, 0._wp), 1._wp)
        a_a = 1._wp - a_j
        rhoe = rho_safe*e

        if (a_j <= sgm_eps) then
            pres = air_gamma*rhoe
            call s_jwl_floor_positive(pres, sgm_eps)
            return
        end if

        if (a_a <= sgm_eps) then
            rho1 = rho_safe
        else
            rho1 = max(Y_safe*rho_safe/a_j, sgm_eps)
        end if

        call s_jwl_pcold(rho1, A, B, R1, R2, omega0, rho0, pref1)
        Kj = 1._wp/max(omega0, sgm_eps)

        pres = (rhoe + a_j*Kj*pref1)/max(a_j*Kj + a_a/max(air_gamma, sgm_eps), sgm_eps)
        call s_jwl_floor_positive(pres, sgm_eps)

    end subroutine s_jwl_pressure_er

    !> Specific internal energy from pressure and density: exact inverse of the isobaric mixture pressure in s_jwl_pressure_er,
    !! keeping primitive<->conservative consistent.
    subroutine s_jwl_energy_pr(rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, e)

        $:GPU_ROUTINE(function_name='s_jwl_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: e
        real(wp)              :: rho_safe, Y_work, a_j, a_a, rho1, pref1, Kj

        rho_safe = max(rho, sgm_eps)
        Y_work = min(max(Y, 0._wp), 1._wp)
        a_j = min(max(alpha_j, 0._wp), 1._wp)
        a_a = 1._wp - a_j

        if (a_j <= sgm_eps) then
            e = pres/max(air_gamma*rho_safe, sgm_eps)
            call s_jwl_floor_nonnegative(e)
            return
        end if

        if (a_a <= sgm_eps) then
            rho1 = rho_safe
        else
            rho1 = max(Y_work*rho_safe/a_j, sgm_eps)
        end if

        call s_jwl_pcold(rho1, A, B, R1, R2, omega0, rho0, pref1)
        Kj = 1._wp/max(omega0, sgm_eps)

        e = (pres*(a_j*Kj + a_a/max(air_gamma, sgm_eps)) - a_j*Kj*pref1)/rho_safe
        ! Sub-cold-curve finite states give e < 0; floor those while preserving NaNs.
        call s_jwl_floor_nonnegative(e)

    end subroutine s_jwl_energy_pr

    !> Kuhl/Khasainov temperature-form mixture pressure: p = Y*(A*exp(-R1*V) + B*exp(-R2*V)) + rho*R_mix*T, with the temperature T
    !! recovered from e. Both phases share the cell density rho; alpha_j is unused.
    subroutine s_jwl_kuhl_pressure_er(rho, e, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, pres)

        $:GPU_ROUTINE(function_name='s_jwl_kuhl_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: pres
        real(wp)              :: rho_safe, Y_safe, V, exp1, exp2, pbase, ecold, cv_mix, R_mix, T

        rho_safe = max(rho, sgm_eps)
        Y_safe = min(max(Y, 0._wp), 1._wp)
        V = rho0/rho_safe
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        pbase = A*exp1 + B*exp2
        ecold = A/(R1*rho0)*exp1 + B/(R2*rho0)*exp2
        cv_mix = max(Y_safe*jwl_cv_prod + (1._wp - Y_safe)*jwl_cv_air, sgm_eps)
        ! Positive because omega0, air_gamma, and both cv values are validated positive, and Y_safe is bounded.
        R_mix = Y_safe*omega0*jwl_cv_prod + (1._wp - Y_safe)*air_gamma*jwl_cv_air
        T = (e - Y_safe*ecold)/cv_mix
        call s_jwl_floor_positive(T, sgm_eps)
        pres = Y_safe*pbase + rho_safe*R_mix*T
        call s_jwl_floor_positive(pres, sgm_eps)

    end subroutine s_jwl_kuhl_pressure_er

    !> Specific internal energy from pressure and density: exact inverse of the Kuhl/Khasainov closure. alpha_j is unused.
    subroutine s_jwl_kuhl_energy_pr(rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, e)

        $:GPU_ROUTINE(function_name='s_jwl_kuhl_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: e
        real(wp)              :: rho_safe, Y_safe, p_safe, V, exp1, exp2, pbase, ecold, cv_mix, R_mix, T

        rho_safe = max(rho, sgm_eps)
        Y_safe = min(max(Y, 0._wp), 1._wp)
        p_safe = pres
        call s_jwl_floor_positive(p_safe, sgm_eps)
        V = rho0/rho_safe
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        pbase = A*exp1 + B*exp2
        ecold = A/(R1*rho0)*exp1 + B/(R2*rho0)*exp2
        cv_mix = max(Y_safe*jwl_cv_prod + (1._wp - Y_safe)*jwl_cv_air, sgm_eps)
        ! Positive because omega0, air_gamma, and both cv values are validated positive, and Y_safe is bounded.
        R_mix = Y_safe*omega0*jwl_cv_prod + (1._wp - Y_safe)*air_gamma*jwl_cv_air
        T = (p_safe - Y_safe*pbase)/max(rho_safe*R_mix, sgm_eps)
        call s_jwl_floor_positive(T, sgm_eps)
        e = Y_safe*ecold + cv_mix*T
        call s_jwl_floor_nonnegative(e)

    end subroutine s_jwl_kuhl_energy_pr

    !> Sound speed for the Kuhl/Khasainov closure, treated as a Mie-Gruneisen EOS with constant-Y Gamma_eff = R_mix/Cv_mix.
    subroutine s_jwl_kuhl_sound_speed_squared(rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, &
        & c2)

        $:GPU_ROUTINE(function_name='s_jwl_kuhl_sound_speed_squared',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: c2
        real(wp)              :: rho_safe, Y_safe, V, exp1, exp2, pbase, ecold, dpbase, decold, cv_mix, R_mix, Gamma_eff, e, dPeff

        rho_safe = max(rho, sgm_eps)
        Y_safe = min(max(Y, 0._wp), 1._wp)
        V = rho0/rho_safe
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        pbase = A*exp1 + B*exp2
        ecold = A/(R1*rho0)*exp1 + B/(R2*rho0)*exp2
        dpbase = (A*R1*rho0*exp1 + B*R2*rho0*exp2)/(rho_safe*rho_safe)
        decold = pbase/(rho_safe*rho_safe)
        cv_mix = max(Y_safe*jwl_cv_prod + (1._wp - Y_safe)*jwl_cv_air, sgm_eps)
        ! Positive because omega0, air_gamma, and both cv values are validated positive, and Y_safe is bounded.
        R_mix = Y_safe*omega0*jwl_cv_prod + (1._wp - Y_safe)*air_gamma*jwl_cv_air
        Gamma_eff = R_mix/cv_mix
        e = (pres - Y_safe*(pbase - Gamma_eff*rho_safe*ecold))/max(Gamma_eff*rho_safe, sgm_eps)
        call s_jwl_floor_nonnegative(e)
        dPeff = Y_safe*dpbase - Gamma_eff*Y_safe*(ecold + rho_safe*decold)
        c2 = dPeff + Gamma_eff*(e + pres/rho_safe)
        call s_jwl_floor_positive(c2, sgm_eps)

    end subroutine s_jwl_kuhl_sound_speed_squared

    !> p-T-equilibrium mixture pressure (Garno 2020, ref [22]): find the products volume fraction alpha_j that balances the two
    !! phase pressures. The shared temperature is closed-form from the energy split, so the residual p_j - p_a is monotone in
    !! alpha_j and solved with the bracket-preserving Illinois method.
    subroutine s_jwl_ptequil_pressure_er(rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, cv_j, cv_a, pres)

        $:GPU_ROUTINE(function_name='s_jwl_ptequil_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, cv_j, cv_a
        real(wp), intent(out) :: pres
        real(wp)              :: rho_s, Y_s, cv_mix, a_lo, a_hi, a_m, rj, ra, V, ecold, T, pj, pa, f_lo, f_m, pcg
        integer               :: it

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        ! Near-pure cells: the equilibrium root sits within ~(1-Y) of the alpha bracket end, where the
        ! solve is ill-conditioned, so below the cutoff evaluate the single-phase EOS directly.
        if (Y_s <= jwl_pure_cutoff) then
            pres = air_gamma*rho_s*e
            call s_jwl_floor_positive(pres, sgm_eps)
            return
        else if (1._wp - Y_s <= jwl_pure_cutoff) then
            call s_jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0, pcg)
            pres = pcg + omega0*rho_s*e
            call s_jwl_floor_positive(pres, sgm_eps)
            return
        end if

        cv_mix = max(Y_s*cv_j + (1._wp - Y_s)*cv_a, sgm_eps)

        a_lo = sgm_eps
        a_hi = 1._wp - sgm_eps
        rj = max(Y_s*rho_s/a_lo, sgm_eps)
        ra = max((1._wp - Y_s)*rho_s/(1._wp - a_lo), sgm_eps)
        V = rho0/rj
        ecold = A/(R1*rho0)*exp(-R1*V) + B/(R2*rho0)*exp(-R2*V)
        T = (e - Y_s*ecold)/cv_mix
        call s_jwl_floor_positive(T, sgm_eps)
        pj = A*exp(-R1*V) + B*exp(-R2*V) + omega0*rj*cv_j*T
        pa = air_gamma*ra*cv_a*T
        f_lo = pj - pa
        pres = 0.5_wp*(pj + pa)

        ! Bisection: a_m stays strictly inside the bracket every step, so the 1/a_m, 1/(1-a_m)
        ! phasic densities never hit their endpoint singularities.
        do it = 1, 60
            a_m = 0.5_wp*(a_lo + a_hi)
            rj = max(Y_s*rho_s/a_m, sgm_eps)
            ra = max((1._wp - Y_s)*rho_s/(1._wp - a_m), sgm_eps)
            V = rho0/rj
            ecold = A/(R1*rho0)*exp(-R1*V) + B/(R2*rho0)*exp(-R2*V)
            T = (e - Y_s*ecold)/cv_mix
            call s_jwl_floor_positive(T, sgm_eps)
            pj = A*exp(-R1*V) + B*exp(-R2*V) + omega0*rj*cv_j*T
            pa = air_gamma*ra*cv_a*T
            f_m = pj - pa
            pres = 0.5_wp*(pj + pa)
            if (abs(f_m) <= jwl_root_rel_tol*max(abs(pres), sgm_eps)) exit
            if (f_lo*f_m > 0._wp) then
                a_lo = a_m
                f_lo = f_m
            else
                a_hi = a_m
            end if
        end do
        call s_jwl_floor_positive(pres, sgm_eps)

    end subroutine s_jwl_ptequil_pressure_er

    !> Specific internal energy for the p-T-equilibrium mixture: inverse of s_jwl_ptequil_pressure_er. With p known the air phase
    !! fixes T in closed form, so the products-pressure residual p_j(alpha_j) - p is solved for alpha_j, then e is assembled.
    subroutine s_jwl_ptequil_energy_pr(rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, cv_j, cv_a, e)

        $:GPU_ROUTINE(function_name='s_jwl_ptequil_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, cv_j, cv_a
        real(wp), intent(out) :: e
        real(wp)              :: rho_s, Y_s, p_s, a_lo, a_hi, a_m, rj, ra, V, ecold, T, pj, pcg, g_lo, g_hi, g_m
        integer               :: it

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)
        p_s = pres
        call s_jwl_floor_positive(p_s, sgm_eps)

        if (Y_s <= jwl_pure_cutoff) then
            e = p_s/max(air_gamma*rho_s, sgm_eps)
            call s_jwl_floor_nonnegative(e)
            return
        else if (1._wp - Y_s <= jwl_pure_cutoff) then
            call s_jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0, pcg)
            e = (p_s - pcg)/max(omega0*rho_s, sgm_eps)
            call s_jwl_floor_nonnegative(e)
            return
        end if

        a_lo = sgm_eps
        a_hi = 1._wp - sgm_eps
        rj = max(Y_s*rho_s/a_lo, sgm_eps)
        ra = max((1._wp - Y_s)*rho_s/(1._wp - a_lo), sgm_eps)
        V = rho0/rj
        T = p_s/max(air_gamma*ra*cv_a, sgm_eps)
        pj = A*exp(-R1*V) + B*exp(-R2*V) + omega0*rj*cv_j*T
        g_lo = pj - p_s

        rj = max(Y_s*rho_s/a_hi, sgm_eps)
        ra = max((1._wp - Y_s)*rho_s/(1._wp - a_hi), sgm_eps)
        V = rho0/rj
        T = p_s/max(air_gamma*ra*cv_a, sgm_eps)
        pj = A*exp(-R1*V) + B*exp(-R2*V) + omega0*rj*cv_j*T
        g_hi = pj - p_s

        ! No sign change across (0,1) means no root: fall back to pure-JWL energy.
        if (g_lo*g_hi >= 0._wp) then
            call s_jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0, pcg)
            e = (p_s - pcg)/max(omega0*rho_s, sgm_eps)
            call s_jwl_floor_nonnegative(e)
            return
        end if

        ! Bisection: a_m stays strictly inside the bracket every step, so the 1/a_m, 1/(1-a_m)
        ! phasic densities never hit their endpoint singularities.
        do it = 1, 60
            a_m = 0.5_wp*(a_lo + a_hi)
            rj = max(Y_s*rho_s/a_m, sgm_eps)
            ra = max((1._wp - Y_s)*rho_s/(1._wp - a_m), sgm_eps)
            V = rho0/rj
            T = p_s/max(air_gamma*ra*cv_a, sgm_eps)
            pj = A*exp(-R1*V) + B*exp(-R2*V) + omega0*rj*cv_j*T
            g_m = pj - p_s
            if (abs(g_m) <= jwl_root_rel_tol*p_s) exit
            if (g_lo*g_m > 0._wp) then
                a_lo = a_m
                g_lo = g_m
            else
                a_hi = a_m
            end if
        end do

        ecold = A/(R1*rho0)*exp(-R1*V) + B/(R2*rho0)*exp(-R2*V)
        e = (Y_s*cv_j + (1._wp - Y_s)*cv_a)*T + Y_s*ecold
        call s_jwl_floor_nonnegative(e)

    end subroutine s_jwl_ptequil_energy_pr

    !> Garno (2020) "Rocflu" single-fluid blended JWL/ideal-gas EOS. The cold-curve weight g_e ramps with mixture energy and the
    !! Grueneisen omega ramps with mixture density (g_rho); the blend applies only between the pure-material cutoffs, where the
    !! unmodified JWL or ideal-gas EOS is used instead.
    subroutine s_jwl_rocflu_pressure_er(rho, e, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, pres)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: pres
        real(wp)              :: rho_s, Y_s, V, g_rho, g_e, om, e0s, cab

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        if (Y_s <= jwl_pure_cutoff) then
            pres = air_gamma*rho_s*e
            call s_jwl_floor_positive(pres, sgm_eps)
            return
        else if (Y_s >= 1._wp - jwl_pure_cutoff) then
            call s_jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0, cab)
            pres = cab + omega0*rho_s*e
            call s_jwl_floor_positive(pres, sgm_eps)
            return
        end if

        V = rho0/rho_s
        g_rho = min(max((rho_s - air_rho0)/max(rho0 - air_rho0, sgm_eps), 0._wp), 1._wp)
        om = air_gamma + (omega0 - air_gamma)*g_rho
        e0s = E0/max(rho0, sgm_eps)
        g_e = min(max((e - air_e0)/max(e0s - air_e0, sgm_eps), 0._wp), 1._wp)
        call s_jwl_pcold(rho_s, A, B, R1, R2, om, rho0, cab)
        pres = g_e*cab + om*rho_s*e
        call s_jwl_floor_positive(pres, sgm_eps)

    end subroutine s_jwl_rocflu_pressure_er

    !> Specific internal energy for the Rocflu blend: inverse of s_jwl_rocflu_pressure_er. The pure-material guards mirror pressure
    !! recovery; in the mixed branch p is linear in e while the energy blend is unsaturated, so we invert linearly then correct if
    !! g_e clamps.
    subroutine s_jwl_rocflu_energy_pr(rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, e)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: e
        real(wp)              :: rho_s, Y_s, V, g_rho, om, e0s, cab, kab, g_e

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        if (Y_s <= jwl_pure_cutoff) then
            e = pres/max(air_gamma*rho_s, sgm_eps)
            call s_jwl_floor_nonnegative(e)
            return
        else if (Y_s >= 1._wp - jwl_pure_cutoff) then
            call s_jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0, cab)
            e = (pres - cab)/max(omega0*rho_s, sgm_eps)
            call s_jwl_floor_nonnegative(e)
            return
        end if

        V = rho0/rho_s
        g_rho = min(max((rho_s - air_rho0)/max(rho0 - air_rho0, sgm_eps), 0._wp), 1._wp)
        om = air_gamma + (omega0 - air_gamma)*g_rho
        e0s = E0/max(rho0, sgm_eps)
        call s_jwl_pcold(rho_s, A, B, R1, R2, om, rho0, cab)
        kab = cab/max(e0s - air_e0, sgm_eps)
        e = (pres + air_e0*kab)/max(kab + om*rho_s, sgm_eps)
        g_e = (e - air_e0)/max(e0s - air_e0, sgm_eps)
        if (g_e < 0._wp) then
            e = pres/max(om*rho_s, sgm_eps)  ! blend saturated to pure air
        else if (g_e > 1._wp) then
            e = (pres - cab)/max(om*rho_s, sgm_eps)  ! blend saturated to full products
        end if
        call s_jwl_floor_nonnegative(e)

    end subroutine s_jwl_rocflu_energy_pr

    !> Isentropic sound-speed squared for the Rocflu blended EOS (Garno 2020). Writing p = (kab + om*rho)*e - air_e0*kab in
    !! Grueneisen form gives c^2 = dp/drho|_e + dp/de|_rho * p/rho^2 directly. Unlike s_jwl_mixture_sound_speed_squared this never
    !! forms phasic densities, so it is safe where alpha_j is a blending marker rather than a volume fraction.
    subroutine s_jwl_rocflu_sound_speed_squared(rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, c2)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_sound_speed_squared',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in) :: rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: c2
        real(wp) :: rho_s, Y_s, V, g_rho, om, e0s, kab_scale, exp1, exp2, cab, kab, dom_drho, dcab_drho, dcab_dom, e, g_e

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        ! Endpoint guards: pure air or pure products bypass the blend entirely.
        if (Y_s <= jwl_pure_cutoff) then
            c2 = (air_gamma + 1._wp)*pres/rho_s
            call s_jwl_floor_positive(c2, sgm_eps)
            return
        else if (Y_s >= 1._wp - jwl_pure_cutoff) then
            call s_jwl_sound_speed_squared(rho_s, pres, A, B, R1, R2, omega0, rho0, c2)
            return
        end if

        ! Blended Gruneisen coefficient om(rho) and cold-curve cab(rho, om).
        V = rho0/rho_s
        g_rho = min(max((rho_s - air_rho0)/max(rho0 - air_rho0, sgm_eps), 0._wp), 1._wp)
        om = air_gamma + (omega0 - air_gamma)*g_rho
        e0s = E0/max(rho0, sgm_eps)
        kab_scale = max(e0s - air_e0, sgm_eps)
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        call s_jwl_pcold(rho_s, A, B, R1, R2, om, rho0, cab)
        kab = cab/kab_scale

        ! Recover e from (rho, pres) using the Rocflu inverse (matches s_jwl_rocflu_energy_pr).
        e = (pres + air_e0*kab)/max(kab + om*rho_s, sgm_eps)
        call s_jwl_floor_nonnegative(e)
        g_e = (e - air_e0)/kab_scale

        ! d(om)/d(rho): non-zero only in the g_rho interior.
        if (g_rho > sgm_eps .and. 1._wp - g_rho > sgm_eps) then
            dom_drho = (omega0 - air_gamma)/max(rho0 - air_rho0, sgm_eps)
        else
            dom_drho = 0._wp
        end if

        ! d(cab)/d(rho): frozen-om part plus om-variation part (dcab_dom * dom_drho).
        call s_jwl_dpcold_drho(rho_s, A, B, R1, R2, om, rho0, dcab_drho)
        dcab_dom = -(A*exp1/(R1*V) + B*exp2/(R2*V))
        dcab_drho = dcab_drho + dcab_dom*dom_drho

        ! c^2 = dp/drho|_e + dp/de|_rho * p/rho^2. Three cases based on g_e saturation.
        if (g_e <= 0._wp) then
            ! g_e clamped at 0: p = om*rho*e, pure Mie-Gruneisen with blended om.
            c2 = (dom_drho*rho_s + om)*e + om*pres/rho_s
        else if (g_e >= 1._wp) then
            ! g_e clamped at 1: p = cab + om*rho*e, JWL form at blended om.
            c2 = dcab_drho + (dom_drho*rho_s + om)*e + om*pres/rho_s
        else
            ! Unsaturated interior: dp/drho|_e = (dcab/drho/kab_scale)*(e-air_e0) + (dom/drho*rho+om)*e
            !                       dp/de|_rho = kab + om*rho
            c2 = (dcab_drho/kab_scale)*(e - air_e0) + (dom_drho*rho_s + om)*e + (kab + om*rho_s)*pres/rho_s**2
        end if
        call s_jwl_floor_positive(c2, sgm_eps)

    end subroutine s_jwl_rocflu_sound_speed_squared

    !> Dispatch pressure-from-energy to the active JWL closure.
    !! @param[in] rho Mixture density.
    !! @param[in] e Mixture specific internal energy.
    !! @param[in] Y Bounded JWL products mass fraction.
    !! @param[in] alpha_j Bounded JWL products volume fraction.
    !! @param[in] jidx Index into the module jwl_*s arrays.
    !! @param[out] pres Mixture pressure.
    subroutine s_jwl_mix_pressure_er(rho, e, Y, alpha_j, jidx, pres)

        $:GPU_ROUTINE(function_name='s_jwl_mix_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, alpha_j
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: pres

        if (Y <= sgm_eps .or. alpha_j <= sgm_eps) then
            pres = jwl_air_gammas(jidx)*max(rho, sgm_eps)*e
            call s_jwl_floor_positive(pres, sgm_eps)
            return
        else if (Y >= 1._wp - jwl_endpoint_tol .and. alpha_j >= 1._wp - jwl_endpoint_tol) then
            call s_jwl_pcold(max(rho, sgm_eps), jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                             & jwl_rho0s(jidx), pres)
            pres = pres + jwl_omegas(jidx)*max(rho, sgm_eps)*e
            call s_jwl_floor_positive(pres, sgm_eps)
            return
        end if

        select case (jwl_mix_type)
        case (jwl_mix_type_kuhl)
            call s_jwl_kuhl_pressure_er(rho, e, Y, alpha_j, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                        & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), &
                                        & jwl_air_rho0s(jidx), jwl_air_gammas(jidx), pres)
        case (jwl_mix_type_ptequil)
            call s_jwl_ptequil_pressure_er(rho, e, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                           & jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, pres)
        case (jwl_mix_type_rocflu)
            call s_jwl_rocflu_pressure_er(rho, e, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                          & jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                          & jwl_air_gammas(jidx), pres)
        case (jwl_mix_type_isobaric)
            call s_jwl_pressure_er(rho, e, Y, alpha_j, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                   & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                   & jwl_air_gammas(jidx), pres)
        case default
            call s_jwl_pressure_er(rho, e, Y, alpha_j, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                   & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                   & jwl_air_gammas(jidx), pres)
        end select

    end subroutine s_jwl_mix_pressure_er

    !> Dispatch energy-from-pressure to the active JWL closure (inverse of s_jwl_mix_pressure_er).
    !! @param[in] rho Mixture density.
    !! @param[in] pres Mixture pressure.
    !! @param[in] Y Bounded JWL products mass fraction.
    !! @param[in] alpha_j Bounded JWL products volume fraction.
    !! @param[in] jidx Index into the module jwl_*s arrays.
    !! @param[out] e Mixture specific internal energy.
    subroutine s_jwl_mix_energy_pr(rho, pres, Y, alpha_j, jidx, e)

        $:GPU_ROUTINE(function_name='s_jwl_mix_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: e

        if (Y <= sgm_eps .or. alpha_j <= sgm_eps) then
            e = pres/max(jwl_air_gammas(jidx)*max(rho, sgm_eps), sgm_eps)
            call s_jwl_floor_nonnegative(e)
            return
        else if (Y >= 1._wp - jwl_endpoint_tol .and. alpha_j >= 1._wp - jwl_endpoint_tol) then
            call s_jwl_pcold(max(rho, sgm_eps), jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                             & jwl_rho0s(jidx), e)
            e = (pres - e)/max(jwl_omegas(jidx)*max(rho, sgm_eps), sgm_eps)
            call s_jwl_floor_nonnegative(e)
            return
        end if

        select case (jwl_mix_type)
        case (jwl_mix_type_kuhl)
            call s_jwl_kuhl_energy_pr(rho, pres, Y, alpha_j, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                      & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                      & jwl_air_gammas(jidx), e)
        case (jwl_mix_type_ptequil)
            call s_jwl_ptequil_energy_pr(rho, pres, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                         & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, e)
        case (jwl_mix_type_rocflu)
            call s_jwl_rocflu_energy_pr(rho, pres, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                        & jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                        & jwl_air_gammas(jidx), e)
        case (jwl_mix_type_isobaric)
            call s_jwl_energy_pr(rho, pres, Y, alpha_j, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                 & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                 & jwl_air_gammas(jidx), e)
        case default
            call s_jwl_energy_pr(rho, pres, Y, alpha_j, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                 & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                 & jwl_air_gammas(jidx), e)
        end select

        ! Each closure already floors its own output (e >= 0), so no post-call mask is applied
        ! here: a NaN from a misbehaving closure propagates instead of being hidden.

    end subroutine s_jwl_mix_energy_pr

    !> Allocate and populate the per-fluid JWL parameter tables from fluid_pp, identify the JWL fluid (jwl_idx), and cache the
    !! products/air specific heats used by the p-T-equilibrium closure. Called from the variables-conversion module initializer.
    impure subroutine s_initialize_jwl_module

        use m_mpi_common, only: s_mpi_abort
        use m_helper_basic, only: f_is_default

        integer :: i, n_jwl

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
                ! Fail fast on the dflt_real sentinel: unset JWL parameters would otherwise
                ! propagate as -1e6 through exp() and produce NaNs deep in the solver.
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
            end if
        end do

        if (n_jwl > 1) then
            call s_mpi_abort('At most one fluid may use eos_jwl; found more than one.')
        end if

        if (jwl_idx > 0 .and. model_eqns /= model_eqns_5eq) then
            call s_mpi_abort('eos_jwl is only supported with model_eqns_5eq.')
        end if

        ! Specific heats for the p-T-equilibrium closure: products from the JWL fluid, air from the first ideal-gas fluid.
        jwl_cv_prod = 0._wp
        jwl_cv_air = 0._wp
        if (jwl_idx > 0) jwl_cv_prod = fluid_pp(jwl_idx)%cv
        do i = 1, num_fluids
            if (fluid_pp(i)%eos == eos_stiffened_gas) then
                jwl_cv_air = fluid_pp(i)%cv
                exit
            end if
        end do

        ! The Kuhl and p-T-equilibrium closures split energy via specific heats; an unset
        ! cv (dflt_real sentinel) or nonpositive cv silently corrupts T and p.
        if (jwl_idx > 0 .and. (jwl_mix_type == jwl_mix_type_kuhl .or. jwl_mix_type == jwl_mix_type_ptequil)) then
            if (f_is_default(jwl_cv_prod) .or. jwl_cv_prod <= 0._wp .or. f_is_default(jwl_cv_air) .or. jwl_cv_air <= 0._wp) then
                call s_mpi_abort('jwl_mix_type_kuhl and jwl_mix_type_ptequil require positive ' &
                                 & // 'fluid_pp%cv for both the JWL fluid and the ideal-gas fluid.')
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
