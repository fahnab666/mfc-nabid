!>
!! @file
!! @brief JWL EOS and mixture closures.

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief JWL equation of state and mixture closures.
!!
!! Closure path: rho_m, e_m, Y(:), component_ids(:) -> p_m, T_m, c_m.
!! Mixture rules: isobaric, Kuhl & Khasainov (2007) piece-wise quadratic,
!! p-T equilibrium, and Rocflu blend. Helpers below own the cold curve,
!! mass fractions, Kuhl coefficients, temperature inversion, and sound speed.
module m_jwl

    use m_global_parameters
    use m_constants, only: jwl_mix_type_isobaric, jwl_mix_type_kuhl, jwl_mix_type_ptequil, jwl_mix_type_rocflu

    implicit none

    private
    public :: s_initialize_jwl_module, s_finalize_jwl_module, s_jwl_mixture_sound_speed_squared, &
        & s_jwl_rocflu_sound_speed_squared, s_jwl_mix_pressure_er, s_jwl_mix_energy_pr, s_jwl_mix_sound_speed, jwl_idx, &
        & s_jwl_mixture_closure

    ! Simulation builds use m_global_parameters tables.
#ifndef MFC_SIMULATION
    real(wp), allocatable, public, dimension(:) :: jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s
    real(wp), allocatable, public, dimension(:) :: jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas
    $:GPU_DECLARE(create='[jwl_As, jwl_Bs, jwl_R1s, jwl_R2s, jwl_omegas, jwl_rho0s, jwl_E0s]')
    $:GPU_DECLARE(create='[jwl_air_e0s, jwl_air_rho0s, jwl_air_gammas]')
#endif

    !> Pure-state cutoff.
    real(wp), parameter :: jwl_pure_cutoff = 1.e-4_wp
    !> Root tolerance.
    real(wp), parameter :: jwl_root_rel_tol = max(1.e-12_wp, 100._wp*epsilon(1._wp))
    !> Endpoint tolerance.
    real(wp), parameter :: jwl_endpoint_tol = max(1.e-12_wp, 100._wp*epsilon(1._wp))
    integer             :: jwl_idx                  !< JWL fluid index.
    real(wp)            :: jwl_cv_prod, jwl_cv_air  !< Products/air specific heats.
    $:GPU_DECLARE(create='[jwl_idx]')
    $:GPU_DECLARE(create='[jwl_cv_prod, jwl_cv_air]')

    ! Kuhl piece-wise constants; component ids match the source tables.
    real(wp), parameter :: kpw_T1 = 2340._wp
    real(wp), parameter :: kpw_T2 = 3700._wp
    real(wp), parameter :: kpw_T3 = 4150._wp
    real(wp), parameter :: kpw_T4 = 4530._wp
    real(wp), parameter :: kpw_cal2j = 4184._wp
    integer, parameter  :: kpw_air_id = 1
    integer, parameter  :: kpw_products_id = 2
    integer, parameter  :: kpw_products_alt_id = 4
    ! R_k = 8314/MW_k.
    real(wp), parameter :: kpw_R(6) = [8314._wp/28.85_wp, 8314._wp/28.76_wp, 8314._wp/27.75_wp, 8314._wp/26.93_wp, &
         & 8314._wp/29.65_wp, 8314._wp/40.78_wp]
    ! cal/g -> J/kg.
    real(wp), parameter :: kpw_a(5, 6) = kpw_cal2j*reshape([2.02768e-5_wp, 1.34322e-4_wp, 7.01281e-5_wp, -1.02084e-4_wp, &
         & 4.04923e-5_wp, 3.31674e-5_wp, 5.97088e-5_wp, 1.9052e-4_wp, 2.28177e-4_wp, 1.78281e-4_wp, 4.745e-5_wp, 4.6038e-4_wp, &
         & 4.9083e-4_wp, -6.1549e-4_wp, -2.8216e-4_wp, 5.3244e-5_wp, 7.9903e-5_wp, 0._wp, 4.5108e-4_wp, 2.578e-3_wp, &
         & 3.5282e-6_wp, 2.5302e-4_wp, -6.1238e-5_wp, -3.9217e-4_wp, 2.7654e-5_wp, 1.76153e-5_wp, 1.49115e-5_wp, 1.13e-3_wp, &
         & 8.26e-3_wp, 5.03544e-5_wp], shape=[5, 6])
    real(wp), parameter :: kpw_b(5, 6) = kpw_cal2j*reshape([0.16498_wp, -0.41045_wp, 0.11507_wp, 1.53731_wp, 0.11381_wp, &
         & 0.20867_wp, 0.0377_wp, -0.89226_wp, -1.20053_wp, -0.774255_wp, 0.1549_wp, -1.7722_wp, -1.841_wp, 7.3463_wp, 3.8022_wp, &
         & 0.17393_wp, 0.035886_wp, 1.80555_wp, -2.7713_wp, -22.917_wp, 0.25361_wp, -0.80169_wp, 1.5345_wp, 4.2413_wp, 0.2432_wp, &
         & 0.20186_wp, 0.2502_wp, -7.95255_wp, -67.29752_wp, -0.07059_wp], shape=[5, 6])
    real(wp), parameter :: kpw_c(5, 6) = kpw_cal2j*reshape([-71.9172_wp, 658.24424_wp, -403.36139_wp, -3340.674_wp, 198.38643_wp, &
         & -1890.164_wp, -1634.868_wp, 20.04935_wp, 651.0422_wp, -248.616_wp, -1555.6_wp, 711.74_wp, 558.87_wp, -18515.0_wp, &
         & -9254.5_wp, -941.33_wp, -760.12_wp, -6211.8_wp, 5014._wp, 52697._wp, -949.3_wp, 168.08_wp, -4178._wp, -9713.6_wp, &
         & -195._wp, -1553.62_wp, -1554.5182_wp, 13553.8_wp, 137084.51_wp, 1216.0279_wp], shape=[5, 6])

contains

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

    !> JWL cold pressure.
    subroutine s_jwl_pcold(rho, A, B, R1, R2, omega, rho0, pcold)

        $:GPU_ROUTINE(function_name='s_jwl_pcold',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, A, B, R1, R2, omega, rho0
        real(wp), intent(out) :: pcold
        real(wp)              :: V

        V = rho0/max(rho, sgm_eps)
        pcold = A*(1._wp - omega/(R1*V))*exp(-R1*V) + B*(1._wp - omega/(R2*V))*exp(-R2*V)

    end subroutine s_jwl_pcold

    !> JWL cold pressure and density derivative.
    subroutine s_jwl_pcold_dpcold(rho, A, B, R1, R2, omega, rho0, pcold, dpcold_drho)

        $:GPU_ROUTINE(function_name='s_jwl_pcold_dpcold',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, A, B, R1, R2, omega, rho0
        real(wp), intent(out) :: pcold, dpcold_drho
        real(wp)              :: V, rho_safe, invV, invV2, exp1, exp2, Vrho

        rho_safe = max(rho, sgm_eps)
        V = rho0/rho_safe
        invV = 1._wp/V
        invV2 = invV*invV
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        Vrho = V/rho_safe

        pcold = A*(1._wp - omega*invV/R1)*exp1 + B*(1._wp - omega*invV/R2)*exp2
        dpcold_drho = A*exp1*Vrho*(R1 - omega*invV - omega*invV2/R1) + B*exp2*Vrho*(R2 - omega*invV - omega*invV2/R2)

    end subroutine s_jwl_pcold_dpcold

    !> Isentropic sound-speed squared for a single-fluid JWL state (Rocflu/Stanley form).
    subroutine s_jwl_sound_speed_squared(rho, pres, A, B, R1, R2, omega, rho0, c2)

        $:GPU_ROUTINE(function_name='s_jwl_sound_speed_squared',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, A, B, R1, R2, omega, rho0
        real(wp), intent(out) :: c2
        real(wp)              :: rho_safe, pcold, dpcold_drho, e

        rho_safe = max(rho, sgm_eps)
        call s_jwl_pcold_dpcold(rho_safe, A, B, R1, R2, omega, rho0, pcold, dpcold_drho)

        e = (pres - pcold)/(omega*rho_safe)
        c2 = dpcold_drho + omega*e + omega*pres/rho_safe
        call s_jwl_floor_positive(c2, sgm_eps)

    end subroutine s_jwl_sound_speed_squared

    !> Frozen mass-weighted JWL/ideal-gas mixture sound speed squared.
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

        c2 = c2_air + Y_safe*(c2_jwl - c2_air)
        call s_jwl_floor_positive(c2, sgm_eps)

    end subroutine s_jwl_mixture_sound_speed_squared

    !> Closed-form isobaric JWL/ideal-gas mixture pressure.
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

    !> Inverse of s_jwl_pressure_er.
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
        call s_jwl_floor_nonnegative(e)

    end subroutine s_jwl_energy_pr

    !> p-T-equilibrium mixture pressure.
    subroutine s_jwl_ptequil_pressure_er(rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, cv_j, cv_a, pres)

        $:GPU_ROUTINE(function_name='s_jwl_ptequil_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, A, B, R1, R2, omega0, rho0, air_gamma, cv_j, cv_a
        real(wp), intent(out) :: pres
        real(wp)              :: rho_s, Y_s, cv_mix, a_lo, a_hi, a_m, rj, ra, V, ecold, T, pj, pa, f_lo, f_m, pcg
        real(wp)              :: exp1, exp2, Aexp1, Bexp2
        integer               :: it

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

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

        cv_mix = max(cv_a + Y_s*(cv_j - cv_a), sgm_eps)

        a_lo = sgm_eps
        a_hi = 1._wp - sgm_eps
        rj = max(Y_s*rho_s/a_lo, sgm_eps)
        ra = max((1._wp - Y_s)*rho_s/(1._wp - a_lo), sgm_eps)
        V = rho0/rj
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        Aexp1 = A*exp1
        Bexp2 = B*exp2
        ecold = Aexp1/(R1*rho0) + Bexp2/(R2*rho0)
        T = (e - Y_s*ecold)/cv_mix
        call s_jwl_floor_positive(T, sgm_eps)
        pj = Aexp1 + Bexp2 + omega0*rj*cv_j*T
        pa = air_gamma*ra*cv_a*T
        f_lo = pj - pa
        pres = 0.5_wp*(pj + pa)

        do it = 1, 60
            a_m = 0.5_wp*(a_lo + a_hi)
            rj = max(Y_s*rho_s/a_m, sgm_eps)
            ra = max((1._wp - Y_s)*rho_s/(1._wp - a_m), sgm_eps)
            V = rho0/rj
            exp1 = exp(-R1*V)
            exp2 = exp(-R2*V)
            Aexp1 = A*exp1
            Bexp2 = B*exp2
            ecold = Aexp1/(R1*rho0) + Bexp2/(R2*rho0)
            T = (e - Y_s*ecold)/cv_mix
            call s_jwl_floor_positive(T, sgm_eps)
            pj = Aexp1 + Bexp2 + omega0*rj*cv_j*T
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

    !> Inverse of s_jwl_ptequil_pressure_er.
    subroutine s_jwl_ptequil_energy_pr(rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, cv_j, cv_a, e)

        $:GPU_ROUTINE(function_name='s_jwl_ptequil_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, A, B, R1, R2, omega0, rho0, air_gamma, cv_j, cv_a
        real(wp), intent(out) :: e
        real(wp)              :: rho_s, Y_s, p_s, a_lo, a_hi, a_m, rj, ra, V, ecold, T, pj, pcg, g_lo, g_hi, g_m
        real(wp)              :: exp1, exp2, Aexp1, Bexp2
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
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        Aexp1 = A*exp1
        Bexp2 = B*exp2
        pj = Aexp1 + Bexp2 + omega0*rj*cv_j*T
        g_lo = pj - p_s

        rj = max(Y_s*rho_s/a_hi, sgm_eps)
        ra = max((1._wp - Y_s)*rho_s/(1._wp - a_hi), sgm_eps)
        V = rho0/rj
        T = p_s/max(air_gamma*ra*cv_a, sgm_eps)
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        Aexp1 = A*exp1
        Bexp2 = B*exp2
        pj = Aexp1 + Bexp2 + omega0*rj*cv_j*T
        g_hi = pj - p_s

        if (g_lo*g_hi >= 0._wp) then
            call s_jwl_pcold(rho_s, A, B, R1, R2, omega0, rho0, pcg)
            e = (p_s - pcg)/max(omega0*rho_s, sgm_eps)
            call s_jwl_floor_nonnegative(e)
            return
        end if

        do it = 1, 60
            a_m = 0.5_wp*(a_lo + a_hi)
            rj = max(Y_s*rho_s/a_m, sgm_eps)
            ra = max((1._wp - Y_s)*rho_s/(1._wp - a_m), sgm_eps)
            V = rho0/rj
            T = p_s/max(air_gamma*ra*cv_a, sgm_eps)
            exp1 = exp(-R1*V)
            exp2 = exp(-R2*V)
            Aexp1 = A*exp1
            Bexp2 = B*exp2
            pj = Aexp1 + Bexp2 + omega0*rj*cv_j*T
            g_m = pj - p_s
            if (abs(g_m) <= jwl_root_rel_tol*p_s) exit
            if (g_lo*g_m > 0._wp) then
                a_lo = a_m
                g_lo = g_m
            else
                a_hi = a_m
            end if
        end do

        ecold = Aexp1/(R1*rho0) + Bexp2/(R2*rho0)
        e = (cv_a + Y_s*(cv_j - cv_a))*T + Y_s*ecold
        call s_jwl_floor_nonnegative(e)

    end subroutine s_jwl_ptequil_energy_pr

    !> Garno/Rocflu single-fluid blended JWL/ideal-gas pressure.
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

    !> Inverse of s_jwl_rocflu_pressure_er.
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
            e = pres/max(om*rho_s, sgm_eps)
        else if (g_e > 1._wp) then
            e = (pres - cab)/max(om*rho_s, sgm_eps)
        end if
        call s_jwl_floor_nonnegative(e)

    end subroutine s_jwl_rocflu_energy_pr

    !> Isentropic sound-speed squared for the Rocflu blended EOS.
    subroutine s_jwl_rocflu_sound_speed_squared(rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma, c2)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_sound_speed_squared',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in) :: rho, pres, Y, A, B, R1, R2, omega0, rho0, E0, air_e0, air_rho0, air_gamma
        real(wp), intent(out) :: c2
        real(wp) :: rho_s, Y_s, V, invV, invV2, Vrho, g_rho, om, e0s, kab_scale, exp1, exp2, cab, kab, dom_drho, dcab_drho, &
             & dcab_dom, e, g_e

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)

        if (Y_s <= jwl_pure_cutoff) then
            c2 = (air_gamma + 1._wp)*pres/rho_s
            call s_jwl_floor_positive(c2, sgm_eps)
            return
        else if (Y_s >= 1._wp - jwl_pure_cutoff) then
            call s_jwl_sound_speed_squared(rho_s, pres, A, B, R1, R2, omega0, rho0, c2)
            return
        end if

        V = rho0/rho_s
        g_rho = min(max((rho_s - air_rho0)/max(rho0 - air_rho0, sgm_eps), 0._wp), 1._wp)
        om = air_gamma + (omega0 - air_gamma)*g_rho
        e0s = E0/max(rho0, sgm_eps)
        kab_scale = max(e0s - air_e0, sgm_eps)
        exp1 = exp(-R1*V)
        exp2 = exp(-R2*V)
        invV = 1._wp/V
        invV2 = invV*invV
        Vrho = V/rho_s
        cab = A*(1._wp - om*invV/R1)*exp1 + B*(1._wp - om*invV/R2)*exp2
        kab = cab/kab_scale

        e = (pres + air_e0*kab)/max(kab + om*rho_s, sgm_eps)
        call s_jwl_floor_nonnegative(e)
        g_e = (e - air_e0)/kab_scale

        if (g_rho > sgm_eps .and. 1._wp - g_rho > sgm_eps) then
            dom_drho = (omega0 - air_gamma)/max(rho0 - air_rho0, sgm_eps)
        else
            dom_drho = 0._wp
        end if

        dcab_drho = A*exp1*Vrho*(R1 - om*invV - om*invV2/R1) + B*exp2*Vrho*(R2 - om*invV - om*invV2/R2)
        dcab_dom = -(A*exp1/(R1*V) + B*exp2/(R2*V))
        dcab_drho = dcab_drho + dcab_dom*dom_drho

        if (g_e <= 0._wp) then
            c2 = (dom_drho*rho_s + om)*e + om*pres/rho_s
        else if (g_e >= 1._wp) then
            c2 = dcab_drho + (dom_drho*rho_s + om)*e + om*pres/rho_s
        else
            c2 = (dcab_drho/kab_scale)*(e - air_e0) + (dom_drho*rho_s + om)*e + (kab + om*rho_s)*pres/(rho_s*rho_s)
        end if
        call s_jwl_floor_positive(c2, sgm_eps)

    end subroutine s_jwl_rocflu_sound_speed_squared

    !> Kuhl piece-wise cold-pressure and thermal-coefficient sums.
    subroutine s_kuhl_pw_cold_sums(rho_s, Y, comp_ids, ncomp, jidx, p_cold_sum, p_thermal_coeff, dp_cold_sum)

        $:GPU_ROUTINE(function_name='s_kuhl_pw_cold_sums',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho_s, Y(ncomp)
        integer, intent(in)   :: comp_ids(ncomp), ncomp, jidx
        real(wp), intent(out) :: p_cold_sum, p_thermal_coeff, dp_cold_sum
        real(wp)              :: rho_k, pcold_k, dpcold_k
        integer               :: k, cid

        p_cold_sum = 0._wp
        p_thermal_coeff = 0._wp
        dp_cold_sum = 0._wp
        do k = 1, ncomp
            cid = comp_ids(k)
            rho_k = Y(k)*rho_s
            p_thermal_coeff = p_thermal_coeff + rho_k*kpw_R(cid)
            if (cid == kpw_products_id .or. cid == kpw_products_alt_id) then
                call s_jwl_pcold_dpcold(rho_k, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                        & jwl_rho0s(jidx), pcold_k, dpcold_k)
                p_cold_sum = p_cold_sum + pcold_k
                dp_cold_sum = dp_cold_sum + Y(k)*Y(k)*dpcold_k
            end if
        end do

    end subroutine s_kuhl_pw_cold_sums

    !> Kuhl & Khasainov mixture coefficients.
    subroutine s_kuhl_pw_mix_coeffs(Y, ncomp, comp_ids, reg, a_m, b_m, c_m, R_m)

        $:GPU_ROUTINE(function_name='s_kuhl_pw_mix_coeffs',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: Y(ncomp)
        integer, intent(in)   :: ncomp, comp_ids(ncomp), reg
        real(wp), intent(out) :: a_m, b_m, c_m, R_m
        integer               :: k, cid

        a_m = 0._wp
        b_m = 0._wp
        c_m = 0._wp
        R_m = 0._wp
        do k = 1, ncomp
            cid = comp_ids(k)
            a_m = a_m + Y(k)*kpw_a(reg, cid)
            b_m = b_m + Y(k)*kpw_b(reg, cid)
            c_m = c_m + Y(k)*kpw_c(reg, cid)
            R_m = R_m + Y(k)*kpw_R(cid)
        end do

    end subroutine s_kuhl_pw_mix_coeffs

    !> Invert the Kuhl piece-wise quadratic caloric EOS.
    subroutine s_kuhl_pw_invert_T(e_m, Y, ncomp, comp_ids, T_m, reg_m, a_m, b_m, c_m, R_m)

        $:GPU_ROUTINE(function_name='s_kuhl_pw_invert_T',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: e_m
        real(wp), intent(in)  :: Y(ncomp)
        integer, intent(in)   :: ncomp, comp_ids(ncomp)
        real(wp), intent(out) :: T_m
        integer, intent(out)  :: reg_m
        real(wp), intent(out) :: a_m, b_m, c_m, R_m
        real(wp)              :: a_c, b_c, c_c, R_c, disc, T_cand, T_lo, T_hi
        integer               :: reg
        real(wp), parameter   :: T_bounds(6) = [300._wp, kpw_T1, kpw_T2, kpw_T3, kpw_T4, 6000._wp]

        T_m = 300._wp
        reg_m = 1
        call s_kuhl_pw_mix_coeffs(Y, ncomp, comp_ids, reg_m, a_m, b_m, c_m, R_m)
        do reg = 1, 5
            if (reg == reg_m) then
                a_c = a_m
                b_c = b_m
                c_c = c_m
                R_c = R_m
            else
                call s_kuhl_pw_mix_coeffs(Y, ncomp, comp_ids, reg, a_c, b_c, c_c, R_c)
            end if
            if (abs(a_c) < sgm_eps) then
                if (abs(b_c) > sgm_eps) then
                    T_cand = (e_m - c_c)/b_c
                else
                    cycle
                end if
            else
                disc = b_c*b_c - 4._wp*a_c*(c_c - e_m)
                if (disc < 0._wp) cycle
                T_cand = (-b_c + sqrt(disc))/(2._wp*a_c)
            end if
            T_lo = T_bounds(reg)
            T_hi = T_bounds(reg + 1)
            ! Region 5 is the last: allow extrapolation above 6000 K.
            if (T_cand >= T_lo .and. (T_cand <= T_hi .or. reg == 5)) then
                T_m = T_cand
                reg_m = reg
                a_m = a_c
                b_m = b_c
                c_m = c_c
                R_m = R_c
                return
            end if
        end do
        T_m = max(T_m, 300._wp)

    end subroutine s_kuhl_pw_invert_T

    !> Kuhl piece-wise quadratic pressure from density and energy.
    subroutine s_kuhl_pw_pressure_er(rho_m, e_m, Y, comp_ids, ncomp, jidx, p_m)

        $:GPU_ROUTINE(function_name='s_kuhl_pw_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho_m, e_m, Y(ncomp)
        integer, intent(in)   :: comp_ids(ncomp), ncomp, jidx
        real(wp), intent(out) :: p_m
        real(wp)              :: rho_s, T_m, p_cold_sum, p_thermal_coeff, dp_cold_sum, a_m, b_m, c_m, R_m
        integer               :: reg

        rho_s = max(rho_m, sgm_eps)
        call s_kuhl_pw_invert_T(e_m, Y, ncomp, comp_ids, T_m, reg, a_m, b_m, c_m, R_m)
        call s_jwl_floor_positive(T_m, 1._wp)
        call s_kuhl_pw_cold_sums(rho_s, Y, comp_ids, ncomp, jidx, p_cold_sum, p_thermal_coeff, dp_cold_sum)

        p_m = p_cold_sum + p_thermal_coeff*T_m
        call s_jwl_floor_positive(p_m, sgm_eps)

    end subroutine s_kuhl_pw_pressure_er

    !> Inverse of s_kuhl_pw_pressure_er.
    subroutine s_kuhl_pw_energy_pr(rho_m, p_m, Y, comp_ids, ncomp, jidx, e_m)

        $:GPU_ROUTINE(function_name='s_kuhl_pw_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho_m, p_m, Y(ncomp)
        integer, intent(in)   :: comp_ids(ncomp), ncomp, jidx
        real(wp), intent(out) :: e_m
        real(wp)              :: rho_s, p_s, T_m, p_cold_sum, p_thermal_coeff, dp_cold_sum, a_m, b_m, c_m, R_m
        integer               :: reg

        rho_s = max(rho_m, sgm_eps)
        p_s = p_m
        call s_jwl_floor_positive(p_s, sgm_eps)
        call s_kuhl_pw_cold_sums(rho_s, Y, comp_ids, ncomp, jidx, p_cold_sum, p_thermal_coeff, dp_cold_sum)

        T_m = (p_s - p_cold_sum)/max(p_thermal_coeff, sgm_eps)
        call s_jwl_floor_positive(T_m, 1._wp)

        if (T_m < kpw_T1) then
            reg = 1
        else if (T_m < kpw_T2) then
            reg = 2
        else if (T_m < kpw_T3) then
            reg = 3
        else if (T_m < kpw_T4) then
            reg = 4
        else
            reg = 5
        end if

        call s_kuhl_pw_mix_coeffs(Y, ncomp, comp_ids, reg, a_m, b_m, c_m, R_m)
        e_m = (a_m*T_m + b_m)*T_m + c_m
        call s_jwl_floor_nonnegative(e_m)

    end subroutine s_kuhl_pw_energy_pr

    !> Clamp component mass fractions to [0,1] and renormalize to sum to one.
    subroutine s_jwl_clamp_mass_fractions(Y_in, ncomp, Y_out)

        $:GPU_ROUTINE(function_name='s_jwl_clamp_mass_fractions',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: Y_in(ncomp)
        integer, intent(in)   :: ncomp
        real(wp), intent(out) :: Y_out(ncomp)
        real(wp)              :: Y_sum
        integer               :: k

        Y_sum = 0._wp
        do k = 1, ncomp
            Y_out(k) = min(max(Y_in(k), 0._wp), 1._wp)
            Y_sum = Y_sum + Y_out(k)
        end do

        if (Y_sum > sgm_eps) then
            do k = 1, ncomp
                Y_out(k) = Y_out(k)/Y_sum
            end do
        end if

    end subroutine s_jwl_clamp_mass_fractions

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

    !> Black-box isobaric closure: (rho, e, Y) -> (p, T, c).
    subroutine s_jwl_isobaric_closure(rho_m, e_m, Y, component_ids, ncomp, jidx, p_m, T_m, c_m)

        $:GPU_ROUTINE(function_name='s_jwl_isobaric_closure',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho_m, e_m, Y(ncomp)
        integer, intent(in)   :: component_ids(ncomp), ncomp, jidx
        real(wp), intent(out) :: p_m, T_m, c_m
        real(wp)              :: rho_s, Y_s, R_m, pcold_k, c2

        rho_s = max(rho_m, sgm_eps)
        call s_jwl_products_mass_fraction(Y, component_ids, ncomp, Y_s)
        R_m = kpw_R(kpw_air_id) + Y_s*(jwl_omegas(jidx)*kpw_R(kpw_products_id) - kpw_R(kpw_air_id))

        call s_jwl_pressure_er(rho_s, e_m, Y_s, Y_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                               & jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), jwl_air_gammas(jidx), p_m)

        call s_jwl_pcold(Y_s*rho_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), jwl_rho0s(jidx), &
                         & pcold_k)
        T_m = (p_m - Y_s*pcold_k)/max(rho_s*R_m, sgm_eps)
        call s_jwl_floor_positive(T_m, sgm_eps)

        call s_jwl_mixture_sound_speed_squared(rho_s, p_m, Y_s, Y_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                               & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), &
                                               & jwl_air_rho0s(jidx), jwl_air_gammas(jidx), c2)
        c_m = sqrt(max(c2, sgm_eps))

    end subroutine s_jwl_isobaric_closure

    !> Black-box Kuhl & Khasainov (2007) piece-wise closure: (rho, e, Y) -> (p, T, c). Returns the additive-pressure sound speed c2
    !! = (p_thermal + Gamma*p)/rho + dp_cold.
    subroutine s_kuhl_pw_closure(rho_m, e_m, Y, component_ids, ncomp, jidx, p_m, T_m, c_m)

        $:GPU_ROUTINE(function_name='s_kuhl_pw_closure',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho_m, e_m, Y(ncomp)
        integer, intent(in)   :: component_ids(ncomp), ncomp, jidx
        real(wp), intent(out) :: p_m, T_m, c_m
        real(wp)              :: Y_work(ncomp), rho_s, a_m, b_m, c_coef, R_m, cv_m, Gamma_m
        real(wp)              :: p_cold_sum, p_thermal, p_thermal_coeff, dp_cold_sum, c2_m
        integer               :: reg

        call s_jwl_clamp_mass_fractions(Y, ncomp, Y_work)
        rho_s = max(rho_m, sgm_eps)

        call s_kuhl_pw_invert_T(e_m, Y_work, ncomp, component_ids, T_m, reg, a_m, b_m, c_coef, R_m)
        call s_jwl_floor_positive(T_m, 1._wp)
        call s_kuhl_pw_cold_sums(rho_s, Y_work, component_ids, ncomp, jidx, p_cold_sum, p_thermal_coeff, dp_cold_sum)

        p_thermal = p_thermal_coeff*T_m
        p_m = p_cold_sum + p_thermal
        call s_jwl_floor_positive(p_m, sgm_eps)

        cv_m = max(2._wp*a_m*T_m + b_m, sgm_eps)
        Gamma_m = R_m/cv_m
        c2_m = (p_thermal + Gamma_m*p_m)/rho_s + dp_cold_sum
        call s_jwl_floor_positive(c2_m, sgm_eps)
        c_m = sqrt(c2_m)

    end subroutine s_kuhl_pw_closure

    !> Black-box p-T equilibrium closure: (rho, e, Y) -> (p, T, c).
    subroutine s_jwl_ptequil_closure(rho_m, e_m, Y, component_ids, ncomp, jidx, p_m, T_m, c_m)

        $:GPU_ROUTINE(function_name='s_jwl_ptequil_closure',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho_m, e_m, Y(ncomp)
        integer, intent(in)   :: component_ids(ncomp), ncomp, jidx
        real(wp), intent(out) :: p_m, T_m, c_m
        real(wp)              :: rho_s, Y_s, V, ecold, cv_m, c2

        rho_s = max(rho_m, sgm_eps)
        call s_jwl_products_mass_fraction(Y, component_ids, ncomp, Y_s)

        call s_jwl_ptequil_pressure_er(rho_s, e_m, Y_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                       & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, p_m)

        V = jwl_rho0s(jidx)/rho_s
        ecold = jwl_As(jidx)/(jwl_R1s(jidx)*jwl_rho0s(jidx))*exp(-jwl_R1s(jidx)*V) + jwl_Bs(jidx)/(jwl_R2s(jidx)*jwl_rho0s(jidx)) &
                       & *exp(-jwl_R2s(jidx)*V)
        cv_m = max(Y_s*jwl_cv_prod + (1._wp - Y_s)*jwl_cv_air, sgm_eps)
        T_m = (e_m - Y_s*ecold)/cv_m
        call s_jwl_floor_positive(T_m, sgm_eps)

        call s_jwl_mixture_sound_speed_squared(rho_s, p_m, Y_s, Y_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                               & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), &
                                               & jwl_air_rho0s(jidx), jwl_air_gammas(jidx), c2)
        c_m = sqrt(max(c2, sgm_eps))

    end subroutine s_jwl_ptequil_closure

    !> Black-box Rocflu blend closure: (rho, e, Y) -> (p, T, c).
    subroutine s_jwl_rocflu_closure(rho_m, e_m, Y, component_ids, ncomp, jidx, p_m, T_m, c_m)

        $:GPU_ROUTINE(function_name='s_jwl_rocflu_closure',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho_m, e_m, Y(ncomp)
        integer, intent(in)   :: component_ids(ncomp), ncomp, jidx
        real(wp), intent(out) :: p_m, T_m, c_m
        real(wp)              :: rho_s, Y_s, R_m, pcold_k, c2

        rho_s = max(rho_m, sgm_eps)
        call s_jwl_products_mass_fraction(Y, component_ids, ncomp, Y_s)
        R_m = kpw_R(kpw_air_id) + Y_s*(jwl_omegas(jidx)*kpw_R(kpw_products_id) - kpw_R(kpw_air_id))

        call s_jwl_rocflu_pressure_er(rho_s, e_m, Y_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                      & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                      & jwl_air_gammas(jidx), p_m)

        call s_jwl_pcold(Y_s*rho_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), jwl_rho0s(jidx), &
                         & pcold_k)
        T_m = (p_m - Y_s*pcold_k)/max(rho_s*R_m, sgm_eps)
        call s_jwl_floor_positive(T_m, sgm_eps)

        call s_jwl_rocflu_sound_speed_squared(rho_s, p_m, Y_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                              & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), &
                                              & jwl_air_rho0s(jidx), jwl_air_gammas(jidx), c2)
        c_m = sqrt(max(c2, sgm_eps))

    end subroutine s_jwl_rocflu_closure

    !> Black-box mixture closure dispatcher: (rho, e, Y, component_ids) -> (p, T, c).
    subroutine s_jwl_mixture_closure(rho_m, e_m, Y, component_ids, ncomp, jidx, p_m, T_m, c_m)

        $:GPU_ROUTINE(function_name='s_jwl_mixture_closure',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho_m, e_m, Y(ncomp)
        integer, intent(in)   :: component_ids(ncomp), ncomp, jidx
        real(wp), intent(out) :: p_m, T_m, c_m

        select case (jwl_mix_type)
        case (jwl_mix_type_kuhl)
            call s_kuhl_pw_closure(rho_m, e_m, Y, component_ids, ncomp, jidx, p_m, T_m, c_m)
        case (jwl_mix_type_ptequil)
            call s_jwl_ptequil_closure(rho_m, e_m, Y, component_ids, ncomp, jidx, p_m, T_m, c_m)
        case (jwl_mix_type_rocflu)
            call s_jwl_rocflu_closure(rho_m, e_m, Y, component_ids, ncomp, jidx, p_m, T_m, c_m)
        case default
            call s_jwl_isobaric_closure(rho_m, e_m, Y, component_ids, ncomp, jidx, p_m, T_m, c_m)
        end select

    end subroutine s_jwl_mixture_closure

    !> Dispatch pressure-from-energy to the active JWL closure.
    subroutine s_jwl_mix_pressure_er(rho, e, Y, alpha_j, jidx, pres)

        $:GPU_ROUTINE(function_name='s_jwl_mix_pressure_er',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, e, Y, alpha_j
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: pres
        real(wp)              :: Y_mix(2), Y_s
        integer               :: component_ids(2)

        select case (jwl_mix_type)
        case (jwl_mix_type_kuhl)
            Y_s = min(max(Y, 0._wp), 1._wp)
            Y_mix(1) = Y_s
            Y_mix(2) = 1._wp - Y_s
            component_ids(1) = kpw_products_id
            component_ids(2) = kpw_air_id
            call s_kuhl_pw_pressure_er(rho, e, Y_mix, component_ids, 2, jidx, pres)
        case (jwl_mix_type_ptequil)
            call s_jwl_ptequil_pressure_er(rho, e, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                           & jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, pres)
        case (jwl_mix_type_rocflu)
            call s_jwl_rocflu_pressure_er(rho, e, Y, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), &
                                          & jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                          & jwl_air_gammas(jidx), pres)
        case default
            call s_jwl_pressure_er(rho, e, Y, alpha_j, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                   & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                   & jwl_air_gammas(jidx), pres)
        end select

    end subroutine s_jwl_mix_pressure_er

    !> Dispatch energy-from-pressure to the active JWL closure.
    subroutine s_jwl_mix_energy_pr(rho, pres, Y, alpha_j, jidx, e)

        $:GPU_ROUTINE(function_name='s_jwl_mix_energy_pr',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: e
        real(wp)              :: rho_s, Y_s, alpha_s, Y_mix(2)
        integer               :: component_ids(2)

        rho_s = max(rho, sgm_eps)
        Y_s = min(max(Y, 0._wp), 1._wp)
        alpha_s = min(max(alpha_j, 0._wp), 1._wp)

        if (Y_s <= sgm_eps .or. alpha_s <= sgm_eps) then
            e = pres/max(jwl_air_gammas(jidx)*rho_s, sgm_eps)
            call s_jwl_floor_nonnegative(e)
            return
        else if (Y_s >= 1._wp - jwl_endpoint_tol .and. alpha_s >= 1._wp - jwl_endpoint_tol) then
            call s_jwl_pcold(rho_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), jwl_omegas(jidx), jwl_rho0s(jidx), e)
            e = (pres - e)/max(jwl_omegas(jidx)*rho_s, sgm_eps)
            call s_jwl_floor_nonnegative(e)
            return
        end if

        select case (jwl_mix_type)
        case (jwl_mix_type_kuhl)
            Y_mix(1) = Y_s
            Y_mix(2) = 1._wp - Y_s
            component_ids(1) = kpw_products_id
            component_ids(2) = kpw_air_id
            call s_kuhl_pw_energy_pr(rho_s, pres, Y_mix, component_ids, 2, jidx, e)
        case (jwl_mix_type_ptequil)
            call s_jwl_ptequil_energy_pr(rho_s, pres, Y_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                         & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_air_gammas(jidx), jwl_cv_prod, jwl_cv_air, e)
        case (jwl_mix_type_rocflu)
            call s_jwl_rocflu_energy_pr(rho_s, pres, Y_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                        & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), &
                                        & jwl_air_rho0s(jidx), jwl_air_gammas(jidx), e)
        case (jwl_mix_type_isobaric)
            call s_jwl_energy_pr(rho_s, pres, Y_s, alpha_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                 & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                 & jwl_air_gammas(jidx), e)
        case default
            call s_jwl_energy_pr(rho_s, pres, Y_s, alpha_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                 & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), jwl_air_rho0s(jidx), &
                                 & jwl_air_gammas(jidx), e)
        end select

    end subroutine s_jwl_mix_energy_pr

    !> Mixture sound speed for the Kuhl closure via the black-box (rho, e, Y) -> c map. Recovers internal energy from pressure, then
    !! evaluates the additive-pressure sound speed consistently with s_jwl_mix_pressure_er.
    subroutine s_jwl_mix_sound_speed(rho, pres, Y, alpha_j, jidx, c)

        $:GPU_ROUTINE(function_name='s_jwl_mix_sound_speed',parallelism='[seq]', cray_noinline=True)

        real(wp), intent(in)  :: rho, pres, Y, alpha_j
        integer, intent(in)   :: jidx
        real(wp), intent(out) :: c
        real(wp)              :: Y_s, e, p_m, T_m, Y_mix(2), c2
        integer               :: component_ids(2)

        Y_s = min(max(Y, 0._wp), 1._wp)

        select case (jwl_mix_type)
        case (jwl_mix_type_kuhl)
            Y_mix(1) = Y_s
            Y_mix(2) = 1._wp - Y_s
            component_ids(1) = kpw_products_id
            component_ids(2) = kpw_air_id
            call s_jwl_mix_energy_pr(rho, pres, Y, alpha_j, jidx, e)
            call s_kuhl_pw_closure(rho, e, Y_mix, component_ids, 2, jidx, p_m, T_m, c)
        case (jwl_mix_type_rocflu)
            call s_jwl_rocflu_sound_speed_squared(rho, pres, Y_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                                  & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), &
                                                  & jwl_air_rho0s(jidx), jwl_air_gammas(jidx), c2)
            c = sqrt(max(c2, sgm_eps))
        case default
            call s_jwl_mixture_sound_speed_squared(rho, pres, Y_s, Y_s, jwl_As(jidx), jwl_Bs(jidx), jwl_R1s(jidx), jwl_R2s(jidx), &
                                                   & jwl_omegas(jidx), jwl_rho0s(jidx), jwl_E0s(jidx), jwl_air_e0s(jidx), &
                                                   & jwl_air_rho0s(jidx), jwl_air_gammas(jidx), c2)
            c = sqrt(max(c2, sgm_eps))
        end select

    end subroutine s_jwl_mix_sound_speed

    !> Initialize JWL parameter tables.
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

        jwl_cv_prod = 0._wp
        jwl_cv_air = 0._wp
        if (jwl_idx > 0) jwl_cv_prod = fluid_pp(jwl_idx)%cv
        do i = 1, num_fluids
            if (fluid_pp(i)%eos == eos_stiffened_gas) then
                jwl_cv_air = fluid_pp(i)%cv
                exit
            end if
        end do

        if (jwl_idx > 0 .and. jwl_mix_type == jwl_mix_type_ptequil) then
            if (f_is_default(jwl_cv_prod) .or. jwl_cv_prod <= 0._wp .or. f_is_default(jwl_cv_air) .or. jwl_cv_air <= 0._wp) then
                call s_mpi_abort('jwl_mix_type_ptequil requires positive ' &
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
