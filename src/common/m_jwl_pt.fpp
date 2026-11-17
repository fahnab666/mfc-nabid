!>
!! @file
!! @brief Scalar JWL/stiffened-gas pressure-temperature equilibrium.

#:include 'macros.fpp'

!> @brief Conservative scalar PT equilibrium for one JWL and one ideal or stiffened constituent.
module m_jwl_pt

    use m_precision_select, only: wp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    private
    public :: jwl_pt_parameters, s_jwl_pt_equilibrium, s_jwl_pt_equilibrium_matrix, s_pt_cw_guess, s_jwl_wall_density

    !> Constituent data in the same energy and temperature convention as m_eos.
    type jwl_pt_parameters
        real(wp) :: a, b, r1, r2, omega, rho0, t0
        real(wp) :: cv_j, cv_a, gamma_a, pi_a, qv_j, qv_a
    end type jwl_pt_parameters

contains

    !> Caloric reference e_c and p_c = rho^2 de_c/drho. JWL's nonzero reference temperature is retained.
    subroutine s_pt_curve(rho, mat, ec, pc, dpc)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: rho
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: ec, pc, dpc
        real(wp)                            :: ea, eb, tref

        ea = mat%a*exp(-mat%r1*mat%rho0/rho)
        eb = mat%b*exp(-mat%r2*mat%rho0/rho)
        tref = 0._wp
        if (mat%t0 > 0._wp) tref = mat%t0*(rho/mat%rho0)**mat%omega
        ec = (ea/mat%r1 + eb/mat%r2)/mat%rho0 - mat%cv_j*tref
        pc = ea + eb - mat%omega*rho*mat%cv_j*tref
        dpc = ((ea*mat%r1 + eb*mat%r2)*mat%rho0/rho)/rho - mat%omega*mat%cv_j*(1._wp + mat%omega)*tref

    end subroutine s_pt_curve

    !> Composition-weighted (CW) law with coefficients blended by the products' heat-capacity share, mapped to a product volume
    !! fraction through the ambient EOS. The ambient cold energy pihat/rho keeps T exact at Y = 0 for a stiffened gas.
    subroutine s_pt_cw_guess(mj, ma, rhoe, mat, alpha, ok)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: mj, ma, rhoe
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: alpha
        logical, intent(out)                :: ok
        real(wp)                            :: rho, e, weight, ec, pc, dpc, omega, pihat, pres, temp, rhoa

        alpha = 0._wp; ok = .false.
        rho = mj + ma
        if (.not. (mj > 0._wp .and. ma > 0._wp .and. rho > 0._wp)) return
        if (.not. (mat%cv_j > 0._wp .and. mat%cv_a > 0._wp .and. mat%gamma_a > 0._wp)) return
        e = (rhoe - mj*mat%qv_j - ma*mat%qv_a)/rho
        weight = mj*mat%cv_j/(mj*mat%cv_j + ma*mat%cv_a)
        omega = mat%gamma_a + weight*(mat%omega - mat%gamma_a)
        pihat = (1._wp - weight)*mat%pi_a
        call s_pt_curve(rho, mat, ec, pc, dpc)
        pres = weight*(pc - omega*rho*ec) + omega*rho*e - (mat%gamma_a + 1._wp)*pihat
        temp = (rho*(e - weight*ec) - pihat)/(mj*mat%cv_j + ma*mat%cv_a)
        if (.not. (ieee_is_finite(pres) .and. ieee_is_finite(temp))) return
        if (.not. (temp > 0._wp .and. pres + mat%pi_a > 0._wp)) return
        rhoa = (pres + mat%pi_a)/(mat%gamma_a*mat%cv_a*temp)
        if (.not. (ieee_is_finite(rhoa) .and. rhoa > 0._wp)) return
        alpha = 1._wp - ma/rhoa
        ok = ieee_is_finite(alpha) .and. alpha > epsilon(1._wp) .and. alpha < 1._wp - epsilon(1._wp)

    end subroutine s_pt_cw_guess

    !> Damped isothermal wall-density inverse with a checked pressure residual.
    subroutine s_jwl_wall_density(rho_from, p_from, p_to, mat, rho_to, ok)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: rho_from, p_from, p_to
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: rho_to
        logical, intent(out)                :: ok
        real(wp)                            :: ec, pc, dpc, temp, thermal, residual, derivative, allowance, step, trial, fc, dc
        integer                             :: it, inner
        logical                             :: accepted

        rho_to = rho_from
        ok = .false.
        if (.not. (rho_from > 0._wp .and. p_from > 0._wp .and. p_to > 0._wp)) return
        if (.not. (mat%omega > 0._wp .and. mat%cv_j > 0._wp .and. mat%t0 >= 0._wp)) return
        if (.not. (ieee_is_finite(rho_from) .and. ieee_is_finite(p_from) .and. ieee_is_finite(p_to))) return
        call s_pt_curve(rho_from, mat, ec, pc, dpc)
        temp = (p_from - pc)/(mat%omega*rho_from*mat%cv_j)
        if (.not. (temp > 0._wp .and. ieee_is_finite(temp))) return

        do it = 1, 80
            call s_pt_curve(rho_to, mat, ec, pc, dpc)
            thermal = mat%omega*rho_to*mat%cv_j*temp
            residual = pc + thermal - p_to
            derivative = dpc + mat%omega*mat%cv_j*temp
            if (.not. (ieee_is_finite(residual) .and. ieee_is_finite(derivative))) return
            if (derivative <= 0._wp) return
            allowance = 16._wp*epsilon(1._wp)*(abs(pc) + abs(thermal) + abs(p_to))
            if (abs(residual) <= allowance) then
                ok = .true.
                return
            end if
            ! Stay positive without flooring density; backtrack until the pressure residual decreases.
            step = max(-0.5_wp*rho_to, min(0.5_wp*rho_to, residual/derivative))
            accepted = .false.
            do inner = 1, 24
                trial = rho_to - step
                call s_pt_curve(trial, mat, ec, pc, dpc)
                fc = pc + mat%omega*trial*mat%cv_j*temp - p_to
                dc = dpc + mat%omega*mat%cv_j*temp
                if (ieee_is_finite(fc) .and. dc > 0._wp .and. abs(fc) < abs(residual)) then
                    accepted = .true.
                    exit
                end if
                step = 0.5_wp*step
            end do
            if (.not. accepted) return
            rho_to = trial
        end do

    end subroutine s_jwl_wall_density

    !> Residual in log(alpha_j/alpha_a). Compute both volumes independently, including near a pure phase.
    subroutine s_pt_trial(x, mj, ma, rhoe, mat, aj, aa, temp, pj, pa, residual, derivative)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: x, mj, ma, rhoe
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: aj, aa, temp, pj, pa, residual, derivative
        real(wp)                            :: z, rj, ra, ec, pc, dpc, cvsum, dtemp

        z = exp(-abs(x))
        if (x >= 0._wp) then
            aj = 1._wp/(1._wp + z); aa = z/(1._wp + z)
        else
            aj = z/(1._wp + z); aa = 1._wp/(1._wp + z)
        end if
        rj = mj/aj; ra = ma/aa
        call s_pt_curve(rj, mat, ec, pc, dpc)
        cvsum = mj*mat%cv_j + ma*mat%cv_a
        temp = (rhoe - mj*(ec + mat%qv_j) - ma*mat%qv_a - mat%pi_a*aa)/cvsum
        dtemp = (pc + mat%pi_a)*aj*aa/cvsum
        pj = pc + mat%omega*rj*mat%cv_j*temp
        pa = mat%gamma_a*ra*mat%cv_a*temp - mat%pi_a
        residual = pj - pa
        derivative = -(dpc + mat%omega*mat%cv_j*temp)*rj*aa + mat%omega*rj*mat%cv_j*dtemp - mat%gamma_a*mat%cv_a*ra*(aj*temp &
                       & + dtemp)

    end subroutine s_pt_trial

    !> Equilibrium sound speed from dv/dp along de = -p dv at fixed constituent masses.
    subroutine s_pt_sound(mj, ma, aj, aa, pres, temp, mat, c2, ok)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: mj, ma, aj, aa, pres, temp
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: c2
        logical, intent(out)                :: ok
        real(wp)                            :: rj, ra, rho, yj, ya, ec, pc, dpc, jp, jt, ap, at, vp, vt, ep, et, dvs

        rho = mj + ma; rj = mj/aj; ra = ma/aa
        yj = mj/rho; ya = ma/rho
        call s_pt_curve(rj, mat, ec, pc, dpc)
        dpc = dpc + mat%omega*mat%cv_j*temp
        ok = dpc > 0._wp
        c2 = 0._wp
        if (.not. ok) return
        jp = 1._wp/dpc; jt = -mat%omega*rj*mat%cv_j/dpc
        ap = 1._wp/(mat%gamma_a*mat%cv_a*temp); at = -ra/temp
        vp = -yj*jp/rj**2 - ya*ap/ra**2
        vt = -yj*jt/rj**2 - ya*at/ra**2
        ep = yj*pc*jp/rj**2 - ya*mat%pi_a*ap/ra**2
        et = yj*mat%cv_j + ya*mat%cv_a + yj*pc*jt/rj**2 - ya*mat%pi_a*at/ra**2
        ok = et + pres*vt > 0._wp
        if (.not. ok) return
        dvs = vp - vt*(ep + pres*vp)/(et + pres*vt)
        ok = dvs < 0._wp
        if (.not. ok) return
        c2 = -1._wp/(rho**2*dvs)
        ok = ieee_is_finite(c2) .and. c2 > 0._wp

    end subroutine s_pt_sound

    !> Accept only a finite, positive, hyperbolic equilibrium with a measured pressure residual.
    subroutine s_pt_accept(x, mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: x, mj, ma, rhoe
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: pres, temp, aj, aa, c2
        logical, intent(out)                :: ok
        real(wp)                            :: pj, pa, residual, derivative, tol, roundoff, scale, allowance

        call s_pt_trial(x, mj, ma, rhoe, mat, aj, aa, temp, pj, pa, residual, derivative)
        tol = max(64._wp*epsilon(1._wp), 1.e-10_wp)
        pres = pj
        ! Subtracting the ambient stiffness has an unavoidable precision-dependent error bound.
        roundoff = 8._wp*epsilon(1._wp)*(abs(pa + mat%pi_a) + abs(mat%pi_a))
        scale = max(abs(pj), abs(pa), tiny(1._wp))
        allowance = tol*scale
        ! Do not turn an unresolved pressure into success through an oversized roundoff allowance.
        if (roundoff <= sqrt(epsilon(1._wp))*scale) allowance = max(allowance, roundoff)
        c2 = 0._wp
        ok = ieee_is_finite(pres) .and. ieee_is_finite(temp) .and. temp > 0._wp .and. pres > 0._wp .and. abs(residual) <= allowance
        if (ok) call s_pt_sound(mj, ma, aj, aa, pres, temp, mat, c2, ok)

    end subroutine s_pt_accept

    !> CW residual gate: if |p_j - p_a| <= 1e-4 max(|p_j|, |p_a|) at the CW split, return the split refined by one Newton step in
    !! rho_j, with p the mean of the constituent pressures, if the refined residual is within tau**2. Since x = ln(m_j) - ln(rho_j -
    !! Y rho), that step is x - ln(1 + R/R_x).
    subroutine s_pt_gate(x, mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: x, mj, ma, rhoe
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: pres, temp, aj, aa, c2
        logical, intent(out)                :: ok
        real(wp), parameter                 :: tau = 1.e-4_wp
        real(wp)                            :: pj, pa, residual, derivative, ratio

        pres = 0._wp; c2 = 0._wp; ok = .false.
        call s_pt_trial(x, mj, ma, rhoe, mat, aj, aa, temp, pj, pa, residual, derivative)
        if (.not. (abs(residual) <= tau*max(abs(pj), abs(pa)) .and. derivative /= 0._wp)) return
        ratio = 1._wp + residual/derivative
        if (.not. (ieee_is_finite(ratio) .and. ratio > 0._wp)) return
        call s_pt_trial(x - log(ratio), mj, ma, rhoe, mat, aj, aa, temp, pj, pa, residual, derivative)
        pres = 0.5_wp*(pj + pa)
        ! A Newton step from a residual of tau leaves O(tau**2); a larger one means the step failed, so the solver takes over.
        ok = ieee_is_finite(pres) .and. ieee_is_finite(temp) .and. temp > 0._wp .and. pres > 0._wp .and. abs(residual) &
                            & <= max(tau*tau, 64._wp*epsilon(1._wp))*max(abs(pj), abs(pa))
        if (ok) call s_pt_sound(mj, ma, aj, aa, pres, temp, mat, c2, ok)

    end subroutine s_pt_gate

    !> Input check and the closed forms for a cell holding one phase; done is false when both phases are present.
    subroutine s_pt_single_phase(mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok, done)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: mj, ma, rhoe
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: pres, temp, aj, aa, c2
        logical, intent(out)                :: ok, done
        real(wp)                            :: ec, pc, dpc, tref, pref, dpref

        pres = 0._wp; temp = 0._wp; aj = 0._wp; aa = 0._wp; c2 = 0._wp; ok = .false.; done = .true.
        if (.not. (mj >= 0._wp .and. ma >= 0._wp .and. mj + ma > 0._wp .and. ieee_is_finite(rhoe))) return
        if (.not. (mat%cv_j > 0._wp .and. mat%cv_a > 0._wp .and. mat%gamma_a > 0._wp)) return
        if (mj == 0._wp) then
            aa = 1._wp
            temp = (rhoe - ma*mat%qv_a - mat%pi_a)/(ma*mat%cv_a)
            pres = mat%gamma_a*ma*mat%cv_a*temp - mat%pi_a
            c2 = (1._wp + mat%gamma_a)*(pres + mat%pi_a)/ma
        else if (ma == 0._wp) then
            aj = 1._wp
            call s_pt_curve(mj, mat, ec, pc, dpc)
            temp = (rhoe/mj - ec - mat%qv_j)/mat%cv_j
            pres = pc + mat%omega*mj*mat%cv_j*temp
            tref = mat%t0*(mj/mat%rho0)**mat%omega
            pref = pc + mat%omega*mj*mat%cv_j*tref
            dpref = dpc + mat%omega*mat%cv_j*(1._wp + mat%omega)*tref
            c2 = dpref + (1._wp + mat%omega)*(pres - pref)/mj
        else
            done = .false.
            return
        end if
        ok = temp > 0._wp .and. pres > 0._wp .and. c2 > 0._wp .and. ieee_is_finite(c2)

    end subroutine s_pt_single_phase

    !> Safeguarded Newton in x from x0: steps limited to 4 and halved until |p_j - p_a| decreases.
    subroutine s_pt_newton(x0, limit, mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: x0, limit, mj, ma, rhoe
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: pres, temp, aj, aa, c2
        logical, intent(out)                :: ok
        real(wp)                            :: x, candidate, dx, pj, pa, residual, derivative, fc, dc
        integer                             :: it, inner

        x = x0
        do it = 1, 24
            call s_pt_accept(x, mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok)
            if (ok) return
            call s_pt_trial(x, mj, ma, rhoe, mat, aj, aa, temp, pj, pa, residual, derivative)
            if (.not. ieee_is_finite(derivative) .or. derivative == 0._wp) return
            dx = max(-4._wp, min(4._wp, residual/derivative))
            do inner = 1, 12
                candidate = max(-limit, min(limit, x - dx))
                call s_pt_trial(candidate, mj, ma, rhoe, mat, aj, aa, temp, pj, pa, fc, dc)
                if (ieee_is_finite(fc) .and. abs(fc) < abs(residual)) exit
                dx = 0.5_wp*dx
            end do
            if (inner > 12) return
            x = candidate
        end do

    end subroutine s_pt_newton

    !> Root in the bracket [left, right] (fl = residual at left) by Newton steps kept inside it, else bisection.
    subroutine s_pt_bracket(left, right, fl, mj, ma, rhoe, mat, x, pres, temp, aj, aa, c2, ok)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(inout)             :: left, right, fl
        real(wp), intent(in)                :: mj, ma, rhoe
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: x, pres, temp, aj, aa, c2
        logical, intent(out)                :: ok
        real(wp)                            :: candidate, dx, pj, pa, fc, dc
        integer                             :: it

        x = 0.5_wp*(left + right)
        do it = 1, 100
            call s_pt_accept(x, mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok)
            if (ok) return
            call s_pt_trial(x, mj, ma, rhoe, mat, aj, aa, temp, pj, pa, fc, dc)
            if (.not. ieee_is_finite(fc)) return
            if ((fl <= 0._wp .and. fc >= 0._wp) .or. (fl >= 0._wp .and. fc <= 0._wp)) then
                right = x
            else
                left = x; fl = fc
            end if
            candidate = 0.5_wp*(left + right)
            if (dc /= 0._wp) then
                dx = x - fc/dc
                if (dx > left .and. dx < right) candidate = dx
            end if
            if (candidate == x) return
            x = candidate
        end do

    end subroutine s_pt_bracket

    !> Scan [-limit, limit] in 128 cells, refine every sign change, and return the accepted root nearest x0.
    subroutine s_pt_scan(x0, limit, mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: x0, limit, mj, ma, rhoe
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: pres, temp, aj, aa, c2
        logical, intent(out)                :: ok
        real(wp)                            :: x, left, right, fl, fprev, previous, trial, pj, pa, residual, derivative, distance
        real(wp)                            :: pbest, tbest, abest, bbest, cbest
        integer                             :: scan
        logical                             :: accepted

        previous = -limit
        call s_pt_trial(previous, mj, ma, rhoe, mat, aj, aa, temp, pj, pa, fprev, derivative)
        distance = huge(1._wp)
        ok = .false.
        do scan = 1, 128
            trial = -limit + 2._wp*limit*real(scan, wp)/128._wp
            call s_pt_trial(trial, mj, ma, rhoe, mat, aj, aa, temp, pj, pa, residual, derivative)
            if (ieee_is_finite(fprev) .and. ieee_is_finite(residual)) then
                if ((fprev <= 0._wp .and. residual >= 0._wp) .or. (fprev >= 0._wp .and. residual <= 0._wp)) then
                    left = previous; right = trial; fl = fprev
                    call s_pt_bracket(left, right, fl, mj, ma, rhoe, mat, x, pres, temp, aj, aa, c2, accepted)
                    if (accepted .and. abs(x - x0) < distance) then
                        distance = abs(x - x0)
                        pbest = pres; tbest = temp; abest = aj; bbest = aa; cbest = c2
                        ok = .true.
                    end if
                end if
            end if
            previous = trial; fprev = residual
        end do
        if (ok) then
            pres = pbest; temp = tbest; aj = abest; aa = bbest; c2 = cbest
        end if

    end subroutine s_pt_scan

    !> Fast PT: keep the existing split if it is in equilibrium; else the CW predictor behind the residual gate; else a safeguarded
    !! Newton in x = ln(alpha_j/alpha_a); else a bracketed scan. Mass and energy are conserved and failure is explicit.
    subroutine s_jwl_pt_equilibrium(mj, ma, rhoe, alpha_guess, mat, pres, temp, aj, aa, c2, ok)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: mj, ma, rhoe, alpha_guess
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: pres, temp, aj, aa, c2
        logical, intent(out)                :: ok
        real(wp)                            :: limit, x0, alpha_cw
        logical                             :: done, cw_ok

        call s_pt_single_phase(mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok, done)
        if (done) return
        x0 = log(mj) - log(ma)
        if (alpha_guess > epsilon(1._wp) .and. alpha_guess < 1._wp - epsilon(1._wp)) then
            x0 = log(alpha_guess) - log(1._wp - alpha_guess)
            call s_pt_accept(x0, mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok)
            if (ok) return
        end if
        call s_pt_cw_guess(mj, ma, rhoe, mat, alpha_cw, cw_ok)
        if (cw_ok) then
            x0 = log(alpha_cw) - log(1._wp - alpha_cw)
            call s_pt_gate(x0, mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok)
            if (ok) return
        end if
        ! The volume coordinate has no subtraction of nearly equal densities and no dimensional floor.
        limit = -2.5_wp*log(epsilon(1._wp))
        x0 = max(-limit, min(limit, x0))
        call s_pt_newton(x0, limit, mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok)
        if (.not. ok) call s_pt_scan(x0, limit, mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok)

    end subroutine s_jwl_pt_equilibrium

    !> Full PT equilibrium in (rho_j, rho_a, T, p) by a scaled 4x4 Newton system: from the given volume fraction, then, if that
    !! fails or none is usable, from the CW predictor (else alpha_j = Y_j).
    subroutine s_jwl_pt_equilibrium_matrix(mj, ma, rhoe, alpha_guess, mat, pres, temp, aj, aa, c2, ok)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: mj, ma, rhoe, alpha_guess
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: pres, temp, aj, aa, c2
        logical, intent(out)                :: ok
        real(wp)                            :: alpha_cw
        logical                             :: done, cw_ok

        call s_pt_single_phase(mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok, done)
        if (done) return
        if (alpha_guess > epsilon(1._wp) .and. alpha_guess < 1._wp - epsilon(1._wp)) then
            call s_pt_matrix_newton(mj, ma, rhoe, alpha_guess, mat, pres, temp, aj, aa, c2, ok)
            if (ok) return
        end if
        call s_pt_cw_guess(mj, ma, rhoe, mat, alpha_cw, cw_ok)
        if (.not. cw_ok) alpha_cw = mj/(mj + ma)
        call s_pt_matrix_newton(mj, ma, rhoe, alpha_cw, mat, pres, temp, aj, aa, c2, ok)

    end subroutine s_jwl_pt_equilibrium_matrix

    !> Damped 4x4 Newton from the split alpha_guess. Every iterate is tested with the Fast PT acceptance, so the pressure tolerance
    !! is relative to p and the returned state conserves energy exactly.
    subroutine s_pt_matrix_newton(mj, ma, rhoe, alpha_guess, mat, pres, temp, aj, aa, c2, ok)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: mj, ma, rhoe, alpha_guess
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: pres, temp, aj, aa, c2
        logical, intent(out)                :: ok
        real(wp)                            :: rho, x(4), f(4), jac(4, 4), step(4), trial(4), ft(4), jt(4, 4)
        real(wp)                            :: tscale, scale, norm0, norm1, damp, aguess, pj, pa, residual, derivative
        integer                             :: it, ls
        logical                             :: solved

        rho = mj + ma
        ! Scales use the energy without qv, so a common shift of the formation energies leaves the iteration unchanged.
        tscale = max(abs(rhoe - mj*mat%qv_j - ma*mat%qv_a)/(mj*mat%cv_j + ma*mat%cv_a), 1._wp)
        scale = max(rho*tscale*max(mat%cv_j, mat%cv_a), 1._wp)
        ! Unknowns rho_k/rho, T/tscale, p/scale; start at the guessed split with T from energy and p = p_j.
        aguess = min(1._wp - epsilon(1._wp), max(epsilon(1._wp), alpha_guess))
        call s_pt_trial(log(aguess) - log(1._wp - aguess), mj, ma, rhoe, mat, aj, aa, temp, pj, pa, residual, derivative)
        x(1) = mj/(rho*aj); x(2) = ma/(rho*aa); x(3) = max(temp/tscale, 1.e-8_wp); x(4) = max(pj/scale, 1.e-8_wp)
        do it = 1, 2000
            call s_pt_accept(log(mj*x(2)) - log(ma*x(1)), mj, ma, rhoe, mat, pres, temp, aj, aa, c2, ok)
            if (ok) return
            call s_pt_matrix_system(x, mj, ma, rhoe, mat, rho, tscale, scale, f, jac)
            norm0 = maxval(abs(f))
            if (.not. ieee_is_finite(norm0)) return
            call s_pt_solve4(jac, -f, step, solved)
            if (.not. solved) return
            damp = 1._wp
            do ls = 1, 24
                trial = x + damp*step
                if (minval(trial(1:3)) > 0._wp) then
                    call s_pt_matrix_system(trial, mj, ma, rhoe, mat, rho, tscale, scale, ft, jt)
                    norm1 = maxval(abs(ft))
                    if (ieee_is_finite(norm1) .and. norm1 < norm0) exit
                end if
                damp = 0.5_wp*damp
            end do
            if (ls > 24) return
            x = trial
        end do

    end subroutine s_pt_matrix_newton

    !> Residual and analytic Jacobian of the four PT equations; pressure and energy rows share the scale ps.
    subroutine s_pt_matrix_system(x, mj, ma, rhoe, mat, rs, ts, ps, f, j)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)                :: x(4), mj, ma, rhoe, rs, ts, ps
        type(jwl_pt_parameters), intent(in) :: mat
        real(wp), intent(out)               :: f(4), j(4, 4)
        real(wp)                            :: rj, ra, t, p, ec, pc, dpc, pj, pa

        rj = x(1)*rs; ra = x(2)*rs; t = x(3)*ts; p = x(4)*ps
        call s_pt_curve(rj, mat, ec, pc, dpc)
        pj = pc + mat%omega*rj*mat%cv_j*t
        pa = mat%gamma_a*ra*mat%cv_a*t - mat%pi_a
        f(1) = (pj - p)/ps; f(2) = (pa - p)/ps
        f(3) = mj/rj + ma/ra - 1._wp
        f(4) = (mj*(ec + mat%qv_j + mat%cv_j*t) + ma*(mat%qv_a + mat%cv_a*t + mat%pi_a/ra) - rhoe)/ps
        j = 0._wp
        j(1, 1) = rs*(dpc + mat%omega*mat%cv_j*t)/ps; j(1, 3) = mat%omega*rj*mat%cv_j*ts/ps; j(1, 4) = -1._wp
        j(2, 2) = rs*mat%gamma_a*mat%cv_a*t/ps; j(2, 3) = mat%gamma_a*ra*mat%cv_a*ts/ps; j(2, 4) = -1._wp
        j(3, 1) = -mj*rs/rj**2; j(3, 2) = -ma*rs/ra**2
        j(4, 1) = mj*pc*rs/rj**2/ps; j(4, 2) = -ma*mat%pi_a*rs/ra**2/ps
        j(4, 3) = (mj*mat%cv_j + ma*mat%cv_a)*ts/ps

    end subroutine s_pt_matrix_system

    !> Gaussian elimination with partial pivoting for the 4x4 Newton update.
    subroutine s_pt_solve4(a, b, x, ok)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in)  :: a(4, 4), b(4)
        real(wp), intent(out) :: x(4)
        logical, intent(out)  :: ok
        real(wp)              :: m(4, 5), scale, tmp
        integer               :: i, j, k, piv

        m(:,1:4) = a; m(:,5) = b; ok = .false.
        do k = 1, 3
            piv = k
            do i = k + 1, 4
                if (abs(m(i, k)) > abs(m(piv, k))) piv = i
            end do
            if (abs(m(piv, k)) <= tiny(1._wp)) return
            if (piv /= k) then
                do j = k, 5
                    tmp = m(k, j); m(k, j) = m(piv, j); m(piv, j) = tmp
                end do
            end if
            do i = k + 1, 4
                scale = m(i, k)/m(k, k); m(i,k:5) = m(i,k:5) - scale*m(k,k:5)
            end do
        end do
        if (abs(m(4, 4)) <= tiny(1._wp)) return
        x(4) = m(4, 5)/m(4, 4)
        do i = 3, 1, -1
            tmp = m(i, 5)
            do j = i + 1, 4
                tmp = tmp - m(i, j)*x(j)
            end do
            if (abs(m(i, i)) <= tiny(1._wp)) return
            x(i) = tmp/m(i, i)
        end do
        ok = all(ieee_is_finite(x))

    end subroutine s_pt_solve4

end module m_jwl_pt
