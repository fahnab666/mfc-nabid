!>
!! @file
!! @brief Contains module m_helper_basic

#:include 'macros.fpp'

!> @brief Basic floating-point utilities: approximate equality, default detection, and coordinate bounds
module m_helper_basic

    use m_derived_types
    use m_precision_select
    use m_constants, only: recon_type_weno, recon_type_muscl

    implicit none

    private
    public :: f_approx_equal, f_approx_in_array, f_is_default, f_all_default, f_is_integer, s_configure_coordinate_bounds, &
        & s_update_cell_bounds, f_mg_reference, f_mg_dref_drho, f_mg_pressure, f_mg_internal_energy, f_mg_temperature, &
        & f_mg_sound_speed_sq

contains

    !> Check if two floating point numbers of wp are within tolerance.
    !! @param tol_input Relative error (default = 1.e-10_wp).
    logical elemental function f_approx_equal(a, b, tol_input) result(res)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in)           :: a, b
        real(wp), optional, intent(in) :: tol_input
        real(wp)                       :: tol

        if (present(tol_input)) then
            tol = tol_input
        else
            if (wp == single_precision) then
                tol = 1.e-6_wp
            else
                tol = 1.e-10_wp
            end if
        end if

        if (a == b) then
            res = .true.
        else if (a == 0._wp .or. b == 0._wp .or. (abs(a) + abs(b) < tiny(a))) then
            res = (abs(a - b) < (tol*tiny(a)))
        else
            res = (abs(a - b)/min(abs(a) + abs(b), huge(a)) < tol)
        end if

    end function f_approx_equal

    !> Check if a wp value approximately matches any element of an array within tolerance.
    !! @param tol_input Relative error (default = 1e-10_wp).
    logical function f_approx_in_array(a, b, tol_input) result(res)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in)           :: a
        real(wp), intent(in)           :: b(:)
        real(wp), optional, intent(in) :: tol_input
        real(wp)                       :: tol
        integer                        :: i

        res = .false.

        if (present(tol_input)) then
            tol = tol_input
        else
            if (wp == single_precision) then
                tol = 1.e-6_wp
            else
                tol = 1.e-10_wp
            end if
        end if

        do i = 1, size(b)
            if (f_approx_equal(a, b(i), tol)) then
                res = .true.
                exit
            end if
        end do

    end function f_approx_in_array

    !> Check if a real(wp) variable is of default value.
    logical elemental function f_is_default(var) result(res)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: var

        res = f_approx_equal(var, dflt_real)

    end function f_is_default

    !> Check if ALL elements of a real(wp) array are of default value.
    logical function f_all_default(var_array) result(res)

        real(wp), intent(in) :: var_array(:)

        res = all(f_is_default(var_array))

    end function f_all_default

    !> Check if a real(wp) variable is an integer.
    logical elemental function f_is_integer(var) result(res)

        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: var

        res = f_approx_equal(var, real(nint(var), wp))

    end function f_is_integer

    !> Compute ghost-cell buffer size and set interior/buffered coordinate index bounds.
    subroutine s_configure_coordinate_bounds(recon_type, weno_polyn, muscl_polyn, igr_order, buff_size, idwint, idwbuff, viscous, &
        & bubbles_lagrange, m, n, p, num_dims, igr, ib, fd_number)

        integer, intent(in)                                :: recon_type, weno_polyn, muscl_polyn
        integer, intent(in)                                :: m, n, p, num_dims, igr_order, fd_number
        integer, intent(inout)                             :: buff_size
        type(int_bounds_info), dimension(3), intent(inout) :: idwint, idwbuff
        logical, intent(in)                                :: viscous, bubbles_lagrange
        logical, intent(in)                                :: igr
        logical, intent(in)                                :: ib

        ! Determine ghost cell buffer size for boundary conditions

        if (igr) then
            buff_size = (igr_order - 1)/2 + 2
        else if (recon_type == recon_type_weno) then
            if (viscous) then
                buff_size = 2*weno_polyn + 2
            else
                buff_size = weno_polyn + 2
            end if
        else if (recon_type == recon_type_muscl) then
            buff_size = muscl_polyn + 2
        end if

        ! Correction for smearing function in the lagrangian subgrid bubble model
        if (bubbles_lagrange) then
            buff_size = max(buff_size + fd_number, mapCells + 1 + fd_number)
        end if

        if (ib) then
            buff_size = max(buff_size, 10)
        end if

        ! Configuring Coordinate Direction Indexes
        idwint(1)%beg = 0; idwint(2)%beg = 0; idwint(3)%beg = 0
        idwint(1)%end = m; idwint(2)%end = n; idwint(3)%end = p

        idwbuff(1)%beg = -buff_size
        if (num_dims > 1) then; idwbuff(2)%beg = -buff_size; else; idwbuff(2)%beg = 0; end if
        if (num_dims > 2) then; idwbuff(3)%beg = -buff_size; else; idwbuff(3)%beg = 0; end if

        idwbuff(1)%end = idwint(1)%end - idwbuff(1)%beg
        idwbuff(2)%end = idwint(2)%end - idwbuff(2)%beg
        idwbuff(3)%end = idwint(3)%end - idwbuff(3)%beg

    end subroutine s_configure_coordinate_bounds

    !> Update the min and max number of cells in each set of axes
    !! @param bounds min and max values to update
    elemental subroutine s_update_cell_bounds(bounds, m, n, p)

        type(cell_num_bounds), intent(out) :: bounds
        integer, intent(in)                :: m, n, p

        bounds%mn_max = max(m, n)
        bounds%np_max = max(n, p)
        bounds%mp_max = max(m, p)
        bounds%mnp_max = max(m, n, p)
        bounds%mn_min = min(m, n)
        bounds%np_min = min(n, p)
        bounds%mp_min = min(m, p)
        bounds%mnp_min = min(m, n, p)

    end subroutine s_update_cell_bounds

    !> Mie-Grueneisen reference pressure curve f(rho) at constant Grueneisen coefficient. Compressed states follow the linear-Us
    !! Hugoniot; expanded states an isentropic tail.
    real(wp) elemental function f_mg_reference(rho, rho0, c0, s, p0, Gamma) result(f_ref)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, rho0, c0, s, p0, Gamma
        real(wp)             :: eta, p_h

        eta = (rho - rho0)/rho
        if (eta >= 0._wp) then
            p_h = p0 + rho0*c0**2*eta/(1._wp - s*eta)**2
            f_ref = p_h*(1._wp - Gamma*eta/(2._wp*(1._wp - eta))) - p0*Gamma*eta/(2._wp*(1._wp - eta))
        else
            f_ref = p0 + c0**2*(rho - rho0)
        end if

    end function f_mg_reference

    !> Density derivative of the Mie-Grueneisen reference curve, df/drho.
    real(wp) elemental function f_mg_dref_drho(rho, rho0, c0, s, p0, Gamma) result(dref)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, rho0, c0, s, p0, Gamma
        real(wp)             :: eta, p_h

        eta = (rho - rho0)/rho
        if (eta >= 0._wp) then
            p_h = p0 + rho0*c0**2*eta/(1._wp - s*eta)**2
            dref = c0**2*(1._wp - eta)*(1._wp - (1._wp + Gamma/2._wp)*eta)*(1._wp + s*eta)/(1._wp - s*eta)**3 - Gamma/(2._wp*rho0) &
                          & *(p_h + p0)
        else
            dref = c0**2
        end if

    end function f_mg_dref_drho

    !> Mie-Grueneisen pressure from density and internal energy per unit mass.
    real(wp) elemental function f_mg_pressure(rho, e, e0, Gamma, f_ref) result(pres)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, e, e0, Gamma, f_ref

        pres = Gamma*rho*(e - e0) + f_ref

    end function f_mg_pressure

    !> Inverse of f_mg_pressure: internal energy per unit mass from pressure.
    real(wp) elemental function f_mg_internal_energy(rho, pres, e0, Gamma, f_ref) result(e)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, pres, e0, Gamma, f_ref

        e = e0 + (pres - f_ref)/(Gamma*rho)

    end function f_mg_internal_energy

    !> Mie-Grueneisen temperature at constant specific heat.
    real(wp) elemental function f_mg_temperature(e, e0, cv, T0) result(T)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: e, e0, cv, T0

        T = T0 + (e - e0)/cv

    end function f_mg_temperature

    !> Single-material Mie-Grueneisen sound speed squared.
    real(wp) elemental function f_mg_sound_speed_sq(rho, e, e0, pres, Gamma, dref) result(c2)
        $:GPU_ROUTINE(parallelism='[seq]')
        real(wp), intent(in) :: rho, e, e0, pres, Gamma, dref

        c2 = Gamma*(e - e0) + dref + Gamma*pres/rho

    end function f_mg_sound_speed_sq

end module m_helper_basic
