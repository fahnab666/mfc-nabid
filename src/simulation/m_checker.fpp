!>
!!@file
!!@brief Contains module m_checker

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Validates simulation input parameters for consistency and supported configurations
module m_checker

    use m_global_parameters
    use m_mpi_proxy
    use m_helper
    use m_helper_basic
    use m_constants, only: recon_type_weno, recon_type_muscl, muscl_order_first_order, eos_jwl, wave_speeds_pressure, &
        & BC_CHAR_SLIP_WALL, BC_CHAR_SUP_OUTFLOW, model_eqns_5eq, bubble_model_gilmore

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file. Used by the simulation stage
    impure subroutine s_check_inputs

        call s_check_inputs_compilers

        if (igr) then
            call s_check_inputs_nvidia_uvm
        else
            if (recon_type == recon_type_weno) then
                call s_check_inputs_weno
            else if (recon_type == recon_type_muscl) then
                call s_check_inputs_muscl
            end if
        end if

        ! Paths that evaluate stiffened-gas gamma/pi_inf mixture relations, which are
        ! meaningless for a JWL fluid: the pressure-based wave-speed shock-Mach
        ! correction, characteristic (CBC) boundary treatments, the alternative
        ! (Wood-like) sound speed, and the elastic pressure recovery.
        if (any(fluid_pp(1:num_fluids)%eos == eos_jwl)) then
            @:PROHIBIT(wave_speeds == wave_speeds_pressure, &
                       & "wave_speeds = 2 (pressure-based) is not supported with eos_jwl; use wave_speeds = 1")
            @:PROHIBIT(any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) <= BC_CHAR_SLIP_WALL .and. (/bc_x%beg, &
                       & bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) >= BC_CHAR_SUP_OUTFLOW), &
                       & "characteristic (CBC) boundary conditions are not supported with eos_jwl")
            @:PROHIBIT(alt_soundspeed, "alt_soundspeed is not supported with eos_jwl")
            @:PROHIBIT(hypoelasticity, "elasticity is not supported with eos_jwl")
            ! These solver paths compute pressure from stiffened-gas relations and
            ! bypass the JWL closure entirely, so JWL cells would be silently wrong.
            @:PROHIBIT(igr, "igr is not supported with eos_jwl (IGR evaluates stiffened-gas pressure on JWL cells)")
            @:PROHIBIT(bubbles_euler .or. mhd .or. chemistry, &
                       & "bubbles_euler, mhd, and chemistry are not supported with eos_jwl (their pressure paths bypass the JWL closure)")

            ! ib is compatible: ghost/fresh cells (including %rxn/%abn) are rebuilt through
            ! the JWL closure in m_ibm, and reaction sources are gated off solid cells.
        else
            @:PROHIBIT(jwl_afterburn .or. prog_burn .or. jwl_reactive, &
                       & "jwl_afterburn, prog_burn, and jwl_reactive require a fluid with eos_jwl")
        end if

        if (jwl_reactive) then
            @:PROHIBIT(riemann_solver /= 2, "jwl_reactive requires riemann_solver = 2 (HLLC)")
            @:PROHIBIT(prog_burn, "jwl_reactive and prog_burn are mutually exclusive detonation models")
            @:PROHIBIT(f_is_default(jwl_G) .or. jwl_G <= 0._wp, "jwl_reactive requires positive jwl_G")
            @:PROHIBIT(f_is_default(jwl_b_exp) .or. jwl_b_exp <= 0._wp, "jwl_reactive requires positive jwl_b_exp")
        end if

        if (jwl_afterburn) then
            @:PROHIBIT(riemann_solver /= 2, "jwl_afterburn requires riemann_solver = 2 (HLLC)")
            @:PROHIBIT(jwl_ab_model /= 1 .and. jwl_ab_model /= 2, "jwl_ab_model must be 1 (mixing-rate) or 2 (Arrhenius)")
            @:PROHIBIT(f_is_default(jwl_q_ab) .or. jwl_q_ab <= 0._wp, "jwl_afterburn requires positive jwl_q_ab")
            @:PROHIBIT(jwl_ab_model == 1 .and. (f_is_default(jwl_ab_tau) .or. jwl_ab_tau <= 0._wp), &
                       & "jwl_ab_model = 1 requires positive jwl_ab_tau")
            @:PROHIBIT(jwl_ab_model == 2 .and. (f_is_default(jwl_ab_A) .or. jwl_ab_A <= 0._wp .or. f_is_default(jwl_ab_theta) &
                       & .or. jwl_ab_theta <= 0._wp), "jwl_ab_model = 2 requires positive jwl_ab_A and jwl_ab_theta")
            ! Afterburn is oxygen combustion with air; a stiffened (liquid) ambient has none.
            @:PROHIBIT(any(fluid_pp(1:num_fluids)%eos /= eos_jwl .and. fluid_pp(1:num_fluids)%pi_inf > 0._wp), &
                       & "jwl_afterburn requires an ideal-gas ambient fluid")
        end if

        if (prog_burn) then
            @:PROHIBIT(f_is_default(pb_D_cj) .or. pb_D_cj <= 0._wp, "prog_burn requires positive pb_D_cj")
            @:PROHIBIT(f_is_default(pb_width) .or. pb_width <= 0._wp, "prog_burn requires positive pb_width")
            ! In 3D cylindrical coordinates z is the azimuthal angle, so the Cartesian
            ! front distance sqrt(dx^2 + dy^2 + dz^2) would mix radians into meters.
            @:PROHIBIT(cyl_coord .and. p > 0, "prog_burn is not supported in 3D cylindrical coordinates")
            ! The burn front advances pb_D_cj*dt per step; if that exceeds the band
            ! width, annuli of cells are never swept and detonation energy is lost.
            @:PROHIBIT(.not. cfl_dt .and. .not. f_is_default(pb_D_cj) .and. .not. f_is_default(pb_width) &
                       & .and. pb_D_cj*dt > pb_width, &
                       & "prog_burn requires pb_D_cj*dt <= pb_width so the front advances at most one band width per step")
        end if

        if (ib .and. chemistry) then
            call s_check_inputs_ib_injection
        end if

    end subroutine s_check_inputs

    !> Checks constraints on compiler options
    impure subroutine s_check_inputs_compilers

#if !defined(MFC_OpenACC) && !(defined(__PGI) || defined(_CRAYFTN))
        @:PROHIBIT(rdma_mpi, "Unsupported value of rdma_mpi for the current compiler")
#endif

    end subroutine s_check_inputs_compilers

    !> Checks constraints on WENO scheme parameters
    impure subroutine s_check_inputs_weno

        character(len=5) :: numStr  !< for int to string conversion

        call s_int_to_str(num_stcls_min*weno_order, numStr)
        ! lint: runtime-check m/n/p are per-rank extents after MPI decomposition, not the case-file values
        @:PROHIBIT(m + 1 < num_stcls_min*weno_order, &
                   & "m must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " // trim(numStr))
        ! lint: runtime-check per-rank n
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*weno_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))
        ! lint: runtime-check per-rank p
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*weno_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))

    end subroutine s_check_inputs_weno

    !> Validate that the grid resolution is sufficient for the MUSCL reconstruction order
    impure subroutine s_check_inputs_muscl

        character(len=5) :: numStr  !< for int to string conversion

        call s_int_to_str(num_stcls_min*muscl_order, numStr)
        ! lint: runtime-check m/n/p are per-rank extents after MPI decomposition, not the case-file values
        @:PROHIBIT(m + 1 < num_stcls_min*muscl_order, &
                   & "m must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " // trim(numStr))
        ! lint: runtime-check per-rank n
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*muscl_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))
        ! lint: runtime-check per-rank p
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*muscl_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))

    end subroutine s_check_inputs_muscl

    !> Validate NVIDIA unified virtual memory configuration parameters
    impure subroutine s_check_inputs_nvidia_uvm

#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu > 3 .or. nv_uvm_igr_temps_on_gpu < 0, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 3]")
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu == 3 .and. igr_iter_solver == 2, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 2] for igr_iter_solver == 2")
#endif

    end subroutine s_check_inputs_nvidia_uvm

    !> Validates that each burning immersed-boundary patch injects a species index within the mechanism. inj_species indexes the
    !! image-point mass-fraction array Ys_IP(1:num_species) in m_ibm; an out-of-range value is an out-of-bounds write (silent
    !! corruption). Only reachable with chemistry.
    impure subroutine s_check_inputs_ib_injection

        integer :: i

        do i = 1, num_ibs
            @:PROHIBIT(patch_ib(i)%inj_species > num_species, &
                       & "patch_ib inj_species must be <= num_species (it indexes the image-point species mass fractions; an out-of-range value writes out of bounds)")
        end do

    end subroutine s_check_inputs_ib_injection

end module m_checker
