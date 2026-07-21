!>
!!@file
!!@brief Contains module m_checker_common

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Shared input validation checks for grid dimensions and AMD GPU compiler limits
module m_checker_common

    use m_global_parameters
    use m_mpi_proxy
    use m_helper_basic
    use m_helper
    use m_constants, only: eos_jwl

    implicit none

    private; public :: s_check_inputs_common

contains

    !> Checks compatibility of parameters in the input file. Used by all three stages
    impure subroutine s_check_inputs_common

#ifndef MFC_SIMULATION
        call s_check_total_cells
#endif
        call s_check_eos
        #:if USING_AMD
            call s_check_amd
        #:endif

    end subroutine s_check_inputs_common

    !> Restrict the per-fluid EOS selector to the currently supported adapters. Only stiffened_gas (non-chemistry) and
    !! ideal_gas_mixture (chemistry) are backed by a thermodynamics backend, and a single run uses one family for every fluid, so
    !! the unimplemented values and any intra-cell EOS mixing are both rejected here.
    impure subroutine s_check_eos

        integer :: i

        do i = 1, num_fluids
            @:PROHIBIT(chemistry .and. fluid_pp(i)%eos /= eos_ideal_gas_mixture, &
                       & "fluid_pp(:)%eos must be 'ideal_gas_mixture' for every fluid when chemistry is enabled")
            @:PROHIBIT(.not. chemistry .and. fluid_pp(i)%eos /= eos_stiffened_gas .and. fluid_pp(i)%eos /= eos_jwl, &
                       & "fluid_pp(:)%eos selector is not supported; only 'stiffened_gas' and 'jwl' are available " &
                       & // "(or 'ideal_gas_mixture' with a chemistry build)")
        end do

        @:PROHIBIT(mixture_closure /= mixture_closure_pressure_equilibrium, &
                   & "mixture_closure selects a reserved closure; only 'pressure_equilibrium' is implemented")

        if (any(fluid_pp(1:num_fluids)%eos == eos_jwl)) call s_check_jwl_inputs

    end subroutine s_check_eos

    !> JWL runs solve the 5-equation model through the generalized Mie-Gruneisen mixture accumulators; restrict to that path and
    !! require a well-posed user-supplied isentrope. Features whose EOS algebra bypasses the accumulators are rejected.
    impure subroutine s_check_jwl_inputs

        integer :: i

        @:PROHIBIT(model_eqns /= model_eqns_5eq, "fluid_pp(:)%eos = 'jwl' requires model_eqns = 2 (5-equation model)")
        @:PROHIBIT(bubbles_euler .or. bubbles_lagrange, "fluid_pp(:)%eos = 'jwl' is not supported with bubbles")
        @:PROHIBIT(relax, "fluid_pp(:)%eos = 'jwl' is not supported with phase change (relax)")
        @:PROHIBIT(hypoelasticity .or. hyperelasticity, "fluid_pp(:)%eos = 'jwl' is not supported with elasticity")
        @:PROHIBIT(mhd, "fluid_pp(:)%eos = 'jwl' is not supported with mhd")
        @:PROHIBIT(surface_tension, "fluid_pp(:)%eos = 'jwl' is not supported with surface_tension")
        @:PROHIBIT(igr, "fluid_pp(:)%eos = 'jwl' is not supported with igr")
        @:PROHIBIT(ib, "fluid_pp(:)%eos = 'jwl' is not supported with immersed boundaries")

        do i = 1, num_fluids
            if (fluid_pp(i)%eos /= eos_jwl) cycle
            @:PROHIBIT(.not. fluid_pp(i)%jwl_A > 0._wp, "fluid_pp(i)%jwl_A must be set and positive for a JWL fluid")
            @:PROHIBIT(.not. fluid_pp(i)%jwl_B > 0._wp, "fluid_pp(i)%jwl_B must be set and positive for a JWL fluid")
            @:PROHIBIT(.not. (fluid_pp(i)%jwl_R2 > 0._wp .and. fluid_pp(i)%jwl_R1 > fluid_pp(i)%jwl_R2), &
                       & "a JWL fluid requires R1 > R2 > 0")
            @:PROHIBIT(.not. fluid_pp(i)%jwl_omega > 0._wp, "fluid_pp(i)%jwl_omega must be set and positive for a JWL fluid")
            @:PROHIBIT(.not. fluid_pp(i)%jwl_rho0 > 0._wp, "fluid_pp(i)%jwl_rho0 must be set and positive for a JWL fluid")
            @:PROHIBIT(abs(fluid_pp(i)%qv) > 0._wp .or. abs(fluid_pp(i)%qvp) > 0._wp, &
                       & "fluid_pp(i)%qv and qvp are unused by a JWL fluid; its reference energy is the principal isentrope")
        end do

    end subroutine s_check_jwl_inputs

#ifndef MFC_SIMULATION
    !> Verify that the total number of grid cells meets the minimum required by the number of dimensions and MPI ranks.
    impure subroutine s_check_total_cells

        character(len=18) :: numStr  !< for int to string conversion
        integer(kind=8)   :: min_cells

        min_cells = int(2, kind=8)**int(min(1, m) + min(1, n) + min(1, p), kind=8)*int(num_procs, kind=8)
        call s_int_to_str(2**(min(1, m) + min(1, n) + min(1, p))*num_procs, numStr)

        @:PROHIBIT(nGlobal < min_cells, &
                   & "Total number of cells must be at least (2^[number of dimensions])*num_procs, " // "which is currently " &
                   & // trim(numStr))

    end subroutine s_check_total_cells
#endif

    !> Check that simulation parameters stay within AMD GPU compiler limits when case optimization is disabled.
    impure subroutine s_check_amd

        #:if not MFC_CASE_OPTIMIZATION
            @:PROHIBIT(num_fluids > 3, "num_fluids <= 3 for AMDFLang when Case optimization is off")
            @:PROHIBIT((bubbles_euler .or. bubbles_lagrange) .and. nb > 3, "nb <= 3 for AMDFLang when Case optimization is off")
            @:PROHIBIT(chemistry .and. num_species > 10, "num_species > 10 for AMDFLang when Case optimization is off")
        #:endif

    end subroutine s_check_amd

#ifndef MFC_POST_PROCESS
#endif
end module m_checker_common
