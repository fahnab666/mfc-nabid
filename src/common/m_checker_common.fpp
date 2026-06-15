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

    implicit none

    private; public :: s_check_inputs_common, wp

contains

    !> Checks compatibility of parameters in the input file. Used by all three stages
    impure subroutine s_check_inputs_common

#ifndef MFC_SIMULATION
        call s_check_total_cells
#endif
        #:if USING_AMD
            call s_check_amd
        #:endif
        call s_check_jwl

    end subroutine s_check_inputs_common

    !> Require every JWL fluid (eos == 2) to set its parameters explicitly; they default to dflt_real, so an unset value is a
    !! missing input rather than a silent magic number.
    impure subroutine s_check_jwl

        integer :: i

        do i = 1, num_fluids
            if (fluid_pp(i)%eos /= 2) cycle
            @:PROHIBIT(f_is_default(fluid_pp(i)%jwl_A) .or. f_is_default(fluid_pp(i)%jwl_B) .or. f_is_default(fluid_pp(i)%jwl_R1) &
                       & .or. f_is_default(fluid_pp(i)%jwl_R2) .or. f_is_default(fluid_pp(i)%jwl_omega) &
                       & .or. f_is_default(fluid_pp(i)%jwl_rho0) .or. f_is_default(fluid_pp(i)%jwl_E0), &
                       & "A JWL fluid (eos = 2) requires jwl_A, jwl_B, jwl_R1, jwl_R2, jwl_omega, jwl_rho0 and jwl_E0")
            @:PROHIBIT(f_is_default(fluid_pp(i)%jwl_air_e0) .or. f_is_default(fluid_pp(i)%jwl_air_rho0) &
                       & .or. f_is_default(fluid_pp(i)%jwl_air_gamma), &
                       & "A JWL fluid (eos = 2) requires jwl_air_e0, jwl_air_rho0 and jwl_air_gamma")
        end do

    end subroutine s_check_jwl

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
