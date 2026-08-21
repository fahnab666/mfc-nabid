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
    impure subroutine s_check_inputs_common(check_total_cells, n_global)

        logical, intent(in)         :: check_total_cells
        integer(kind=8), intent(in) :: n_global

        if (check_total_cells) call s_check_total_cells(n_global)
        call s_derive_jwl_energy_defaults
        #:if USING_AMD
            call s_check_amd
        #:endif

    end subroutine s_check_inputs_common

    !> Fills in the JWL fluid's jwl_E0 or jwl_Q from the other plus jwl_rho0, whichever the case file omitted. All JWL input
    !! constraints, including the E0/Q/rho0 consistency check when the case file sets all three, live in
    !! toolchain/mfc/case_validator.py (check_stiffened_eos).
    impure subroutine s_derive_jwl_energy_defaults

        integer :: i

        do i = 1, num_fluids
            if (fluid_pp(i)%eos == eos_jwl .and. .not. f_is_default(fluid_pp(i)%jwl_rho0)) then
                if (f_is_default(fluid_pp(i)%jwl_E0) .and. .not. f_is_default(fluid_pp(i)%jwl_Q)) then
                    fluid_pp(i)%jwl_E0 = fluid_pp(i)%jwl_rho0*fluid_pp(i)%jwl_Q
                else if (.not. f_is_default(fluid_pp(i)%jwl_E0) .and. f_is_default(fluid_pp(i)%jwl_Q)) then
                    fluid_pp(i)%jwl_Q = fluid_pp(i)%jwl_E0/fluid_pp(i)%jwl_rho0
                end if
            end if
        end do

    end subroutine s_derive_jwl_energy_defaults

    !> Verify that the total number of grid cells meets the minimum required by the number of dimensions and MPI ranks.
    impure subroutine s_check_total_cells(n_global)

        character(len=18)           :: numStr  !< for int to string conversion
        integer(kind=8)             :: min_cells
        integer(kind=8), intent(in) :: n_global

        min_cells = int(2, kind=8)**int(min(1, m) + min(1, n) + min(1, p), kind=8)*int(num_procs, kind=8)
        call s_int_to_str(2**(min(1, m) + min(1, n) + min(1, p))*num_procs, numStr)

        @:PROHIBIT(n_global < min_cells, &
                   & "Total number of cells must be at least (2^[number of dimensions])*num_procs, " // "which is currently " &
                   & // trim(numStr))

    end subroutine s_check_total_cells

    !> Check that simulation parameters stay within AMD GPU compiler limits when case optimization is disabled.
    impure subroutine s_check_amd

        #:if not MFC_CASE_OPTIMIZATION
            @:PROHIBIT(num_fluids > 3, "num_fluids <= 3 for AMDFLang when Case optimization is off")
            @:PROHIBIT((bubbles_euler .or. bubbles_lagrange) .and. nb > 3, "nb <= 3 for AMDFLang when Case optimization is off")
            @:PROHIBIT(chemistry .and. num_species > 10, "num_species > 10 for AMDFLang when Case optimization is off")
        #:endif

    end subroutine s_check_amd

end module m_checker_common
