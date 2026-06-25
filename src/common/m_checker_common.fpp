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
    use m_constants, only: eos_stiffened_gas, eos_jwl, model_eqns_5eq, jwl_mix_type_rocflu

    implicit none

    private; public :: s_check_inputs_common, wp

contains

    !> Checks compatibility of parameters in the input file. Used by all three stages
    impure subroutine s_check_inputs_common

#ifndef MFC_SIMULATION
        call s_check_total_cells
#endif
        call s_check_jwl_inputs
        #:if USING_AMD
            call s_check_amd
        #:endif

    end subroutine s_check_inputs_common

    !> Validate JWL EOS selector combinations and required material parameters.
    impure subroutine s_check_jwl_inputs

        integer           :: i, n_jwl, jwl_fluid, air_fluid
        character(len=10) :: iStr

        n_jwl = 0
        jwl_fluid = 0
        air_fluid = 0

        do i = 1, num_fluids
            call s_int_to_str(i, iStr)

            @:PROHIBIT(fluid_pp(i)%eos /= eos_stiffened_gas .and. fluid_pp(i)%eos /= eos_jwl, &
                       & "fluid_pp(" // trim(iStr) // ")%eos must be eos_stiffened_gas (1) or eos_jwl (2)")

            if (fluid_pp(i)%eos == eos_jwl) then
                n_jwl = n_jwl + 1
                jwl_fluid = i

                @:PROHIBIT(f_is_default(fluid_pp(i)%jwl_A) .or. f_is_default(fluid_pp(i)%jwl_B) &
                           & .or. f_is_default(fluid_pp(i)%jwl_R1) .or. f_is_default(fluid_pp(i)%jwl_R2) &
                           & .or. f_is_default(fluid_pp(i)%jwl_omega) .or. f_is_default(fluid_pp(i)%jwl_rho0) &
                           & .or. f_is_default(fluid_pp(i)%jwl_E0), &
                           & "fluid_pp(" // trim(iStr) // ")%eos = eos_jwl requires jwl_A, jwl_B, jwl_R1, " &
                           & // "jwl_R2, jwl_omega, jwl_rho0, and jwl_E0")
                @:PROHIBIT(f_is_default(fluid_pp(i)%jwl_air_e0) .or. f_is_default(fluid_pp(i)%jwl_air_rho0) &
                           & .or. f_is_default(fluid_pp(i)%jwl_air_gamma), &
                           & "fluid_pp(" // trim(iStr) // ")%eos = eos_jwl requires jwl_air_e0, " &
                           & // "jwl_air_rho0, and jwl_air_gamma")
                @:PROHIBIT(fluid_pp(i)%jwl_R1 <= 0._wp .or. fluid_pp(i)%jwl_R2 <= 0._wp .or. fluid_pp(i)%jwl_omega <= 0._wp &
                           & .or. fluid_pp(i)%jwl_rho0 <= 0._wp .or. fluid_pp(i)%jwl_air_rho0 <= 0._wp &
                           & .or. fluid_pp(i)%jwl_air_gamma <= 0._wp, &
                           & "JWL parameters jwl_R1, jwl_R2, jwl_omega, jwl_rho0, " &
                           & // "jwl_air_rho0, and jwl_air_gamma must be positive")
            else if (air_fluid == 0) then
                air_fluid = i
            end if
        end do

        @:PROHIBIT(n_jwl > 1, "At most one fluid may use eos_jwl")
        @:PROHIBIT(jwl_fluid > 0 .and. model_eqns /= model_eqns_5eq, "JWL EOS is only supported with model_eqns_5eq")

        @:PROHIBIT(jwl_fluid > 0 .and. jwl_mix_type /= jwl_mix_type_rocflu, &
                   & "Only jwl_mix_type = 3 (Rocflu) is supported; modes 0, 1, and 2 have been removed")

        if (jwl_fluid > 0 .and. jwl_mix_type == jwl_mix_type_rocflu) then
            @:PROHIBIT(num_fluids > 1 .and. air_fluid == 0, "Rocflu closure requires one non-JWL ideal-gas fluid")
            @:PROHIBIT(f_is_default(fluid_pp(jwl_fluid)%cv) .or. fluid_pp(jwl_fluid)%cv <= 0._wp, &
                       & "Rocflu closure requires positive fluid_pp%cv for the JWL fluid")
            @:PROHIBIT(air_fluid > 0 .and. (f_is_default(fluid_pp(air_fluid)%cv) .or. fluid_pp(air_fluid)%cv <= 0._wp), &
                       & "Rocflu closure requires positive fluid_pp%cv for the non-JWL air fluid")
            @:PROHIBIT(fluid_pp(jwl_fluid)%jwl_rho0 <= fluid_pp(jwl_fluid)%jwl_air_rho0 &
                       & .or. fluid_pp(jwl_fluid)%jwl_E0/fluid_pp(jwl_fluid)%jwl_rho0 <= fluid_pp(jwl_fluid)%jwl_air_e0, &
                       & "Rocflu closure requires increasing air-to-products reference density and energy")
        end if

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
