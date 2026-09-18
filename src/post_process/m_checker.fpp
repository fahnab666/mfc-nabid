!>
!!@file
!!@brief Contains module m_checker

#:include 'macros.fpp'

!> @brief Validates post-process input parameters and output format consistency
module m_checker

    use m_global_parameters
    use m_mpi_proxy
    use m_helper_basic
    use m_helper
    use m_constants, only: eos_stiffened_gas, eos_ideal_gas

    implicit none

    private; public :: s_check_inputs, s_check_inputs_fft

contains

    !> Checks compatibility of parameters in the input file. Used by the post_process stage
    impure subroutine s_check_inputs

        if (lso_filter_wrt .or. lso_stat_wrt .or. lso_pp_filter .or. lso_closure_wrt) then
            @:PROHIBIT(lso_stat_wrt .and. .not. parallel_io, "LSO statistical output requires parallel_io")
            @:PROHIBIT(lso_stat_wrt .and. .not. lso_filter_wrt, "lso_stat_wrt requires lso_filter_wrt")
            @:PROHIBIT(lso_pp_filter .and. .not. lso_filter_wrt, "lso_pp_filter requires lso_filter_wrt")
            @:PROHIBIT(lso_closure_wrt .and. .not. lso_stat_wrt, "lso_closure_wrt requires lso_stat_wrt")
        end if

        if (lso_closure_wrt) then
            @:PROHIBIT(num_fluids /= 1, "LSO closures currently require num_fluids = 1")
            @:PROHIBIT(chemistry, "LSO closures do not support chemistry")
            @:PROHIBIT(fluid_pp(1)%eos /= eos_stiffened_gas .and. fluid_pp(1)%eos /= eos_ideal_gas, &
                       & "LSO closures require a calorically perfect ideal or stiffened gas")
            @:PROHIBIT(lso_R_gas <= 0._wp, "LSO closures require lso_R_gas > 0")
        end if

    end subroutine s_check_inputs

    !> Checks constraints on fft_wrt
    impure subroutine s_check_inputs_fft

        integer :: num_procs_y, num_procs_z

        @:PROHIBIT(fft_wrt .and. MOD(n_glb+1,n+1) /= 0, "FFT WRT requires n_glb to be divisible by num_procs_y")
        @:PROHIBIT(fft_wrt .and. MOD(p_glb+1,p+1) /= 0, "FFT WRT requires p_glb to be divisible by num_procs_z")
        num_procs_y = (n_glb + 1)/(n + 1)
        num_procs_z = (p_glb + 1)/(p + 1)
        @:PROHIBIT(fft_wrt .and. MOD(m_glb+1,num_procs_y) /= 0, "FFT WRT requires m_glb to be divisible by num_procs_y")
        @:PROHIBIT(fft_wrt .and. MOD(n_glb+1,num_procs_z) /= 0, "FFT WRT requires n_glb to be divisible by num_procs_z")

    end subroutine s_check_inputs_fft

end module m_checker
