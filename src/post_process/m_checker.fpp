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
    use m_constants, only: eos_jwl

    implicit none

    private; public :: s_check_inputs, s_check_inputs_fft

contains

    !> Checks compatibility of parameters in the input file. Used by the post_process stage
    impure subroutine s_check_inputs

        ! The sim_data energy diagnostics build the mixture gamma/pi_inf from the stiffened
        ! gammas/pi_infs arrays, which a JWL fluid leaves at their dflt_real placeholders, so the
        ! derived sound speed and Mach output would be garbage. Prohibit the combination.
        @:PROHIBIT(sim_data .and. any(fluid_pp(1:num_fluids)%eos == eos_jwl), &
                   & "sim_data is not supported with fluid_pp(:)%eos = 'jwl'")

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
