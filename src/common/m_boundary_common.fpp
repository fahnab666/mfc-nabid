!>
!! @file
!! Contains module m_boundary_common

!> @brief Noncharacteristic and processor boundary condition application for ghost cells and buffer regions
#:include 'case.fpp'
#:include 'macros.fpp'

module m_boundary_common

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_constants
    use m_boundary_primitives
    use m_boundary_io
#ifndef MFC_PRE_PROCESS
    use m_mpi_common, only: s_mpi_sendrecv_grid_variables_buffers
#endif

    implicit none

    private; public :: s_initialize_boundary_common_module, s_populate_variables_buffers, s_populate_beta_buffers, &
        & s_populate_capillary_buffers, s_populate_F_igr_buffers, s_populate_grid_variables_buffers, &
        & s_finalize_boundary_common_module

    type(int_bounds_info), dimension(3) :: beta_bc_bounds
    $:GPU_DECLARE(create='[beta_bc_bounds]')

    public :: bc_buffers

#ifdef MFC_MPI
    public :: MPI_BC_TYPE_TYPE, MPI_BC_BUFFER_TYPE
#endif

contains

    !> Allocates and sets up boundary condition buffer arrays for all coordinate directions.
    impure subroutine s_initialize_boundary_common_module()

        integer :: i, j, sys_size_alloc

        @:ALLOCATE(bc_buffers(1:3, 1:2))

        if (bc_io) then
            sys_size_alloc = sys_size
            if (chemistry) sys_size_alloc = sys_size + 1

            @:ALLOCATE(bc_buffers(1, 1)%sf(1:sys_size_alloc, 0:n, 0:p))
            @:ALLOCATE(bc_buffers(1, 2)%sf(1:sys_size_alloc, 0:n, 0:p))
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                if (n > 0) then
                    @:ALLOCATE(bc_buffers(2,1)%sf(-buff_size:m+buff_size,1:sys_size_alloc,0:p))
                    @:ALLOCATE(bc_buffers(2,2)%sf(-buff_size:m+buff_size,1:sys_size_alloc,0:p))
                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                        if (p > 0) then
                            @:ALLOCATE(bc_buffers(3,1)%sf(-buff_size:m+buff_size,-buff_size:n+buff_size,1:sys_size_alloc))
                            @:ALLOCATE(bc_buffers(3,2)%sf(-buff_size:m+buff_size,-buff_size:n+buff_size,1:sys_size_alloc))
                        end if
                    #:endif
                end if
            #:endif
            do i = 1, num_dims
                do j = 1, 2
                    @:ACC_SETUP_SFs(bc_buffers(i,j))
                end do
            end do
        end if

        if (bubbles_lagrange .or. particles_lagrange) then
            beta_bc_bounds(1)%beg = -mapcells - 1
            beta_bc_bounds(1)%end = m + mapcells + 1
            ! n > 0 always for bubbles_lagrange
            beta_bc_bounds(2)%beg = -mapcells - 1
            beta_bc_bounds(2)%end = n + mapcells + 1
            if (p == 0) then
                beta_bc_bounds(3)%beg = 0
                beta_bc_bounds(3)%end = 0
            else
                beta_bc_bounds(3)%beg = -mapcells - 1
                beta_bc_bounds(3)%end = p + mapcells + 1
            end if
        end if
        $:GPU_UPDATE(device='[beta_bc_bounds]')

    end subroutine s_initialize_boundary_common_module

    !> Populate the buffers of the primitive variables based on the selected boundary conditions.
    impure subroutine s_populate_variables_buffers(bc_type, q_prim_vf, pb_in, mv_in, q_T_sf)

        type(scalar_field), dimension(sys_size), intent(inout)                                               :: q_prim_vf
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        type(integer_field), dimension(1:num_dims,1:2), intent(in)                                           :: bc_type
        type(scalar_field), optional, intent(inout)                                                          :: q_T_sf

        call s_populate_bc_direction(1, -1, bc_x, bc_type(1, 1), q_prim_vf, pb_in, mv_in, q_T_sf)
        call s_populate_bc_direction(1, 1, bc_x, bc_type(1, 2), q_prim_vf, pb_in, mv_in, q_T_sf)

        if (n == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            call s_populate_bc_direction(2, -1, bc_y, bc_type(2, 1), q_prim_vf, pb_in, mv_in, q_T_sf)
            call s_populate_bc_direction(2, 1, bc_y, bc_type(2, 2), q_prim_vf, pb_in, mv_in, q_T_sf)
        #:endif

        if (p == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            call s_populate_bc_direction(3, -1, bc_z, bc_type(3, 1), q_prim_vf, pb_in, mv_in, q_T_sf)
            call s_populate_bc_direction(3, 1, bc_z, bc_type(3, 2), q_prim_vf, pb_in, mv_in, q_T_sf)
        #:endif
        ! END: Population of Buffers in z-direction

    end subroutine s_populate_variables_buffers

    !> Populate the variable buffers along one direction and location, via MPI exchange for processor boundaries or by dispatching
    !! the per-cell BC routines over the boundary face.
    impure subroutine s_populate_bc_direction(bc_dir, bc_loc, bc_bounds, bc_type_edge, q_prim_vf, pb_in, mv_in, q_T_sf)

        integer, intent(in) :: bc_dir, bc_loc
        type(int_bounds_info), intent(in) :: bc_bounds
        type(integer_field), intent(in) :: bc_type_edge
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        real(stp), optional, dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:,1:), intent(inout) :: pb_in, mv_in
        type(scalar_field), optional, intent(inout) :: q_T_sf
        integer :: bc_edge, k_beg, k_end, l_beg, l_end
        integer :: bc_code, k, l

        if (bc_loc == -1) then
            bc_edge = bc_bounds%beg
        else
            bc_edge = bc_bounds%end
        end if

        ! BC type codes defined in m_constants.fpp; non-negative values are MPI boundaries
        if (bc_edge >= 0) then
            call s_mpi_sendrecv_variables_buffers(q_prim_vf, bc_dir, bc_loc, sys_size, pb_in, mv_in, q_T_sf)
            return
        end if

        if (bc_dir == 1) then
            k_beg = 0; k_end = n; l_beg = 0; l_end = p
        else if (bc_dir == 2) then
            k_beg = -buff_size; k_end = m + buff_size; l_beg = 0; l_end = p
        else
            k_beg = -buff_size; k_end = m + buff_size; l_beg = -buff_size; l_end = n + buff_size
        end if

        $:GPU_PARALLEL_LOOP(private='[l, k, bc_code]', collapse=2)
        do l = l_beg, l_end
            do k = k_beg, k_end
                if (bc_dir == 1) then
                    bc_code = int(bc_type_edge%sf(0, k, l))
                else if (bc_dir == 2) then
                    bc_code = int(bc_type_edge%sf(k, 0, l))
                else
                    bc_code = int(bc_type_edge%sf(k, l, 0))
                end if

                select case (bc_code)
                case (BC_CHAR_SUP_OUTFLOW:BC_GHOST_EXTRAP)
                    call s_ghost_cell_extrapolation(q_prim_vf, bc_dir, bc_loc, k, l, q_T_sf)
                case (BC_AXIS)
                    if (bc_dir == 2 .and. bc_loc == -1) call s_axis(q_prim_vf, pb_in, mv_in, k, l)
                case (BC_REFLECTIVE)
                    call s_symmetry(q_prim_vf, bc_dir, bc_loc, k, l, pb_in, mv_in, q_T_sf)
                case (BC_PERIODIC)
                    call s_periodic(q_prim_vf, bc_dir, bc_loc, k, l, pb_in, mv_in, q_T_sf)
                case (BC_SLIP_WALL)
                    call s_slip_wall(q_prim_vf, bc_dir, bc_loc, k, l, q_T_sf)
                case (BC_NO_SLIP_WALL)
                    call s_no_slip_wall(q_prim_vf, bc_dir, bc_loc, k, l, q_T_sf)
                case (BC_DIRICHLET)
                    call s_dirichlet(q_prim_vf, bc_dir, bc_loc, k, l, q_T_sf)
                end select

                if (qbmm .and. (.not. polytropic) .and. present(pb_in) .and. present(mv_in) .and. (bc_code <= BC_GHOST_EXTRAP) &
                    & .and. .not. (bc_dir == 2 .and. bc_loc == -1 .and. bc_code == BC_AXIS)) then
                    call s_qbmm_extrapolation(bc_dir, bc_loc, k, l, pb_in, mv_in)
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_populate_bc_direction

    impure subroutine s_populate_beta_buffers(q_beta, kahan_comp, bc_type, nvar, vars_comm)

        type(scalar_field), dimension(:), intent(inout)            :: q_beta
        type(scalar_field), dimension(:), intent(inout)            :: kahan_comp
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        integer, intent(in)                                        :: nvar
        integer, dimension(:), intent(in)                          :: vars_comm
        integer                                                    :: k, l

        !> x-direction

        if (bc_x%beg < 0) then
            $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2, copyin = '[vars_comm]')
            do l = beta_bc_bounds(3)%beg, beta_bc_bounds(3)%end
                do k = beta_bc_bounds(2)%beg, beta_bc_bounds(2)%end
                    select case (bc_x%beg)
                    case (BC_PERIODIC)
                        call s_beta_periodic(q_beta, kahan_comp, 1, -1, k, l, nvar, vars_comm)
                    case (BC_REFLECTIVE)
                        call s_beta_reflective(q_beta, kahan_comp, 1, -1, k, l, nvar, vars_comm)
                    case default
                    end select
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        if (bc_x%beg >= 0 .or. bc_x%end >= 0) then
            call s_mpi_reduce_beta_variables_buffers(q_beta, kahan_comp, 1, -1, nvar, vars_comm)
        end if

        if (bc_x%end < 0) then
            $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2, copyin = '[vars_comm]')
            do l = beta_bc_bounds(3)%beg, beta_bc_bounds(3)%end
                do k = beta_bc_bounds(2)%beg, beta_bc_bounds(2)%end
                    select case (bc_x%end)
                    case (BC_PERIODIC)
                        call s_beta_periodic(q_beta, kahan_comp, 1, 1, k, l, nvar, vars_comm)
                    case (BC_REFLECTIVE)
                        call s_beta_reflective(q_beta, kahan_comp, 1, 1, k, l, nvar, vars_comm)
                    case default
                    end select
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        if (bc_x%beg >= 0 .or. bc_x%end >= 0) then
            call s_mpi_reduce_beta_variables_buffers(q_beta, kahan_comp, 1, 1, nvar, vars_comm)
        end if

        !> y-direction
        if (bc_y%beg < 0) then
            $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2, copyin = '[vars_comm]')
            do l = beta_bc_bounds(3)%beg, beta_bc_bounds(3)%end
                do k = beta_bc_bounds(1)%beg, beta_bc_bounds(1)%end
                    select case (bc_y%beg)
                    case (BC_PERIODIC)
                        call s_beta_periodic(q_beta, kahan_comp, 2, -1, k, l, nvar, vars_comm)
                    case (BC_REFLECTIVE)
                        call s_beta_reflective(q_beta, kahan_comp, 2, -1, k, l, nvar, vars_comm)
                    case default
                    end select
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        if (bc_y%beg >= 0 .or. bc_y%end >= 0) then
            call s_mpi_reduce_beta_variables_buffers(q_beta, kahan_comp, 2, -1, nvar, vars_comm)
        end if

        if (bc_y%end < 0) then
            $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2, copyin = '[vars_comm]')
            do l = beta_bc_bounds(3)%beg, beta_bc_bounds(3)%end
                do k = beta_bc_bounds(1)%beg, beta_bc_bounds(1)%end
                    select case (bc_y%end)
                    case (BC_PERIODIC)
                        call s_beta_periodic(q_beta, kahan_comp, 2, 1, k, l, nvar, vars_comm)
                    case (BC_REFLECTIVE)
                        call s_beta_reflective(q_beta, kahan_comp, 2, 1, k, l, nvar, vars_comm)
                    case default
                    end select
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if
        if (bc_y%beg >= 0 .or. bc_y%end >= 0) then
            call s_mpi_reduce_beta_variables_buffers(q_beta, kahan_comp, 2, 1, nvar, vars_comm)
        end if

        if (num_dims == 2) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            !> z-direction
            if (bc_z%beg < 0) then
                $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2, copyin = '[vars_comm]')
                do l = beta_bc_bounds(2)%beg, beta_bc_bounds(2)%end
                    do k = beta_bc_bounds(1)%beg, beta_bc_bounds(1)%end
                        select case (bc_type(3, 1)%sf(k, l, 0))
                        case (BC_PERIODIC)
                            call s_beta_periodic(q_beta, kahan_comp, 3, -1, k, l, nvar, vars_comm)
                        case (BC_REFLECTIVE)
                            call s_beta_reflective(q_beta, kahan_comp, 3, -1, k, l, nvar, vars_comm)
                        case default
                        end select
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            if (bc_z%beg >= 0 .or. bc_z%end >= 0) then
                call s_mpi_reduce_beta_variables_buffers(q_beta, kahan_comp, 3, -1, nvar, vars_comm)
            end if

            if (bc_z%end < 0) then
                $:GPU_PARALLEL_LOOP(private='[l, k]', collapse=2, copyin = '[vars_comm]')
                do l = beta_bc_bounds(2)%beg, beta_bc_bounds(2)%end
                    do k = beta_bc_bounds(1)%beg, beta_bc_bounds(1)%end
                        select case (bc_type(3, 2)%sf(k, l, 0))
                        case (BC_PERIODIC)
                            call s_beta_periodic(q_beta, kahan_comp, 3, 1, k, l, nvar, vars_comm)
                        case (BC_REFLECTIVE)
                            call s_beta_reflective(q_beta, kahan_comp, 3, 1, k, l, nvar, vars_comm)
                        case default
                        end select
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
            if (bc_z%beg >= 0 .or. bc_z%end >= 0) then
                call s_mpi_reduce_beta_variables_buffers(q_beta, kahan_comp, 3, 1, nvar, vars_comm)
            end if
        #:endif

    end subroutine s_populate_beta_buffers

    subroutine s_beta_periodic(q_beta, kahan_comp, bc_dir, bc_loc, k, l, nvar, vars_comm)

        $:GPU_ROUTINE(function_name='s_beta_periodic', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(1:), intent(inout) :: q_beta
        type(scalar_field), dimension(1:), intent(inout) :: kahan_comp
        integer, intent(in)                              :: bc_dir, bc_loc
        integer, intent(in)                              :: k, l
        integer, intent(in)                              :: nvar
        integer, dimension(:), intent(in)                :: vars_comm
        integer                                          :: j, i
        real(wp)                                         :: y_kahan, t_kahan

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc_x%beg
                do i = 1, nvar
                    do j = -mapCells - 1, mapCells
                        ! Kahan-compensated addition of ghost to interior
                        y_kahan = real(q_beta(vars_comm(i))%sf(m + j + 1, k, l), &
                                       & kind=wp) + kahan_comp(vars_comm(i))%sf(m + j + 1, k, l) - kahan_comp(vars_comm(i))%sf(j, &
                                       & k, l)
                        t_kahan = real(q_beta(vars_comm(i))%sf(j, k, l), kind=wp) + y_kahan
                        kahan_comp(vars_comm(i))%sf(j, k, l) = (t_kahan - q_beta(vars_comm(i))%sf(j, k, l)) - y_kahan
                        q_beta(vars_comm(i))%sf(j, k, l) = t_kahan
                    end do
                end do
            else  !< bc_x%end
                do i = 1, nvar
                    do j = -mapcells, mapcells + 1
                        q_beta(vars_comm(i))%sf(m + j, k, l) = q_beta(vars_comm(i))%sf(j - 1, k, l)
                        kahan_comp(vars_comm(i))%sf(m + j, k, l) = kahan_comp(vars_comm(i))%sf(j - 1, k, l)
                    end do
                end do
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, nvar
                    do j = -mapcells - 1, mapcells
                        y_kahan = real(q_beta(vars_comm(i))%sf(k, n + j + 1, l), kind=wp) + kahan_comp(vars_comm(i))%sf(k, &
                                       & n + j + 1, l) - kahan_comp(vars_comm(i))%sf(k, j, l)
                        t_kahan = real(q_beta(vars_comm(i))%sf(k, j, l), kind=wp) + y_kahan
                        kahan_comp(vars_comm(i))%sf(k, j, l) = (t_kahan - q_beta(vars_comm(i))%sf(k, j, l)) - y_kahan
                        q_beta(vars_comm(i))%sf(k, j, l) = t_kahan
                    end do
                end do
            else  !< bc_y%end
                do i = 1, nvar
                    do j = -mapcells, mapcells + 1
                        q_beta(vars_comm(i))%sf(k, n + j, l) = q_beta(vars_comm(i))%sf(k, j - 1, l)
                        kahan_comp(vars_comm(i))%sf(k, n + j, l) = kahan_comp(vars_comm(i))%sf(k, j - 1, l)
                    end do
                end do
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, nvar
                    do j = -mapcells - 1, mapcells
                        y_kahan = real(q_beta(vars_comm(i))%sf(k, l, p + j + 1), kind=wp) + kahan_comp(vars_comm(i))%sf(k, l, &
                                       & p + j + 1) - kahan_comp(vars_comm(i))%sf(k, l, j)
                        t_kahan = real(q_beta(vars_comm(i))%sf(k, l, j), kind=wp) + y_kahan
                        kahan_comp(vars_comm(i))%sf(k, l, j) = (t_kahan - q_beta(vars_comm(i))%sf(k, l, j)) - y_kahan
                        q_beta(vars_comm(i))%sf(k, l, j) = t_kahan
                    end do
                end do
            else  !< bc_z%end
                do i = 1, nvar
                    do j = -mapcells, mapcells + 1
                        q_beta(vars_comm(i))%sf(k, l, p + j) = q_beta(vars_comm(i))%sf(k, l, j - 1)
                        kahan_comp(vars_comm(i))%sf(k, l, p + j) = kahan_comp(vars_comm(i))%sf(k, l, j - 1)
                    end do
                end do
            end if
        end if

    end subroutine s_beta_periodic

    subroutine s_beta_extrapolation(q_beta, bc_dir, bc_loc, k, l, nvar, vars_comm)

        $:GPU_ROUTINE(function_name='s_beta_extrapolation', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(1:), intent(inout) :: q_beta
        integer, intent(in)                              :: bc_dir, bc_loc
        integer, intent(in)                              :: k, l
        integer, intent(in)                              :: nvar
        integer, dimension(:), intent(in)                :: vars_comm
        integer                                          :: j, i

        ! Set beta in buffer regions equal to zero

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  ! bc_x%beg
                do i = 1, nvar
                    do j = 1, buff_size
                        q_beta(vars_comm(i))%sf(-j, k, l) = 0._wp
                    end do
                end do
            else  !< bc_x%end
                do i = 1, nvar
                    do j = 1, buff_size
                        q_beta(vars_comm(i))%sf(m + j, k, l) = 0._wp
                    end do
                end do
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, nvar
                    do j = 1, buff_size
                        q_beta(vars_comm(i))%sf(k, -j, l) = 0._wp
                    end do
                end do
            else  !< bc_y%end
                do i = 1, nvar
                    do j = 1, buff_size
                        q_beta(vars_comm(i))%sf(k, n + j, l) = 0._wp
                    end do
                end do
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, nvar
                    do j = 1, buff_size
                        q_beta(vars_comm(i))%sf(k, l, -j) = 0._wp
                    end do
                end do
            else  !< bc_z%end
                do i = 1, nvar
                    do j = 1, buff_size
                        q_beta(vars_comm(i))%sf(k, l, p + j) = 0._wp
                    end do
                end do
            end if
        end if

    end subroutine s_beta_extrapolation

    subroutine s_beta_reflective(q_beta, kahan_comp, bc_dir, bc_loc, k, l, nvar, vars_comm)

        $:GPU_ROUTINE(function_name='s_beta_reflective', parallelism='[seq]', cray_inline=True)
        type(scalar_field), dimension(1:), intent(inout) :: q_beta
        type(scalar_field), dimension(1:), intent(inout) :: kahan_comp
        integer, intent(in)                              :: bc_dir, bc_loc
        integer, intent(in)                              :: k, l
        integer, intent(in)                              :: nvar
        integer, dimension(:), intent(in)                :: vars_comm
        integer                                          :: j, i
        real(wp)                                         :: y_kahan, t_kahan

        ! Reflective BC for void fraction: 1) Fold ghost-cell contributions back onto their mirror interior cells (Kahan) 2) Set
        ! ghost cells = mirror of (now-folded) interior values

        if (bc_dir == 1) then  !< x-direction
            if (bc_loc == -1) then  !< bc_x%beg
                do i = 1, nvar
                    do j = 1, mapCells + 1
                        y_kahan = real(q_beta(vars_comm(i))%sf(-j, k, l), kind=wp) + kahan_comp(vars_comm(i))%sf(-j, k, &
                                       & l) - kahan_comp(vars_comm(i))%sf(j - 1, k, l)
                        t_kahan = real(q_beta(vars_comm(i))%sf(j - 1, k, l), kind=wp) + y_kahan
                        kahan_comp(vars_comm(i))%sf(j - 1, k, l) = (t_kahan - q_beta(vars_comm(i))%sf(j - 1, k, l)) - y_kahan
                        q_beta(vars_comm(i))%sf(j - 1, k, l) = t_kahan
                    end do
                    do j = 1, mapCells + 1
                        q_beta(vars_comm(i))%sf(-j, k, l) = q_beta(vars_comm(i))%sf(j - 1, k, l)
                        kahan_comp(vars_comm(i))%sf(-j, k, l) = kahan_comp(vars_comm(i))%sf(j - 1, k, l)
                    end do
                end do
            else  !< bc_x%end
                do i = 1, nvar
                    do j = 1, mapCells + 1
                        y_kahan = real(q_beta(vars_comm(i))%sf(m + j, k, l), kind=wp) + kahan_comp(vars_comm(i))%sf(m + j, k, &
                                       & l) - kahan_comp(vars_comm(i))%sf(m - (j - 1), k, l)
                        t_kahan = real(q_beta(vars_comm(i))%sf(m - (j - 1), k, l), kind=wp) + y_kahan
                        kahan_comp(vars_comm(i))%sf(m - (j - 1), k, l) = (t_kahan - q_beta(vars_comm(i))%sf(m - (j - 1), k, &
                                   & l)) - y_kahan
                        q_beta(vars_comm(i))%sf(m - (j - 1), k, l) = t_kahan
                    end do
                    do j = 1, mapCells + 1
                        q_beta(vars_comm(i))%sf(m + j, k, l) = q_beta(vars_comm(i))%sf(m - (j - 1), k, l)
                        kahan_comp(vars_comm(i))%sf(m + j, k, l) = kahan_comp(vars_comm(i))%sf(m - (j - 1), k, l)
                    end do
                end do
            end if
        else if (bc_dir == 2) then  !< y-direction
            if (bc_loc == -1) then  !< bc_y%beg
                do i = 1, nvar
                    do j = 1, mapCells + 1
                        y_kahan = real(q_beta(vars_comm(i))%sf(k, -j, l), kind=wp) + kahan_comp(vars_comm(i))%sf(k, -j, &
                                       & l) - kahan_comp(vars_comm(i))%sf(k, j - 1, l)
                        t_kahan = real(q_beta(vars_comm(i))%sf(k, j - 1, l), kind=wp) + y_kahan
                        kahan_comp(vars_comm(i))%sf(k, j - 1, l) = (t_kahan - q_beta(vars_comm(i))%sf(k, j - 1, l)) - y_kahan
                        q_beta(vars_comm(i))%sf(k, j - 1, l) = t_kahan
                    end do
                    do j = 1, mapCells + 1
                        q_beta(vars_comm(i))%sf(k, -j, l) = q_beta(vars_comm(i))%sf(k, j - 1, l)
                        kahan_comp(vars_comm(i))%sf(k, -j, l) = kahan_comp(vars_comm(i))%sf(k, j - 1, l)
                    end do
                end do
            else  !< bc_y%end
                do i = 1, nvar
                    do j = 1, mapCells + 1
                        y_kahan = real(q_beta(vars_comm(i))%sf(k, n + j, l), kind=wp) + kahan_comp(vars_comm(i))%sf(k, n + j, &
                                       & l) - kahan_comp(vars_comm(i))%sf(k, n - (j - 1), l)
                        t_kahan = real(q_beta(vars_comm(i))%sf(k, n - (j - 1), l), kind=wp) + y_kahan
                        kahan_comp(vars_comm(i))%sf(k, n - (j - 1), l) = (t_kahan - q_beta(vars_comm(i))%sf(k, n - (j - 1), &
                                   & l)) - y_kahan
                        q_beta(vars_comm(i))%sf(k, n - (j - 1), l) = t_kahan
                    end do
                    do j = 1, mapCells + 1
                        q_beta(vars_comm(i))%sf(k, n + j, l) = q_beta(vars_comm(i))%sf(k, n - (j - 1), l)
                        kahan_comp(vars_comm(i))%sf(k, n + j, l) = kahan_comp(vars_comm(i))%sf(k, n - (j - 1), l)
                    end do
                end do
            end if
        else if (bc_dir == 3) then  !< z-direction
            if (bc_loc == -1) then  !< bc_z%beg
                do i = 1, nvar
                    do j = 1, mapCells + 1
                        y_kahan = real(q_beta(vars_comm(i))%sf(k, l, -j), kind=wp) + kahan_comp(vars_comm(i))%sf(k, l, &
                                       & -j) - kahan_comp(vars_comm(i))%sf(k, l, j - 1)
                        t_kahan = real(q_beta(vars_comm(i))%sf(k, l, j - 1), kind=wp) + y_kahan
                        kahan_comp(vars_comm(i))%sf(k, l, j - 1) = (t_kahan - q_beta(vars_comm(i))%sf(k, l, j - 1)) - y_kahan
                        q_beta(vars_comm(i))%sf(k, l, j - 1) = t_kahan
                    end do
                    do j = 1, mapCells + 1
                        q_beta(vars_comm(i))%sf(k, l, -j) = q_beta(vars_comm(i))%sf(k, l, j - 1)
                        kahan_comp(vars_comm(i))%sf(k, l, -j) = kahan_comp(vars_comm(i))%sf(k, l, j - 1)
                    end do
                end do
            else  !< bc_z%end
                do i = 1, nvar
                    do j = 1, mapCells + 1
                        y_kahan = real(q_beta(vars_comm(i))%sf(k, l, p + j), kind=wp) + kahan_comp(vars_comm(i))%sf(k, l, &
                                       & p + j) - kahan_comp(vars_comm(i))%sf(k, l, p - (j - 1))
                        t_kahan = real(q_beta(vars_comm(i))%sf(k, l, p - (j - 1)), kind=wp) + y_kahan
                        kahan_comp(vars_comm(i))%sf(k, l, p - (j - 1)) = (t_kahan - q_beta(vars_comm(i))%sf(k, l, &
                                   & p - (j - 1))) - y_kahan
                        q_beta(vars_comm(i))%sf(k, l, p - (j - 1)) = t_kahan
                    end do
                    do j = 1, mapCells + 1
                        q_beta(vars_comm(i))%sf(k, l, p + j) = q_beta(vars_comm(i))%sf(k, l, p - (j - 1))
                        kahan_comp(vars_comm(i))%sf(k, l, p + j) = kahan_comp(vars_comm(i))%sf(k, l, p - (j - 1))
                    end do
                end do
            end if
        end if

    end subroutine s_beta_reflective

    !> Populate ghost cell buffers for the color function and its divergence used in capillary surface tension.
    impure subroutine s_populate_capillary_buffers(c_divs, bc_type, bc)

        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(bc_xyz_info), intent(in)                              :: bc

        call s_populate_capillary_bc_direction(1, -1, bc%x, bc_type(1, 1), c_divs)
        call s_populate_capillary_bc_direction(1, 1, bc%x, bc_type(1, 2), c_divs)

        if (n == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            call s_populate_capillary_bc_direction(2, -1, bc%y, bc_type(2, 1), c_divs)
            call s_populate_capillary_bc_direction(2, 1, bc%y, bc_type(2, 2), c_divs)
        #:endif

        if (p == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            call s_populate_capillary_bc_direction(3, -1, bc%z, bc_type(3, 1), c_divs)
            call s_populate_capillary_bc_direction(3, 1, bc%z, bc_type(3, 2), c_divs)
        #:endif

    end subroutine s_populate_capillary_buffers

    !> Populate ghost cell buffers for one capillary BC direction and location, via MPI exchange for processor boundaries or by
    !! dispatching the per-cell capillary BC routines over the boundary face.
    impure subroutine s_populate_capillary_bc_direction(bc_dir, bc_loc, bc_bounds, bc_type_edge, c_divs)

        integer, intent(in)                                        :: bc_dir, bc_loc
        type(int_bounds_info), intent(in)                          :: bc_bounds
        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs
        type(integer_field), intent(in)                            :: bc_type_edge
        integer                                                    :: bc_edge, k_beg, k_end, l_beg, l_end, k, l, bc_code

        if (bc_loc == -1) then
            bc_edge = bc_bounds%beg
        else
            bc_edge = bc_bounds%end
        end if

        if (bc_edge >= 0) then
            call s_mpi_sendrecv_variables_buffers(c_divs, bc_dir, bc_loc, num_dims + 1)
            return
        end if

        if (bc_dir == 1) then
            k_beg = 0; k_end = n; l_beg = 0; l_end = p
        else if (bc_dir == 2) then
            k_beg = -buff_size; k_end = m + buff_size; l_beg = 0; l_end = p
        else
            k_beg = -buff_size; k_end = m + buff_size; l_beg = -buff_size; l_end = n + buff_size
        end if

        $:GPU_PARALLEL_LOOP(private='[l, k, bc_code]', collapse=2)
        do l = l_beg, l_end
            do k = k_beg, k_end
                if (bc_dir == 1) then
                    bc_code = int(bc_type_edge%sf(0, k, l))
                else if (bc_dir == 2) then
                    bc_code = int(bc_type_edge%sf(k, 0, l))
                else
                    bc_code = int(bc_type_edge%sf(k, l, 0))
                end if

                select case (bc_code)
                case (BC_PERIODIC)
                    call s_color_function_periodic(c_divs, bc_dir, bc_loc, k, l)
                case (BC_REFLECTIVE)
                    call s_color_function_reflective(c_divs, bc_dir, bc_loc, k, l)
                case default
                    call s_color_function_ghost_cell_extrapolation(c_divs, bc_dir, bc_loc, k, l)
                end select
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_populate_capillary_bc_direction

    !> Populates ghost cell buffers for the Jacobian scalar field used in the IGR elliptic solver.
    impure subroutine s_populate_F_igr_buffers(bc_type, jac_sf)

        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), dimension(1:), intent(inout)           :: jac_sf

        call s_populate_F_igr_bc_direction(1, -1, bc_x, bc_type(1, 1), jac_sf)
        call s_populate_F_igr_bc_direction(1, 1, bc_x, bc_type(1, 2), jac_sf)

        if (n == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            call s_populate_F_igr_bc_direction(2, -1, bc_y, bc_type(2, 1), jac_sf)
            call s_populate_F_igr_bc_direction(2, 1, bc_y, bc_type(2, 2), jac_sf)
        #:endif

        if (p == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            call s_populate_F_igr_bc_direction(3, -1, bc_z, bc_type(3, 1), jac_sf)
            call s_populate_F_igr_bc_direction(3, 1, bc_z, bc_type(3, 2), jac_sf)
        #:endif

    end subroutine s_populate_F_igr_buffers

    !> Populate ghost cell buffers for one IGR Jacobian BC direction and location, via MPI exchange for processor boundaries or by
    !! dispatching the per-cell IGR Jacobian BC routines over the boundary face.
    impure subroutine s_populate_F_igr_bc_direction(bc_dir, bc_loc, bc_bounds, bc_type_edge, jac_sf)

        integer, intent(in)                              :: bc_dir, bc_loc
        type(int_bounds_info), intent(in)                :: bc_bounds
        type(integer_field), intent(in)                  :: bc_type_edge
        type(scalar_field), dimension(1:), intent(inout) :: jac_sf
        integer                                          :: bc_edge, k_beg, k_end, l_beg, l_end, k, l, j, bc_code

        if (bc_loc == -1) then
            bc_edge = bc_bounds%beg
        else
            bc_edge = bc_bounds%end
        end if

        if (bc_edge >= 0) then
            call s_mpi_sendrecv_variables_buffers(jac_sf, bc_dir, bc_loc, 1)
            return
        end if

        if (bc_dir == 1) then
            k_beg = 0; k_end = n; l_beg = 0; l_end = p
        else if (bc_dir == 2) then
            k_beg = idwbuff(1)%beg; k_end = idwbuff(1)%end; l_beg = 0; l_end = p
        else
            k_beg = idwbuff(1)%beg; k_end = idwbuff(1)%end; l_beg = idwbuff(2)%beg; l_end = idwbuff(2)%end
        end if

        $:GPU_PARALLEL_LOOP(private='[l, k, bc_code]', collapse=2)
        do l = l_beg, l_end
            do k = k_beg, k_end
                if (bc_dir == 1) then
                    bc_code = int(bc_type_edge%sf(0, k, l))
                else if (bc_dir == 2) then
                    bc_code = int(bc_type_edge%sf(k, 0, l))
                else
                    bc_code = int(bc_type_edge%sf(k, l, 0))
                end if

                select case (bc_code)
                case (BC_PERIODIC)
                    call s_F_igr_periodic(jac_sf, bc_dir, bc_loc, k, l)
                case (BC_REFLECTIVE)
                    call s_F_igr_reflective(jac_sf, bc_dir, bc_loc, k, l)
                case default
                    call s_F_igr_ghost_cell_extrapolation(jac_sf, bc_dir, bc_loc, k, l)
                end select
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_populate_F_igr_bc_direction

    !> The purpose of this subroutine is to populate the buffers of the grid variables, which are constituted of the cell- boundary
    !! locations and cell-width distributions, based on the boundary conditions.
    subroutine s_populate_grid_variables_buffers

#ifdef MFC_SIMULATION
        ! Simulation allocates the cell-boundary arrays with buff_size ghost layers, so the
        ! ghost extrapolation extends over the full buffer. In post-process the module-level
        ! offset_x/y/z (which the x_cb/y_cb/z_cb allocations are sized to) are used instead.
        type(int_bounds_info) :: offset_x, offset_y, offset_z
        offset_x%beg = buff_size; offset_x%end = buff_size
        offset_y%beg = buff_size; offset_y%end = buff_size
        offset_z%beg = buff_size; offset_z%end = buff_size

#ifdef MFC_MPI
        ! Populate global domain boundaries with stretched grids
        call s_mpi_allreduce_min(x_cb(-1), glb_bounds(1)%beg)
        call s_mpi_allreduce_max(x_cb(m), glb_bounds(1)%end)

        if (n > 0) then
            call s_mpi_allreduce_min(y_cb(-1), glb_bounds(2)%beg)
            call s_mpi_allreduce_max(y_cb(n), glb_bounds(2)%end)
            if (p > 0) then
                call s_mpi_allreduce_min(z_cb(-1), glb_bounds(3)%beg)
                call s_mpi_allreduce_max(z_cb(p), glb_bounds(3)%end)
            end if
        end if
#else
        glb_bounds(1)%beg = x_cb(-1); glb_bounds(1)%end = x_cb(m)
        if (n > 0) then
            glb_bounds(2)%beg = y_cb(-1); glb_bounds(2)%end = y_cb(n)
            if (p > 0) then
                glb_bounds(3)%beg = z_cb(-1); glb_bounds(3)%end = z_cb(p)
            end if
        end if
#endif
        $:GPU_UPDATE(device='[glb_bounds]')
#endif

#ifndef MFC_PRE_PROCESS
        call s_populate_grid_bc_direction(1, -1, bc_x, offset_x)
        call s_populate_grid_bc_direction(1, 1, bc_x, offset_x)

        if (n == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
            call s_populate_grid_bc_direction(2, -1, bc_y, offset_y)
            call s_populate_grid_bc_direction(2, 1, bc_y, offset_y)
        #:endif

        if (p == 0) return

        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            call s_populate_grid_bc_direction(3, -1, bc_z, offset_z)
            call s_populate_grid_bc_direction(3, 1, bc_z, offset_z)
        #:endif
#endif

    end subroutine s_populate_grid_variables_buffers

#ifndef MFC_PRE_PROCESS
    !> Populate grid variable buffers (cell widths and centers) for one direction and location.
    subroutine s_populate_grid_bc_direction(bc_dir, bc_loc, bc_bounds, offset_dir)

        integer, intent(in)               :: bc_dir, bc_loc
        type(int_bounds_info), intent(in) :: bc_bounds, offset_dir
        integer                           :: bc_edge

        if (bc_loc == -1) then
            bc_edge = bc_bounds%beg
        else
            bc_edge = bc_bounds%end
        end if

        if (bc_edge >= 0) then
            call s_mpi_sendrecv_grid_variables_buffers(bc_dir, bc_loc, offset_dir)
            return
        end if

        select case (bc_edge)
        case (BC_PERIODIC)
            call s_grid_periodic_bc(bc_dir, bc_loc, offset_dir)
        case (BC_REFLECTIVE)
            call s_grid_reflective_bc(bc_dir, bc_loc, offset_dir)
        case (BC_AXIS)
            if (bc_dir == 2) call s_grid_axis_bc(bc_loc, offset_dir)
        case default
            call s_grid_ghost_cell_extrapolation_bc(bc_dir, bc_loc, offset_dir)
        end select

    end subroutine s_populate_grid_bc_direction
#endif

    !> Deallocate boundary condition buffer arrays allocated during module initialization.
    subroutine s_finalize_boundary_common_module()

        if (bc_io) then
            @:DEALLOCATE(bc_buffers(1, 1)%sf)
            @:DEALLOCATE(bc_buffers(1, 2)%sf)
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                if (n > 0) then
                    @:DEALLOCATE(bc_buffers(2, 1)%sf)
                    @:DEALLOCATE(bc_buffers(2, 2)%sf)
                    #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                        if (p > 0) then
                            @:DEALLOCATE(bc_buffers(3, 1)%sf)
                            @:DEALLOCATE(bc_buffers(3, 2)%sf)
                        end if
                    #:endif
                end if
            #:endif
        end if
        @:DEALLOCATE(bc_buffers)

    end subroutine s_finalize_boundary_common_module

end module m_boundary_common
