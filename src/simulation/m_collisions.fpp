!>
!! @file
!! @brief Contains module m_collisions

#:include 'macros.fpp'

!> @brief Geometric sphere/circle contacts and soft-sphere force evaluation
module m_collisions

    use m_derived_types      !< Definitions of the derived types
    use m_global_parameters  !< Definitions of the global parameters
    use m_helper, only: s_cross_product
    use m_helper_basic       !< Functions to compare floating point numbers
    use m_constants
    use m_mpi_proxy

    implicit none

    private; public :: s_apply_collision_forces, s_initialize_collisions_module, s_finalize_collisions_module, &
        & f_local_rank_owns_location, f_neighborhood_ranks_own_location, ib_gbl_idx_lookup
    ! Wall overlap distances for computing contacts
    real(wp), allocatable, dimension(:,:) :: wall_overlap_distances
    real(wp)                              :: spring_stiffness, damping_parameter
    $:GPU_DECLARE(create='[spring_stiffness, damping_parameter]')
    $:GPU_DECLARE(create='[wall_overlap_distances]')

    integer, dimension(:), allocatable :: ib_gbl_idx_lookup
    $:GPU_DECLARE(create='[ib_gbl_idx_lookup]')

contains

    subroutine s_initialize_collisions_module()

        real(wp) :: e

        e = coefficient_of_restitution
        damping_parameter = -2._wp*log(e)/collision_time
        spring_stiffness = (pi**2 + log(e)**2)/(collision_time**2)
        $:GPU_UPDATE(device='[damping_parameter, spring_stiffness]')

        @:ALLOCATE(wall_overlap_distances(num_local_ibs_max*27, 6))

        wall_overlap_distances = 0
        $:GPU_UPDATE(device='[wall_overlap_distances]')
        $:GPU_UPDATE(device='[ib_coefficient_of_friction]')

    end subroutine s_initialize_collisions_module

    subroutine s_apply_collision_forces(forces, torques)

        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques

        if (collision_model == 0) return

        ! A rank may legitimately have no IBs in its current neighborhood.
        ! In that case forces/torques are zero-size automatic arrays, whose
        ! OpenACC device address is null. Do not launch collision kernels that
        ! reference them.
        if (num_ibs == 0) return

        ! Distances used in the particle-wall force calculation.
        call s_detect_wall_collisions()

        select case (collision_model)
        case (1)  ! soft sphere model
            call s_apply_wall_collision_forces_soft_sphere(forces, torques)
            call s_apply_ib_collision_forces_soft_sphere(forces, torques)
        end select

    end subroutine s_apply_collision_forces

    !> @brief applies collision forces to IBs assuming a soft-sphere collision model (all IBs are circles or spheres)
    subroutine s_apply_ib_collision_forces_soft_sphere(forces, torques)

        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques
        integer :: pid1, pid2, l
        real(wp) :: overlap_distance, distance, period
        real(wp), dimension(3) :: normal_vector, centroid_1, centroid_2
        real(wp), dimension(3) :: normal_velocity, tangential_vector, normal_force, tangential_force, torque, radial_vector, &
             & rotation_velocity, vel1, vel2
        real(wp) :: k, eta, effective_mass  ! the spring stiffness and damping coefficient and mass of a specific interaction

        ! Geometric broad phase over this rank's IB neighborhood. Every global
        ! pair is evaluated once, on the owner of its smaller global ID, even
        ! when marker overwrite hides a particle in a multiple contact.
        ! This allocation-free local pair scan is independent of fluid resolution.

        $:GPU_PARALLEL_LOOP(private='[pid1, pid2, l, centroid_1, centroid_2, normal_vector, overlap_distance, distance, period, &
                            & effective_mass, k, eta, normal_velocity, tangential_vector, normal_force, tangential_force, torque, &
                            & radial_vector, rotation_velocity, vel1, vel2]', copy='[forces, torques]', collapse=2)
        do pid1 = 1, num_ibs
            do pid2 = 1, num_ibs
                if (patch_ib(pid1)%gbl_patch_id >= patch_ib(pid2)%gbl_patch_id) cycle
                centroid_1 = [patch_ib(pid1)%x_centroid, patch_ib(pid1)%y_centroid, patch_ib(pid1)%z_centroid]
                centroid_2 = [patch_ib(pid2)%x_centroid, patch_ib(pid2)%y_centroid, patch_ib(pid2)%z_centroid]
                ! Wrap the owner position too: RK stages can cross a periodic face before the end-of-step ownership handoff.
                #:for X, ID in [('x', 1), ('y', 2), ('z', 3)]
                    if (num_dims >= ${ID}$ .and. ib_bc_${X}$%beg == BC_PERIODIC) then
                        period = glb_bounds(${ID}$)%end - glb_bounds(${ID}$)%beg
                        centroid_1(${ID}$) = glb_bounds(${ID}$)%beg + modulo(centroid_1(${ID}$) - glb_bounds(${ID}$)%beg, period)
                    end if
                #:endfor
                if (.not. f_local_rank_owns_location(centroid_1)) cycle
                normal_vector = centroid_2 - centroid_1
                #:for X, ID in [('x', 1), ('y', 2), ('z', 3)]
                    if (num_dims >= ${ID}$ .and. ib_bc_${X}$%beg == BC_PERIODIC) then
                        period = glb_bounds(${ID}$)%end - glb_bounds(${ID}$)%beg
                        normal_vector(${ID}$) = normal_vector(${ID}$) - anint(normal_vector(${ID}$)/period)*period
                    end if
                #:endfor
                if (num_dims == 2) normal_vector(3) = 0._wp
                distance = norm2(normal_vector)
                overlap_distance = patch_ib(pid1)%radius + patch_ib(pid2)%radius - distance
                if (overlap_distance > 0._wp) then  ! if the two patches are close enough to collide
                    ! A coincident pair has no geometric normal. Use relative motion
                    ! (or a deterministic axis at rest) rather than dividing by zero.
                    if (distance > 0._wp) then
                        normal_vector = normal_vector/distance
                    else
                        normal_vector = patch_ib(pid1)%vel - patch_ib(pid2)%vel
                        if (norm2(normal_vector) > 0._wp) then
                            normal_vector = normal_vector/norm2(normal_vector)
                        else
                            normal_vector = [1._wp, 0._wp, 0._wp]
                        end if
                    end if
                    ! compute constants of the collision
                    effective_mass = 1.0_wp/((1.0_wp/patch_ib(pid1)%mass) + (1._wp/(patch_ib(pid2)%mass)))
                    k = spring_stiffness*effective_mass
                    eta = damping_parameter*effective_mass

                    ! Contact-point velocities include translation and rotation.
                    radial_vector = normal_vector*(patch_ib(pid1)%radius - 0.5_wp*overlap_distance)
                    call s_cross_product(patch_ib(pid1)%angular_vel, radial_vector, rotation_velocity)
                    vel1 = patch_ib(pid1)%vel + rotation_velocity
                    radial_vector = normal_vector*(-1.0_wp)*(patch_ib(pid2)%radius - 0.5_wp*overlap_distance)
                    call s_cross_product(patch_ib(pid2)%angular_vel, radial_vector, rotation_velocity)
                    vel2 = patch_ib(pid2)%vel + rotation_velocity

                    normal_velocity = dot_product(vel1 - vel2, normal_vector)*normal_vector
                    tangential_vector = (vel1 - vel2) - normal_velocity
                    if (.not. f_approx_equal(norm2(tangential_vector), &
                        & 0._wp)) tangential_vector = tangential_vector/norm2(tangential_vector)

                    ! compute force and torque
                    normal_force = -k*overlap_distance*normal_vector - eta*normal_velocity
                    tangential_force = -ib_coefficient_of_friction*norm2(normal_force)*tangential_vector
                    call s_cross_product(normal_vector*patch_ib(pid1)%radius, tangential_force, torque)

                    do l = 1, 3
                        ! update the first IB
                        $:GPU_ATOMIC(atomic='update')
                        forces(pid1, l) = forces(pid1, l) + (normal_force(l) + tangential_force(l))
                        $:GPU_ATOMIC(atomic='update')
                        torques(pid1, l) = torques(pid1, l) + torque(l)

                        ! apply equal and opposite force/torque to second IB
                        $:GPU_ATOMIC(atomic='update')
                        forces(pid2, l) = forces(pid2, l) - (normal_force(l) + tangential_force(l))
                        $:GPU_ATOMIC(atomic='update')
                        torques(pid2, l) = torques(pid2, l) + torque(l)*patch_ib(pid2)%radius/patch_ib(pid1)%radius
                    end do
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_ib_collision_forces_soft_sphere

    !> @brief applies collision forces to IBs assuming a soft-sphere collision model (all IBs are circles or spheres)
    subroutine s_apply_wall_collision_forces_soft_sphere(forces, torques)

        real(wp), dimension(num_ibs, 3), intent(inout) :: forces, torques
        integer :: patch_id, i, l
        real(wp), dimension(3) :: normal_force, tangential_force, normal_vector, normal_velocity, tangential_vector, &
             & collision_location, torque, radial_vector, rotation_velocity, velocity
        real(wp) :: k, eta  ! the spring stiffness and damping coefficient for a specific IB

        $:GPU_PARALLEL_LOOP(private='[patch_id, i, l, collision_location, normal_vector, k, eta, normal_velocity, &
                            & tangential_vector, normal_force, tangential_force, torque, radial_vector, rotation_velocity, &
                            & velocity]', copy='[forces, torques]', collapse=2)
        do patch_id = 1, num_ibs
            do i = 1, num_dims*2
                ! only compute force contributions if there was an overlap
                if (f_approx_equal(wall_overlap_distances(patch_id, i), 0._wp)) cycle

                select case (i)
                case (1)  ! x domain left
                    normal_vector = [-1._wp, 0._wp, 0._wp]
                case (2)  ! x domain right
                    normal_vector = [1._wp, 0._wp, 0._wp]
                case (3)  ! y domain bottom
                    normal_vector = [0._wp, -1._wp, 0._wp]
                case (4)  ! y domain top
                    normal_vector = [0._wp, 1._wp, 0._wp]
                case (5)  ! z domain back
                    normal_vector = [0._wp, 0._wp, -1._wp]
                case (6)  ! z domain front
                    normal_vector = [0._wp, 0._wp, 1._wp]
                end select

                ! ensure the local rank owns that collision before proceeding
                collision_location = [patch_ib(patch_id)%x_centroid, patch_ib(patch_id)%y_centroid, 0._wp]
                if (num_dims == 3) collision_location(3) = patch_ib(patch_id)%z_centroid
                if (f_local_rank_owns_location(collision_location)) then
                    k = spring_stiffness*patch_ib(patch_id)%mass
                    eta = damping_parameter*patch_ib(patch_id)%mass

                    ! get the vector that points from the centroid to the point of collision
                    radial_vector = normal_vector*(patch_ib(patch_id)%radius - wall_overlap_distances(patch_id, i))
                    ! convert the angular velocity to linear velocity
                    call s_cross_product(patch_ib(patch_id)%angular_vel, radial_vector, rotation_velocity)
                    velocity = patch_ib(patch_id)%vel + rotation_velocity

                    ! standard soft-sphere collision  with the wall
                    normal_velocity = dot_product(velocity, normal_vector)*normal_vector
                    tangential_vector = velocity - normal_velocity
                    if (.not. f_approx_equal(norm2(tangential_vector), &
                        & 0._wp)) tangential_vector = tangential_vector/norm2(tangential_vector)
                    normal_force = -k*wall_overlap_distances(patch_id, i)*normal_vector - eta*normal_velocity
                    tangential_force = -ib_coefficient_of_friction*norm2(normal_force)*tangential_vector
                    call s_cross_product(normal_vector*patch_ib(patch_id)%radius, tangential_force, torque)

                    do l = 1, 3
                        $:GPU_ATOMIC(atomic='update')
                        forces(patch_id, l) = forces(patch_id, l) + (normal_force(l) + tangential_force(l))
                        $:GPU_ATOMIC(atomic='update')
                        torques(patch_id, l) = torques(patch_id, l) + torque(l)
                    end do
                end if
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_apply_wall_collision_forces_soft_sphere

    !> @brief uses boundary conditions and particle locations to check for wall conditions
    subroutine s_detect_wall_collisions()

        integer  :: patch_id
        real(wp) :: edge_location, overlap_distance

        ! iterate over all ghost points to detect the one that is most-overlapping in each direction

        $:GPU_PARALLEL_LOOP(private='[patch_id, edge_location, overlap_distance]')
        do patch_id = 1, num_ibs
            #:for X, DIR, IDX in [('x', 1, 1), ('y', 2, 3), ('z', 3, 5)]
                ! check if the boundaries are either of the two conditions we should compute collisions with
                if (ib_bc_${X}$%beg == BC_SLIP_WALL .or. ib_bc_${X}$%beg == BC_NO_SLIP_WALL) then
                    ! get the location of the true IB surface towards the domain boundary
                    edge_location = patch_ib(patch_id)%${X}$_centroid - patch_ib(patch_id)%radius
                    ! check if that edge actually extends out of the comutational domain
                    if (edge_location < glb_bounds(${DIR}$)%beg) then
                        ! the distance that the IB extends out of the domain
                        overlap_distance = glb_bounds(${DIR}$)%beg - edge_location
                    else
                        overlap_distance = 0._wp
                    end if
                    wall_overlap_distances(patch_id, ${IDX}$) = overlap_distance
                end if

                if (ib_bc_${X}$%end == BC_SLIP_WALL .or. ib_bc_${X}$%end == BC_NO_SLIP_WALL) then
                    edge_location = patch_ib(patch_id)%${X}$_centroid + patch_ib(patch_id)%radius
                    if (edge_location > glb_bounds(${DIR}$)%end) then
                        overlap_distance = edge_location - glb_bounds(${DIR}$)%end
                    else
                        overlap_distance = 0._wp
                    end if
                    wall_overlap_distances(patch_id, ${IDX}$ + 1) = overlap_distance
                end if
            #:endfor
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_detect_wall_collisions

    !> @brief Check whether this rank owns a contact location.
    function f_local_rank_owns_location(location) result(owns_collision)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), dimension(3), intent(in) :: location
        logical                            :: owns_collision
        real(wp), dimension(3)             :: projected_location

        owns_collision = .true.

#ifdef MFC_MPI
        if (num_procs > 1) then
            projected_location(:) = location(:)

            ! catch the edge case where th collision lies just outside the computational domain
            #:for X, ID, DIM in [('x', 1, 'm'), ('y', 2, 'n'), ('z', 3, 'p')]
                if (num_dims >= ${ID}$) then
                    if (ib_bc_${X}$%beg /= BC_PERIODIC) then
                        ! if it is outside the domain in one direction, project it somewhere inside so at least one rank owns it
                        if (location(${ID}$) < glb_bounds(${ID}$)%beg) then
                            projected_location(${ID}$) = glb_bounds(${ID}$)%beg
                        else if (glb_bounds(${ID}$)%end < location(${ID}$)) then
                            projected_location(${ID}$) = glb_bounds(${ID}$)%end - 1.0e-10_wp
                        end if
                    end if
                    owns_collision = owns_collision .and. ${X}$_cb(-1) <= projected_location(${ID}$) &
                        & .and. projected_location(${ID}$) < ${X}$_cb(${DIM}$)
                end if
            #:endfor
        end if
#endif

    end function f_local_rank_owns_location

    !> @brief function checks if this local MPI processor owns this specific collision
    function f_neighborhood_ranks_own_location(location) result(owns_collision)

        real(wp), dimension(3), intent(in) :: location
        logical                            :: owns_collision, periodic_owner
        real(wp)                           :: temp_neighbor_domain
        integer                            :: i

        owns_collision = .true.

#ifdef MFC_MPI
        if (num_procs > 2) then
            ! catch the edge case where th collision lies just outside the computational domain
            owns_collision = .true.
            #:for X, ID in [('x', 1), ('y', 2,), ('z', 3,)]
                if (num_dims >= ${ID}$) then
                    if (ib_bc_${X}$%beg == BC_PERIODIC .and. neighbor_domain_${X}$%beg >= neighbor_domain_${X}$%end) then
                        ! project right side to the left
                        temp_neighbor_domain = neighbor_domain_${X}$%end + (glb_bounds(${ID}$)%end - glb_bounds(${ID}$)%beg)
                        periodic_owner = neighbor_domain_${X}$%beg <= location(${ID}$) .and. location(${ID}$) < temp_neighbor_domain
                        ! project the left side to the right
                        temp_neighbor_domain = neighbor_domain_${X}$%beg - (glb_bounds(${ID}$)%end - glb_bounds(${ID}$)%beg)
                        periodic_owner = periodic_owner .or. (temp_neighbor_domain <= location(${ID}$) .and. location(${ID}$) &
                                                              & < neighbor_domain_${X}$%end)

                        owns_collision = owns_collision .and. periodic_owner
                    else
                        owns_collision = owns_collision .and. neighbor_domain_${X}$%beg <= location(${ID}$) .and. location(${ID}$) &
                            & < neighbor_domain_${X}$%end
                    end if
                end if
            #:endfor
        end if
#endif

    end function f_neighborhood_ranks_own_location

    subroutine s_finalize_collisions_module()

        @:DEALLOCATE(wall_overlap_distances)

    end subroutine s_finalize_collisions_module

end module m_collisions
