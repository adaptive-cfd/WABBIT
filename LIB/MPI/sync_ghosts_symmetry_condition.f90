subroutine sync_ghosts_symmetry_condition( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    integer :: k, hvy_id, lgt_id, ineighbor, ijk(1:2,1:3), ix, iy ,iz, Bs(1:3), g, a, j, neqn
    integer(kind=2) :: n_domain(1:3)
    integer(kind=ik) :: neighborhoods(1:4)

    if (.not.allocated(params%symmetry_vector_component)) then
        call abort(11351, "WABBIT internal error: Required arrays not allocated?!")
    endif

    if (.not. ghost_nodes_module_ready) call abort(210129, "not initialized?!")

    Bs = params%Bs
    g = params%n_ghosts
    neqn = params%n_eqn

    ! First step: faces (3D) or edges (2D)
    ! This enables us to mirror corners from the ghost nodes layer, and not from
    ! the blocks interior. Then, their treatment is identical for symmetric/periodic and symmetric/symmetric cases

    !  ____  ____
    ! |___ \|  _ \
    !   __) | | | |
    !  / __/| |_| |
    ! |_____|____/
    !
    if (params%dim == 2) then
        ! loop is parallel (hvy_n), this is a local operation without communication overhead
        do k = 1, hvy_n
            hvy_id = hvy_active(k)
            call hvy_id_to_lgt_id(lgt_id, hvy_id, params%rank, params%number_blocks)

            call get_adjacent_boundary_surface_normal( lgt_block(lgt_id, 1:lgt_block(lgt_id,params%max_treelevel+IDX_MESH_LVL)), &
            params%domain_size, params%Bs, params%dim, n_domain )

            ! this block is not concerned (its an interior block)
            if (n_domain(1) == 0 .and. n_domain(2) == 0) cycle


            if (n_domain(1) == +1 .and. params%symmetry_BC(1)) then
                ineighbor = 1
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,inverse_neighbor(ineighbor,2)) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ! bounds of recver patch (the face we will fill)
                ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, EXCLUDE_REDUNDANT, RECVER)

                do j = 1, neqn
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                    hvy_block( Bs(1)+g-1:Bs(1):-1, ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="x") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

            if (n_domain(1) == -1 .and. params%symmetry_BC(1)) then
                ineighbor = 3
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,inverse_neighbor(ineighbor,2)) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, EXCLUDE_REDUNDANT, RECVER)

                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( g+2+g-1:g+2:-1, ijk(1,2):ijk(2,2) , ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="x") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

            if (n_domain(2) == +1 .and. params%symmetry_BC(2)) then
                ineighbor = 4
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,inverse_neighbor(ineighbor,2)) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, EXCLUDE_REDUNDANT, RECVER)

                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( ijk(1,1):ijk(2,1), Bs(2)+g-1:Bs(2):-1 , ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="y") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

            if (n_domain(2) == -1 .and. params%symmetry_BC(2)) then
                ineighbor = 2
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,inverse_neighbor(ineighbor,2)) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, EXCLUDE_REDUNDANT, RECVER)

                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( ijk(1,1):ijk(2,1), g+2+g-1:g+2:-1 , ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="y") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

            if (n_domain(1) == +1 .and. params%symmetry_BC(1)) then
                ineighbor = 6
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,inverse_neighbor(ineighbor,2)) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, EXCLUDE_REDUNDANT, RECVER)

                ! this means we mirror in x direction (i.e. w.r.t. y-axis)
                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( Bs(1)+g-1:Bs(1):-1, ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="x") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo

                ineighbor = 5
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,inverse_neighbor(ineighbor,2)) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, EXCLUDE_REDUNDANT, RECVER)

                ! this means we mirror in x direction (i.e. w.r.t. y-axis)
                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( Bs(1)+g-1:Bs(1):-1, ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="x") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

            if (n_domain(1) == -1 .and. params%symmetry_BC(1)) then
                ineighbor = 7
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,inverse_neighbor(ineighbor,2)) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, EXCLUDE_REDUNDANT, RECVER)

                ! this means we mirror in x direction (i.e. w.r.t. y-axis)
                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( g+2+g-1:g+2:-1, ijk(1,2):ijk(2,2) , ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="x") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo

                ineighbor = 8
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,inverse_neighbor(ineighbor,2)) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, EXCLUDE_REDUNDANT, RECVER)

                ! this means we mirror in x direction (i.e. w.r.t. y-axis)
                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( g+2+g-1:g+2:-1, ijk(1,2):ijk(2,2) , ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="x") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

            if (n_domain(2) == +1 .and. params%symmetry_BC(2)) then
                ineighbor = 8
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,inverse_neighbor(ineighbor,2)) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, EXCLUDE_REDUNDANT, RECVER)

                ! this means we mirror in y direction (i.e. w.r.t. x-axis)
                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( ijk(1,1):ijk(2,1), Bs(2)+g-1:Bs(2):-1 , ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="y") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo

                ineighbor = 6
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,inverse_neighbor(ineighbor,2)) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, EXCLUDE_REDUNDANT, RECVER)

                ! this means we mirror in y direction (i.e. w.r.t. x-axis)
                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( ijk(1,1):ijk(2,1), Bs(2)+g-1:Bs(2):-1 , ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="y") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

            if (n_domain(2) == -1 .and. params%symmetry_BC(2)) then
                ineighbor = 7
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,inverse_neighbor(ineighbor,2)) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, EXCLUDE_REDUNDANT, RECVER)

                ! this means we mirror in y direction (i.e. w.r.t. x-axis)
                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( ijk(1,1):ijk(2,1), g+2+g-1:g+2:-1 , ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="y") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo


                ineighbor = 5
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,inverse_neighbor(ineighbor,2)) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, EXCLUDE_REDUNDANT, RECVER)

                ! this means we mirror in y direction (i.e. w.r.t. x-axis)
                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( ijk(1,1):ijk(2,1), g+2+g-1:g+2:-1 , ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="y") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, ineighbor, 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

        enddo

    !  _____ ____
    ! |___ /|  _ \
    !   |_ \| | | |
    !  ___) | |_| |
    ! |____/|____/
    else
        ! loop is parallel (hvy_n), this is a local operation without communication overhead
        do k = 1, hvy_n
            hvy_id = hvy_active(k)
            call hvy_id_to_lgt_id(lgt_id, hvy_id, params%rank, params%number_blocks)

            call get_adjacent_boundary_surface_normal( lgt_block(lgt_id, 1:lgt_block(lgt_id,params%max_treelevel+IDX_MESH_LVL)), &
            params%domain_size, params%Bs, params%dim, n_domain )

            ! this block is not concerned (its an interior block)
            if (n_domain(1) == 0 .and. n_domain(2) == 0 .and. n_domain(3) == 0) cycle

            !-------------------------------------------------------------------
            !---faces---
            !-------------------------------------------------------------------
            if (n_domain(1) == +1 .and. params%symmetry_BC(1)) then
                ineighbor = 3
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ! bounds of recver patch (the face we will fill)
                ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)

                do j = 1, neqn
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                    hvy_block( Bs(1)+g-1:Bs(1):-1, ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="x") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

            if (n_domain(1) == -1 .and. params%symmetry_BC(1)) then
                ineighbor = 5
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)

                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( g+2+g-1:g+2:-1, ijk(1,2):ijk(2,2) , ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="x") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

            if (n_domain(2) == +1 .and. params%symmetry_BC(2)) then
                ineighbor = 4
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)

                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( ijk(1,1):ijk(2,1), Bs(2)+g-1:Bs(2):-1 , ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="y") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

            if (n_domain(2) == -1 .and. params%symmetry_BC(2)) then
                ineighbor = 2
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)

                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( ijk(1,1):ijk(2,1), g+2+g-1:g+2:-1 , ijk(1,3):ijk(2,3), :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="y") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

            if (n_domain(3) == +1 .and. params%symmetry_BC(3)) then
                ineighbor = 1
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)

                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), Bs(3)+g-1:Bs(3):-1, :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="z") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

            if (n_domain(3) == -1 .and. params%symmetry_BC(3)) then
                ineighbor = 6
                ! if the neighbor is not -1 (no neighbor), peridic sync' may not be disabled and overwrite the redundant point.
                if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry BC error: go play outside. Take the kids. They like it.")
                ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)

                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), g+2+g-1:g+2:-1,  :, hvy_id)

                do j = 1, neqn
                    ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                    if (params%symmetry_vector_component(j)=="z") then
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                        -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                        ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                        hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                    endif
                enddo
            endif

            !-------------------------------------------------------------------
            !--edges---
            !-------------------------------------------------------------------
            if (n_domain(1) == +1 .and. params%symmetry_BC(1)) then
                neighborhoods = (/8,17,12,15/)
                do a = 1, 4
                    ineighbor = neighborhoods(a)
                    if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry edge BC error: go play outside. Take the kids. They like it.")
                    ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                    hvy_block( Bs(1)+g-1:Bs(1):-1, ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id)

                    do j = 1, neqn
                        ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                        if (params%symmetry_vector_component(j)=="x") then
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                            -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                            ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                        endif
                    enddo
                enddo
            endif

            if (n_domain(1) == -1 .and. params%symmetry_BC(1)) then
                neighborhoods = (/10,14,16,18/)
                do a = 1, 4
                    ineighbor = neighborhoods(a)
                    if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry edge BC error: go play outside. Take the kids. They like it.")
                    ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                    hvy_block( g+2+g-1:g+2:-1, ijk(1,2):ijk(2,2) , ijk(1,3):ijk(2,3), :, hvy_id)

                    do j = 1, neqn
                        ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                        if (params%symmetry_vector_component(j)=="x") then
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                            -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                            ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                        endif
                    enddo
                enddo
            endif

            if (n_domain(2) == +1 .and. params%symmetry_BC(2)) then
                neighborhoods = (/9,13,17,18/)
                do a = 1, 4
                    ineighbor = neighborhoods(a)
                    if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry edge BC error: go play outside. Take the kids. They like it.")
                    ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                    hvy_block( ijk(1,1):ijk(2,1), Bs(2)+g-1:Bs(2):-1 , ijk(1,3):ijk(2,3), :, hvy_id)

                    do j = 1, neqn
                        ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                        if (params%symmetry_vector_component(j)=="y") then
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                            -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                            ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                        endif
                    enddo
                enddo
            endif

            if (n_domain(2) == -1 .and. params%symmetry_BC(2)) then
                neighborhoods = (/7,15,11,16/)
                do a = 1, 4
                    ineighbor = neighborhoods(a)
                    if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry edge BC error: go play outside. Take the kids. They like it.")
                    ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                    hvy_block( ijk(1,1):ijk(2,1), g+2+g-1:g+2:-1 , ijk(1,3):ijk(2,3), :, hvy_id)

                    do j = 1, neqn
                        ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                        if (params%symmetry_vector_component(j)=="y") then
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                            -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                            ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                        endif
                    enddo
                enddo
            endif

            if (n_domain(3) == +1 .and. params%symmetry_BC(3)) then
                neighborhoods = (/7,8,9,10/)
                do a = 1, 4
                    ineighbor = neighborhoods(a)
                    if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry edge BC error: go play outside. Take the kids. They like it.")
                    ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), Bs(3)+g-1:Bs(3):-1, :, hvy_id)

                    do j = 1, neqn
                        ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                        if (params%symmetry_vector_component(j)=="z") then
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                            -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                            ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                        endif
                    enddo
                enddo
            endif

            if (n_domain(3) == -1 .and. params%symmetry_BC(3)) then
                neighborhoods = (/11,12,13,14/)
                do a = 1, 4
                    ineighbor = neighborhoods(a)
                    if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry edge BC error: go play outside. Take the kids. They like it.")
                    ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), g+2+g-1:g+2:-1,  :, hvy_id)

                    do j = 1, neqn
                        ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                        if (params%symmetry_vector_component(j)=="z") then
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                            -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                            ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                        endif
                    enddo
                enddo
            endif

            !-------------------------------------------------------------------
            !--corners---
            !-------------------------------------------------------------------
            if (n_domain(1) == +1 .and. params%symmetry_BC(1)) then
                neighborhoods = (/19,20,23,24/)
                do a = 1, 4
                    ineighbor = neighborhoods(a)
                    if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry edge BC error: go play outside. Take the kids. They like it.")
                    ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                    hvy_block( Bs(1)+g-1:Bs(1):-1, ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id)

                    do j = 1, neqn
                        ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                        if (params%symmetry_vector_component(j)=="x") then
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                            -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                            ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                        endif
                    enddo
                enddo
            endif

            if (n_domain(1) == -1 .and. params%symmetry_BC(1)) then
                neighborhoods = (/21,22,25,26/)
                do a = 1, 4
                    ineighbor = neighborhoods(a)
                    if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry edge BC error: go play outside. Take the kids. They like it.")
                    ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                    hvy_block( g+2+g-1:g+2:-1, ijk(1,2):ijk(2,2) , ijk(1,3):ijk(2,3), :, hvy_id)

                    do j = 1, neqn
                        ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                        if (params%symmetry_vector_component(j)=="x") then
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                            -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                            ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                        endif
                    enddo
                enddo
            endif

            if (n_domain(2) == +1 .and. params%symmetry_BC(2)) then
                neighborhoods = (/20,21,24,25/)
                do a = 1, 4
                    ineighbor = neighborhoods(a)
                    if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry edge BC error: go play outside. Take the kids. They like it.")
                    ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                    hvy_block( ijk(1,1):ijk(2,1), Bs(2)+g-1:Bs(2):-1 , ijk(1,3):ijk(2,3), :, hvy_id)

                    do j = 1, neqn
                        ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                        if (params%symmetry_vector_component(j)=="y") then
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                            -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                            ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                        endif
                    enddo
                enddo
            endif

            if (n_domain(2) == -1 .and. params%symmetry_BC(2)) then
                neighborhoods = (/22,19,26,23/)
                do a = 1, 4
                    ineighbor = neighborhoods(a)
                    if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry edge BC error: go play outside. Take the kids. They like it.")
                    ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                    hvy_block( ijk(1,1):ijk(2,1), g+2+g-1:g+2:-1 , ijk(1,3):ijk(2,3), :, hvy_id)

                    do j = 1, neqn
                        ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                        if (params%symmetry_vector_component(j)=="y") then
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                            -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                            ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                        endif
                    enddo
                enddo
            endif

            if (n_domain(3) == +1 .and. params%symmetry_BC(3)) then
                neighborhoods = (/19,20,21,22/)
                do a = 1, 4
                    ineighbor = neighborhoods(a)
                    if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry edge BC error: go play outside. Take the kids. They like it.")
                    ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), Bs(3)+g-1:Bs(3):-1, :, hvy_id)

                    do j = 1, neqn
                        ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                        if (params%symmetry_vector_component(j)=="z") then
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                            -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                            ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                        endif
                    enddo
                enddo
            endif

            if (n_domain(3) == -1 .and. params%symmetry_BC(3)) then
                neighborhoods = (/26,23,24,25/)
                do a = 1, 4
                    ineighbor = neighborhoods(a)
                    if (hvy_neighbor(hvy_id,ineighbor) /= -1) call abort(ineighbor,"Symmetry edge BC error: go play outside. Take the kids. They like it.")
                    ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, EXCLUDE_REDUNDANT, RECVER)
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_id) = &
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), g+2+g-1:g+2:-1,  :, hvy_id)

                    do j = 1, neqn
                        ! if this is the velocity component normal to the symmetry plane, invert sign and enforce Dirichlet BC on the line of symmetry
                        if (params%symmetry_vector_component(j)=="z") then
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = &
                            -1.0_rk*hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id)

                            ijk = ijkGhosts(1:2, 1:3, inverse_neighbor(ineighbor,3), 0, ONLY_REDUNDANT, RECVER)
                            hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), j, hvy_id) = 0.0_rk
                        endif
                    enddo
                enddo
            endif
        enddo
    endif
end subroutine
