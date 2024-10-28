subroutine componentWiseNorm_tree(params, hvy_block, tree_ID, which_norm, norm)
    implicit none

    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)             !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID
    character(len=*), intent(in)        :: which_norm                           !> which norm to use ? "L2", "Linfty"
    real(kind=rk), intent(inout)        :: norm(:)                              !> the computed norm for each component of the vector
    real(kind=rk)                       :: x0(1:3), dx(1:3)
    integer(kind=ik) :: k, hvy_id, n_eqn, Bs(1:3), g(1:3), p, mpierr, lgt_id, D

    ! note: if norm and hvy_block components are of different size, we use the smaller one.
    n_eqn = min( size(norm, 1), size(hvy_block,4) )
    Bs = params%Bs
    g = 0
    g(1:params%dim) = params%g
    D = params%dim

    select case (which_norm)
    case ("L2", "H1")
        norm = 0.0_rk
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

            ! only sum up over leafs, we sadly don't have access to function block_is_leaf
            if (any(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)) /= -1)) cycle

            call get_block_spacing_origin_b( get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2)), params%domain_size, &
                params%Bs, x0, dx, dim=params%dim, level=lgt_block(lgt_id, IDX_MESH_LVL), max_level=params%Jmax)

            ! JB ToDo: This for now goes against the idea of physics types, if we decide type by type the norm, then this has to move into the physics module
            if (params%physics_type == "ACM_new") then
                ! velocity treated together
                norm(1) = norm(1) + product(dx(1:params%dim))*sum( hvy_block(g(1)+1:Bs(1)+g(1), g(2)+1:Bs(2)+g(2), g(3)+1:Bs(3)+g(3), 1:params%dim, hvy_id )**2 )
                norm(1:params%dim) = norm(1)  ! transcribe to other velocity components as well
                ! rest - should be only pressure?
                do p = params%dim+1, n_eqn
                    norm(p) = norm(p) + product(dx(1:params%dim))*sum( hvy_block(g(1)+1:Bs(1)+g(1), g(2)+1:Bs(2)+g(2), g(3)+1:Bs(3)+g(3), p, hvy_id )**2 )
                enddo
            else
                do p = 1, n_eqn
                    norm(p) = norm(p) + product(dx(1:params%dim))*sum( hvy_block(g(1)+1:Bs(1)+g(1), g(2)+1:Bs(2)+g(2), g(3)+1:Bs(3)+g(3), p, hvy_id )**2 )
                enddo
            endif
        enddo

        call MPI_ALLREDUCE(MPI_IN_PLACE, norm, n_eqn, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        ! norm = sqrt( norm )
        norm = sqrt( norm )
    

    case ("threshold-image-denoise")
        norm = 0.0_rk
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

            ! only sum up over leafs, we sadly don't have access to function block_is_leaf
            if (any(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)) /= -1)) cycle
            
            call get_block_spacing_origin_b( get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2)), params%domain_size, &
                params%Bs, x0, dx, dim=params%dim, level=lgt_block(lgt_id, IDX_MESH_LVL), max_level=params%Jmax)

            do p = 1, n_eqn
                norm(p) = norm(p) + product(dx(1:params%dim))*sum( hvy_block(g(1)+1:Bs(1)+g(1), g(2)+1:Bs(2)+g(2), g(3)+1:Bs(3)+g(3), p, hvy_id )**2 )
            enddo
        enddo

        call MPI_ALLREDUCE(MPI_IN_PLACE, norm, n_eqn, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        ! we have to normalize by the area as this is area weighted
        do p = 1, n_eqn
            norm(p) = norm(p) / product(params%domain_size(1:params%dim))
        enddo
    
    case ("threshold-cvs")
        norm = 0.0_rk

        ! Variant 1: compute after enstrophy

        ! Variant 2: compute after energy - kinetic energy and pressure energy
        norm = 0.0_rk
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

            ! only sum up over leafs, we sadly don't have access to function block_is_leaf
            if (any(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)) /= -1)) cycle

            call get_block_spacing_origin_b( get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2)), params%domain_size, &
                params%Bs, x0, dx, dim=params%dim, level=lgt_block(lgt_id, IDX_MESH_LVL), max_level=params%Jmax)
            
            ! JB ToDo: This for now goes against the idea of physics types, if we decide type by type the norm, then this has to move into the physics module
            if (params%physics_type == "ACM_new") then
                ! velocity treated together
                norm(1) = norm(1) + product(dx(1:params%dim))*sum( hvy_block(g(1)+1:Bs(1)+g(1), g(2)+1:Bs(2)+g(2), g(3)+1:Bs(3)+g(3), 1:params%dim, hvy_id )**2 )
                norm(1:params%dim) = norm(1)  ! transcribe to other velocity components as well
                ! rest - should be only pressure?
                do p = params%dim+1, n_eqn
                    norm(p) = norm(p) + product(dx(1:params%dim))*sum( hvy_block(g(1)+1:Bs(1)+g(1), g(2)+1:Bs(2)+g(2), g(3)+1:Bs(3)+g(3), p, hvy_id )**2 )
                enddo
            else
                do p = 1, n_eqn
                    norm(p) = norm(p) + product(dx(1:params%dim))*sum( hvy_block(g(1)+1:Bs(1)+g(1), g(2)+1:Bs(2)+g(2), g(3)+1:Bs(3)+g(3), p, hvy_id )**2 )
                enddo
            endif
        enddo

        call MPI_ALLREDUCE(MPI_IN_PLACE, norm, n_eqn, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        ! we have to normalize by the area as this is area weighted
        do p = 1, n_eqn
            norm(p) = norm(p) / product(params%domain_size(1:params%dim))
        enddo

    case ("Linfty")
        norm = -1.0_rk
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)

            ! only sum up over leafs, we sadly don't have access to function block_is_leaf
            if (any(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)) /= -1)) cycle

            do p = 1, n_eqn
                norm(p) = max(norm(p), abs(maxval( hvy_block(g(1)+1:Bs(1)+g(1), g(2)+1:Bs(2)+g(2), g(3)+1:Bs(3)+g(3), p, hvy_id ))) )
            enddo
        enddo

        call MPI_ALLREDUCE(MPI_IN_PLACE, norm, n_eqn, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

    case default
        write(*,'(A)') which_norm
        call abort(20030201, "The tree norm you desire is not implemented. How dare you.")

    end select

end subroutine
