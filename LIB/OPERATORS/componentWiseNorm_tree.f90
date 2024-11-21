subroutine componentWiseNorm_tree(params, hvy_block, tree_ID, which_norm, norm, norm_case, n_val)
    implicit none

    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)             !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID
    character(len=*), intent(in)        :: which_norm                           !> which norm to use ? "L2", "Linfty"
    real(kind=rk), intent(inout)        :: norm(:)                              !> the computed norm for each component of the vector
    !> String representing to choose special norm cases, can be tree (default), level or ref and additionally
    !! if we add "SC" to the norm case logic, we treat the input as wavelet decomposed and only compute the norms over the SCs
    character(len=*), intent(in), optional  :: norm_case
    !> Additional value to be considered for norm logic, can be level or refinement status to which should be synced, used if sync case includes ref or level
    integer(kind=ik), intent(in), optional  :: n_val

    real(kind=rk)                       :: x0(1:3), dx(1:3), volume
    integer(kind=ik) :: k, hvy_id, n_eqn, Bs(1:3), g(1:3), p, p_norm, l, mpierr, lgt_id, D, level_me, ref_me, norm_case_id
    logical          :: SC_only
    integer, dimension(:), allocatable  :: mask_i, mask

    ! note: if norm and hvy_block components are of different size, we use the smaller one.
    n_eqn = min( size(norm, 1), size(hvy_block,4) )
    Bs = params%Bs
    g = 0
    g(1:params%dim) = params%g
    D = params%dim

    ! for threshold state vector components we do some slicing magic
    ! fortran does not allow logical mask slicing, so we build a mask of indices to slice with
    if (allocated(mask_i)) then
        if (size(mask_i) < n_eqn) deallocate(mask_i)
    endif
    if (.not. allocated(mask_i)) allocate(mask_i(1:n_eqn))
    do l = 1, n_eqn
        mask_i(l) = l
    end do

    if (.not. present(norm_case)) then
        norm_case_id = 1
    else
        ! treat the norm computations as for SC
        ! with this, we will skip every second point in each direction and multiply the L2norms for the fields by 2^d, as the cells are larger
        if (index(norm_case, "SC") > 0) then
            norm_case_id = 2
        else 
            norm_case_id = 1
        endif
    
        ! now lets treat the special restrictions, set to the second digit
        if (index(norm_case, "level") > 0) norm_case_id = norm_case_id + 10*1
        if (index(norm_case, "ref") > 0) norm_case_id = norm_case_id + 10*2
    endif

    select case (which_norm)
    case ("L2", "H1")
        norm = 0.0_rk
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            level_me = lgt_block( lgt_ID, IDX_MESH_LVL )
            ref_me   = lgt_block( lgt_ID, IDX_REFINE_STS)
            ! skip some blocks in case we want to compute the norm over specific parts
            if (norm_case_ID/10 == 1) then
                if (level_me /= n_val) cycle
            elseif (norm_case_ID/10 == 2) then
                if (ref_me /= n_val) cycle
            endif

            ! only sum up over leafs, we sadly don't have access to function block_is_leaf
            if (any(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)) /= -1)) cycle

            call get_block_spacing_origin_b( get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2)), params%domain_size, &
                params%Bs, x0, dx, dim=params%dim, level=lgt_block(lgt_id, IDX_MESH_LVL), max_level=params%Jmax)

            ! compute norm for thresholding components that are treated on their own
            do p = 1, n_eqn
                if (params%threshold_state_vector_component(p) == 1) then
                    norm(p) = norm(p) + product(dx(1:params%dim))*mod(norm_case_id,10)**params%dim* &
                    sum( hvy_block(g(1)+1:Bs(1)+g(1):mod(norm_case_id,10), g(2)+1:Bs(2)+g(2):mod(norm_case_id,10), g(3)+1:Bs(3)+g(3):mod(norm_case_id,10), p, hvy_id )**2 )
                endif
            enddo

            ! some thresholding components want to be treated together, as for example the velocity, so we compute a norm-equivalent if wanted
            ! every value >1 in thresholding_component are treated as belonging together
            ! so we could set 2 2 3 3 and actually have those treated as independent vector fields, where we compute the norm
            ! maybe a bit overkill because usually we only have the velocity, but why not have the capacity that we can build up on
            do l = 2, maxval(params%threshold_state_vector_component(:))
                ! convert logical mask to index mask
                mask = pack(mask_i, params%threshold_state_vector_component(:)==l)
                norm(mask) = norm(mask) + product(dx(1:params%dim))*mod(norm_case_id,10)**params%dim* &
                sum( hvy_block(g(1)+1:Bs(1)+g(1):mod(norm_case_id,10), g(2)+1:Bs(2)+g(2):mod(norm_case_id,10), g(3)+1:Bs(3)+g(3):mod(norm_case_id,10), mask, hvy_id )**2 )
            enddo

            ! if (params%physics_type == "ACM_new") then
            !     ! velocity treated together
            !     norm(1) = norm(1) + product(dx(1:params%dim))*mod(norm_case_id,10)**params%dim* &
            !         sum( hvy_block(g(1)+1:Bs(1)+g(1):mod(norm_case_id,10), g(2)+1:Bs(2)+g(2):mod(norm_case_id,10), g(3)+1:Bs(3)+g(3):mod(norm_case_id,10), 1:params%dim, hvy_id )**2 )
            !     norm(1:params%dim) = norm(1)  ! transcribe to other velocity components as well
            !     ! rest - should be only pressure?
            !     do p = params%dim+1, n_eqn
            !         norm(p) = norm(p) + product(dx(1:params%dim))*mod(norm_case_id,10)**params%dim* &
            !         sum( hvy_block(g(1)+1:Bs(1)+g(1):mod(norm_case_id,10), g(2)+1:Bs(2)+g(2):mod(norm_case_id,10), g(3)+1:Bs(3)+g(3):mod(norm_case_id,10), p, hvy_id )**2 )
            !     enddo
            ! else
            !     do p = 1, n_eqn
            !         norm(p) = norm(p) + product(dx(1:params%dim))*mod(norm_case_id,10)**params%dim* &
            !         sum( hvy_block(g(1)+1:Bs(1)+g(1):mod(norm_case_id,10), g(2)+1:Bs(2)+g(2):mod(norm_case_id,10), g(3)+1:Bs(3)+g(3):mod(norm_case_id,10), p, hvy_id )**2 )
            !     enddo
            ! endif
        enddo

        call MPI_ALLREDUCE(MPI_IN_PLACE, norm, n_eqn, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        ! norm = sqrt( norm )
        norm = sqrt( norm )
    

    case ("threshold-image-denoise")
        norm = 0.0_rk
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            level_me = lgt_block( lgt_ID, IDX_MESH_LVL )
            ref_me   = lgt_block( lgt_ID, IDX_REFINE_STS)
            ! skip some blocks in case we want to compute the norm over specific parts
            if (norm_case_ID/10 == 1) then
                if (level_me /= n_val) cycle
            elseif (norm_case_ID/10 == 2) then
                if (ref_me /= n_val) cycle
            endif

            ! only sum up over leafs, we sadly don't have access to function block_is_leaf
            if (any(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)) /= -1)) cycle
            
            call get_block_spacing_origin_b( get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2)), params%domain_size, &
                params%Bs, x0, dx, dim=params%dim, level=lgt_block(lgt_id, IDX_MESH_LVL), max_level=params%Jmax)

            ! compute integral over the volume
            volume = volume + product(dx(1:params%dim))*mod(norm_case_id,10)**params%dim

            do p = 1, n_eqn
                norm(p) = norm(p) + product(dx(1:params%dim))*mod(norm_case_id,10)**params%dim* &
                sum( hvy_block(g(1)+1:Bs(1)+g(1):mod(norm_case_id,10), g(2)+1:Bs(2)+g(2):mod(norm_case_id,10), g(3)+1:Bs(3)+g(3):mod(norm_case_id,10), p, hvy_id )**2 )
            enddo
        enddo

        call MPI_ALLREDUCE(MPI_IN_PLACE, norm, n_eqn, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        ! we have to normalize by the volume as this is volume weighted
        do p = 1, n_eqn
            ! we integrate over the full leaf layer, this is the domain size so we do not have to call MPI_ALLREDUCE for the volume integral
            if (norm_case_ID/10 == 0) then
                norm(p) = norm(p) / product(params%domain_size(1:params%dim))
            ! we compute only a part of it, so we compute the actual integral of the volume
            else
                call MPI_ALLREDUCE(MPI_IN_PLACE, volume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
                norm(p) = norm(p) / volume
            endif
        enddo
    
    case ("threshold-cvs")
        norm = 0.0_rk

        ! Variant 1: compute after enstrophy

        ! Variant 2: compute after energy - kinetic energy and pressure energy
        norm = 0.0_rk
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            level_me = lgt_block( lgt_ID, IDX_MESH_LVL )
            ref_me   = lgt_block( lgt_ID, IDX_REFINE_STS)
            ! skip some blocks in case we want to compute the norm over specific parts
            if (norm_case_ID/10 == 1) then
                if (level_me /= n_val) cycle
            elseif (norm_case_ID/10 == 2) then
                if (ref_me /= n_val) cycle
            endif

            ! only sum up over leafs, we sadly don't have access to function block_is_leaf
            if (any(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)) /= -1)) cycle

            call get_block_spacing_origin_b( get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2)), params%domain_size, &
                params%Bs, x0, dx, dim=params%dim, level=lgt_block(lgt_id, IDX_MESH_LVL), max_level=params%Jmax)

            ! compute integral over the volume
            volume = volume + product(dx(1:params%dim))*mod(norm_case_id,10)**params%dim

            ! compute norm for thresholding components that are treated on their own
            do p = 1, n_eqn
                if (params%threshold_state_vector_component(p) == 1) then
                    norm(p) = norm(p) + product(dx(1:params%dim))*mod(norm_case_id,10)**params%dim* &
                    sum( hvy_block(g(1)+1:Bs(1)+g(1):mod(norm_case_id,10), g(2)+1:Bs(2)+g(2):mod(norm_case_id,10), g(3)+1:Bs(3)+g(3):mod(norm_case_id,10), p, hvy_id )**2 )
                endif
            enddo

            ! some thresholding components want to be treated together, as for example the velocity, so we compute a norm-equivalent if wanted
            ! every value >1 in thresholding_component are treated as belonging together
            ! so we could set 2 2 3 3 and actually have those treated as independent vector fields, where we compute the norm
            ! maybe a bit overkill because usually we only have the velocity, but why not have the capacity that we can build up on
            do l = 2, maxval(params%threshold_state_vector_component(:))
                ! convert logical mask to index mask
                mask = pack(mask_i, params%threshold_state_vector_component(:)==l)
                norm(mask) = norm(mask) + product(dx(1:params%dim))*mod(norm_case_id,10)**params%dim* &
                sum( hvy_block(g(1)+1:Bs(1)+g(1):mod(norm_case_id,10), g(2)+1:Bs(2)+g(2):mod(norm_case_id,10), g(3)+1:Bs(3)+g(3):mod(norm_case_id,10), mask, hvy_id )**2 )
            enddo
            
            ! ! JB ToDo: This for now goes against the idea of physics types, if we decide type by type the norm, then this has to move into the physics module
            ! if (params%physics_type == "ACM_new") then
            !     ! velocity treated together
            !     norm(1) = norm(1) + product(dx(1:params%dim))*mod(norm_case_id,10)**params%dim* &
            !         sum( hvy_block(g(1)+1:Bs(1)+g(1):mod(norm_case_id,10), g(2)+1:Bs(2)+g(2):mod(norm_case_id,10), g(3)+1:Bs(3)+g(3):mod(norm_case_id,10), 1:params%dim, hvy_id )**2 )
            !     norm(1:params%dim) = norm(1)  ! transcribe to other velocity components as well
            !     ! rest - should be only pressure?
            !     do p = params%dim+1, n_eqn
            !         norm(p) = norm(p) + product(dx(1:params%dim))*mod(norm_case_id,10)**params%dim* &
            !             sum( hvy_block(g(1)+1:Bs(1)+g(1):mod(norm_case_id,10), g(2)+1:Bs(2)+g(2):mod(norm_case_id,10), g(3)+1:Bs(3)+g(3):mod(norm_case_id,10), p, hvy_id )**2 )
            !     enddo
            ! else
            !     do p = 1, n_eqn
            !         norm(p) = norm(p) + product(dx(1:params%dim))*mod(norm_case_id,10)**params%dim* &
            !             sum( hvy_block(g(1)+1:Bs(1)+g(1):mod(norm_case_id,10), g(2)+1:Bs(2)+g(2):mod(norm_case_id,10), g(3)+1:Bs(3)+g(3):mod(norm_case_id,10), p, hvy_id )**2 )
            !     enddo
            ! endif
        enddo

        call MPI_ALLREDUCE(MPI_IN_PLACE, norm, n_eqn, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        ! we have to normalize by the volume as this is volume weighted
        do p = 1, n_eqn
            ! we integrate over the full leaf layer, this is the domain size so we do not have to call MPI_ALLREDUCE for the volume integral
            if (norm_case_ID/10 == 0) then
                norm(p) = norm(p) / product(params%domain_size(1:params%dim))
            ! we compute only a part of it, so we compute the actual integral of the volume
            else
                call MPI_ALLREDUCE(MPI_IN_PLACE, volume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
                norm(p) = norm(p) / volume
            endif
        enddo

    case ("Linfty")
        norm = -1.0_rk
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            level_me = lgt_block( lgt_ID, IDX_MESH_LVL )
            ref_me   = lgt_block( lgt_ID, IDX_REFINE_STS)
            ! skip some blocks in case we want to compute the norm over specific parts
            if (norm_case_ID/10 == 1) then
                if (level_me /= n_val) cycle
            elseif (norm_case_ID/10 == 2) then
                if (ref_me /= n_val) cycle
            endif

            ! only sum up over leafs, we sadly don't have access to function block_is_leaf
            if (any(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)) /= -1)) cycle

            do p = 1, n_eqn
                if (params%threshold_state_vector_component(p) == 1) then
                    norm(p) = max(norm(p), maxval( abs(hvy_block(g(1)+1:Bs(1)+g(1):mod(norm_case_id,10), g(2)+1:Bs(2)+g(2):mod(norm_case_id,10), g(3)+1:Bs(3)+g(3):mod(norm_case_id,10), p, hvy_id ))) )
                endif
            enddo

            ! some thresholding components want to be treated together, as for example the velocity, so we compute a norm-equivalent if wanted
            ! every value >1 in thresholding_component are treated as belonging together
            ! so we could set 2 2 3 3 and actually have those treated as independent vector fields, where we compute the norm
            ! maybe a bit overkill because usually we only have the velocity, but why not have the capacity that we can build up on
            do l = 2, maxval(params%threshold_state_vector_component(:))
                ! convert logical mask to index mask
                mask = pack(mask_i, params%threshold_state_vector_component(:)==l)
                norm(mask) = max(norm(mask), maxval(sum( hvy_block(g(1)+1:Bs(1)+g(1):mod(norm_case_id,10), g(2)+1:Bs(2)+g(2):mod(norm_case_id,10), g(3)+1:Bs(3)+g(3):mod(norm_case_id,10), mask, hvy_id )**2, dim=4)) )
            enddo
        enddo

        call MPI_ALLREDUCE(MPI_IN_PLACE, norm, n_eqn, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

    case default
        write(*,'(A)') which_norm
        call abort(20030201, "The tree norm you desire is not implemented. How dare you.")

    end select

end subroutine
