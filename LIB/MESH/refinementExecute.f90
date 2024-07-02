!> \brief Refine block by interpolating it to the next finer level and creating
!!  four(2D) or eight(3D) new blocks, each carrying a part of the interpolated data.
!> Interpolated data is set in array, and last daughter overwrites motherblock
subroutine refineBlock(params, hvy_block, hvyID, tree_ID)
    implicit none

    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)             !> heavy data array - block data
    integer(kind=ik), intent(in)        :: hvyID                                !> Block to be reined
    integer(kind=ik), intent(in)        :: tree_ID

    integer(kind=ik)                    :: k, N, dF, k_daughter                 ! loop variables
    integer(kind=ik)                    :: rank, dim
    integer(kind=ik)                    :: g                                    ! grid parameter
    integer(kind=ik), dimension(3)      :: Bs
    real(kind=rk), allocatable, save    :: new_data(:,:,:,:), data_predict_fine(:,:,:)  ! data fields for interpolation
    integer(kind=ik)                    :: lgt_free_id, free_heavy_id, lgt_id
    integer(kind=tsize)                 :: treecode
    integer(kind=ik)                    :: level

    dim = params%dim
    N = params%number_blocks
    rank = params%rank
    Bs = params%Bs
    g  = params%g

    ! data fields for interpolation
    ! coarse: current data, fine: new (refine) data, new_data: gather all refined data for all data fields
    ! NOTE: the predictor for the refinement acts on the extended blocks i.e. it
    ! includes the ghost nodes layer. Therefore, you MUST call sync_ghosts before this routine.
    ! The datafield for prediction is one level up, i.e. it contains Bs+g + (Bs+2g-1) points
    if (dim == 2) then
        if (.not. allocated(data_predict_fine)) allocate( data_predict_fine( 2*(Bs(1)+2*g)-1, 2*(Bs(2)+2*g)-1, 1) )
    else
        if (.not. allocated(data_predict_fine)) allocate( data_predict_fine( 2*(Bs(1)+2*g)-1, 2*(Bs(2)+2*g)-1, 2*(Bs(3)+2*g)-1) )
    endif
    ! the new_data field holds the interior part of the new, refined block (which
    ! will become four/eight blocks), without the ghost nodes.
    ! uniqueGrid modification
    if (allocated(new_data)) then
        if (size(new_data, 4) < size(hvy_block,4)) deallocate(new_data)
    endif
    if (dim == 2) then
        if (.not. allocated(new_data)) allocate( new_data(2*Bs(1), 2*Bs(2), 1, size(hvy_block,4)) )
    else
        if (.not. allocated(new_data)) allocate( new_data(2*Bs(1), 2*Bs(2), 2*Bs(3), size(hvy_block,4)) )
    endif

    ! extract treecode and mesh level
    call hvy2lgt(lgt_id, hvyID, rank, N)
    treecode = get_tc(lgt_block( lgt_id, IDX_TC_1 : IDX_TC_2 ))
    level    = lgt_block( lgt_id, IDX_MESH_LVL )

    ! ------------------------------------------------------------------------------------------------------
    ! first: interpolate block data
    ! loop over all components
    do dF = 1, size(hvy_block,4)
        ! NOTE: the refinement interpolation acts on the entire block including ghost nodes.
        ! interpolate data, then save new data, but cut ghost nodes.
        ! uniqueGrid modification
        if (dim == 2) then
            call prediction_2D(hvy_block(:, :, 1, dF, hvyID), data_predict_fine(:,:,1), params%order_predictor)
            new_data(:,:,1,dF) = data_predict_fine( 2*g+1:2*Bs(1)+2*g, 2*g+1:2*Bs(2)+2*g, 1 )
        else
            call prediction_3D(hvy_block(:, :, :, dF, hvyID), data_predict_fine, params%order_predictor)
            new_data(:,:,:,dF) = data_predict_fine( 2*g+1:2*Bs(1)+2*g, 2*g+1:2*Bs(2)+2*g, 2*g+1:2*Bs(3)+2*g )
        endif
    end do

    ! ------------------------------------------------------------------------------------------------------
    ! second: split new data and write into new blocks
    !--------------------------
    ! create all four daughters
    do k_daughter = 0,2**(params%dim) -1
        ! find a free light id on this rank
        if (k_daughter < 2**(params%dim) -1) then
            call get_free_local_light_id( params, rank, lgt_free_id, message="refinement_execute")
            call lgt2hvy( free_heavy_id, lgt_free_id, rank, N )
        ! last daughter overwrites mother
        else
            free_heavy_id = hvyID
            call hvy2lgt( lgt_free_id, free_heavy_id, rank, N )
        endif

        treecode = tc_set_digit_at_level_b(treecode, k_daughter, dim=params%dim, level=level+1, max_level=params%Jmax)

        ! init array - needed to change values if never adressed
        lgt_block( lgt_free_id, : ) = -1

        ! write new light data
        ! new treecode
        call set_tc(lgt_block( lgt_free_id, IDX_TC_1 : IDX_TC_2 ), treecode)
        ! new level + 1
        lgt_block( lgt_free_id, IDX_MESH_LVL ) = level+1
        ! new blocks have refinement_status==0 (STAY)
        lgt_block( lgt_free_id, idx_refine_sts ) = 0
        ! the tree_ID is the same as the one of the mother block
        lgt_block( lgt_free_id, IDX_TREE_ID ) = tree_ID

        ! save interpolated data, loop over all datafields - select data from correct quadrant
        ! k=0 -> X_l, Y_l; k=1 -> X_l, Y_r; k=2 -> X_r, Y_l; k=3 -> X_r, Y_r;
        do dF = 1, size(hvy_block,4)
            if (dim == 2) then
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, 1, dF, free_heavy_id ) =  new_data( &
                        (k_daughter/2)*Bs(1) +1 : Bs(1)*(1+      (k_daughter/2)), &
                    modulo(k_daughter,2)*Bs(2) +1 : Bs(2)*(1+modulo(k_daughter,2)), 1, dF)
            else
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, free_heavy_id ) =  new_data( &
                    modulo(k_daughter/2,2)*Bs(1) +1 : Bs(1)*(1+modulo(k_daughter/2,2)), &
                    modulo(k_daughter  ,2)*Bs(2) +1 : Bs(2)*(1+modulo(k_daughter  ,2)), &
                        (k_daughter/4  )*Bs(3) +1 : Bs(3)*(1+      (k_daughter/4  )), dF)
            endif
        end do
    end do

end subroutine refineBlock



!> \brief Refine block directly to spaghetti decomposition
!!    SC is copied to correct quadrant / octant
!!    WC is set to 0
!> All daughter blocks are created and at last mother is set as inactive
!  Module MPI could be used in order to do the copying more easily, as we have the indices there in ijkPatches
subroutine refineBlock2SpaghettiWD(params, hvy_block, hvyID, tree_ID)
    implicit none

    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)             !> heavy data array - block data
    integer(kind=ik), intent(in)        :: hvyID                                !> Block to be reined
    integer(kind=ik), intent(in)        :: tree_ID

    integer(kind=ik)                    :: k, N, dF, k_daughter                 ! loop variables
    integer(kind=ik)                    :: rank, dim, g, data_bounds(2,3)
    integer(kind=ik), dimension(3)      :: Bs
    integer(kind=ik)                    :: lgt_free_id, free_heavy_id, lgt_id
    integer(kind=tsize)                 :: treecode
    integer(kind=ik)                    :: level, lgtID
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: tmp_wd  ! used for data juggling

    if (allocated(tmp_wd)) then
        if (size(tmp_wd, 4) < size(hvy_block, 4)) deallocate(tmp_wd)
    endif
    if (.not. allocated(tmp_wd)) allocate(tmp_wd(1:size(hvy_block, 1), 1:size(hvy_block, 2), 1:size(hvy_block, 3), 1:size(hvy_block, 4)) )


    dim = params%dim
    N = params%number_blocks
    rank = params%rank
    Bs = params%Bs
    g  = params%g
    call hvy2lgt(lgtID, hvyID, rank, N)

    ! extract treecode and mesh level
    treecode = get_tc(lgt_block( lgt_id, IDX_TC_1 : IDX_TC_2 ))
    level    = lgt_block( lgt_id, IDX_MESH_LVL )

    ! ------------------------------------------------------------------------------------------------------
    ! second: split new data and write into new blocks
    !--------------------------
    ! create all four daughters
    do k_daughter = 0,2**(params%dim) -1
        ! find a free light id on this rank
        if (k_daughter < 2**(params%dim) -1) then
            call get_free_local_light_id( params, rank, lgt_free_id, message="refinement_execute")
            call lgt2hvy( free_heavy_id, lgt_free_id, rank, N )
        endif

        treecode = tc_set_digit_at_level_b(treecode, k_daughter, dim=params%dim, level=level+1, max_level=params%Jmax)

        ! init array - needed to change values if never adressed
        lgt_block( lgt_free_id, : ) = -1

        ! write new light data
        ! new treecode
        call set_tc(lgt_block( lgt_free_id, IDX_TC_1 : IDX_TC_2 ), treecode)
        ! new level + 1
        lgt_block( lgt_free_id, IDX_MESH_LVL ) = level+1
        ! new blocks have refinement_status==0 (STAY)
        lgt_block( lgt_free_id, idx_refine_sts ) = 0
        ! the tree_ID is the same as the one of the mother block
        lgt_block( lgt_free_id, IDX_TREE_ID ) = tree_ID

        ! save interpolated data, loop over all datafields - select data from correct quadrant
        ! k=0 -> X_l, Y_l; k=1 -> X_l, Y_r; k=2 -> X_r, Y_l; k=3 -> X_r, Y_r;
        do dF = 1, size(hvy_block,4)
            ! Init values as 0 to set all WC to 0
            tmp_wd(:, :, :, df) = 0.0_rk

            ! compute bounds from which quadrant / octant to copy
            data_bounds = -1
            data_bounds(1, 1) = g+1       + Bs(1)/2 * modulo(k_daughter/2, 2)
            data_bounds(2, 1) = g+Bs(1)/2 + Bs(1)/2 * modulo(k_daughter/2, 2)
            data_bounds(1, 2) = g+1       + Bs(2)/2 * modulo(k_daughter  , 2)
            data_bounds(2, 2) = g+Bs(2)/2 + Bs(2)/2 * modulo(k_daughter  , 2)
            data_bounds(1, 3) = g+1       + Bs(3)/2 * modulo(k_daughter/4, 2)
            data_bounds(2, 3) = g+Bs(3)/2 + Bs(3)/2 * modulo(k_daughter/4, 2)

            ! now copy values from this quadrant / octant to SC area in Mallat format
            if (dim == 2) then
                tmp_wd(ceiling(g/2.0)+1:ceiling(g/2.0)+Bs(1)/2, ceiling(g/2.0)+1:ceiling(g/2.0)+Bs(2)/2, 1, df) = &
                    hvy_block(data_bounds(1, 1):data_bounds(2, 1), data_bounds(1, 2):data_bounds(2, 2), 1, df, hvyID)
            else
                tmp_wd(ceiling(g/2.0)+1:ceiling(g/2.0)+Bs(1)/2, ceiling(g/2.0)+1:ceiling(g/2.0)+Bs(2)/2, &
                    ceiling(g/2.0)+1:ceiling(g/2.0)+Bs(3)/2, df) = &
                    hvy_block(data_bounds(1, 1):data_bounds(2, 1), data_bounds(1, 2):data_bounds(2, 2), &
                    data_bounds(1, 3):data_bounds(2, 3), df, hvyID)
            endif

        end do
        ! transform to spaghetti form as this is what we usually work with
        call Mallat2Spaghetti_block(params, tmp_wd(:, :, :, 1:size(hvy_block,4)), hvy_block(:, :, :, 1:size(hvy_block,4), free_heavy_id))
    end do

    ! delete mother
    lgt_block(lgtID, :) = -1
    lgt_block(lgtID, IDX_REFINE_STS) = 0

end subroutine refineBlock2SpaghettiWD



!  ToDo: Merge with other function and use refineBlock
!> \brief Refine mesh (2D version). All cpu loop over their heavy data and check if the refinement
!! flag +1 is set on the block. If so, we take this block, interpolate it to the next finer
!! level and create four(2D) or eight(3D) new blocks, each carrying a part of the interpolated data.
!! As all CPU first work individually, the light data array is synced afterwards.
!
!> \note The interpolation (or prediction) operator here is applied to a block INCLUDING
!! any ghost nodes. You must sync first.
!
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays \n
! ********************************************************************************************
subroutine refinementExecute2D_tree( params, hvy_block, tree_ID )

    implicit none

    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :)                !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID

    integer(kind=ik)                    :: k, N, dF, k_daughter                 ! loop variables
    integer(kind=ik)                    :: rank                                 ! process rank
    integer(kind=ik)                    :: g                                    ! grid parameter
    integer(kind=ik), dimension(3)      :: Bs
    real(kind=rk), allocatable, save    :: new_data(:,:,:), data_predict_fine(:,:)  ! data fields for interpolation
    integer(kind=ik)                    :: lgt_free_id, free_heavy_id, lgt_id
    integer(kind=tsize)                 :: treecode
    integer(kind=ik)                    :: level

    N = params%number_blocks
    rank = params%rank
    Bs = params%Bs
    g  = params%g

    ! data fields for interpolation
    ! coarse: current data, fine: new (refine) data, new_data: gather all refined data for all data fields
    ! NOTE: the predictor for the refinement acts on the extended blocks i.e. it
    ! includes the ghost nodes layer. Therefore, you MUST call sync_ghosts before this routine.
    ! The datafield for prediction is one level up, i.e. it contains Bs+g + (Bs+2g-1) points
    if (.not. allocated(data_predict_fine)) allocate( data_predict_fine( 2*(Bs(1)+2*g)-1, 2*(Bs(2)+2*g)-1 ) )
    ! the new_data field holds the interior part of the new, refined block (which
    ! will become four/eight blocks), without the ghost nodes.
    ! uniqueGrid modification
    if (allocated(new_data)) then
        if (size(new_data, 3) < size(hvy_block,3)) deallocate(new_data)
    endif
    if (.not. allocated(new_data)) allocate( new_data(2*Bs(1), 2*Bs(2), size(hvy_block,3)) )


    ! every proc loop over its active heavy data array
    do k = 1, hvy_n(tree_ID)

        ! calculate light id
        call hvy2lgt( lgt_id, hvy_active(k, tree_ID), rank, N )

        ! block wants to refine
        if ( (lgt_block( lgt_id, idx_refine_sts) == +1) ) then

            ! extract treecode and mesh level
            treecode = get_tc(lgt_block( lgt_id, IDX_TC_1 : IDX_TC_2 ))
            level    = lgt_block( lgt_id, IDX_MESH_LVL )

            ! ------------------------------------------------------------------------------------------------------
            ! first: interpolate block data
            ! loop over all components
            do dF = 1, size(hvy_block,3)
                ! NOTE: the refinement interpolation acts on the entire block including ghost nodes.
                ! interpolate data
                call prediction_2D(hvy_block(:, :, dF, hvy_active(k, tree_ID)), data_predict_fine, params%order_predictor)
                ! save new data, but cut ghost nodes.
                ! uniqueGrid modification
                new_data(:,:,dF) = data_predict_fine( 2*g+1:2*Bs(1)+2*g, 2*g+1:2*Bs(2)+2*g )
            end do

            ! ------------------------------------------------------------------------------------------------------
            ! second: split new data and write into new blocks
            !--------------------------
            ! create all four daughters
            do k_daughter = 0,3
                ! find a free light id on this rank
                if (k_daughter < 3) then
                    call get_free_local_light_id( params, rank, lgt_free_id, message="refinement_execute")
                    call lgt2hvy( free_heavy_id, lgt_free_id, rank, N )
                ! last daughter overwrites mother
                else
                    free_heavy_id = hvy_active(k, tree_ID)
                    call hvy2lgt( lgt_free_id, free_heavy_id, rank, N )
                endif

                treecode = tc_set_digit_at_level_b(treecode, k_daughter, dim=params%dim, level=level+1, max_level=params%Jmax)

                ! init array - needed to change values if never adressed
                lgt_block( lgt_free_id, : ) = -1

                ! write new light data
                ! new treecode
                call set_tc(lgt_block( lgt_free_id, IDX_TC_1 : IDX_TC_2 ), treecode)
                ! new level + 1
                lgt_block( lgt_free_id, IDX_MESH_LVL ) = level+1
                ! new blocks have refinement_status==0 (STAY)
                lgt_block( lgt_free_id, idx_refine_sts ) = 0
                ! the tree_ID is the same as the one of the mother block
                lgt_block( lgt_free_id, IDX_TREE_ID ) = tree_ID

                ! save interpolated data, loop over all datafields - select data from correct quadrant
                ! k=0 -> X_l, Y_l; k=1 -> X_l, Y_r; k=2 -> X_r, Y_l; k=3 -> X_r, Y_r;
                do dF = 1, size(hvy_block,3)
                    hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, dF, free_heavy_id ) =  new_data( &
                              (k_daughter/2)*Bs(1) +1 : Bs(1)*(1+      (k_daughter/2)), &
                        modulo(k_daughter,2)*Bs(2) +1 : Bs(2)*(1+modulo(k_daughter,2)), dF)
                end do
            end do
        end if

    end do

    ! synchronize light data
    call synchronize_lgt_data( params, refinement_status_only=.false. )

end subroutine refinementExecute2D_tree



!  ToDo: Merge with other function and use refineBlock
!> \brief Refine mesh (3D version). All cpu loop over their heavy data and check if the refinement
!! flag +1 is set on the block. If so, we take this block, interpolate it to the next finer
!! level and create four(2D) or eight(3D) new blocks, each carrying a part of the interpolated data.
!! As all CPU first work individually, the light data array is synced afterwards.
!
!> \note The interpolation (or prediction) operator here is applied to a block INCLUDING
!! any ghost nodes. You must sync first.
!
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays \n
! ********************************************************************************************
subroutine refinementExecute3D_tree( params, hvy_block, tree_ID )

    implicit none

    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)             !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID

    integer(kind=ik)                    :: k, N, dF, k_daughter                 ! loop variables
    integer(kind=ik)                    :: rank                                 ! process rank
    integer(kind=ik)                    :: g                                    ! grid parameter
    integer(kind=ik), dimension(3)      :: Bs
    real(kind=rk), allocatable, save    :: new_data(:,:,:,:), data_predict_fine(:,:,:)  ! data fields for interpolation
    integer(kind=ik)                    :: lgt_free_id, free_heavy_id, lgt_id 
    integer(kind=tsize)                 :: treecode
    integer(kind=ik)                    :: level

    N = params%number_blocks
    rank = params%rank
    Bs = params%Bs
    g  = params%g

    ! NOTE: the predictor for the refinement acts on the extended blocks i.e. it
    ! includes the ghost nodes layer. Therefore, you MUST call sync_ghosts before this routine.
    ! The datafield for prediction is one level up, i.e. it contains Bs+g + (Bs+2g-1) points
    if (.not.allocated(data_predict_fine)) allocate( data_predict_fine(2*(Bs(1)+2*g)-1, 2*(Bs(2)+2*g)-1, 2*(Bs(3)+2*g)-1) )
    ! the new_data field holds the interior part of the new, refined block (which
    ! will become four/eight blocks), without the ghost nodes.
    ! uniqueGrid modification
    if (allocated(new_data)) then
        if (size(new_data, 4) < size(hvy_block,4)) deallocate(new_data)
    endif
    if (.not.allocated(new_data)) allocate( new_data(2*Bs(1), 2*Bs(2), 2*Bs(3), size(hvy_block,4)) )


    ! every proc loop over its active heavy data array
    do k = 1, hvy_n(tree_ID)

        ! light data id
        call hvy2lgt( lgt_id, hvy_active(k, tree_ID), rank, N )

        ! block wants to refine
        if ( (lgt_block( lgt_id, idx_refine_sts) == +1) ) then

            ! treecode and mesh level
            treecode = get_tc(lgt_block( lgt_id,IDX_TC_1 : IDX_TC_2 ))
            level    = lgt_block( lgt_id, IDX_MESH_LVL )

            ! ------------------------------------------------------------------------------------------------------
            ! first: interpolate block data
            ! loop over all components
            do dF = 1, size(hvy_block,4)
                ! NOTE: the refinement interpolation acts on the entire block including ghost nodes.
                ! interpolate data
                call prediction_3D(hvy_block(:, :, :, dF, hvy_active(k, tree_ID)), data_predict_fine, params%order_predictor)
                ! save new data, but cut ghost nodes.
                new_data(:, :, :, dF) = data_predict_fine(2*g+1:2*g+2*Bs(1), 2*g+1:2*g+2*Bs(2), 2*g+1:2*g+2*Bs(3))
            end do

            ! ------------------------------------------------------------------------------------------------------
            ! second: split new data and write into new blocks
            !--------------------------
            ! create all eight daughters
            do k_daughter = 0,7
                ! find a free light id on this rank
                if (k_daughter < 7) then
                    call get_free_local_light_id( params, rank, lgt_free_id, message="refinement_execute")
                    call lgt2hvy( free_heavy_id, lgt_free_id, rank, N )
                ! last daughter overwrites mother
                else
                    free_heavy_id = hvy_active(k, tree_ID)
                    call hvy2lgt( lgt_free_id, free_heavy_id, rank, N )
                endif

                treecode = tc_set_digit_at_level_b(treecode, k_daughter, dim=params%dim, level=level+1, max_level=params%Jmax)

                ! init array - needed to change values if never adressed
                lgt_block( lgt_free_id, : ) = -1

                ! write new light data
                ! new treecode
                call set_tc(lgt_block( lgt_free_id, IDX_TC_1 : IDX_TC_2 ), treecode)
                ! new level + 1
                lgt_block( lgt_free_id, IDX_MESH_LVL ) = level+1
                ! new blocks have refinement_status==0 (STAY)
                lgt_block( lgt_free_id, idx_refine_sts ) = 0
                ! the tree_ID is the same as the one of the mother block
                lgt_block( lgt_free_id, IDX_TREE_ID ) = tree_ID

                ! save interpolated data, loop over all datafields - select data from correct octant
                ! k=0 -> X_l, Y_l, Z_l;     k=1 -> X_l, Y_r, Z_l;
                ! k=2 -> X_r, Y_l, Z_l;     k=3 -> X_r, Y_r, Z_l;
                ! k=4 -> X_l, Y_l, Z_r;     k=5 -> X_l, Y_r, Z_r;
                ! k=6 -> X_r, Y_l, Z_r;     k=7 -> X_r, Y_r, Z_r;
                do dF = 1, size(hvy_block,4)
                    hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, free_heavy_id ) =  new_data( &
                        modulo(k_daughter/2,2)*Bs(1) +1 : Bs(1)*(1+modulo(k_daughter/2,2)), &
                        modulo(k_daughter  ,2)*Bs(2) +1 : Bs(2)*(1+modulo(k_daughter  ,2)), &
                              (k_daughter/4  )*Bs(3) +1 : Bs(3)*(1+      (k_daughter/4  )),dF)
                end do
            end do

        end if

    end do

    ! synchronize light data
    call synchronize_lgt_data( params, refinement_status_only=.false. )

end subroutine refinementExecute3D_tree



!> \brief Refine mesh. All cpu loop over their heavy data and check if the refinement
!! flag +1 is set on the block. If so, we either take this block, interpolate it to the next finer
!! level and create four(2D) or eight(3D) new blocks, each carrying a part of the interpolated data,
!! or copy them to Mallat decomposed blocks if they are on the input level
!! As all CPU first work individually, the light data array is synced afterwards.
!> This function is used when parts of the blocks are already waveled decomposed in adapt_tree,
!! so for block on the given level it refines directly to mallat wavelet decomposition (and then spaghetti),
!! for other blocks it does normal refinement
!
!> \note The interpolation (or prediction) operator here is applied to a block INCLUDING
!! any ghost nodes. You must sync first. This does not apply to blocks which will be in MallatWD though
!
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays \n
! ********************************************************************************************
subroutine refinementExecute_lvl2SpaghettiWD( params, hvy_block, tree_ID, level )

    implicit none

    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)             !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID                              !> Tree to refine
    integer(kind=ik), intent(in)        :: level                                !> Level, where blocks will be refined to MallatWD

    integer(kind=ik)                    :: k, N, dF, k_daughter                 ! loop variables
    integer(kind=ik)                    :: rank                                 ! process rank
    integer(kind=ik)                    :: g                                    ! grid parameter
    integer(kind=ik), dimension(3)      :: Bs
    integer(kind=ik)                    :: lgt_free_id, free_heavy_id, lgt_id 
    integer(kind=tsize)                 :: treecode
    integer(kind=ik)                    :: level_me, hvyID

    N = params%number_blocks
    rank = params%rank
    Bs = params%Bs
    g  = params%g

    ! every proc loop over its active heavy data array
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k, tree_ID)

        ! light data id
        call hvy2lgt( lgt_id, hvyID, rank, N )

        ! block wants to refine
        if ( (lgt_block( lgt_id, idx_refine_sts) == +1) ) then
            level_me = lgt_block( lgt_id, IDX_MESH_LVL)
            
            ! Daughters created as Mallat WDed blocks
            if (level_me == level) then
                call refineBlock2SpaghettiWD(params, hvy_block, hvyID, tree_ID)
            ! Daughters refined normally
            else
                call refineBlock(params, hvy_block, hvyID, tree_ID)
            endif

        end if
    end do

    ! synchronize light data
    call synchronize_lgt_data( params, refinement_status_only=.false. )

end subroutine refinementExecute_lvl2SpaghettiWD