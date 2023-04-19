!> \brief Refine mesh (2D version). All cpu loop over their heavy data and check if the refinement
!! flag +1 is set on the block. If so, we take this block, interpolate it to the next finer
!! level and create four new blocks, each carrying a part of the interpolated data.
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

    integer(kind=ik)                    :: k, N, dF                             ! loop variables
    integer(kind=ik)                    :: rank                                 ! process rank
    integer(kind=ik)                    :: g                                    ! grid parameter
    integer(kind=ik), dimension(3)      :: Bs
    real(kind=rk), allocatable, save    :: new_data(:,:,:), data_predict_fine(:,:)  ! data fields for interpolation
    integer(kind=ik)                    :: lgt_free_id, free_heavy_id, lgt_id
    integer(kind=ik)                    :: treecode(params%Jmax)
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
    ! will become four blocks), without the ghost nodes.
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
        if ( (lgt_block( lgt_id, params%Jmax + idx_refine_sts) == +1) ) then

            ! extract treecode and mesh level
            treecode = lgt_block( lgt_id, 1:params%Jmax )
            level    = lgt_block( lgt_id, params%Jmax + IDX_MESH_LVL )
            ! ------------------------------------------------------------------------------------------------------
            ! first: interpolate block data
            ! loop over all data fields
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
            ! first new block
            ! find a free light id on this rank
            call get_free_local_light_id( params, rank, lgt_free_id)
            call lgt2hvy( free_heavy_id, lgt_free_id, rank, N )

            ! write new light data
            ! old treecode
            lgt_block( lgt_free_id, 1:params%Jmax ) = treecode
            ! new treecode one level up - "0" block
            lgt_block( lgt_free_id, level+1 ) = 0
            ! new level + 1
            lgt_block( lgt_free_id, params%Jmax + IDX_MESH_LVL ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%Jmax + idx_refine_sts ) = 0
            ! the tree_ID is the same as the one of the mother block
            lgt_block( lgt_free_id, params%Jmax + IDX_TREE_ID ) = tree_ID


            ! save interpolated data, loop over all datafields
            do dF = 1, size(hvy_block,3)
                ! hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, dF, free_heavy_id ) = new_data(1:Bs(1), 1:Bs(2), dF)
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, dF, free_heavy_id ) = new_data(1:Bs(1), 1:Bs(2), dF)
            end do

            !--------------------------
            ! second new block
            ! find a free light id on this rank
            call get_free_local_light_id( params, rank, lgt_free_id)
            call lgt2hvy( free_heavy_id, lgt_free_id, rank, N )

            ! write new light data
            ! old treecode
            lgt_block( lgt_free_id, 1:params%Jmax ) = treecode
            ! new treecode one level up - "1" block
            lgt_block( lgt_free_id, level+1 ) = 1
            ! new level + 1
            lgt_block( lgt_free_id, params%Jmax + IDX_MESH_LVL ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%Jmax + idx_refine_sts ) = 0
            ! the tree_ID is the same as the one of the mother block
            lgt_block( lgt_free_id, params%Jmax + IDX_TREE_ID ) = tree_ID


            ! save interpolated data, loop over all datafields
            do dF = 1, size(hvy_block,3)
                ! hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, dF, free_heavy_id ) = new_data(1:Bs(1), Bs(2):2*Bs(2)-1, dF)
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, dF, free_heavy_id ) = new_data(1:Bs(1), Bs(2)+1:2*Bs(2), dF)
            end do

            !--------------------------
            ! third new block
            ! find a free light id on this rank
            call get_free_local_light_id( params, rank, lgt_free_id)
            call lgt2hvy( free_heavy_id, lgt_free_id, rank, N )

            ! write new light data
            ! old treecode
            lgt_block( lgt_free_id, 1:params%Jmax ) = treecode
            ! new treecode one level up - "1" block
            lgt_block( lgt_free_id, level+1 ) = 2
            ! new level + 1
            lgt_block( lgt_free_id, params%Jmax + IDX_MESH_LVL ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%Jmax + idx_refine_sts ) = 0
            ! the tree_ID is the same as the one of the mother block
            lgt_block( lgt_free_id, params%Jmax + IDX_TREE_ID ) = tree_ID


            ! save interpolated data, loop over all datafields
            do dF = 1, size(hvy_block,3)
                ! hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, dF, free_heavy_id ) = new_data(Bs(1):2*Bs(1)-1, 1:Bs(2), dF)
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, dF, free_heavy_id ) = new_data(Bs(1)+1:2*Bs(1), 1:Bs(2), dF)
            end do

            !--------------------------
            ! fourth new block
            ! write data on current heavy id
            free_heavy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_free_id, free_heavy_id, rank, N )

            ! write new light data
            ! old treecode
            lgt_block( lgt_free_id, 1:params%Jmax ) = treecode
            ! new treecode one level up - "1" block
            lgt_block( lgt_free_id, level+1 ) = 3
            ! new level + 1
            lgt_block( lgt_free_id, params%Jmax + IDX_MESH_LVL ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%Jmax + idx_refine_sts ) = 0
            ! the tree_ID is the same as the one of the mother block
            lgt_block( lgt_free_id, params%Jmax + IDX_TREE_ID ) = tree_ID


            ! save interpolated data, loop over all datafields
            do dF = 1, size(hvy_block,3)
                ! hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, dF, free_heavy_id ) = new_data(Bs(1):2*Bs(1)-1, Bs(2):2*Bs(2)-1, dF)
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, dF, free_heavy_id ) = new_data(Bs(1)+1:2*Bs(1), Bs(2)+1:2*Bs(2), dF)
            end do

        end if

    end do

    ! synchronize light data
    call synchronize_lgt_data( params, refinement_status_only=.false. )

end subroutine refinementExecute2D_tree
