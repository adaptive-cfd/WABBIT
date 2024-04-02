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

                treecode = tc_set_level_b(treecode, k_daughter, dim=params%dim, level=level+1, max_level=params%Jmax)

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
                ! as this block is upsampled, it has zero details (no need to compute them usinng FWT,
                ! even if coarseExtension modifies the block, it would just reset some WC coeffs to zero...)
                hvy_details(:, free_heavy_id) = 0.0_rk

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
