
!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name refinement_execute_3D.f90
!> \version 0.5
!> \author msr, engels
!
!> \brief Refine mesh (3D version). All cpu loop over their heavy data and check if the refinement
!! flag +1 is set on the block. If so, we take this block, interpolate it to the next finer
!! level and create 8 new blocks, each carrying a part of the interpolated data.
!! As all CPU first work individually, the light data array is synced afterwards.
!
!> \note The interpolation (or prediction) operator here is applied to a block INCLUDING
!! any ghost nodes. You must sync first.
!
!> \details
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays \n
!!
!!
!! = log ======================================================================================
!! \n
!! 03/02/17 - create
!
! ********************************************************************************************

subroutine refinement_execute_3D( params, lgt_block, hvy_block, hvy_active, hvy_n )

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! loop variables
    integer(kind=ik)                    :: k, N, dF, lgt_id
    ! process rank
    integer(kind=ik)                    :: rank
    ! grid parameter
    integer(kind=ik)                    :: g
    integer(kind=ik), dimension(3)      :: Bs
    ! data fields for interpolation
    real(kind=rk), allocatable, save    :: new_data(:,:,:,:), data_predict_fine(:,:,:)
    ! free light/heavy data id
    integer(kind=ik)                    :: lgt_free_id, free_heavy_id
    ! treecode varaible
    integer(kind=ik)                    :: treecode(params%max_treelevel)
    ! mesh level
    integer(kind=ik)                    :: level, tree_id

    N = params%number_blocks
    ! set MPI parameter
    rank = params%rank
    ! grid parameter
    Bs = params%Bs
    g  = params%n_ghosts

    ! NOTE: the predictor for the refinement acts on the extended blocks i.e. it
    ! includes the ghost nodes layer. Therefore, you MUST call sync_ghosts before this routine.
    ! The datafield for prediction is one level up, i.e. it contains Bs+g + (Bs+2g-1) points
    if (.not.allocated(data_predict_fine)) allocate( data_predict_fine(2*(Bs(1)+2*g)-1, 2*(Bs(2)+2*g)-1, 2*(Bs(3)+2*g)-1) )
    ! the new_data field holds the interior part of the new, refined block (which
    ! will become four blocks), without the ghost nodes.
    if (.not.allocated(new_data)) allocate( new_data(2*Bs(1)-1, 2*Bs(2)-1, 2*Bs(3)-1, N_MAX_COMPONENTS) )


    ! every proc loop over its active heavy data array
    do k = 1, hvy_n

        ! light data id
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )

        ! block wants to refine
        if ( (lgt_block( lgt_id, params%max_treelevel + idx_refine_sts) == +1) ) then

            ! treecode and mesh level
            treecode = lgt_block( lgt_id, 1:params%max_treelevel )
            level    = lgt_block( lgt_id, params%max_treelevel + IDX_MESH_LVL )
            tree_id  = lgt_block( lgt_id, params%max_treelevel + IDX_TREE_ID )

            ! ------------------------------------------------------------------------------------------------------
            ! first: interpolate block data
            ! loop over all data fields
            do dF = 1, size(hvy_block,4)
                ! NOTE: the refinement interpolation acts on the entire block including ghost nodes.
                ! interpolate data
                call prediction_3D(hvy_block(:, :, :, dF, hvy_active(k)), data_predict_fine, params%order_predictor)
                ! save new data, but cut ghost nodes.
                new_data(:, :, :, dF) = data_predict_fine(2*g+1:2*g+1+2*Bs(1)-2, 2*g+1:2*g+1+2*Bs(2)-2, 2*g+1:2*g+1+2*Bs(3)-2)
            end do

            ! ------------------------------------------------------------------------------------------------------
            ! second: split new data and write into new blocks
            !--------------------------
            ! first new block
            ! find a free light id on this rank
            call get_free_local_light_id( params, rank, lgt_block, lgt_free_id)
            call lgt_id_to_hvy_id( free_heavy_id, lgt_free_id, rank, N )
            ! write new light data
            ! old treecode
            lgt_block( lgt_free_id, 1:params%max_treelevel ) = treecode(:)
            ! new treecode one level up - "0" block
            lgt_block( lgt_free_id, level+1 )                = 0
            ! new level + 1
            lgt_block( lgt_free_id, params%max_treelevel + IDX_MESH_LVL ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0
            ! the tree_id is the same as the one of the mother block
            lgt_block( lgt_free_id, params%max_treelevel + IDX_TREE_ID ) = tree_id


            ! save interpolated data, loop over all datafields
            do dF = 1, size(hvy_block,4)
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, free_heavy_id ) = new_data(1:Bs(1), 1:Bs(2), 1:Bs(3), dF)
            end do

            !--------------------------
            ! second new block
            ! find a free light id on this rank
            call get_free_local_light_id( params, rank, lgt_block, lgt_free_id)
            call lgt_id_to_hvy_id( free_heavy_id, lgt_free_id, rank, N )
            ! write new light data
            ! old treecode
            lgt_block( lgt_free_id, 1:params%max_treelevel ) = treecode(:)
            ! new treecode one level up - "0" block
            lgt_block( lgt_free_id, level+1 )                = 1
            ! new level + 1
            lgt_block( lgt_free_id, params%max_treelevel + IDX_MESH_LVL ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0
            ! the tree_id is the same as the one of the mother block
            lgt_block( lgt_free_id, params%max_treelevel + IDX_TREE_ID ) = tree_id


            ! save interpolated data, loop over all datafields
            do dF = 1, size(hvy_block,4)
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, free_heavy_id ) = new_data(1:Bs(1), Bs(2):2*Bs(2)-1, 1:Bs(3), dF)
            end do

            !--------------------------
            ! third new block
            ! find a free light id on this rank
            call get_free_local_light_id( params, rank, lgt_block, lgt_free_id)
            call lgt_id_to_hvy_id( free_heavy_id, lgt_free_id, rank, N )
            ! write new light data
            ! old treecode
            lgt_block( lgt_free_id, 1:params%max_treelevel ) = treecode(:)
            ! new treecode one level up - "0" block
            lgt_block( lgt_free_id, level+1 )                = 2
            ! new level + 1
            lgt_block( lgt_free_id, params%max_treelevel + IDX_MESH_LVL ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0
            ! the tree_id is the same as the one of the mother block
            lgt_block( lgt_free_id, params%max_treelevel + IDX_TREE_ID ) = tree_id


            ! save interpolated data, loop over all datafields
            do dF = 1, size(hvy_block,4)
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, free_heavy_id ) = new_data(Bs(1):2*Bs(1)-1, 1:Bs(2), 1:Bs(3), dF)
            end do

            !--------------------------
            ! fourth new block
            ! find a free light id on this rank
            call get_free_local_light_id( params, rank, lgt_block, lgt_free_id)
            call lgt_id_to_hvy_id( free_heavy_id, lgt_free_id, rank, N )
            ! write new light data
            ! old treecode
            lgt_block( lgt_free_id, 1:params%max_treelevel ) = treecode(:)
            ! new treecode one level up - "0" block
            lgt_block( lgt_free_id, level+1 )                = 3
            ! new level + 1
            lgt_block( lgt_free_id, params%max_treelevel + IDX_MESH_LVL ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0
            ! the tree_id is the same as the one of the mother block
            lgt_block( lgt_free_id, params%max_treelevel + IDX_TREE_ID ) = tree_id


            ! save interpolated data, loop over all datafields
            do dF = 1, size(hvy_block,4)
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, free_heavy_id ) = new_data(Bs(1):2*Bs(1)-1, Bs(2):2*Bs(2)-1, 1:Bs(3), dF)
            end do

            !--------------------------
            ! fifth new block
            ! find a free light id on this rank
            call get_free_local_light_id( params, rank, lgt_block, lgt_free_id)
            call lgt_id_to_hvy_id( free_heavy_id, lgt_free_id, rank, N )
            ! write new light data
            ! old treecode
            lgt_block( lgt_free_id, 1:params%max_treelevel ) = treecode(:)
            ! new treecode one level up - "0" block
            lgt_block( lgt_free_id, level+1 )                = 4
            ! new level + 1
            lgt_block( lgt_free_id, params%max_treelevel + IDX_MESH_LVL ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0
            ! the tree_id is the same as the one of the mother block
            lgt_block( lgt_free_id, params%max_treelevel + IDX_TREE_ID ) = tree_id


            ! save interpolated data, loop over all datafields
            do dF = 1, size(hvy_block,4)
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, free_heavy_id ) = new_data(1:Bs(1), 1:Bs(2), Bs(3):2*Bs(3)-1, dF)
            end do

            !--------------------------
            ! sixth new block
            ! find a free light id on this rank
            call get_free_local_light_id( params, rank, lgt_block, lgt_free_id)
            call lgt_id_to_hvy_id( free_heavy_id, lgt_free_id, rank, N )
            ! write new light data
            ! old treecode
            lgt_block( lgt_free_id, 1:params%max_treelevel ) = treecode(:)
            ! new treecode one level up - "0" block
            lgt_block( lgt_free_id, level+1 )                = 5
            ! new level + 1
            lgt_block( lgt_free_id, params%max_treelevel + IDX_MESH_LVL ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0
            ! the tree_id is the same as the one of the mother block
            lgt_block( lgt_free_id, params%max_treelevel + IDX_TREE_ID ) = tree_id


            ! save interpolated data, loop over all datafields
            do dF = 1, size(hvy_block,4)
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, free_heavy_id ) = new_data(1:Bs(1), Bs(2):2*Bs(2)-1, Bs(3):2*Bs(3)-1, dF)
            end do

            !--------------------------
            ! seventh new block
            ! find a free light id on this rank
            call get_free_local_light_id( params, rank, lgt_block, lgt_free_id)
            call lgt_id_to_hvy_id( free_heavy_id, lgt_free_id, rank, N )
            ! write new light data
            ! old treecode
            lgt_block( lgt_free_id, 1:params%max_treelevel ) = treecode(:)
            ! new treecode one level up - "0" block
            lgt_block( lgt_free_id, level+1 )                = 6
            ! new level + 1
            lgt_block( lgt_free_id, params%max_treelevel + IDX_MESH_LVL ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0
            ! the tree_id is the same as the one of the mother block
            lgt_block( lgt_free_id, params%max_treelevel + IDX_TREE_ID ) = tree_id


            ! save interpolated data, loop over all datafields
            do dF = 1, size(hvy_block,4)
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, free_heavy_id ) = new_data(Bs(1):2*Bs(1)-1, 1:Bs(2), Bs(3):2*Bs(3)-1, dF)
            end do

            !--------------------------
            ! eigth new block
            ! write data on current heavy id
            free_heavy_id = hvy_active(k)
            call hvy_id_to_lgt_id( lgt_free_id, free_heavy_id, rank, N )

            ! write new light data
            ! old treecode
            lgt_block( lgt_free_id, 1:params%max_treelevel ) = treecode(:)
            ! new treecode one level up - "1" block
            lgt_block( lgt_free_id, level+1 )                = 7
            ! new level + 1
            lgt_block( lgt_free_id, params%max_treelevel + IDX_MESH_LVL ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0
            ! the tree_id is the same as the one of the mother block
            lgt_block( lgt_free_id, params%max_treelevel + IDX_TREE_ID ) = tree_id


            ! save interpolated data, loop over all datafields
            do dF = 1, size(hvy_block,4)
                hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, free_heavy_id ) = new_data(Bs(1):2*Bs(1)-1, Bs(2):2*Bs(2)-1, Bs(3):2*Bs(3)-1, dF)
            end do

        end if

    end do

    ! synchronize light data
    call synchronize_lgt_data( params, lgt_block, refinement_status_only=.false. )

end subroutine refinement_execute_3D
