
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

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

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
    integer(kind=ik)                    :: Bs, g
    ! data fields for interpolation
    real(kind=rk), allocatable, save    :: new_data(:,:,:,:), data_predict_coarse(:,:,:), data_predict_fine(:,:,:)
    ! free light/heavy data id
    integer(kind=ik)                    :: lgt_free_id, free_heavy_id
    ! treecode varaible
    integer(kind=ik)                    :: treecode(params%max_treelevel)
    ! mesh level
    integer(kind=ik)                    :: level


!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks
    ! set MPI parameter
    rank = params%rank
    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! NOTE: the predictor for the refinement acts on the extended blocks i.e. it
    ! includes the ghost nodes layer. Therefore, you MUST call sync_ghosts before this routine.
    ! The datafield for prediction is one level up, i.e. it contains Bs+g + (Bs+2g-1) points
    if (.not.allocated(data_predict_fine)) allocate( data_predict_fine(2*(Bs+2*g)-1, 2*(Bs+2*g)-1, 2*(Bs+2*g)-1) )
    ! the coarse field has the same size as the block.
    if (.not.allocated(data_predict_coarse)) allocate( data_predict_coarse(Bs+2*g, Bs+2*g, Bs+2*g) )
    ! the new_data field holds the interior part of the new, refined block (which
    ! will become four blocks), without the ghost nodes.
    if (.not.allocated(new_data)) allocate( new_data(2*Bs-1, 2*Bs-1, 2*Bs-1, params%number_data_fields) )


!---------------------------------------------------------------------------------------------
! main body

    ! every proc loop over his active heavy data array
    do k = 1, hvy_n

        ! light data id
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )

        ! block wants to refine
        if ( (lgt_block( lgt_id, params%max_treelevel + idx_refine_sts) == +1) ) then

            ! treecode and mesh level
            treecode = lgt_block( lgt_id, 1:params%max_treelevel )
            level    = lgt_block( lgt_id, params%max_treelevel + idx_mesh_lvl )

            ! ------------------------------------------------------------------------------------------------------
            ! first: interpolate block data
            ! loop over all data fields
            do dF = 1, params%number_data_fields
                data_predict_coarse = hvy_block(:, :, :, dF, hvy_active(k) )
                ! interpolate data
                call prediction_3D(data_predict_coarse, data_predict_fine, params%order_predictor)
                ! save new data, but cut ghost nodes.
                new_data(:, :, :, dF) = data_predict_fine(2*g+1:2*g+1+2*Bs-2, 2*g+1:2*g+1+2*Bs-2, 2*g+1:2*g+1+2*Bs-2)
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
            lgt_block( lgt_free_id, params%max_treelevel + idx_mesh_lvl ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(1:Bs, 1:Bs, 1:Bs, dF)
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
            lgt_block( lgt_free_id, params%max_treelevel + idx_mesh_lvl ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(1:Bs, Bs:2*Bs-1, 1:Bs, dF)
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
            lgt_block( lgt_free_id, params%max_treelevel + idx_mesh_lvl ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(Bs:2*Bs-1, 1:Bs, 1:Bs, dF)
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
            lgt_block( lgt_free_id, params%max_treelevel + idx_mesh_lvl ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(Bs:2*Bs-1, Bs:2*Bs-1, 1:Bs, dF)
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
            lgt_block( lgt_free_id, params%max_treelevel + idx_mesh_lvl ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(1:Bs, 1:Bs, Bs:2*Bs-1, dF)
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
            lgt_block( lgt_free_id, params%max_treelevel + idx_mesh_lvl ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(1:Bs, Bs:2*Bs-1, Bs:2*Bs-1, dF)
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
            lgt_block( lgt_free_id, params%max_treelevel + idx_mesh_lvl ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(Bs:2*Bs-1, 1:Bs, Bs:2*Bs-1, dF)
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
            lgt_block( lgt_free_id, params%max_treelevel + idx_mesh_lvl ) = level+1
            ! reset refinement status
            lgt_block( lgt_free_id, params%max_treelevel + idx_refine_sts ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(Bs:2*Bs-1, Bs:2*Bs-1, Bs:2*Bs-1, dF)
            end do

        end if

    end do

    ! synchronize light data
    call synchronize_lgt_data( params, lgt_block, refinement_status_only=.false. )

end subroutine refinement_execute_3D
