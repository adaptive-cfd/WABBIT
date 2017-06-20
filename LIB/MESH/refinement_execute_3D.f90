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
!! level and create four new blocks, each carrying a part of the interpolated data.
!! As all CPU first work individually, the light data array is synced.
!
!> \note The interpolation (or prediction) operator here is applied to a block EXCLUDING
!! any ghost nodes.
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
    ! light id start
    integer(kind=ik)                    :: my_light_start

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! grid parameter
    integer(kind=ik)                    :: Bs, g
    ! data fields for interpolation
    real(kind=rk), allocatable          :: new_data(:,:,:,:), data_predict_coarse(:,:,:), data_predict_fine(:,:,:)


    ! free light/heavy data id
    integer(kind=ik)                    :: free_light_id, free_heavy_id

    ! treecode varaible
    integer(kind=ik)                    :: treecode(params%max_treelevel)
    ! mesh level
    integer(kind=ik)                    :: level

    ! light data list for working
    integer(kind=ik)                    :: my_lgt_block( size(lgt_block, 1), params%max_treelevel+2)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks

    ! set MPI parameter
    rank         = params%rank

    ! light data start line
    my_light_start = rank*N

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! data fields for interpolation
    ! coarse: current data, fine: new (refine) data, new_data: gather all refined data for all data fields
    allocate( data_predict_fine(2*Bs-1, 2*Bs-1, 2*Bs-1)  )
    allocate( data_predict_coarse(Bs, Bs, Bs)  )
    allocate( new_data(2*Bs-1, 2*Bs-1, 2*Bs-1, params%number_data_fields)  )
    ! set light data list for working, only light data coresponding to proc are not zero
    my_lgt_block = 0
    my_lgt_block( my_light_start+1: my_light_start+N, :) = lgt_block( my_light_start+1: my_light_start+N, :)

!---------------------------------------------------------------------------------------------
! main body

    ! every proc loop over his active heavy data array
    do k = 1, hvy_n

        ! light data id
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )

        ! block wants to refine
        if ( (my_lgt_block( lgt_id, params%max_treelevel+2) == 1) ) then

            ! treecode and mesh level
            treecode = my_lgt_block( lgt_id, 1:params%max_treelevel )
            level    = my_lgt_block( lgt_id, params%max_treelevel+1 )

            ! ------------------------------------------------------------------------------------------------------
            ! first: interpolate block data
            ! loop over all data fields
            do dF = 1, params%number_data_fields
                ! reset data
                data_predict_coarse = hvy_block(g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, hvy_active(k) )
                !data_predict_fine   = 9.0e9_rk
                ! interpolate data
                call prediction_3D(data_predict_coarse, data_predict_fine, params%order_predictor)
                ! save new data
                new_data(:,:,:,dF) = data_predict_fine
            end do

            ! ------------------------------------------------------------------------------------------------------
            ! second: split new data and write into new blocks
            !--------------------------
            ! first new block
            ! find free heavy id, use free light id subroutine with reduced light data list for this
            call get_free_light_id( free_heavy_id, my_lgt_block( my_light_start+1 : my_light_start+N , 1 ), N )
            ! calculate light id
            call hvy_id_to_lgt_id( free_light_id, free_heavy_id, rank, N )

            ! write new light data
            ! old treecode
            my_lgt_block( free_light_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "0" block
            my_lgt_block( free_light_id, level+1 )                = 0
            ! new level + 1
            my_lgt_block( free_light_id, params%max_treelevel+1 ) = level+1
            ! reset refinement status
            my_lgt_block( free_light_id, params%max_treelevel+2 ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(1:Bs, 1:Bs, 1:Bs, dF)
            end do

            !--------------------------
            ! second new block
            ! find free heavy id, use free light id subroutine with reduced light data list for this
            call get_free_light_id( free_heavy_id, my_lgt_block( my_light_start+1 : my_light_start+N , 1 ), N )
            ! calculate light id
            call hvy_id_to_lgt_id( free_light_id, free_heavy_id, rank, N )

            ! write new light data
            ! old treecode
            my_lgt_block( free_light_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "0" block
            my_lgt_block( free_light_id, level+1 )                = 1
            ! new level + 1
            my_lgt_block( free_light_id, params%max_treelevel+1 ) = level+1
            ! reset refinement status
            my_lgt_block( free_light_id, params%max_treelevel+2 ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(1:Bs, Bs:2*Bs-1, 1:Bs, dF)
            end do

            !--------------------------
            ! third new block
            ! find free heavy id, use free light id subroutine with reduced light data list for this
            call get_free_light_id( free_heavy_id, my_lgt_block( my_light_start+1 : my_light_start+N , 1 ), N )
            ! calculate light id
            call hvy_id_to_lgt_id( free_light_id, free_heavy_id, rank, N )

            ! write new light data
            ! old treecode
            my_lgt_block( free_light_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "0" block
            my_lgt_block( free_light_id, level+1 )                = 2
            ! new level + 1
            my_lgt_block( free_light_id, params%max_treelevel+1 ) = level+1
            ! reset refinement status
            my_lgt_block( free_light_id, params%max_treelevel+2 ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(Bs:2*Bs-1, 1:Bs, 1:Bs, dF)
            end do

            !--------------------------
            ! fourth new block
            ! find free heavy id, use free light id subroutine with reduced light data list for this
            call get_free_light_id( free_heavy_id, my_lgt_block( my_light_start+1 : my_light_start+N , 1 ), N )
            ! calculate light id
            call hvy_id_to_lgt_id( free_light_id, free_heavy_id, rank, N )

            ! write new light data
            ! old treecode
            my_lgt_block( free_light_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "0" block
            my_lgt_block( free_light_id, level+1 )                = 3
            ! new level + 1
            my_lgt_block( free_light_id, params%max_treelevel+1 ) = level+1
            ! reset refinement status
            my_lgt_block( free_light_id, params%max_treelevel+2 ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(Bs:2*Bs-1, Bs:2*Bs-1, 1:Bs, dF)
            end do

            !--------------------------
            ! fifth new block
            ! find free heavy id, use free light id subroutine with reduced light data list for this
            call get_free_light_id( free_heavy_id, my_lgt_block( my_light_start+1 : my_light_start+N , 1 ), N )
            ! calculate light id
            call hvy_id_to_lgt_id( free_light_id, free_heavy_id, rank, N )

            ! write new light data
            ! old treecode
            my_lgt_block( free_light_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "0" block
            my_lgt_block( free_light_id, level+1 )                = 4
            ! new level + 1
            my_lgt_block( free_light_id, params%max_treelevel+1 ) = level+1
            ! reset refinement status
            my_lgt_block( free_light_id, params%max_treelevel+2 ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(1:Bs, 1:Bs, Bs:2*Bs-1, dF)
            end do

            !--------------------------
            ! sixth new block
            ! find free heavy id, use free light id subroutine with reduced light data list for this
            call get_free_light_id( free_heavy_id, my_lgt_block( my_light_start+1 : my_light_start+N , 1 ), N )
            ! calculate light id
            call hvy_id_to_lgt_id( free_light_id, free_heavy_id, rank, N )

            ! write new light data
            ! old treecode
            my_lgt_block( free_light_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "0" block
            my_lgt_block( free_light_id, level+1 )                = 5
            ! new level + 1
            my_lgt_block( free_light_id, params%max_treelevel+1 ) = level+1
            ! reset refinement status
            my_lgt_block( free_light_id, params%max_treelevel+2 ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(1:Bs, Bs:2*Bs-1, Bs:2*Bs-1, dF)
            end do

            !--------------------------
            ! seventh new block
            ! find free heavy id, use free light id subroutine with reduced light data list for this
            call get_free_light_id( free_heavy_id, my_lgt_block( my_light_start+1 : my_light_start+N , 1 ), N )
            ! calculate light id
            call hvy_id_to_lgt_id( free_light_id, free_heavy_id, rank, N )

            ! write new light data
            ! old treecode
            my_lgt_block( free_light_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "0" block
            my_lgt_block( free_light_id, level+1 )                = 6
            ! new level + 1
            my_lgt_block( free_light_id, params%max_treelevel+1 ) = level+1
            ! reset refinement status
            my_lgt_block( free_light_id, params%max_treelevel+2 ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(Bs:2*Bs-1, 1:Bs, Bs:2*Bs-1, dF)
            end do

            !--------------------------
            ! eigth new block
            ! write data on current heavy id
            free_heavy_id = hvy_active(k)
            ! calculate light id
            call hvy_id_to_lgt_id( free_light_id, free_heavy_id, rank, N )

            ! write new light data
            ! old treecode
            my_lgt_block( free_light_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "1" block
            my_lgt_block( free_light_id, level+1 )                = 7
            ! new level + 1
            my_lgt_block( free_light_id, params%max_treelevel+1 ) = level+1
            ! reset refinement status
            my_lgt_block( free_light_id, params%max_treelevel+2 ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(Bs:2*Bs-1, Bs:2*Bs-1, Bs:2*Bs-1, dF)
            end do

        end if

    end do

    ! synchronize light data
    lgt_block = 0
    call MPI_Allreduce(my_lgt_block, lgt_block, size(lgt_block,1)*size(lgt_block,2), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! clean up
    deallocate( data_predict_fine  )
    deallocate( data_predict_coarse  )
    deallocate( new_data  )

end subroutine refinement_execute_3D
