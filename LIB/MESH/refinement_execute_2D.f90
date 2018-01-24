!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name refinement_execute_2D.f90
!> \version 0.5
!> \author msr
!
!> \brief Refine mesh (2D version). All cpu loop over their heavy data and check if the refinement
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
!! 08/11/16 - switch to v0.4, split old interpolate_mesh subroutine into two refine/coarsen
!!            subroutines
!! 09/06/17 - speed up light data synchronization, uses 8bit integers and send/receive changed data only
!!
! ********************************************************************************************

subroutine refinement_execute_2D( params, lgt_block, hvy_block, hvy_active, hvy_n )

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
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! loop variables
    integer(kind=ik)                    :: k, N, dF
    ! light id start
    integer(kind=ik)                    :: my_light_start

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! grid parameter
    integer(kind=ik)                    :: Bs, g
    ! data fields for interpolation
    real(kind=rk), allocatable          :: new_data(:,:,:), data_predict_coarse(:,:), data_predict_fine(:,:)
    ! free light/heavy data id
    integer(kind=ik)                    :: free_heavy_id

    ! treecode varaible
    integer(kind=1)                    :: treecode(params%max_treelevel)
    ! mesh level
    integer(kind=1)                    :: level

    ! light data list for working
    integer(kind=1), allocatable       :: my_lgt_block(:,:)

    ! send/receive buffer for data synchronization
    integer(kind=1), allocatable        :: my_lgt_block_send_buffer(:,:), my_lgt_block_receive_buffer(:,:)

    ! maximum heavy id, use to synchronize reduced light data array, sum of heavy blocks, start of send buffer
    integer(kind=ik)                    :: heavy_max, block_sum, buffer_start
    ! list of max heavy ids, use to build send/receive buffer
    integer(kind=ik)                    :: proc_heavy_max(params%number_procs), my_proc_heavy_max(params%number_procs)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks

    ! set MPI parameter
    rank  = params%rank

    ! light data start line
    my_light_start = rank*N

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! data fields for interpolation
    ! coarse: current data, fine: new (refine) data, new_data: gather all refined data for all data fields
    allocate( data_predict_fine(Bs+2*g+(Bs+2*g-1),Bs+2*g+(Bs+2*g-1)) )
    ! allocate( data_predict_fine(2*Bs-1, 2*Bs-1) )
    ! allocate( data_predict_coarse(Bs, Bs) )
    allocate( data_predict_coarse(Bs+2*g, Bs+2*g) )
    allocate( new_data(2*Bs-1, 2*Bs-1, params%number_data_fields) )

    ! allocate lgt data working array
    allocate( my_lgt_block(N, params%max_treelevel+2 ) )
    ! set light data list for working, only light data coresponding to proc
    my_lgt_block = int(lgt_block( my_light_start+1: my_light_start+N, :),kind=1)

    ! reset max heavy id
    heavy_max = 0

!---------------------------------------------------------------------------------------------
! main body

    ! every proc loop over his active heavy data array
    ! note: in light data working array heavy id = block id
    do k = 1, hvy_n

        ! save heavy id, if new block id is larger than old one
        heavy_max = max(heavy_max, hvy_active(k))

        ! block wants to refine
        if ( (my_lgt_block( hvy_active(k), params%max_treelevel+2) == 1) ) then

            ! treecode and mesh level
            treecode = my_lgt_block( hvy_active(k), 1:params%max_treelevel )
            level    = my_lgt_block( hvy_active(k), params%max_treelevel+1 )

            ! ------------------------------------------------------------------------------------------------------
            ! first: interpolate block data
            ! loop over all data fields
            do dF = 1, params%number_data_fields
                ! NOTE: the refinement interpolation acts on the blocks interior
                ! nodes and ignores ghost nodes.
                ! data_predict_coarse = hvy_block(g+1:Bs+g, g+1:Bs+g, dF, hvy_active(k) )
                data_predict_coarse = hvy_block(:, :, dF, hvy_active(k) )
                ! reset data
                data_predict_fine   = 9.0e9_rk
                ! interpolate data
                call prediction_2D(data_predict_coarse, data_predict_fine, params%order_predictor)
                ! save new data
                new_data(:,:,dF) = data_predict_fine(2*g+1:Bs+2*g+(Bs+2*g-1)-2*g,2*g+1:Bs+2*g+(Bs+2*g-1)-2*g)
            end do

            ! ------------------------------------------------------------------------------------------------------
            ! second: split new data and write into new blocks
            !--------------------------
            ! first new block
            ! find free heavy id, use free light id subroutine with reduced light data list for this
            call get_free_light_id( free_heavy_id, int(my_lgt_block( : , 1 ), kind=ik ), N )

            ! write new light data
            ! old treecode
            my_lgt_block( free_heavy_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "0" block
            my_lgt_block( free_heavy_id, level+1 )                = 0
            ! new level + 1
            my_lgt_block( free_heavy_id, params%max_treelevel+1 ) = level+1_1
            ! reset refinement status
            my_lgt_block( free_heavy_id, params%max_treelevel+2 ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(1:Bs, 1:Bs, dF)
            end do

            !--------------------------
            ! second new block
            ! find free heavy id, use free light id subroutine with reduced light data list for this
            call get_free_light_id( free_heavy_id, int(my_lgt_block( : , 1 ), kind=ik ), N )

            ! write new light data
            ! old treecode
            my_lgt_block( free_heavy_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "1" block
            my_lgt_block( free_heavy_id, level+1 )                = 1
            ! new level + 1
            my_lgt_block( free_heavy_id, params%max_treelevel+1 ) = level+1_1
            ! reset refinement status
            my_lgt_block( free_heavy_id, params%max_treelevel+2 ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(1:Bs, Bs:2*Bs-1, dF)
            end do

            !--------------------------
            ! third new block
            ! find free heavy id, use free light id subroutine with reduced light data list for this
            call get_free_light_id( free_heavy_id, int(my_lgt_block( : , 1 ), kind=ik ), N )

            ! write new light data
            ! old treecode
            my_lgt_block( free_heavy_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "1" block
            my_lgt_block( free_heavy_id, level+1 )                = 2
            ! new level + 1
            my_lgt_block( free_heavy_id, params%max_treelevel+1 ) = level+1_1
            ! reset refinement status
            my_lgt_block( free_heavy_id, params%max_treelevel+2 ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(Bs:2*Bs-1, 1:Bs, dF)
            end do

            ! save heavy id, if new block id is larger than old one
            ! note: heavy id of third block is always larger than id from block 1,2
            heavy_max = max(heavy_max, free_heavy_id)

            !--------------------------
            ! fourth new block
            ! write data on current heavy id
            free_heavy_id = hvy_active(k)

            ! write new light data
            ! old treecode
            my_lgt_block( free_heavy_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "1" block
            my_lgt_block( free_heavy_id, level+1 )                = 3
            ! new level + 1
            my_lgt_block( free_heavy_id, params%max_treelevel+1 ) = level+1_1
            ! reset refinement status
            my_lgt_block( free_heavy_id, params%max_treelevel+2 ) = 0

            ! save interpolated data, loop over all datafields
            do dF = 1, params%number_data_fields
                hvy_block( g+1:Bs+g, g+1:Bs+g, dF, free_heavy_id ) = new_data(Bs:2*Bs-1, Bs:2*Bs-1, dF)
            end do

        end if

    end do

    ! set array for max heavy ids
    my_proc_heavy_max = 0
    my_proc_heavy_max(rank+1) = heavy_max

    ! synchronize array
    call MPI_Allreduce(my_proc_heavy_max, proc_heavy_max, size(proc_heavy_max,1), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! for readability, calc sum of all max heavy ids
    block_sum = sum(proc_heavy_max)

    ! now we can allocate send/receive buffer arrays
    allocate( my_lgt_block_send_buffer( block_sum, size(lgt_block,2) ), my_lgt_block_receive_buffer( block_sum, size(lgt_block,2) ) )

    ! reset send buffer
    my_lgt_block_send_buffer = 0
    buffer_start = sum(proc_heavy_max(1:rank))
    my_lgt_block_send_buffer( buffer_start+1 : buffer_start+heavy_max, : ) = my_lgt_block( 1 : heavy_max, :)

    ! synchronize light data
    call MPI_Allreduce(my_lgt_block_send_buffer, my_lgt_block_receive_buffer, block_sum*size(lgt_block,2), MPI_INTEGER1, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! write synchronized light data
    ! loop over number of procs and reset lgt_block array
    do k = 1, params%number_procs
        ! proc k-1 has send data
        if ( proc_heavy_max(k) /= 0 ) then
            ! write received light data
            lgt_block( (k-1)*N+1 : (k-1)*N + proc_heavy_max(k), : ) =  my_lgt_block_receive_buffer( sum(proc_heavy_max(1:k-1))+1 : sum(proc_heavy_max(1:k-1))+proc_heavy_max(k), : )
        else
            ! nothing to do
        end if
    end do

    ! clean up
    deallocate( data_predict_fine )
    deallocate( data_predict_coarse )
    deallocate( new_data )
    deallocate( my_lgt_block_send_buffer, my_lgt_block_receive_buffer, my_lgt_block )

end subroutine refinement_execute_2D
