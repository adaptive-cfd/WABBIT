! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: refine_mesh.f90
! version: 0.4
! author: msr
!
! refine the mesh:
! every proc work on his own heavy data and change the corresponding light data
! after this: synchronize light data array
!
! input:    - params, light and heavy data
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4, split old interpolate_mesh subroutine into two refine/coarsen
!            subroutines
! ********************************************************************************************

subroutine refine_mesh( params, block_list, block_data )

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! interpolation routines
    use module_interpolation

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)

    ! loop variables
    integer(kind=ik)                    :: k, N, dF, i
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
    ! new coordinates vectors
    real(kind=rk), allocatable          :: new_coord_x(:), new_coord_y(:)

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! free light/heavy data id
    integer(kind=ik)                    :: free_light_id, free_heavy_id

    ! treecode varaible
    integer(kind=ik)                    :: treecode(params%max_treelevel)
    ! mesh level
    integer(kind=ik)                    :: level

    ! light data list for working
    integer(kind=ik)                    :: my_block_list( size(block_list, 1), params%max_treelevel+2)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! light data start line
    my_light_start = rank*N

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! data fields for interpolation
    ! coarse: current data, fine: new (refine) data, new_data: gather all refined data for all data fields
    allocate( data_predict_fine(2*Bs-1, 2*Bs-1), stat=allocate_error )
    allocate( data_predict_coarse(Bs, Bs), stat=allocate_error )
    allocate( new_data(2*Bs-1, 2*Bs-1, params%number_data_fields), stat=allocate_error )

    ! new coordinates vectors
    allocate( new_coord_x(Bs), stat=allocate_error )
    allocate( new_coord_y(Bs), stat=allocate_error )

    ! set light data list for working, only light data coresponding to proc are not zero
    my_block_list = 0
    my_block_list( my_light_start+1: my_light_start+N, :) = block_list( my_light_start+1: my_light_start+N, :)

!---------------------------------------------------------------------------------------------
! main body

    ! every proc loop over the light data array corresponding to his heavy data
    do k = 1, N

        ! block is active and wants to refine
        if ( (my_block_list( my_light_start + k, 1) /= -1) .and. (my_block_list( my_light_start + k, params%max_treelevel+2) == 1) ) then

            ! treecode and mesh level
            treecode = my_block_list( my_light_start + k, 1:params%max_treelevel )
            level    = my_block_list( my_light_start + k, params%max_treelevel+1 )

            ! ------------------------------------------------------------------------------------------------------
            ! first: interpolate block data
            ! loop over all data fields
            do dF = 2, params%number_data_fields+1
                ! reset data
                data_predict_coarse = block_data(k, g+1:Bs*g, g+1:Bs*g, dF )
                data_predict_fine   = 9.0e9_rk
                ! interpolate data
                call prediction_2D(data_predict_coarse, data_predict_fine, params%order_predictor)
                ! save new data
                new_data(:,:,dF-1) = data_predict_fine
            end do

            ! ------------------------------------------------------------------------------------------------------
            ! second: split new data and write into new blocks
            !--------------------------
            ! first new block
            ! find free light id, work only on corresponding light data, so returned id is the heavy id
            call get_free_light_id( free_heavy_id, my_block_list( my_light_start+1 : my_light_start+1+N , 1 ), N )
            ! calculate light id
            free_light_id = my_light_start + free_heavy_id

            ! write new light data
            ! old treecode
            my_block_list( free_light_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "0" block
            my_block_list( free_light_id, level+1 )                = 0
            ! new level + 1
            my_block_list( free_light_id, params%max_treelevel+1 ) = level+1
            ! reset refinement status
            my_block_list( free_light_id, params%max_treelevel+2 ) = 0

            ! interpolate new coordinates
            new_coord_x(1:Bs:2) = block_data( k, 1, 1:(Bs-1)/2+1, 1 )
            new_coord_y(1:Bs:2) = block_data( k, 2, 1:(Bs-1)/2+1, 1 )
            do i = 2, Bs, 2
                new_coord_x(i)  = ( new_coord_x(i-1) + new_coord_x(i+1) ) / 2.0_rk
                new_coord_y(i)  = ( new_coord_y(i-1) + new_coord_y(i+1) ) / 2.0_rk
            end do

            ! save coordinates
            block_data( free_heavy_id, 1, 1:Bs, 1 ) = new_coord_x
            block_data( free_heavy_id, 2, 1:Bs, 1 ) = new_coord_y

            ! save interpolated data, loop over all datafields
            do dF = 2, params%number_data_fields+1
                block_data( free_heavy_id, g+1:Bs+g, g+1:Bs+g, dF ) = new_data(1:Bs, 1:Bs, dF-1)
            end do

            !--------------------------
            ! second new block
            ! find free light id, work only on corresponding light data, so returned id is the heavy id
            call get_free_light_id( free_heavy_id, my_block_list( my_light_start+1 : my_light_start+1+N , 1 ), N )
            ! calculate light id
            free_light_id = my_light_start + free_heavy_id

            ! write new light data
            ! old treecode
            my_block_list( free_light_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "1" block
            my_block_list( free_light_id, level+1 )                = 1
            ! new level + 1
            my_block_list( free_light_id, params%max_treelevel+1 ) = level+1
            ! reset refinement status
            my_block_list( free_light_id, params%max_treelevel+2 ) = 0

            ! interpolate new coordinates
            new_coord_x(1:Bs:2) = block_data( k, 1, (Bs-1)/2+1:Bs, 1 )
            new_coord_y(1:Bs:2) = block_data( k, 2, 1:(Bs-1)/2+1, 1 )
            do i = 2, Bs, 2
                new_coord_x(i)  = ( new_coord_x(i-1) + new_coord_x(i+1) ) / 2.0_rk
                new_coord_y(i)  = ( new_coord_y(i-1) + new_coord_y(i+1) ) / 2.0_rk
            end do

            ! save coordinates
            block_data( free_heavy_id, 1, 1:Bs, 1 ) = new_coord_x
            block_data( free_heavy_id, 2, 1:Bs, 1 ) = new_coord_y

            ! save interpolated data, loop over all datafields
            do dF = 2, params%number_data_fields+1
                block_data( free_heavy_id, g+1:Bs+g, g+1:Bs+g, dF ) = new_data(1:Bs, Bs:2*Bs-1, dF-1)
            end do

            !--------------------------
            ! third new block
            ! find free light id, work only on corresponding light data, so returned id is the heavy id
            call get_free_light_id( free_heavy_id, my_block_list( my_light_start+1 : my_light_start+1+N , 1 ), N )
            ! calculate light id
            free_light_id = my_light_start + free_heavy_id

            ! write new light data
            ! old treecode
            my_block_list( free_light_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "1" block
            my_block_list( free_light_id, level+1 )                = 2
            ! new level + 1
            my_block_list( free_light_id, params%max_treelevel+1 ) = level+1
            ! reset refinement status
            my_block_list( free_light_id, params%max_treelevel+2 ) = 0

            ! interpolate new coordinates
            new_coord_x(1:Bs:2) = block_data( k, 1, 1:(Bs-1)/2+1, 1 )
            new_coord_y(1:Bs:2) = block_data( k, 2, (Bs-1)/2+1:Bs, 1 )
            do i = 2, Bs, 2
                new_coord_x(i)  = ( new_coord_x(i-1) + new_coord_x(i+1) ) / 2.0_rk
                new_coord_y(i)  = ( new_coord_y(i-1) + new_coord_y(i+1) ) / 2.0_rk
            end do

            ! save coordinates
            block_data( free_heavy_id, 1, 1:Bs, 1 ) = new_coord_x
            block_data( free_heavy_id, 2, 1:Bs, 1 ) = new_coord_y

            ! save interpolated data, loop over all datafields
            do dF = 2, params%number_data_fields+1
                block_data( free_heavy_id, g+1:Bs+g, g+1:Bs+g, dF ) = new_data(Bs:2*Bs-1, 1:Bs, dF-1)
            end do

            !--------------------------
            ! fourth new block
            ! write data on current heavy id
            free_heavy_id = k
            ! calculate light id
            free_light_id = my_light_start + free_heavy_id

            ! write new light data
            ! old treecode
            my_block_list( free_light_id, 1:params%max_treelevel ) = treecode
            ! new treecode one level up - "1" block
            my_block_list( free_light_id, level+1 )                = 3
            ! new level + 1
            my_block_list( free_light_id, params%max_treelevel+1 ) = level+1
            ! reset refinement status
            my_block_list( free_light_id, params%max_treelevel+2 ) = 0

            ! interpolate new coordinates
            new_coord_x(1:Bs:2) = block_data( k, 1, (Bs-1)/2+1:Bs, 1 )
            new_coord_y(1:Bs:2) = block_data( k, 2, (Bs-1)/2+1:Bs, 1 )
            do i = 2, Bs, 2
                new_coord_x(i)  = ( new_coord_x(i-1) + new_coord_x(i+1) ) / 2.0_rk
                new_coord_y(i)  = ( new_coord_y(i-1) + new_coord_y(i+1) ) / 2.0_rk
            end do

            ! save coordinates
            block_data( free_heavy_id, 1, 1:Bs, 1 ) = new_coord_x
            block_data( free_heavy_id, 2, 1:Bs, 1 ) = new_coord_y

            ! save interpolated data, loop over all datafields
            do dF = 2, params%number_data_fields+1
                block_data( free_heavy_id, g+1:Bs+g, g+1:Bs+g, dF ) = new_data(Bs:2*Bs-1, Bs:2*Bs-1, dF-1)
            end do

        end if

    end do

    ! synchronize light data
    block_list = 0
    call MPI_Allreduce(my_block_list, block_list, size(block_list,1)*size(block_list,2), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! clean up
    deallocate( data_predict_fine, stat=allocate_error )
    deallocate( data_predict_coarse, stat=allocate_error )
    deallocate( new_data, stat=allocate_error )
    deallocate( new_coord_x, stat=allocate_error )
    deallocate( new_coord_y, stat=allocate_error )

end subroutine refine_mesh
