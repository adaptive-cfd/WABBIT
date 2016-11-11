! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: threshold_block.f90
! version: 0.4
! author: msr
!
! thresholding for all blocks
!
! input:    - params, light and heavy data, neighbor list
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 10/11/16 - switch to v0.4
! ********************************************************************************************

subroutine threshold_block( params, block_list, block_data, neighbor_list )

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
    ! neighbor list
    integer(kind=ik), intent(in)        :: neighbor_list(:)

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! loop parameter
    integer(kind=ik)                    :: k, N, dF, i, j

    ! detail
    real(kind=rk)                       :: detail

    ! grid parameter
    integer(kind=ik)                    :: Bs, g

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! interpolation fields
    real(kind=rk), allocatable          :: u1(:,:), u2(:,:), u3(:,:)

    ! light data list for working
    integer(kind=ik)                    :: my_block_list( size(block_list, 1), params%max_treelevel+2)

!---------------------------------------------------------------------------------------------
! interfaces

    interface
        subroutine synchronize_ghosts( params, block_list, block_data, neighbor_list )
            use module_params
            type (type_params), intent(in)              :: params
            integer(kind=ik), intent(in)                :: block_list(:, :)
            real(kind=rk), intent(inout)                :: block_data(:, :, :, :)
            integer(kind=ik), intent(in)                :: neighbor_list(:)
        end subroutine synchronize_ghosts

        subroutine respect_min_max_treelevel( block_list, max_level, min_level)
            use module_params
            integer(kind=ik), intent(inout)  :: block_list(:, :)
            integer(kind=ik), intent(in)     :: max_level, min_level
        end subroutine respect_min_max_treelevel

    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! allocate interpolation fields
    allocate( u1(1:Bs+2*g,1:Bs+2*g), stat=allocate_error )
    allocate( u2(1:Bs+2*g,1:Bs+2*g), stat=allocate_error )
    ! coarsen field are half block size + 1/2
    allocate( u3( 1:(Bs+1)/2 + g , 1:(Bs+1)/2 + g), stat=allocate_error )

    ! set light data list for working, only light data coresponding to proc are not zero
    my_block_list = 0
    my_block_list( rank*N+1: rank*N+N, :) = block_list( rank*N+1: rank*N+N, :)

!---------------------------------------------------------------------------------------------
! main body

    ! ------------------------------------------------------------------------------------
    ! first: synchronize ghost nodes - thresholding on block with ghost nodes
    ! synchronize ghostnodes
    call synchronize_ghosts( params, block_list, block_data, neighbor_list )

    ! ------------------------------------------------------------------------------------
    ! second: clear old refinement status
    ! set status "no refine/coarse" for all active blocks
    do k = 1, size(block_list, 1)
        if ( block_list(k, 1) /= -1 ) then
            block_list(k, params%max_treelevel+2 ) = 0
        end if
    end do

    ! ------------------------------------------------------------------------------------
    ! third: calculate detail and set new refinement status
    ! loop over all heavy data, note: light data need to synchronize after this step
    do k = 1, N

        ! reset detail
        detail = 0.0_rk

        ! block is active
        if ( block_list(rank*N + k , 1) /= -1 ) then

            ! loop over all datafields
            do dF = 2, params%number_data_fields+1

                ! reset interpolation fields
                u1        = block_data( :, :, dF, k)
                u2        = 0.0_rk
                u3        = 0.0_rk

                ! wavelet indicator
                call restriction_2D( u1, u3 )  ! fine, coarse
                call prediction_2D ( u3, u2, params%order_predictor )  ! coarse, fine
                ! calculate deatil
                do i = 1, Bs+2*g
                    do j = 1, Bs+2*g
                        detail = max( detail, sqrt( (u1(i,j)-u2(i,j)) * ( u1(i,j)-u2(i,j)) ) )
                    end do
                end do

            end do

            ! threshold
            if (detail < params%eps) then
                ! coarsen block, -1
                my_block_list( rank*N + k, params%max_treelevel+2 ) = -1
            end if

        end if

    end do

    ! ------------------------------------------------------------------------------------
    ! fourth: synchronize light data
    block_list = 0
    call MPI_Allreduce(my_block_list, block_list, size(block_list,1)*size(block_list,2), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! ------------------------------------------------------------------------------------
    ! fifth: check if block has reached maximal level
    call respect_min_max_treelevel( block_list, params%max_treelevel, params%min_treelevel )

    ! clean up
    deallocate( u1, stat=allocate_error )
    deallocate( u2, stat=allocate_error )
    deallocate( u3, stat=allocate_error )

end subroutine threshold_block
