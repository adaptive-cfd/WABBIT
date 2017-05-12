!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name threshold_block.f90
!> \version 0.4
!> \author msr
!
!> \brief thresholding for all blocks
!
!> 
!! The block thresholding is done with the restriction/prediction operators acting on the
!! entire block, INCLUDING GHOST NODES. Ghost node syncing is performed here. \n
!! 
!!
!! input:    - params, light and heavy data, neighbor list \n
!! output:   - light and heavy data arrays \n
!!
!!
!! = log ======================================================================================
!! \n
!! 10/11/16 - switch to v0.4
! ********************************************************************************************
!> \image html threshold.png width=400

subroutine threshold_block( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

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
    !> neighbor list
    integer(kind=ik), intent(in)        :: hvy_neighbor(:, :)

    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! loop parameter
    integer(kind=ik)                    :: k, dF, i, j, l, N, lgt_id

    ! detail
    real(kind=rk)                       :: detail

    ! grid parameter
    integer(kind=ik)                    :: Bs, g

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! interpolation fields
    real(kind=rk), allocatable          :: u1(:,:,:), u2(:,:,:), u3(:,:,:)

    ! light data list for working
    integer(kind=ik)                    :: my_lgt_block( size(lgt_block, 1), params%max_treelevel+2)

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! block number
    N = params%number_blocks

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! set MPI parameter
    rank         = params%rank

    ! allocate interpolation fields
    allocate( u1( 1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g ), stat=allocate_error )
    call check_allocation(allocate_error)
    allocate( u2( 1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g ), stat=allocate_error )
    call check_allocation(allocate_error)
    ! coarsened field is half block size + 1/2
    allocate( u3( 1:(Bs+1)/2 + g , 1:(Bs+1)/2 + g, 1:(Bs+1)/2 + g), stat=allocate_error )
    call check_allocation(allocate_error)

    ! set light data list for working, only light data coresponding to proc are not zero
    my_lgt_block = 0
    my_lgt_block( rank*N+1: rank*N+N, :) = lgt_block( rank*N+1: rank*N+N, :)

!---------------------------------------------------------------------------------------------
! main body

    ! ------------------------------------------------------------------------------------
    ! first: synchronize ghost nodes - thresholding on block with ghost nodes

    ! synchronize ghostnodes
    call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    ! start time
    sub_t0 = MPI_Wtime()

    ! ------------------------------------------------------------------------------------
    ! second: clear old refinement status
    ! set status "no refine/coarse" for all active blocks
    do k = 1, lgt_n
        lgt_block( lgt_active(k), params%max_treelevel+2 ) = 0
    end do

    ! ------------------------------------------------------------------------------------
    ! third: calculate detail and set new refinement status
    ! loop over all active heavy data, note: light data need to synchronize after this step
    do k = 1, hvy_n

        ! calculate light id
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )

        ! reset detail
        detail = 0.0_rk

        ! loop over all datafields
        do dF = 2, params%number_data_fields+1

            ! reset interpolation fields
            u1        = hvy_block( :, :, :, dF, hvy_active(k) )
            !u2        = 0.0_rk
            !u3        = 0.0_rk

            ! wavelet indicator
            if ( params%threeD_case ) then

                ! 3D:
                call restriction_3D( u1, u3 )  ! fine, coarse
                call prediction_3D ( u3, u2, params%order_predictor )  ! coarse, fine

                ! calculate deatil
                do i = 1, Bs+2*g
                    do j = 1, Bs+2*g
                        do l = 1, Bs+2*g
                            detail = max( detail, sqrt( (u1(i,j,l)-u2(i,j,l)) * ( u1(i,j,l)-u2(i,j,l)) ) )
                        end do
                    end do
                end do

            else

                ! 2D:
                call restriction_2D( u1(:,:,1), u3(:,:,1) )  ! fine, coarse
                call prediction_2D ( u3(:,:,1), u2(:,:,1), params%order_predictor )  ! coarse, fine

                ! calculate deatil
                do i = 1, Bs+2*g
                    do j = 1, Bs+2*g
                        detail = max( detail, sqrt( (u1(i,j,1)-u2(i,j,1)) * ( u1(i,j,1)-u2(i,j,1)) ) )
                    end do
                end do

            end if

        end do

        ! threshold
        if (detail < params%eps) then
            ! coarsen block, -1
            my_lgt_block( lgt_id, params%max_treelevel+2 ) = -1
        end if

    end do

    ! ------------------------------------------------------------------------------------
    ! fourth: synchronize light data
    lgt_block = 0
    call MPI_Allreduce(my_lgt_block, lgt_block, size(lgt_block,1)*size(lgt_block,2), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! ------------------------------------------------------------------------------------
    ! fifth: check if block has reached maximal level
    call respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n )

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "treshold_block (w/o ghost synch.)" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "treshold_block (w/o ghost synch.)"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

    ! clean up
    deallocate( u1, stat=allocate_error )
    deallocate( u2, stat=allocate_error )
    deallocate( u3, stat=allocate_error )

end subroutine threshold_block
