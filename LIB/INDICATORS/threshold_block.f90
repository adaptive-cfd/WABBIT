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
!! entire block \n
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
!> \image html threshold.svg width=400

subroutine threshold_block( params, lgt_block, hvy_block, hvy_active, hvy_n)

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
    ! interpolation fields
    real(kind=rk), allocatable          :: u1(:,:,:), u2(:,:,:), u3(:,:,:)
    ! light data (refinement status column) list for working
    integer(kind=1), allocatable       :: my_refinement_status(:)

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1, time_sum

    ! send/receive buffer for data synchronization
    integer(kind=1), allocatable        :: my_lgt_block_send_buffer(:), my_lgt_block_receive_buffer(:)

    ! maximum heavy id, use to synchronize reduced light data array, sum of heavy blocks, start of send buffer
    integer(kind=ik)                    :: heavy_max, block_sum, buffer_start
    ! list of max heavy ids, use to build send/receive buffer
    integer(kind=ik)                    :: proc_heavy_max(params%number_procs), my_proc_heavy_max(params%number_procs)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! start time
    sub_t0 = MPI_Wtime()

    time_sum = 0.0_rk

    ! block number
    N = params%number_blocks

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! set MPI parameter
    rank         = params%rank

    ! allocate interpolation fields
    allocate( u1( 1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g ) )
    allocate( u2( 1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g ) )
    ! coarsened field is half block size + 1/2
    allocate( u3( 1:(Bs+1)/2 + g , 1:(Bs+1)/2 + g, 1:(Bs+1)/2 + g) )

    allocate( my_refinement_status( hvy_n ) )

    my_refinement_status = 0

!---------------------------------------------------------------------------------------------
! main body

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
          if ( params%threeD_case ) then
            ! ********** 3D **********
            ! copy block data to array u1
            u1(:,:,:) = hvy_block( :, :, :, dF, hvy_active(k) )
            ! now, coarsen array u1 (restriction)
            call restriction_3D( u1, u3 )  ! fine, coarse
            ! then, re-interpolate to the initial level (prediciton)
            call prediction_3D ( u3, u2, params%order_predictor )  ! coarse, fine

            ! Calculate detail by comparing u1 (original data) and u2 (result of predict(restrict(u1)))
            ! NOTE: the error (or detail) is evaluated on the entire block, INCLUDING the ghost nodes layer
            do i = 1, Bs+2*g
              do j = 1, Bs+2*g
                do l = 1, Bs+2*g
                  detail = max( detail, sqrt( (u1(i,j,l)-u2(i,j,l)) * ( u1(i,j,l)-u2(i,j,l)) ) )
                end do
              end do
            end do

          else
            ! ********** 2D **********
            ! copy block data to array u1
            u1(:,:,1) = hvy_block( :, :, 1, dF, hvy_active(k) )
            ! now, coarsen array u1 (restriction)
            call restriction_2D( u1(:,:,1), u3(:,:,1) )  ! fine, coarse
            ! then, re-interpolate to the initial level (prediciton)
            call prediction_2D ( u3(:,:,1), u2(:,:,1), params%order_predictor )  ! coarse, fine

            ! Calculate detail by comparing u1 (original data) and u2 (result of predict(restrict(u1)))
            ! NOTE: the error (or detail) is evaluated on the entire block, INCLUDING the ghost nodes layer
            do i = 1, Bs+2*g
              do j = 1, Bs+2*g
                detail = max( detail, sqrt( (u1(i,j,1)-u2(i,j,1)) * ( u1(i,j,1)-u2(i,j,1)) ) )
              end do
            end do

          end if
      end do

      ! evaluate criterion: if this blocks detail is smaller than the prescribed precision,
      ! the block is tagged as "wants to coarsen" by setting the tag -1
      ! note gradedness and completeness may prevent it from actually going through with that
      if (detail < params%eps) then
        ! coarsen block, -1
        my_refinement_status( k ) = -1
      end if

    end do

    ! ------------------------------------------------------------------------------------
    ! fourth: synchronize light data
    ! set array for max heavy ids
    my_proc_heavy_max = 0
    heavy_max = maxval(hvy_active)
    if (heavy_max == -1) heavy_max = 0
    my_proc_heavy_max(rank+1) = heavy_max

    ! synchronize array
    call MPI_Allreduce(my_proc_heavy_max, proc_heavy_max, size(proc_heavy_max,1), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! for readability, calc sum of all max heavy ids
    block_sum = sum(proc_heavy_max)

    ! now we can allocate send/receive buffer arrays
    allocate( my_lgt_block_send_buffer( block_sum ), my_lgt_block_receive_buffer( block_sum ) )

    ! reset send buffer
    my_lgt_block_send_buffer = 0
    buffer_start = sum(proc_heavy_max(1:rank))
    do k = 1, hvy_n
        my_lgt_block_send_buffer( buffer_start + hvy_active(k) ) = my_refinement_status(k)
    end do

    ! synchronize light data
    call MPI_Allreduce(my_lgt_block_send_buffer, my_lgt_block_receive_buffer, block_sum, MPI_INTEGER1, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! write synchronized light data
    ! loop over number of procs and reset lgt_block array
    do k = 1, params%number_procs
        ! proc k-1 has send data
        if ( proc_heavy_max(k) /= 0 ) then
            ! write received light data
            lgt_block( (k-1)*N+1 : (k-1)*N + proc_heavy_max(k), params%max_treelevel+2 ) =  my_lgt_block_receive_buffer( sum(proc_heavy_max(1:k-1))+1 : sum(proc_heavy_max(1:k-1))+proc_heavy_max(k) )
        else
            ! nothing to do
        end if
    end do

    ! clean up
    deallocate( u1, u2, u3, my_refinement_status )
    deallocate( my_lgt_block_send_buffer, my_lgt_block_receive_buffer )

    ! end time
    sub_t1 = MPI_Wtime()
    time_sum = time_sum + (sub_t1 - sub_t0)
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "threshold_block (w/o ghost synch.)" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "threshold_block (w/o ghost synch.)"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + time_sum
    end if

end subroutine threshold_block
