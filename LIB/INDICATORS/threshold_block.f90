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
!> \image html threshold.svg width=400

subroutine threshold_block( params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_active, hvy_n)

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
    ! interpolation fields
    real(kind=rk), allocatable          :: u1(:,:,:), u2(:,:,:), u3(:,:,:)
    ! light data (refinement status column) list for working
    integer(kind=ik), allocatable       :: my_refinement_status(:)

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1, time_sum

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

    allocate( my_refinement_status( size(lgt_block, 1)) )

    my_refinement_status = 0
    my_refinement_status( rank*N+1: rank*N+N ) = lgt_block( rank*N+1: rank*N+N, params%max_treelevel+2)

!---------------------------------------------------------------------------------------------
! main body

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
        my_refinement_status( lgt_id ) = -1
      end if

    end do

    ! ------------------------------------------------------------------------------------
    ! fourth: synchronize light data
    lgt_block(:,params%max_treelevel+2) = 0
    call MPI_Allreduce(my_refinement_status, lgt_block(:,params%max_treelevel+2), size(lgt_block,1), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! clean up
    deallocate( u1, u2, u3, my_refinement_status )

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
