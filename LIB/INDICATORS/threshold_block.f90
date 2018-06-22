!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name threshold_block.f90
!> \version 0.6
!> \author msr, engels
!
!> \brief thresholding for a single block
!
!>
!! The block thresholding is done with the restriction/prediction operators acting on the
!! entire block \n
!!
!!
!!
!!
!! = log ======================================================================================
!! \n
!! 10/11/16 - switch to v0.4
!! 08/09/17 - add linear wavelet filtering - discarding all details, if block level == Jmax
!! 29/05/18 - now consider only one block, not the entire grid. old name was confusing.
!!
! ********************************************************************************************
!> \image html threshold.svg width=400

subroutine threshold_block( params, block_data, thresholding_component, refinement_status, norm )

    !---------------------------------------------------------------------------------------------
    ! modules

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data - this routine is called on one block only, not on the entire grid. hence th 4D array.
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    !> it can be useful not to consider all components for threhsolding here.
    !! e.g. to work only on the pressure or vorticity.
    logical, intent(in)                 :: thresholding_component(:)
    !> main output of this routine is the new satus
    integer(kind=ik), intent(out)       :: refinement_status
    !
    real(kind=rk), intent(inout)        :: norm(1:params%number_data_fields)

    ! loop parameter
    integer(kind=ik)                    :: dF, i, j, l
    ! detail
    real(kind=rk)                       :: detail(1:params%number_data_fields)
    ! grid parameter
    integer(kind=ik)                    :: Bs, g

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0

    !---------------------------------------------------------------------------------------------
    ! variables initialization

    ! start time
    t0 = MPI_Wtime()
    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes


    ! reset detail
    detail = 0.0_rk

    ! loop over all datafields
    do dF = 1, params%number_data_fields
        ! is this component of the block used for thresholding or not?
        if (thresholding_component(dF)) then
            if (abs(norm(dF))<1.e-10_rk) norm(dF) = 1.0_rk ! avoid division by zero

            if ( params%threeD_case ) then
                ! ********** 3D **********
                ! allocate interpolation fields
                if (.not.allocated(u1)) allocate( u1( 1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g ) )
                if (.not.allocated(u2)) allocate( u2( 1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g ) )
                ! coarsened field is half block size + 1/2
                if (.not.allocated(u3)) allocate( u3( 1:(Bs+1)/2 + g , 1:(Bs+1)/2 + g, 1:(Bs+1)/2 + g) )

                ! copy block data to array u1
                u1(:,:,:) = block_data( :, :, :, dF )
                ! now, coarsen array u1 (restriction)
                call restriction_3D( u1, u3 )  ! fine, coarse
                ! then, re-interpolate to the initial level (prediciton)
                call prediction_3D ( u3, u2, params%order_predictor )  ! coarse, fine

                ! Calculate detail by comparing u1 (original data) and u2 (result of predict(restrict(u1)))
                ! NOTE: the error (or detail) is evaluated on the entire block, INCLUDING the ghost nodes layer
                do i = 1, Bs+2*g
                    do j = 1, Bs+2*g
                        do l = 1, Bs+2*g
                            detail(dF) = max( detail(dF), abs(u1(i,j,l)-u2(i,j,l)) / norm(dF) )
                        end do
                    end do
                end do
            else
                ! ********** 2D **********
                ! allocate interpolation fields
                if (.not.allocated(u1)) allocate( u1( 1:Bs+2*g, 1:Bs+2*g, 1 ) )
                if (.not.allocated(u2)) allocate( u2( 1:Bs+2*g, 1:Bs+2*g, 1 ) )
                ! coarsened field is half block size + 1/2
                if (.not.allocated(u3)) allocate( u3( 1:(Bs+1)/2 + g , 1:(Bs+1)/2 + g, 1) )

                ! copy block data to array u1
                u1(:,:,1) = block_data( :, :, 1, dF )
                ! now, coarsen array u1 (restriction)
                call restriction_2D( u1(:,:,1), u3(:,:,1) )  ! fine, coarse
                ! then, re-interpolate to the initial level (prediciton)
                call prediction_2D ( u3(:,:,1), u2(:,:,1), params%order_predictor )  ! coarse, fine

                ! Calculate detail by comparing u1 (original data) and u2 (result of predict(restrict(u1)))
                ! NOTE: the error (or detail) is evaluated on the entire block, INCLUDING the ghost nodes layer
                do i = 1, Bs+2*g
                    do j = 1, Bs+2*g
                        detail(dF) = max( detail(dF), abs(u1(i,j,1)-u2(i,j,1)) / norm(dF) )
                    end do
                end do

            end if
        end if
    end do


    ! evaluate criterion: if this blocks detail is smaller than the prescribed precision,
    ! the block is tagged as "wants to coarsen" by setting the tag -1
    ! note gradedness and completeness may prevent it from actually going through with that
    if ( maxval(detail) < params%eps) then
        ! coarsen block, -1
        refinement_status = -1
    end if


    ! clean up
    deallocate( u1, u2, u3 )

    ! timings
    call toc( params, "threshold_block (w/o ghost synch.)", MPI_Wtime() - t0 )
end subroutine threshold_block
