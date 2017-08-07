!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name wavelet_filter.f90
!> \version 0.5
!> \author msr
!
!> \brief wavelet filter subroutine
!
!>
!! input:    - params, heavy data,  \n
!! output:   - heavy data  \n
!!
!!
!! = log ======================================================================================
!! \n
!! 24/07/17 - create
! ********************************************************************************************

subroutine wavelet_filter( params, block_data)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :)

    ! loop parameter
    integer(kind=ik)                    :: i, j, l, N
    ! detail
    real(kind=rk)                       :: detail
    ! grid parameter
    integer(kind=ik)                    :: Bs, g
    ! interpolation fields
    real(kind=rk), allocatable          :: u1(:,:,:), u2(:,:,:), u3(:,:,:)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization
!
    ! block number
    N = params%number_blocks

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! allocate interpolation fields
    allocate( u1( 1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g ) )
    allocate( u2( 1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g ) )
    ! coarsened field is half block size + 1/2
    allocate( u3( 1:(Bs+1)/2 + g , 1:(Bs+1)/2 + g, 1:(Bs+1)/2 + g) )

!---------------------------------------------------------------------------------------------
! main body

    ! reset detail
    detail = 0.0_rk

    if ( params%threeD_case ) then
        ! ********** 3D **********
        ! copy block data to array u1
        u1(:,:,:) = block_data( :, :, : )
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

        ! evaluate criterion: if this blocks detail is smaller than the prescribed precision,
        ! the block should be filtered, overwrite block data with predicted data
        if (detail < params%eps) then
            ! wavelet filtering
            !block_data(:,:,:) = u2(:,:,:)
            ! note: do not filter redundant nodes, to avoid instabilities
            block_data(g+2:Bs+g-1,g+2:Bs+g-1,g+2:Bs+g-1) = u2(g+2:Bs+g-1,g+2:Bs+g-1,g+2:Bs+g-1)
        end if

    else
        ! ********** 2D **********
        ! copy block data to array u1
         u1(:,:,1) = block_data( :, :, 1 )
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

        ! evaluate criterion: if this blocks detail is smaller than the prescribed precision,
        ! the block should be filtered, overwrite block data with predicted data
        if (detail < params%eps) then
            ! wavelet filtering
            !block_data(:,:,1) = u2(:,:,1)
            ! note: do not filter redundant nodes, to avoid instabilities
            block_data(g+2:Bs+g-1,g+2:Bs+g-1,1) = u2(g+2:Bs+g-1,g+2:Bs+g-1,1)
        end if

    end if

    ! clean up
    deallocate( u1, u2, u3 )

end subroutine wavelet_filter
