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

subroutine threshold_block( params, block_data, thresholding_component, refinement_status, norm, level, eps )

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data - this routine is called on one block only, not on the entire grid. hence th 4D array.
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    !> it can be useful not to consider all components for thresholding here.
    !! e.g. to work only on the pressure or vorticity.
    logical, intent(in)                 :: thresholding_component(:)
    !> main output of this routine is the new satus
    integer(kind=ik), intent(out)       :: refinement_status
    ! If we use L2 or H1 normalization, the threshold eps is level-dependent, hence
    ! we pass the level to this routine
    integer(kind=ik), intent(in)        :: level
    !
    real(kind=rk), intent(inout)        :: norm( size(block_data,4) )
    ! if different from the default eps (params%eps), you can pass a different value here. This is optional
    ! and used for example when thresholding the mask function.
    real(kind=rk), intent(in), optional :: eps

    integer(kind=ik)                    :: dF, i, j, l
    real(kind=rk)                       :: detail( size(block_data,4) )
    integer(kind=ik)                    :: g
    integer(kind=ik), dimension(3)      :: Bs
    real(kind=rk)                       :: t0, eps2

    t0 = MPI_Wtime()
    Bs = params%Bs
    g  = params%n_ghosts

    ! reset detail
    detail = -1.0_rk

    ! --------------------------- 3D -------------------------------------------

    if (params%dim == 3) then

        ! allocate interpolation fields
        if (.not.allocated(u2)) allocate( u2( 1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g ) )
        ! coarsened field is half block size + 1/2
        if (.not.allocated(u3)) allocate( u3( 1:(Bs(1)+1)/2 + g , 1:(Bs(2)+1)/2 + g, 1:(Bs(3)+1)/2 + g) )

        ! loop over all datafields
        do dF = 1, size(block_data,4)
            ! is this component of the block used for thresholding or not?
            if (thresholding_component(dF)) then
                if (params%harten_multiresolution) then
                    ! coarsen block data (restriction)
                    call restriction_3D( block_data( :, :, :, dF ), u3 )  ! fine, coarse
                    ! then, re-interpolate to the initial level (prediciton)
                    call prediction_3D ( u3, u2, params%order_predictor )  ! coarse, fine

                    ! Calculate detail by comparing u1 (original data) and u2 (result of predict(restrict(u1)))
                    ! NOTE: the detail is evaluated on the entire block, INCLUDING the ghost nodes layer
                    do i = 1, Bs(1)+2*g
                        do j = 1, Bs(2)+2*g
                            do l = 1, Bs(3)+2*g
                                detail(dF) = max( detail(dF), abs(block_data(i,j,l,dF)-u2(i,j,l)) / norm(dF) )
                            end do
                        end do
                    end do
                else
                    ! apply a smoothing filter (low-pass, in wavelet terminology h_tilde)
                    call restriction_prefilter_3D( block_data( :, :, :, dF ), u2(:,:,:), params%wavelet )
                    ! now, coarsen block data (restriction)
                    call restriction_3D( u2, u3 )  ! fine, coarse
                    ! then, re-interpolate to the initial level (prediciton)
                    call prediction_3D ( u3, u2, params%order_predictor )  ! coarse, fine

                    ! Calculate detail by comparing u1 (original data) and u2 (result of predict(restrict(u1)))
                    ! NOTE: we EXCLUDE ghost nodes
                    do i = g+1, Bs(1)+g
                        do j = g+1, Bs(2)+g
                            do l = g+1, Bs(3)+g
                                detail(dF) = max( detail(dF), abs(block_data(i,j,l,dF)-u2(i,j,l)) / norm(dF) )
                            end do
                        end do
                    end do
                end if

            end if
        end do
    end if

    ! --------------------------- 2D -------------------------------------------

    if (params%dim == 2 ) then
        ! allocate interpolation fields
        if (.not.allocated(u2)) allocate( u2( 1:Bs(1)+2*g, 1:Bs(2)+2*g, 1 ) )
        ! coarsened field is half block size + 1/2
        if (.not.allocated(u3)) allocate( u3( 1:(Bs(1)+1)/2 + g , 1:(Bs(2)+1)/2 + g, 1) )

        ! loop over all datafields
        do dF = 1, size(block_data,4)
            ! is this component of the block used for thresholding or not?
            if (thresholding_component(dF)) then
                ! Harten multiresolution or biorthogonal?
                if (params%harten_multiresolution) then
                    ! coarsen block data (restriction)
                    call restriction_2D( block_data( :, :, 1, dF ), u3(:,:,1) )  ! fine, coarse
                    ! then, re-interpolate to the initial level (prediciton)
                    call prediction_2D ( u3(:,:,1), u2(:,:,1), params%order_predictor )  ! coarse, fine

                    ! Calculate detail by comparing u1 (original data) and u2 (result of predict(restrict(u1)))
                    ! NOTE: the error (or detail) is evaluated on the entire block, INCLUDING the ghost nodes layer
                    do i = 1, Bs(1)+2*g
                        do j = 1, Bs(2)+2*g
                            detail(dF) = max( detail(dF), abs(block_data(i,j,1,dF)-u2(i,j,1)) / norm(dF) )
                        end do
                    end do
                else
                    ! apply a smoothing filter (low-pass, in wavelet terminology h_tilde)
                    call restriction_prefilter_2D(block_data( :, :, 1, dF ), u2(:,:,1), params%wavelet)
                    ! now, coarsen block data (restriction)
                    call restriction_2D( u2(:,:,1), u3(:,:,1) )  ! fine, coarse
                    ! then, re-interpolate to the initial level (prediciton)
                    call prediction_2D ( u3(:,:,1), u2(:,:,1), params%order_predictor )  ! coarse, fine

                    ! Calculate detail by comparing u1 (original data) and u2 (result of predict(restrict(u1)))
                    do i = g+1, Bs(1)+g
                        do j = g+1, Bs(2)+g
                            detail(dF) = max( detail(dF), abs(block_data(i,j,1,dF)-u2(i,j,1)) / norm(dF) )
                        end do
                    end do
                end if
            end if
        end do
    end if

    ! default threshlding level is the one in the parameter struct
    eps2 = params%eps
    ! but if we pass another one, use that.
    if (present(eps)) eps2 = eps

    select case(params%eps_norm)
    case ("Linfty")
        ! do nothing, our wavelets are normalized in L_infty norm by default, hence
        ! a simple threshold controls this norm
        eps2 = eps2

    case ("L2")
        ! If we want to control the L2 norm (with wavelets that are normalized in Linfty norm)
        ! we have to have a level-dependent threshold
        eps2 = eps2 * ( 2**(-level*params%dim*0.5_rk) )
        
    case ("H1")
        ! H1 norm mimicks filtering of vorticity
        eps2 = eps2 * ( 2**(-level*(params%dim+2.0_rk)*0.5_rk) )

    case default
        call abort(20022811, "ERROR:threshold_block.f90:Unknown wavelet normalization!")

    end select

    ! evaluate criterion: if this blocks detail is smaller than the prescribed precision,
    ! the block is tagged as "wants to coarsen" by setting the tag -1
    ! note gradedness and completeness may prevent it from actually going through with that
    if ( maxval(detail) < eps2) then
        ! coarsen block, -1
        refinement_status = -1
    else
        refinement_status = 0
    end if


    ! timings
    call toc( "threshold_block (w/o ghost synch.)", MPI_Wtime() - t0 )
end subroutine threshold_block
