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
    integer(kind=ik) :: n_coarse(1:3), n_fine(1:3), idim
    real(kind=rk), SAVE, allocatable :: ufine(:,:,:), ucoarse(:,:,:)

    t0 = MPI_Wtime()
    Bs = params%Bs
    g  = params%n_ghosts
    detail = -1.0_rk

    do idim = 1, params%dim
        if ( modulo(Bs(idim),2)==0) then
            ! Bs even
            n_coarse(idim) = (Bs(idim) + 2*g) / 2
            n_fine(idim) = Bs(idim) + 2*g - 1
        else
            ! Bs odd
            n_coarse(idim) = (Bs(idim) + 2*g + 1)/2
            n_fine(idim) = Bs(idim) + 2*g
        endif
    enddo

    ! --------------------------- 3D -------------------------------------------

    if (params%dim == 3) then
        if (.not.allocated(ufine)) allocate( ufine( 1:n_fine(1), 1:n_fine(2), 1:n_fine(3) ) )
        if (.not.allocated(ucoarse)) allocate( ucoarse( 1:n_coarse(1) , 1:n_coarse(2), 1:n_coarse(3)) )

        ! loop over all datafields
        do dF = 1, size(block_data,4)
            ! is this component of the block used for thresholding or not?
            if (thresholding_component(dF)) then
                if (params%harten_multiresolution) then
                    ! coarsen block data (restriction)
                    call restriction_3D( block_data(1:n_fine(1), 1:n_fine(2), 1:n_fine(3), dF ), ucoarse )  ! fine, coarse
                    ! then, re-interpolate to the initial level (prediciton)
                    call prediction_3D ( ucoarse, ufine, params%order_predictor )  ! coarse, fine

                    ! Calculate detail by comparing u1 (original data) and ufine (result of predict(restrict(u1)))
                    ! NOTE: the detail is evaluated on the entire block, INCLUDING the ghost nodes layer
                    ! NOTE: if the block size is EVEN, the last detail in the block cannot be computed (this is a general restriction
                    ! of the predict / restrict cycle and not a bug) therefore the loops run only until n_fine(i)
                    do l = 1, n_fine(3)
                        do j = 1, n_fine(2)
                            do i = 1, n_fine(1)
                                detail(dF) = max( detail(dF), abs(block_data(i,j,l,dF)-ufine(i,j,l)) / norm(dF) )
                            end do
                        end do
                    end do
                else
                    ! apply a smoothing filter (low-pass, in wavelet terminology h_tilde)
                    call restriction_prefilter_3D( block_data(1:n_fine(1), 1:n_fine(2), 1:n_fine(3), dF), ufine(:,:,:), params%wavelet )
                    ! now, coarsen block data (restriction)
                    call restriction_3D( ufine, ucoarse )  ! fine, coarse
                    ! then, re-interpolate to the initial level (prediciton)
                    call prediction_3D ( ucoarse, ufine, params%order_predictor )  ! coarse, fine

                    ! Calculate detail by comparing u1 (original data) and ufine (result of predict(restrict(u1)))
                    ! NOTE: we EXCLUDE ghost nodes
                    do l = g+1, Bs(3)+g
                        do j = g+1, Bs(2)+g
                            do i = g+1, Bs(1)+g
                                detail(dF) = max( detail(dF), abs(block_data(i,j,l,dF)-ufine(i,j,l)) / norm(dF) )
                            end do
                        end do
                    end do
                end if

            end if
        end do
    end if

    ! --------------------------- 2D -------------------------------------------

    if (params%dim == 2 ) then
        if (.not.allocated(ufine)) allocate( ufine(1:n_fine(1), 1:n_fine(2), 1) )
        if (.not.allocated(ucoarse)) allocate( ucoarse(1:n_coarse(1), 1:n_coarse(2), 1) )

        ! loop over all datafields
        do dF = 1, size(block_data,4)
            ! is this component of the block used for thresholding or not?
            if (thresholding_component(dF)) then
                ! Harten multiresolution or biorthogonal?
                if (params%harten_multiresolution) then
                    ! coarsen block data (restriction)
                    call restriction_2D( block_data(1:n_fine(1), 1:n_fine(2), 1, dF), ucoarse(:,:,1) )  ! fine, coarse
                    ! then, re-interpolate to the initial level (prediciton)
                    call prediction_2D ( ucoarse(:,:,1), ufine(:,:,1), params%order_predictor )  ! coarse, fine

                    ! Calculate detail by comparing u1 (original data) and ufine (result of predict(restrict(u1)))
                    ! NOTE: the error (or detail) is evaluated on the entire block, INCLUDING the ghost nodes layer
                    ! NOTE: if the block size is EVEN, the last detail in the block cannot be computed (this is a general restriction
                    ! of the predict / restrict cycle and not a bug) therefore the loops run only until n_fine(i)
                    do j = 1, n_fine(2)
                        do i = 1, n_fine(1)
                            detail(dF) = max( detail(dF), abs(block_data(i,j,1,dF)-ufine(i,j,1)) / norm(dF) )
                        end do
                    end do
                else
                    ! apply a smoothing filter (low-pass, in wavelet terminology h_tilde)
                    call restriction_prefilter_2D(block_data(1:n_fine(1), 1:n_fine(2), 1, dF), ufine(:,:,1), params%wavelet)
                    ! now, coarsen block data (restriction)
                    call restriction_2D( ufine(:,:,1), ucoarse(:,:,1) )  ! fine, coarse
                    ! then, re-interpolate to the initial level (prediciton)
                    call prediction_2D ( ucoarse(:,:,1), ufine(:,:,1), params%order_predictor )  ! coarse, fine

                    ! Calculate detail by comparing u1 (original data) and ufine (result of predict(restrict(u1)))
                    do j = g+1, Bs(2)+g
                        do i = g+1, Bs(1)+g
                            detail(dF) = max( detail(dF), abs(block_data(i,j,1,dF)-ufine(i,j,1)) / norm(dF) )
                        end do
                    end do
                end if
            end if
        end do
    end if




! ich habe die wavelet normalization ausgebruetet und aufgeschrieben.
! ich schicke dir die notizen gleich (photos).
!
! also wir brauchen einen scale(level)- dependent threshold, d.h. \epsilon_j
! zudem ist dieser abhaengig von der raum dimension d.
!
! Fuer die L^2 normalisierung (mit wavelets welche in der L^\infty norm normalisiert sind) haben wir
!
! \epsilon_j = 2^{-jd/2} \epsilon
!
! d.h. der threshold wird kleiner auf kleinen skalen.
!
! Fuer die vorticity (anstatt der velocity) kommt nochmal ein faktor 2^{-j} dazu, d.h.
!
! \epsilon_j = 2^{-j(d+2)/2} \epsilon
!
! Zum testen waere es gut in 1d oder 2d zu pruefen, ob die L^2 norm von u - u_\epsilon
! linear mit epsilon abnimmt, das gleiche koennte man auch fuer H^1 (philipp koennte dies doch mal ausprobieren?).
!
! fuer CVS brauchen wir dann noch \epsilon was von Z (der enstrophy) und der feinsten
! aufloesung abhaengt. fuer L^2 normalisierte wavelets ist
! der threshold:
!
! \epsilon = \sqrt{2/3 \sigma^2 \ln N}
!
! wobei \sigma^2 die varianz (= 2 Z) der incoh. vorticity ist.
! typischerweise erhaelt man diese mit 1-3 iterationen.
! als ersten schritt koennen wir einfach Z der totalen stroemung nehmen.
! N ist die maximale aufloesung, typicherweise 2^{d J}.
!

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
        eps2 = eps2 * ( 2.0_rk**(-dble(level*params%dim)/2.0_rk) )

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
