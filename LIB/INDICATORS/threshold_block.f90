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
!! 08/09/17 - add linear wavelet filtering - discarding all details, if block level == Jmax
!!
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
    integer(kind=ik)                    :: Bs, g, R
    ! interpolation fields
    real(kind=rk), allocatable          :: u1(:,:,:), u2(:,:,:), u3(:,:,:)

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0
    !---------------------------------------------------------------------------------------------
    ! interfaces

    !---------------------------------------------------------------------------------------------
    ! variables initialization

    ! start time
    t0 = MPI_Wtime()

    ! block number
    N = params%number_blocks
    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes
    ! index of refinement status
    R  = params%max_treelevel+2

    ! set MPI parameter
    rank          = params%rank
    ! allocate interpolation fields
    allocate( u1( 1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g ) )
    allocate( u2( 1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g ) )
    ! coarsened field is half block size + 1/2
    allocate( u3( 1:(Bs+1)/2 + g , 1:(Bs+1)/2 + g, 1:(Bs+1)/2 + g) )

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

        if ( lgt_block( lgt_id, params%max_treelevel+1) == params%max_treelevel ) then
            ! email, mario, 19.04.2018:
            !> \todo Once the ghost nodes are fixed, we should remove the unnecessary removal of details at the finest level Jmax
            ! Wenn ein Block (A) auf dem maxlevel ist und ein Nachbarblock (B)
            ! nicht, dann kommt es zu einem Fehler in den redundanten Knoten. Der
            ! Nachbarblock (B) wird vor dem Zeitschritt auf maxlevel verfeinert und
            ! interpoliert dafür neue Punkte. (A) verändert sich nicht. Nun gibt es
            ! aber bei jedem zweiten redundanten Punkt, gemeint sind die gemeinsamen
            ! Punkten zwischen (A) und (B) Unterschiede. Der aktuelle
            ! Korrekturmechanismus (fein überschreibt grob) greift hier aber nicht,
            ! wenn beide Blöcke nach dem Zeitschritt auf maxlevel bleiben.
            !
            ! Daher wird aktuell die rechte Seite auf maxlevel (wenn notwendig)
            ! ausgewertet, dieses level aber auf jeden Fall wieder um 1 reduziert,
            ! so dass vor dem nächsten Zeitschritt alle Blöcke immer ein level nach
            ! oben gehen können. Daher taucht in den Bildschirmausgaben,
            ! gespeicherten Daten, ... auch maximal maxlevel-1 auf.
            lgt_block( lgt_id, R ) = -1
        else
            ! loop over all datafields
            do dF = 1, params%number_data_fields
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
                                detail = max( detail, abs(u1(i,j,l)-u2(i,j,l)) )
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
                            detail = max( detail, abs(u1(i,j,1)-u2(i,j,1)) )
                        end do
                    end do

                end if
            end do

            ! evaluate criterion: if this blocks detail is smaller than the prescribed precision,
            ! the block is tagged as "wants to coarsen" by setting the tag -1
            ! note gradedness and completeness may prevent it from actually going through with that
            if (detail < params%eps) then
                ! coarsen block, -1
                lgt_block( lgt_id, R ) = -1
            end if
        end if
    end do

    ! ------------------------------------------------------------------------------------
    ! synchronize light data
    call synchronize_lgt_data( params, lgt_block, refinement_status_only=.true. )

    ! clean up
    deallocate( u1, u2, u3 )

    ! timings
    call toc( params, "threshold_block (w/o ghost synch.)", MPI_Wtime() - t0 )
end subroutine threshold_block
