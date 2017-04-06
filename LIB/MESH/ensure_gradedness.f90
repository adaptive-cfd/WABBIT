


! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: ensure_gradedness.f90
! version: 0.5
! author: msr, engels
!
! This routine is called after all blocks have been tagged whether to refine or coarsen or stay.
! It now goes through the list of blocks and looks for refinement or coarsening states that would
! result in an non-graded mesh. These mistakes are corrected, their status -1 or 0 is overwritten.
! The status +1 is always conserved (recall to call respect_min_max_treelevel before).
!
! Sine 04/2017, the new code checks all blocks that want to coarsen or remain, NOT the ones that
! want to refine, as was done in prototypes. The reason is MPI: I cannot easily set the flags of
! my neighbors, as they might reside on another proc.
!
! input:    - light data, neighbor list, list of active blocks(light data)
! output:   - light data array
!
! = log ======================================================================================
!
! 10/11/16 - switch to v0.4
! 23/11/16 - rework complete subroutine: use list of active blocks, procs works now on light data
! 03/02/17 - insert neighbor_num variable to use subroutine for 2D and 3D data
! 05/04/17 - Improvement: Ensure a graded mesh in any case, not only in the coarsen states (which was done before)
!
! ********************************************************************************************

subroutine ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    ! neighbor list
    integer(kind=ik), intent(in)        :: hvy_neighbor(:, :)

    ! active_block_list (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! loop variables
    integer(kind=ik)                    :: k, i, N, mylevel, neighbor_level, counter, hvy_id, neighbor_status, max_treelevel, proc_id, lgt_id

    ! status of grid changing
    logical                             :: grid_changed

    ! refinement status change
    integer(kind=ik)                    :: refine_change( size(lgt_block,1) ), my_refine_change( size(lgt_block,1) )

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

    ! number of neighbor relations
    ! 2D: 16, 3D: 74
    integer(kind=ik)                    :: neighbor_num

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks
    max_treelevel = params%max_treelevel

    ! set MPI parameter
    rank         = params%rank

    if ( params%threeD_case ) then
        ! 3D:
        neighbor_num = 74
    else
        ! 2D:
        neighbor_num = 16
    end if

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! we repeat the ensure_gradedness procedure until this flag is .false. since as long
    ! as the grid changes due to gradedness requirements, we have to check it again
    grid_changed = .true. ! set true to tigger the loop
    counter = 0

    do while ( grid_changed )
        ! we hope not to set the flag to .true. again in this iteration
        grid_changed    = .false.

        ! -------------------------------------------------------------------------------------
        ! first: every proc loop over the light data and calculate the refinement status change
        my_refine_change = -999999
        refine_change = -99999

        do k = 1, lgt_n

            ! calculate proc rank respondsible for current block
            call lgt_id_to_proc_rank( proc_id, lgt_active(k), N )
            ! proc is responsible for current block
            if ( proc_id == rank ) then
                ! Get this blocks heavy id
                call lgt_id_to_hvy_id( hvy_id, lgt_active(k), rank, N )

                !-----------------------------------------------------------------------
                ! This block wants to coarsen
                !-----------------------------------------------------------------------
                if ( lgt_block( lgt_active(k) , max_treelevel+2 ) == -1) then
                ! loop over all neighbors
                do i = 1, neighbor_num
                    ! neighbor exists ? If not, this is a bad error
                    if ( hvy_neighbor( hvy_id, i ) /= -1 ) then

                        ! check neighbor treelevel
                        mylevel         = lgt_block( lgt_active(k), max_treelevel+1 )
                        neighbor_level  = lgt_block( hvy_neighbor( hvy_id, i ) , max_treelevel+1 )
                        neighbor_status = lgt_block( hvy_neighbor( hvy_id, i ) , max_treelevel+2 )

                        if (mylevel == neighbor_level) then
                              ! neighbor on same level
                              ! block can not coarsen, if neighbor wants to refine
                              if ( neighbor_status == -1 ) then
                                  ! neighbor wants to coarsen, as do I, we're on the same level -> ok
                              elseif ( neighbor_status == 0 ) then
                                  ! neighbor wants to stay, I want to coarsen, we're on the same level -> ok
                              elseif ( neighbor_status == 1 ) then
                                  ! neighbor wants to refine, I want to coarsen, we're on the same level -> NOT OK
                                  ! I have at least to stay on my level.
                                  ! Note we cannot simply set 0 as we could accidentally overwrite a refinement flag
                                  my_refine_change( lgt_active(k) ) = max( 0, my_refine_change( lgt_active(k) ) )
                              end if
                        elseif (mylevel - neighbor_level == 1) then
                              ! neighbor on lower level
                              if ( neighbor_status == -1 ) then
                                  ! neighbor wants to coarsen, as do I, he is one level coarser, -> ok
                              elseif ( neighbor_status == 0 ) then
                                  ! neighbor wants to stay, I want to coarsen, he is one level coarser, -> ok
                              elseif ( neighbor_status == 1 ) then
                                  ! neighbor wants to refine, I want to coarsen,  he is one level coarser, -> ok
                              end if
                        elseif (neighbor_level - mylevel == 1) then
                              ! neighbor on higher level
                              ! neighbor wants to refine, ...
                              if ( neighbor_status == +1) then
                                  ! ... so I also have to refine (not only can I NOT coarsen, I actually
                                  ! have to refine!)
                                  my_refine_change( lgt_active(k) ) = +1
                              elseif ( neighbor_status == 0) then
                                  ! neighbor wants to stay and I want to coarsen, but
                                  ! I cannot do that (there would be two levels between us)
                                  ! Note we cannot simply set 0 as we could accidentally overwrite a refinement flag
                                  my_refine_change( lgt_active(k) ) = max( 0, my_refine_change( lgt_active(k) ) )
                              elseif ( neighbor_status == -1) then
                                  ! neighbor wants to coarsen, which is what I want too,
                                  ! so we both would just go up one level together - that's fine
                                  ! FIXME: I have no idea why the following line is required.
                                  my_refine_change( lgt_active(k) ) = max( 0, my_refine_change( lgt_active(k) ) )
                              end if
                        else
                          call error_msg("ERROR: ensure_gradedness: my neighbor does not seem to have -1,0,+1 level diff!")
                        end if
                    end if ! if neighbor exists
                    end do ! loop over neighbors

            !-----------------------------------------------------------------------
            ! this block wants to stay on his level
            !-----------------------------------------------------------------------
            elseif (lgt_block( lgt_active(k) , max_treelevel+2 ) == 0) then
                    ! loop over all neighbors
                    do i = 1, neighbor_num
                      ! neighbor exists ? If not, this is a bad error
                      if ( hvy_neighbor( hvy_id, i ) /= -1 ) then
                            mylevel     = lgt_block( lgt_active(k), max_treelevel+1 )
                            neighbor_level = lgt_block( hvy_neighbor( hvy_id, i ) , max_treelevel+1 )
                            neighbor_status = lgt_block( hvy_neighbor( hvy_id, i ) , max_treelevel+2 )

                            if (mylevel == neighbor_level) then
                              ! me and my neighbor are on the same level
                              ! As I'd wish to stay where I am, my neighbor is free to go -1,0,+1
                            elseif (mylevel - neighbor_level == 1) then
                              ! my neighbor is one level coarser
                              ! My neighbor can stay or refine, but not coarsen. This case is however handled above (coarsening inhibited)
                            elseif (neighbor_level - mylevel == 1) then
                              ! my neighbor is one level finer
                              if (neighbor_status == +1) then
                                ! neighbor refines (and we cannot inhibt that) so I HAVE TO do so as well
                                my_refine_change( lgt_active(k) ) = +1
                              end if
                            else
                              call error_msg("ERROR: ensure_gradedness: my neighbor does not seem to have -1,0,+1 level diff!")
                            end if
                      end if ! if neighbor exists
                      end do
                end if ! my refinement status
            end if ! mpi responsibile for block
        end do ! loop over blocks

        ! second: synchronize local refinement changes.
        ! Now all procs have looked at their blocks, and they may have modified their blocks
        ! (remove coarsen states, force refinement through neighbors). None of the procs had to
        ! touch the neighboring blocks, there can be no MPI conflicts.
        ! So we can simply snynchronize and know what changes have to be made
        call MPI_Allreduce(my_refine_change, refine_change, size(lgt_block,1), MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ierr)

        ! -------------------------------------------------------------------------------------
        ! third: change light data and set grid_changed status
        ! loop over active blocks
        do k = 1, lgt_n
            lgt_id = lgt_active(k)
            ! refinement status changed
            if ( refine_change(lgt_id) > -999 ) then
                ! change light data
                lgt_block( lgt_id, max_treelevel+2 ) = max( lgt_block( lgt_id, max_treelevel+2 ), refine_change(lgt_id) )
                ! set grid status since if we changed something, we have to repeat the entire process
                grid_changed = .true.
            end if
        end do


        ! avoid infinite loops
        counter = counter + 1
        if (counter == 10) call error_msg("ERROR: unable to build a graded mesh")
    end do ! end do of repeat procedure until grid_changed==.false.




    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "ensure_gradedness" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "ensure_gradedness"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine ensure_gradedness
