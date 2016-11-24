! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: ensure_gradedness.f90
! version: 0.4
! author: msr, engels
!
! check the gradedness after new refinement status
!
! input:    - light data, neighbor list, list of active blocks(light data)
! output:   - light data array
!
! = log ======================================================================================
!
! 10/11/16 - switch to v0.4
! 23/11/16 - rework complete subroutine: use list of active blocks, procs works now on light data
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
    integer(kind=ik)                    :: k, i, N, mylevel, neighbor_level, counter, hvy_id, neighbor_status, max_treelevel

    ! status of grid changing
    logical                             :: grid_changed

    ! refinement status change
    integer(kind=ik)                    :: refine_change( size(lgt_block,1) ), my_refine_change( size(lgt_block,1) )

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks
    max_treelevel = params%max_treelevel

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

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
        refine_change    = 0
        my_refine_change = 0

        do k = 1, lgt_n

            ! proc is responsible for current block
            if ( ((lgt_active(k)-1) / N) == rank ) then

                !-----------------------------------------------------------------------
                ! block wants to coarsen, refinement only in safety zone step
                !-----------------------------------------------------------------------
                if ( lgt_block( lgt_active(k) , max_treelevel+2 ) == -1) then

                    ! heavy id
                    call lgt_id_to_hvy_id( hvy_id, lgt_active(k), rank, N )

                    ! loop over all neighbors
                    do i = 1, 16
                        ! neighbor exists
                        if ( hvy_neighbor( hvy_id, i ) /= -1 ) then

                            ! check neighbor treelevel
                            mylevel         = lgt_block( lgt_active(k), max_treelevel+1 )
                            neighbor_level  = lgt_block( hvy_neighbor( hvy_id, i ) , max_treelevel+1 )
                            neighbor_status = lgt_block( hvy_neighbor( hvy_id, i ) , max_treelevel+2 )

                            if (mylevel == neighbor_level) then
                                ! neighbor on same level
                                ! block can not coarsen, if neighbor wants to refine
                                if ( neighbor_status == 1) then
                                    ! block can not coarsen, change refinement with +1
                                    my_refine_change( lgt_active(k) ) = 1
                                end if

                            elseif (mylevel - neighbor_level == 1) then
                                ! neighbor on lower level
                                ! nothing to do

                            elseif (neighbor_level - mylevel == 1) then
                                ! neighbor on higher level
                                ! neighbor want to refine, ...
                                if ( neighbor_status == 1) then
                                    ! error case, neighbor can not have +1 refinement status
                                    print*, "ERROR: gradedness, neighbor wants to refine"
                                    stop

                                elseif ( neighbor_status == 0) then
                                    ! neighbor stays on his level,
                                    ! so block can not coarsen
                                    ! change refinement with +1
                                    my_refine_change( lgt_active(k) ) = 1

                                else
                                    ! neighbor want to coarsen -> nothing to do

                                end if

                            else
                              ! error case
                              print*, "ERROR: block has wrong number of neighbors!"
                              stop

                            end if

                        end if
                    end do

                end if

            endif

        end do

        ! -------------------------------------------------------------------------------------
        ! second: synchronize refinement change
        call MPI_Allreduce(my_refine_change, refine_change, size(lgt_block,1), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

        ! -------------------------------------------------------------------------------------
        ! third: change light data and set grid_changed status

        ! loop over active blocks
        do k = 1, lgt_n

            ! refinement status changed
            if ( refine_change(k) > 0 ) then
                ! change light data
                lgt_block( lgt_active(k), max_treelevel+2 ) = lgt_block( lgt_active(k), max_treelevel+2 ) + refine_change(lgt_active(k))
                ! set grid status
                grid_changed = .true.
            end if

        end do

        ! debug error
        counter = counter + 1
        if (counter == 10) then
            print*, "ERROR: unable to build a graded mesh"
            stop
        end if

    ! end do of repeat procedure until grid_changed==.false.
    end do

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find first free line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "ensure_gradedness"
        debug%comp_time(rank+1, k) = sub_t1 - sub_t0
    end if

end subroutine ensure_gradedness
