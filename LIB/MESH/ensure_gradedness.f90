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

subroutine ensure_gradedness( block_list, neighbor_list, N, max_treelevel, lgt_active, n_lgt )

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! light data array
    integer(kind=ik), intent(inout)     :: block_list(:, :)
    ! neighbor list
    integer(kind=ik), intent(in)        :: neighbor_list(:)
    ! number of blocks
    integer(kind=ik), intent(in)        :: N
    ! max treelevel
    integer(kind=ik), intent(in)        :: max_treelevel
    ! active_block_list (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(in)        :: n_lgt

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! loop variables
    integer(kind=ik)                    :: k, i, mylevel, neighbor_level, counter, heavy_id, neighbor_status

    ! status of grid changing
    logical                             :: grid_changed

    ! refinement status change
    integer(kind=ik)                    :: refine_change( n_lgt ), my_refine_change( n_lgt )

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

!---------------------------------------------------------------------------------------------
! main body

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

        do k = 1, n_lgt

            ! proc is responsible for current block
            if ( ((lgt_active(k)-1) / N) == rank ) then

                !-----------------------------------------------------------------------
                ! block wants to coarsen, refinement only in safety zone step
                !-----------------------------------------------------------------------
                if ( block_list( lgt_active(k) , max_treelevel+2 ) == -1) then

                    ! heavy id
                    heavy_id = lgt_active(k) - rank*N

                    ! loop over all neighbors
                    do i = 1, 16
                        ! neighbor exists
                        if ( neighbor_list( (heavy_id - 1)*16 + i ) /= -1 ) then

                            ! check neighbor treelevel
                            mylevel         = block_list( lgt_active(k), max_treelevel+1 )
                            neighbor_level  = block_list( neighbor_list( (heavy_id - 1)*16 + i ) , max_treelevel+1 )
                            neighbor_status = block_list( neighbor_list( (heavy_id - 1)*16 + i ) , max_treelevel+2 )

                            if (mylevel == neighbor_level) then
                                ! neighbor on same level
                                ! block can not coarsen, if neighbor wants to refine
                                if ( neighbor_status == 1) then
                                    ! block can not coarsen, change refinement with +1
                                    my_refine_change(k) = 1
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
                                    my_refine_change(k) = 1

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
        call MPI_Allreduce(my_refine_change, refine_change, n_lgt, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

        ! -------------------------------------------------------------------------------------
        ! third: change light data and set grid_changed status

        ! loop over active blocks
        do k = 1, n_lgt

            ! refinement status changed
            if ( refine_change(k) > 0 ) then
                ! change light data
                block_list( lgt_active(k), max_treelevel+2 ) = block_list( lgt_active(k), max_treelevel+2 ) + refine_change(k)
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

end subroutine ensure_gradedness
