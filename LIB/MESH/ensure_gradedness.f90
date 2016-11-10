! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: ensure_gradedness.f90
! version: 0.4
! author: msr, engels
!
! check the gradedness after new refinement status
!
! input:    - light data, neighbor list
! output:   - light data array
!
! = log ======================================================================================
!
! 10/11/16 - switch to v0.4
! ********************************************************************************************

subroutine ensure_gradedness( block_list, neighbor_list, N, max_treelevel )

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

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! loop variables
    integer(kind=ik)                    :: k, i, mylevel, neighbor_level, counter
    ! status of grid changing
    logical                             :: grid_changed

    ! light data list for working
    integer(kind=ik)                    :: my_block_list( size(block_list, 1), max_treelevel+2)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! set light data list for working, only light data coresponding to proc are not zero
    my_block_list = 0
    my_block_list( rank*N+1: rank*N+N, :) = block_list( rank*N+1: rank*N+N, :)

!---------------------------------------------------------------------------------------------
! main body

    ! we repeat the ensure_gradedness procedure until this flag is .false. since as long
    ! as the grid changes due to gradedness requirements, we have to check it again
    grid_changed = .true. ! set true to tigger the loop
    counter = 0

    do while ( grid_changed )
        ! we hope not to set the flag to .true. again in this iteration
        grid_changed = .false.

        ! loop over heavy data
        do k = 1, N
            ! block is active
            if ( block_list(rank*N + k , 1) /= -1 ) then

                !-----------------------------------------------------------------------
                ! block wants to coarsen, refinement only in safety zone step
                !-----------------------------------------------------------------------
                if ( block_list( rank*N + k , max_treelevel+2 ) == -1) then

                    ! loop over all neighbors
                    do i = 1, 16
                        ! neighbor exists
                        if ( neighbor_list( (k - 1)*16 + i ) /= -1 ) then

                            ! check neighbor treelevel
                            mylevel         = block_list( rank*N + k, max_treelevel+1 )
                            neighbor_level  = block_list( neighbor_list( (k - 1)*16 + i ) , max_treelevel+1 )

                            if (mylevel == neighbor_level) then
                              ! neighbor on same level
                              ! block can not coarsen, if neighbor wants to refine
                              if ( block_list( neighbor_list( (k - 1)*16 + i ) , max_treelevel+2 ) == 1) then
                                my_block_list( rank*N + k, max_treelevel+2 ) = max(0, block_list( rank*N + k , max_treelevel+2 ) )
                                grid_changed = .true. ! we changed something in the grid
                              end if

                            elseif (mylevel - neighbor_level == 1) then
                              ! neighbor on lower level
                              ! nothing to do

                            elseif (neighbor_level - mylevel == 1) then
                              ! neighbor on higher level
                              ! neighbor want to refine, ...
                              if ( block_list( neighbor_list( (k - 1)*16 + i ) , max_treelevel+2 ) == 1) then
                                  ! ... so I also have to refine
                                  my_block_list( rank*N + k, max_treelevel+2 ) = 1
                                  grid_changed = .true. ! we changed something in the grid
                              else
                                  ! neighbor stays on his level or want to coarsen,
                                  ! so block can not coarsen
                                  ! TODO: block can coarsen, if neighbor also want to coarsen
                                  my_block_list( rank*N + k, max_treelevel+2 ) = max(0, my_block_list( rank*N + k, max_treelevel+2 ) )
                                  grid_changed = .true. ! we changed something in the grid
                              end if

                            else
                              ! error case
                              print*, "ERROR: block has wrong number of neighbors!"
                              stop

                            end if

                        end if
                    end do

                end if

            end if

        end do

    counter = counter + 1
    if (counter == 99) then
        print*, "ERROR: unable to build a graded mesh"
        stop
    end if

    ! end do of repeat procedure until grid_changed==.false.
    end do

    ! synchronize light data
    block_list = 0
    call MPI_Allreduce(my_block_list, block_list, size(block_list,1)*size(block_list,2), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

end subroutine ensure_gradedness
