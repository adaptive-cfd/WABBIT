! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: update_neighbors.f90
! version: 0.4
! author: msr
!
! update neighbor relations with light data, store neighbors in neighbor list (heavy data)
!
! input:    - light data array
!           - number of blocks per proc
! output:   - neighbor list array
!
! = log ======================================================================================
!
! 07/11/16 - switch to v0.4
! ********************************************************************************************

subroutine update_neighbors(block_list, neighbor_list, N, max_treelevel)

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! light data array
    integer(kind=ik), intent(in)        :: block_list(:, :)
    ! heavy data array - neifghbor data
    integer(kind=ik), intent(out)       :: neighbor_list(:)
    ! number of blocks per proc
    integer(kind=ik), intent(in)        :: N
    ! max treelevel
    integer(kind=ik), intent(in)        :: max_treelevel

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! loop variable
    integer(kind=ik)                    :: k

!---------------------------------------------------------------------------------------------
! interfaces

    interface
        subroutine find_neighbor_edge(heavy_id, light_id, block_list, max_treelevel, dir, neighbor_list)
            use module_params
            integer(kind=ik), intent(in)                :: heavy_id
            integer(kind=ik), intent(in)                :: light_id
            integer(kind=ik), intent(in)                :: max_treelevel
            integer(kind=ik), intent(in)                :: block_list(:, :)
            character(len=3), intent(in)                :: dir
            integer(kind=ik), intent(out)               :: neighbor_list(:)
        end subroutine find_neighbor_edge

        subroutine find_neighbor_corner(heavy_id, light_id, block_list, max_treelevel, dir, neighbor_list)
            use module_params
            integer(kind=ik), intent(in)                :: heavy_id
            integer(kind=ik), intent(in)                :: light_id
            integer(kind=ik), intent(in)                :: max_treelevel
            integer(kind=ik), intent(in)                :: block_list(:, :)
            character(len=3), intent(in)                :: dir
            integer(kind=ik), intent(out)               :: neighbor_list(:)
        end subroutine find_neighbor_corner

    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset neighbor list
    neighbor_list = -1

!---------------------------------------------------------------------------------------------
! main body

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! special case:
    ! if there is only one block => all neighbors are this block
    ! one block criteria: size of block_list should be one!
    if ( size(block_list,1) == 1 ) then
        neighbor_list(1:8) = 1
    end if

    ! loop over all heavy data blocks
    do k = 1, N

        ! block is active
        if ( block_list( rank*N + k , 1) /= -1 ) then

            ! find edge neighbors
            ! north
            call find_neighbor_edge( k, rank*N + k, block_list, max_treelevel, '__N', neighbor_list)
            ! east
            call find_neighbor_edge( k, rank*N + k, block_list, max_treelevel, '__E', neighbor_list)
            ! south
            call find_neighbor_edge( k, rank*N + k, block_list, max_treelevel, '__S', neighbor_list)
            ! west
            call find_neighbor_edge( k, rank*N + k, block_list, max_treelevel, '__W', neighbor_list)

            ! find corner neighbor
            ! northeast
            call find_neighbor_corner( k, rank*N + k, block_list, max_treelevel, '_NE', neighbor_list)
            ! northwest
            call find_neighbor_corner( k, rank*N + k, block_list, max_treelevel, '_NW', neighbor_list)
            ! southeast
            call find_neighbor_corner( k, rank*N + k, block_list, max_treelevel, '_SE', neighbor_list)
            ! southwest
            call find_neighbor_corner( k, rank*N + k, block_list, max_treelevel, '_SW', neighbor_list)

        end if

    end do

    ! barrier
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

end subroutine update_neighbors