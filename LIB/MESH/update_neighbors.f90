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
!           - params struct
! output:   - neighbor list array
!
! = log ======================================================================================
!
! 07/11/16 - switch to v0.4
! ********************************************************************************************

subroutine update_neighbors(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    ! heavy data array - neifghbor data
    integer(kind=ik), intent(out)       :: hvy_neighbor(:)

    ! list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    ! list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    ! number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! number of blocks per proc
    integer(kind=ik)                    :: N
    ! max treelevel
    integer(kind=ik)                    :: max_treelevel

    ! MPI error variable
    integer(kind=ik)                    :: ierr,a
    ! process rank
    integer(kind=ik)                    :: rank

    ! loop variable
    integer(kind=ik)                    :: k, block_number

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset neighbor list
    hvy_neighbor = -1

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! special case:
    ! if there is only one block => all neighbors are this block
    ! one block criteria: size of block_list should be one!
    if ( size(lgt_block,1) == 1 ) then
        hvy_neighbor(1:8) = 1
    end if

    ! loop over all heavy data blocks
    do k = 1, N

!        ! block is active
!        if ( block_list( rank*N + k , 1) /= -1 ) then
!
!            ! find edge neighbors
!            ! north
!            call find_neighbor_edge( k, rank*N + k, block_list, max_treelevel, '__N', neighbor_list, active)
!            ! east
!            call find_neighbor_edge( k, rank*N + k, block_list, max_treelevel, '__E', neighbor_list, active)
!            ! south
!            call find_neighbor_edge( k, rank*N + k, block_list, max_treelevel, '__S', neighbor_list, active)
!            ! west
!            call find_neighbor_edge( k, rank*N + k, block_list, max_treelevel, '__W', neighbor_list, active)
!
!            ! find corner neighbor
!            ! northeast
!            call find_neighbor_corner( k, rank*N + k, block_list, max_treelevel, '_NE', neighbor_list, active)
!            ! northwest
!            call find_neighbor_corner( k, rank*N + k, block_list, max_treelevel, '_NW', neighbor_list, active)
!            ! southeast
!            call find_neighbor_corner( k, rank*N + k, block_list, max_treelevel, '_SE', neighbor_list, active)
!            ! southwest
!            call find_neighbor_corner( k, rank*N + k, block_list, max_treelevel, '_SW', neighbor_list, active)
!
!        end if

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
        debug%name_comp_time(k) = "update_neighbors"
        debug%comp_time(rank+1, k) = sub_t1 - sub_t0
    end if

end subroutine update_neighbors
