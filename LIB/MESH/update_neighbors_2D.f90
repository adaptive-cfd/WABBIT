!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name update_neighbors_2D.f90
!> \version 0.4
!> \author msr
!
!> \brief update neighbor relations with light data, store neighbors in neighbor list (heavy data)
!!        2D version
!
!> \details
!! input:
!!           - light data array
!!           - params struct
!!
!! output:
!!           - neighbor list array
!!
!! = log ======================================================================================
!! \n
!! 07/11/16 - switch to v0.4 \n
!! 27/01/17 - use process rank from params struct
!
! ********************************************************************************************
!> \image html update_neighbors.png "Update Neighbors" width=400

subroutine update_neighbors_2D(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist,  hvy_active, hvy_n)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(out)       :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(in)     :: lgt_sortednumlist(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! number of blocks per proc
    integer(kind=ik)                    :: N
    ! max treelevel
    integer(kind=ik)                    :: max_treelevel

    ! process rank
    integer(kind=ik)                    :: rank

    ! loop variable
    integer(kind=ik)                    :: k, lgt_id, i,j

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1


!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank         = params%rank

    ! reset neighbor list
    hvy_neighbor = -1

    ! number of blocks per proc
    N = params%number_blocks

    ! max treelevel
    max_treelevel = params%max_treelevel

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! special case:
    ! if there is only one block => all neighbors are this block
    ! one block criteria: size of block_list should be one!
    !> \todo Wouldn't lgt_n==1 be a better criterion for this exception?
    if ( size(lgt_block,1) == 1 ) then
        hvy_neighbor(1,1:8) = 1
    end if

    ! loop over active heavy data blocks
    do k = 1, hvy_n

        ! light id
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )

        ! north
        call find_neighbor_edge_2D( hvy_active(k), lgt_id, lgt_block, max_treelevel, '__N', hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist)
        ! east
        call find_neighbor_edge_2D( hvy_active(k), lgt_id, lgt_block, max_treelevel, '__E', hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist)
        ! south
        call find_neighbor_edge_2D( hvy_active(k), lgt_id, lgt_block, max_treelevel, '__S', hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist)
        ! west
        call find_neighbor_edge_2D( hvy_active(k), lgt_id, lgt_block, max_treelevel, '__W', hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist)

        ! northeast
        call find_neighbor_corner_2D( hvy_active(k), lgt_id, lgt_block, max_treelevel, '_NE', hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist)
        ! northwest
        call find_neighbor_corner_2D( hvy_active(k), lgt_id, lgt_block, max_treelevel, '_NW', hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist)
        ! southeast
        call find_neighbor_corner_2D( hvy_active(k), lgt_id, lgt_block, max_treelevel, '_SE', hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist)
        ! southwest
        call find_neighbor_corner_2D( hvy_active(k), lgt_id, lgt_block, max_treelevel, '_SW', hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist)

    end do

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "update_neighbors" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "update_neighbors"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine update_neighbors_2D
