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
!> \image html update_neighbors2D.svg "Update Neighbors" width=400

subroutine update_neighbors_2D(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist,  &
    hvy_active, hvy_n, error, skip_diagonal_neighbors)
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
    logical, intent(inout)              :: error
    logical, intent(in), optional       :: skip_diagonal_neighbors ! currently not working (Thomas, 02-2021) [unused]

    ! number of blocks per proc
    integer(kind=ik)                    :: N
    ! max treelevel
    integer(kind=ik)                    :: max_treelevel

    ! process rank
    integer(kind=ik)                    :: rank

    ! loop variable
    integer(kind=ik)                    :: k, lgt_id
    integer(kind=2) :: n_domain(1:3)

    rank = params%rank
    max_treelevel = params%max_treelevel
    ! number of blocks per proc
    N = params%number_blocks

    ! special case:
    ! if there is only one block => all neighbors are this block
    ! one block criteria: lgt_n should be one!
    if ( lgt_n == 1 ) then
        hvy_neighbor(1,1:8) = lgt_active(1)
        return
    end if

    ! loop over active heavy data blocks
    do k = 1, hvy_n
        ! delete existing neighbors (they are assumed to be outdated when calling this routine)
        hvy_neighbor(hvy_active(k), : ) = -1

        ! light id
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )

        call get_adjacent_boundary_surface_normal( lgt_block(lgt_id, 1:lgt_block(lgt_id,params%max_treelevel+IDX_MESH_LVL)), &
        params%domain_size, params%Bs, params%dim, n_domain )

        ! north
        call find_neighbor_edge_2D( params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '__N', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
        ! east
        call find_neighbor_edge_2D( params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '__E', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
        ! south
        call find_neighbor_edge_2D( params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '__S', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
        ! west
        call find_neighbor_edge_2D( params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '__W', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)

        if (PRESENT(skip_diagonal_neighbors)) then
            if (skip_diagonal_neighbors) cycle
        endif

        ! northeast
        call find_neighbor_corner_2D( params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_NE', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
        ! northwest
        call find_neighbor_corner_2D( params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_NW', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
        ! southeast
        call find_neighbor_corner_2D( params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_SE', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
        ! southwest
        call find_neighbor_corner_2D( params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_SW', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)

    end do

end subroutine update_neighbors_2D
