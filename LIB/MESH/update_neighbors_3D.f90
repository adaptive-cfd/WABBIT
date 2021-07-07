!>\brief update neighbor relations with light data, store neighbors in neighbor list (heavy data)
!!       3D case
!! input:    - light data array
!!           - params struct
!! output:   - neighbor list array
!  --------------------------------------------------------------------------------------------
!> neighbor codes: \n
!  ---------------
!> for imagination:
!!                   - 6-sided dice with '1'-side on top, '6'-side on bottom, '2'-side in front
!!                   - edge: boundary between two sides - use sides numbers for coding
!!                   - corner: between three sides - so use all three sides numbers
!!                   - block on higher/lower level: block shares face/edge and one unique corner,
!!                     so use this corner code in second part of neighbor code
!!
!! faces:  '__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___' \n
!! edges:  '_12/___', '_13/___', '_14/___', '_15/___' \n
!!         '_62/___', '_63/___', '_64/___', '_65/___' \n
!!         '_23/___', '_25/___', '_43/___', '_45/___' \n
!! corner: '123/___', '134/___', '145/___', '152/___' \n
!!         '623/___', '634/___', '645/___', '652/___' \n
!! \n
!! complete neighbor code array, 74 possible neighbor relations \n
!! neighbors = (/'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___', '_12/___', '_13/___', '_14/___', '_15/___',
!!               '_62/___', '_63/___', '_64/___', '_65/___', '_23/___', '_25/___', '_43/___', '_45/___', '123/___', '134/___',
!!               '145/___', '152/___', '623/___', '634/___', '645/___', '652/___', '__1/123', '__1/134', '__1/145', '__1/152',
!!               '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634',
!!               '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652',
!!               '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652',
!!               '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652',
!!               '_43/134', '_43/634', '_45/145', '_45/645' /) \n
! ********************************************************************************************
!> \image html update_neighbors3D.svg "Update Neighbors" width=400

subroutine update_neighbors_3D(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_active, hvy_n, error, skip_diagonal_neighbors)

    implicit none

    type (type_params), intent(in)      :: params                   !> user defined parameter structure
    integer(kind=ik), intent(in)        :: lgt_block(:, :)          !> light data array
    integer(kind=ik), intent(out)       :: hvy_neighbor(:,:)        !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: lgt_active(:)            !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n                    !> number of active blocks (light data)
    integer(kind=tsize), intent(in)     :: lgt_sortednumlist(:,:)   !> sorted list of numerical treecodes, used for block finding
    integer(kind=ik), intent(in)        :: hvy_active(:)            !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n                    !> number of active blocks (heavy data)
    logical, intent(inout)              :: error
    logical, intent(in), optional       :: skip_diagonal_neighbors  ! currently not working (Thomas, 02-2021) [unused]
    integer(kind=ik)                    :: N                        ! number of blocks per proc
    integer(kind=ik)                    :: max_treelevel
    integer(kind=ik)                    :: rank                     ! process rank
    integer(kind=ik)                    :: k, lgt_id                ! loop variable
    integer(kind=2) :: n_domain(1:3)

    rank = params%rank

    ! number of blocks per proc
    N = params%number_blocks
    max_treelevel = params%max_treelevel

    ! special case:
    ! if there is only one block => all neighbors are this block
    ! one block criteria: lgt_n should be one!
    if ( lgt_n == 1 ) then
        hvy_neighbor(1,1:26) = lgt_active(1)
    end if

    ! loop only if more than one block is active
    if ( lgt_n /= 1 ) then
        ! loop over active heavy data blocks
        do k = 1, hvy_n
            ! delete existing neighbors (they are assumed to be outdated when calling this routine)
            hvy_neighbor(hvy_active(k), : ) = -1

            ! light id
            call hvy2lgt( lgt_id, hvy_active(k), rank, N )

            call get_adjacent_boundary_surface_normal( lgt_block(lgt_id, 1:lgt_block(lgt_id,params%max_treelevel+IDX_MESH_LVL)), &
            params%domain_size, params%Bs, params%dim, n_domain )

            ! faces
            ! '1'-side
            call find_neighbor_face_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '__1/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '2'-side
            call find_neighbor_face_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '__2/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '3'-side
            call find_neighbor_face_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '__3/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '4'-side
            call find_neighbor_face_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '__4/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '5'-side
            call find_neighbor_face_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '__5/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '6'-side
            call find_neighbor_face_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '__6/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)

            if (PRESENT(skip_diagonal_neighbors)) then
                if (skip_diagonal_neighbors) cycle
            endif

            ! edges
            ! '12'-edge
            call find_neighbor_edge_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_12/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '13'-edge
            call find_neighbor_edge_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_13/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '14'-edge
            call find_neighbor_edge_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_14/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '15'-edge
            call find_neighbor_edge_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_15/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '62'-edge
            call find_neighbor_edge_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_62/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '63'-edge
            call find_neighbor_edge_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_63/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '64'-edge
            call find_neighbor_edge_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_64/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '65'-edge
            call find_neighbor_edge_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_65/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '23'-edge
            call find_neighbor_edge_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_23/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '25'-edge
            call find_neighbor_edge_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_25/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '43'-edge
            call find_neighbor_edge_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_43/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '45'-edge
            call find_neighbor_edge_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '_45/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)

            ! corners
            ! '123'-corner
            call find_neighbor_corner_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '123/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '134'-corner
            call find_neighbor_corner_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '134/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '145'-corner
            call find_neighbor_corner_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '145/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '152'-corner
            call find_neighbor_corner_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '152/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '623'-corner
            call find_neighbor_corner_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '623/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '634'-corner
            call find_neighbor_corner_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '634/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '645'-corner
            call find_neighbor_corner_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '645/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            ! '652'-corner
            call find_neighbor_corner_3D(params, hvy_active(k), lgt_id, lgt_block, max_treelevel, '652/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
        end do
    end if

end subroutine update_neighbors_3D
