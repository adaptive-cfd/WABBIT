!> \brief functional wrapper for the 2D and 3D version of update neighbors.
!! input:     - light data array
!!            - params struct
!! output:    - neighbor list array
! ********************************************************************************************
subroutine updateNeighbors_tree(params, tree_ID, search_overlapping, verbose_check)

    implicit none
    type (type_params), intent(in)      :: params               !< user defined parameter structure
    integer(kind=ik), intent(in)        :: tree_ID
    logical, intent(in)                 :: search_overlapping   !< for CVS multiple neighbors can coexist, so we search all of them
    logical, optional, intent(in)  :: verbose_check       !< Output verbose flag


    logical                             :: error = .false.
    integer(kind=ik)                    :: mpierror, k, N, lgt_ID, hvy_ID
    real(kind=rk)                       :: x0(1:3), dx(1:3)
    integer(kind=2)                     :: n_domain(1:3)


    N = params%number_blocks

    ! special case:
    ! if there is only one block => all neighbors are this block
    ! one block criteria: lgt_n should be one!
    if ( lgt_n(tree_ID) == 1 ) then
        hvy_neighbor(1,1:26) = lgt_active(1, tree_ID)
        return
    end if

    ! loop over active heavy data blocks
    do k = 1, hvy_n(tree_ID)
        ! delete existing neighbors (they are assumed to be outdated when calling this routine)
        hvy_neighbor(hvy_active(k, tree_ID), :) = -1
        hvy_ID = hvy_active(k, tree_ID)

        call hvy2lgt( lgt_ID, hvy_ID, params%rank, N )

        call get_adjacent_boundary_surface_normal( params, lgt_ID, n_domain )

        !> direction for neighbor search - number where each digit represents a cardinal direction XYZ
        !> 9 -> lower direction (-) - we would love to use -1 here but it's obviously not possible so we settled for 9 instead (kind of like 10-1)
        !> 1 -> upper direction (+)
        !> 0 -> no direction (but finer neighbors can vary in + and - and coarser neighbors have one configuration of those two) 

        ! faces, 2D edges
        call find_neighbor(params, hvy_ID, lgt_ID, 900, error, n_domain, search_overlapping)  ! -x
        call find_neighbor(params, hvy_ID, lgt_ID, 100, error, n_domain, search_overlapping)  ! +x
        call find_neighbor(params, hvy_ID, lgt_ID, 090, error, n_domain, search_overlapping)  ! -y
        call find_neighbor(params, hvy_ID, lgt_ID, 010, error, n_domain, search_overlapping)  ! +y
        if (params%dim == 3) then
            call find_neighbor(params, hvy_ID, lgt_ID, 009, error, n_domain, search_overlapping)  ! -z
            call find_neighbor(params, hvy_ID, lgt_ID, 001, error, n_domain, search_overlapping)  ! +z
        endif
        ! edges, 2D corners
        call find_neighbor(params, hvy_ID, lgt_ID, 990, error, n_domain, search_overlapping)  ! -x-y
        call find_neighbor(params, hvy_ID, lgt_ID, 910, error, n_domain, search_overlapping)  ! -x+y
        call find_neighbor(params, hvy_ID, lgt_ID, 190, error, n_domain, search_overlapping)  ! +x-y
        call find_neighbor(params, hvy_ID, lgt_ID, 110, error, n_domain, search_overlapping)  ! +x+x
        if (params%dim == 3) then
            call find_neighbor(params, hvy_ID, lgt_ID, 909, error, n_domain, search_overlapping)  ! -x-z
            call find_neighbor(params, hvy_ID, lgt_ID, 901, error, n_domain, search_overlapping)  ! -x+z
            call find_neighbor(params, hvy_ID, lgt_ID, 109, error, n_domain, search_overlapping)  ! +x-z
            call find_neighbor(params, hvy_ID, lgt_ID, 101, error, n_domain, search_overlapping)  ! +x+z
            call find_neighbor(params, hvy_ID, lgt_ID, 099, error, n_domain, search_overlapping)  ! -y-z
            call find_neighbor(params, hvy_ID, lgt_ID, 091, error, n_domain, search_overlapping)  ! -y+z
            call find_neighbor(params, hvy_ID, lgt_ID, 019, error, n_domain, search_overlapping)  ! +y-z
            call find_neighbor(params, hvy_ID, lgt_ID, 011, error, n_domain, search_overlapping)  ! +y+z
            ! ecorners
            call find_neighbor(params, hvy_ID, lgt_ID, 999, error, n_domain, search_overlapping)  ! -x-y-z
            call find_neighbor(params, hvy_ID, lgt_ID, 199, error, n_domain, search_overlapping)  ! +x-y-z
            call find_neighbor(params, hvy_ID, lgt_ID, 919, error, n_domain, search_overlapping)  ! -x+y-z
            call find_neighbor(params, hvy_ID, lgt_ID, 119, error, n_domain, search_overlapping)  ! +x+y-z
            call find_neighbor(params, hvy_ID, lgt_ID, 991, error, n_domain, search_overlapping)  ! -x-y+z
            call find_neighbor(params, hvy_ID, lgt_ID, 191, error, n_domain, search_overlapping)  ! +x-y+z
            call find_neighbor(params, hvy_ID, lgt_ID, 911, error, n_domain, search_overlapping)  ! -x+y+z
            call find_neighbor(params, hvy_ID, lgt_ID, 111, error, n_domain, search_overlapping)  ! +x+y+z
        endif
    end do

    ! error is only set to true inside the functions and never to false, so if one encounters problems then this condition is active
    if (error) then
        call abort(71737, "Grid error: we did not find a neighboring block. Either your data is corrupt (reading from file) or you found a bug (grid generation failed).")
    endif

    ! Is there any non-periodic boundary ?
    if ( .not. All(params%periodic_BC) ) then
        call remove_nonperiodic_neighbors(params, tree_ID, verbose_check=.true.)
    endif

end subroutine updateNeighbors_tree
