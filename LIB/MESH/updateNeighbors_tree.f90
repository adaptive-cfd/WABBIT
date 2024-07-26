!> \brief functional wrapper for the 2D and 3D version of update neighbors.
!! input:     - light data array
!!            - params struct
!! output:    - neighbor list array
! ********************************************************************************************
subroutine updateNeighbors_tree(params, tree_ID)

    implicit none
    type (type_params), intent(in)      :: params                   !> user defined parameter structure
    integer(kind=ik), intent(in)        :: tree_ID
    logical                             :: error = .false.
    integer(kind=ik)                    :: mpierror, k, N, Jmax, lgtID, hvyID
    real(kind=rk)                       :: x0(1:3), dx(1:3)
    integer(kind=2)                     :: n_domain(1:3)


    N = params%number_blocks
    Jmax = params%Jmax

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
        hvyID = hvy_active(k, tree_ID)

        call hvy2lgt( lgtID, hvyID, params%rank, N )

        call get_adjacent_boundary_surface_normal( params, lgtID, n_domain )

        !> direction for neighbor search - number where each digit represents a cardinal direction XYZ
        !> 9 -> lower direction (-)
        !> 1 -> upper direction (+)
        !> 0 -> no direction (but finer neighbors can vary in + and - and coarser neighbors have one configuration of those two) 

        ! faces, 2D edges
        call find_neighbor(params, hvyID, lgtID, Jmax, 900, error, n_domain)  ! -x
        call find_neighbor(params, hvyID, lgtID, Jmax, 100, error, n_domain)  ! +x
        call find_neighbor(params, hvyID, lgtID, Jmax, 090, error, n_domain)  ! -y
        call find_neighbor(params, hvyID, lgtID, Jmax, 010, error, n_domain)  ! +y
        if (params%dim == 3) then
            call find_neighbor(params, hvyID, lgtID, Jmax, 009, error, n_domain)  ! -z
            call find_neighbor(params, hvyID, lgtID, Jmax, 001, error, n_domain)  ! +z
        endif
        ! edges, 2D corners
        call find_neighbor(params, hvyID, lgtID, Jmax, 990, error, n_domain)  ! -x-y
        call find_neighbor(params, hvyID, lgtID, Jmax, 910, error, n_domain)  ! -x+y
        call find_neighbor(params, hvyID, lgtID, Jmax, 190, error, n_domain)  ! +x-y
        call find_neighbor(params, hvyID, lgtID, Jmax, 110, error, n_domain)  ! +x+x
        if (params%dim == 3) then
            call find_neighbor(params, hvyID, lgtID, Jmax, 909, error, n_domain)  ! -x-z
            call find_neighbor(params, hvyID, lgtID, Jmax, 901, error, n_domain)  ! -x+z
            call find_neighbor(params, hvyID, lgtID, Jmax, 109, error, n_domain)  ! +x-z
            call find_neighbor(params, hvyID, lgtID, Jmax, 101, error, n_domain)  ! +x+z
            call find_neighbor(params, hvyID, lgtID, Jmax, 099, error, n_domain)  ! -y-z
            call find_neighbor(params, hvyID, lgtID, Jmax, 091, error, n_domain)  ! -y+z
            call find_neighbor(params, hvyID, lgtID, Jmax, 019, error, n_domain)  ! +y-z
            call find_neighbor(params, hvyID, lgtID, Jmax, 011, error, n_domain)  ! +y+z
            ! ecorners
            call find_neighbor(params, hvyID, lgtID, Jmax, 999, error, n_domain)  ! -x-y-z
            call find_neighbor(params, hvyID, lgtID, Jmax, 199, error, n_domain)  ! +x-y-z
            call find_neighbor(params, hvyID, lgtID, Jmax, 919, error, n_domain)  ! -x+y-z
            call find_neighbor(params, hvyID, lgtID, Jmax, 119, error, n_domain)  ! +x+y-z
            call find_neighbor(params, hvyID, lgtID, Jmax, 991, error, n_domain)  ! -x-y+z
            call find_neighbor(params, hvyID, lgtID, Jmax, 191, error, n_domain)  ! +x-y+z
            call find_neighbor(params, hvyID, lgtID, Jmax, 911, error, n_domain)  ! -x+y+z
            call find_neighbor(params, hvyID, lgtID, Jmax, 111, error, n_domain)  ! +x+y+z
        endif
    end do

    if (error) then
        call abort(71737, "Grid error: we did not find a neighboring block. Either your data is corrupt (reading from file) or you found a bug (grid generation failed).")
    endif

    ! Is there any non-periodic boundary ?
    if ( .not. All(params%periodic_BC) ) then
        call remove_nonperiodic_neighbors(params, tree_ID)
    endif

end subroutine updateNeighbors_tree
