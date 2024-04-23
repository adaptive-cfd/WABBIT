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


    if ( params%dim == 3 ) then
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !               3D
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

            !> direction for neighbor search - number where each digit represents a cardinal direction
            !> 652 -> first 6 (bottom, z-1), then 5 (north, x-1) then 2 (front, y-1) 
            ! faces
            call find_neighbor(params, hvyID, lgtID, Jmax, 1, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 2, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 3, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 4, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 5, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 6, error, n_domain)

            ! edges
            call find_neighbor(params, hvyID, lgtID, Jmax, 12, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 13, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 14, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 15, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 62, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 63, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 64, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 65, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 23, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 25, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 43, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 45, error, n_domain)

            ! corners
            call find_neighbor(params, hvyID, lgtID, Jmax, 123, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 134, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 145, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 152, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 623, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 634, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 645, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, 652, error, n_domain)
        end do

    else
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !               2D
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! special case:
        ! if there is only one block => all neighbors are this block
        ! one block criteria: lgt_n should be one!
        if ( lgt_n(tree_ID) == 1 ) then
            hvy_neighbor(1,1:8) = lgt_active(1, tree_ID)
            return
        end if

        ! loop over active heavy data blocks
        do k = 1, hvy_n(tree_ID)
            ! delete existing neighbors (they are assumed to be outdated when calling this routine)
            hvy_neighbor(hvy_active(k, tree_ID), :) = -1
            hvyID = hvy_active(k, tree_ID)

            call hvy2lgt( lgtID, hvyID, params%rank, N )

            call get_adjacent_boundary_surface_normal( params, lgtID, n_domain )

            !> direction for neighbor search - number where each digit represents a cardinal direction
            !> 652 -> first 6 (bottom, z-1), then 5 (north, x-1) then 2 (front, y-1) 
            call find_neighbor( params, hvyID, lgtID, Jmax, 5, error, n_domain)  ! '__N'
            call find_neighbor( params, hvyID, lgtID, Jmax, 4, error, n_domain)  ! '__E'
            call find_neighbor( params, hvyID, lgtID, Jmax, 3, error, n_domain)  ! '__S'
            call find_neighbor( params, hvyID, lgtID, Jmax, 2, error, n_domain)  ! '__W'

            call find_neighbor( params, hvyID, lgtID, Jmax, 45, error, n_domain)  ! '_NE'
            call find_neighbor( params, hvyID, lgtID, Jmax, 25, error, n_domain)  ! '_NW'
            call find_neighbor( params, hvyID, lgtID, Jmax, 43, error, n_domain)  ! '_SE'
            call find_neighbor( params, hvyID, lgtID, Jmax, 23, error, n_domain)  ! '_SW'
        end do

    end if

    if (error) then
        call abort(71737, "Grid error: we did not find a neighboring block. Either your data is corrupt (reading from file) or you found a bug (grid generation failed).")
    endif

    ! Is there any non-periodic boundary ?
    if ( .not. All(params%periodic_BC) ) then
        call remove_nonperiodic_neighbors(params, tree_ID)
    endif

end subroutine updateNeighbors_tree
