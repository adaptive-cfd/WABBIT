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

            call get_adjacent_boundary_surface_normal( lgt_block(lgtID, 1:lgt_block(lgtID,params%Jmax+IDX_MESH_LVL)), &
            params%domain_size, params%Bs, params%dim, n_domain )

            ! faces
            call find_neighbor(params, hvyID, lgtID, Jmax, '__1/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '__2/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '__3/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '__4/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '__5/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '__6/___', error, n_domain)

            ! edges
            call find_neighbor(params, hvyID, lgtID, Jmax, '_12/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '_13/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '_14/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '_15/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '_62/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '_63/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '_64/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '_65/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '_23/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '_25/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '_43/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '_45/___', error, n_domain)

            ! corners
            call find_neighbor(params, hvyID, lgtID, Jmax, '123/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '134/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '145/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '152/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '623/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '634/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '645/___', error, n_domain)
            call find_neighbor(params, hvyID, lgtID, Jmax, '652/___', error, n_domain)
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

            call get_adjacent_boundary_surface_normal( lgt_block(lgtID, 1:lgt_block(lgtID,params%Jmax+IDX_MESH_LVL)), &
            params%domain_size, params%Bs, params%dim, n_domain )

            call find_neighbor( params, hvyID, lgtID, Jmax, '__N', error, n_domain)
            call find_neighbor( params, hvyID, lgtID, Jmax, '__E', error, n_domain)
            call find_neighbor( params, hvyID, lgtID, Jmax, '__S', error, n_domain)
            call find_neighbor( params, hvyID, lgtID, Jmax, '__W', error, n_domain)

            call find_neighbor( params, hvyID, lgtID, Jmax, '_NE', error, n_domain)
            call find_neighbor( params, hvyID, lgtID, Jmax, '_NW', error, n_domain)
            call find_neighbor( params, hvyID, lgtID, Jmax, '_SE', error, n_domain)
            call find_neighbor( params, hvyID, lgtID, Jmax, '_SW', error, n_domain)
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
