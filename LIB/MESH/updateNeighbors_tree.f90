!> \brief functional wrapper for the 2D and 3D version of update neighbors.
!! input:     - light data array
!!            - params struct
!! output:    - neighbor list array
! ********************************************************************************************
subroutine updateNeighbors_tree(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, tree_ID)

    implicit none
    type (type_params), intent(in)      :: params                   !> user defined parameter structure
    integer(kind=ik), intent(in)        :: lgt_block(:, :)          !> light data array
    integer(kind=ik), intent(out)       :: hvy_neighbor(:,:)        !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: lgt_active(:,:)          !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n(:)                 !> number of active blocks (light data)
    integer(kind=tsize), intent(in)     :: lgt_sortednumlist(:,:,:) !> sorted list of numerical treecodes, used for block finding
    integer(kind=ik), intent(in)        :: hvy_active(:,:)          !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n(:)                 !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: tree_ID
    logical                             :: error = .false.
    integer(kind=ik)                    :: mpierror, k, N, Jmax, lgtID, hvyID
    real(kind=rk)                       :: x0(1:3), dx(1:3)
    integer(kind=2)                     :: n_domain(1:3)


    N = params%number_blocks
    Jmax = params%max_treelevel


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

            call get_adjacent_boundary_surface_normal( lgt_block(lgtID, 1:lgt_block(lgtID,params%max_treelevel+IDX_MESH_LVL)), &
            params%domain_size, params%Bs, params%dim, n_domain )

            ! faces
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '__1/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '__2/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '__3/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '__4/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '__5/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '__6/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)

            ! edges
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '_12/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '_13/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '_14/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '_15/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '_62/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '_63/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '_64/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '_65/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '_23/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '_25/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '_43/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '_45/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)

            ! corners
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '123/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '134/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '145/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '152/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '623/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '634/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '645/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor(params, hvyID, lgtID, lgt_block, Jmax, '652/___', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
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

            call get_adjacent_boundary_surface_normal( lgt_block(lgtID, 1:lgt_block(lgtID,params%max_treelevel+IDX_MESH_LVL)), &
            params%domain_size, params%Bs, params%dim, n_domain )

            call find_neighbor( params, hvyID, lgtID, lgt_block, Jmax, '__N', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor( params, hvyID, lgtID, lgt_block, Jmax, '__E', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor( params, hvyID, lgtID, lgt_block, Jmax, '__S', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor( params, hvyID, lgtID, lgt_block, Jmax, '__W', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)

            call find_neighbor( params, hvyID, lgtID, lgt_block, Jmax, '_NE', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor( params, hvyID, lgtID, lgt_block, Jmax, '_NW', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor( params, hvyID, lgtID, lgt_block, Jmax, '_SE', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
            call find_neighbor( params, hvyID, lgtID, lgt_block, Jmax, '_SW', hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)
        end do

    end if

    if (error) then
        call abort(71737, "Grid error: we did not find a neighboring block. Either your data is corrupt (reading from file) or you found a bug (grid generation failed).")
    endif

    ! Is there any non-periodic boundary ?
    if ( .not. All(params%periodic_BC) ) then
        call remove_nonperiodic_neighbors(params, lgt_block, hvy_neighbor, hvy_active, hvy_n, tree_ID)
    endif

end subroutine updateNeighbors_tree
