!> If some BC are non-periodic, we have to remove their neighborhood relationship
!! by setting it to -1. Otherwise, the sync_ghosts routine will overwrite the redundant points,
!! which is an error in this case.
subroutine remove_nonperiodic_neighbors(params, tree_ID, verbose_check)
    implicit none

    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)        :: tree_ID
    logical, optional, intent(in)       :: verbose_check       !< Output verbose flag


    integer(kind=ik) :: k, hvy_id, lgt_id, i_n
    logical          :: remove
    integer(kind=2)  :: n_domain(1:3)

    do k = 1, hvy_n(tree_ID)
        ! the block we're looking at ...
        hvy_id = hvy_active(k, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )

        ! is this an interior block ? (note: this is completely equivalent to checking if a neighborhood crosses the periodic domain,
        ! because only in this case, n_domain has nonzero value)
        if (maxval(abs(n_domain(1:params%dim)))==0) cycle

        do i_n = 1, size(hvy_neighbor,2)
            ! Always reset variables that we set each loop, this is not the first time we forgot it
            remove = .False.

            ! if this neighborhood is unused (there is no neighbor in this direction): skip it
            if (hvy_neighbor(hvy_id, i_n) < 0) cycle

            ! this checks in -x, +x, -y, +y, -z and +z direction if this is non-periodic and should be removed
            if (any(mod(i_n, 56) == (/  1, 2, 3, 4,  25,26,29,30,33,34,37,38,  49,51,53,55 /)) .and. .not.params%periodic_BC(1) .and. n_domain(1)==-1) remove = .True.
            if (any(mod(i_n, 56) == (/  5, 6, 7, 8,  27,28,31,32,35,36,39,40,  50,52,54,56 /)) .and. .not.params%periodic_BC(1) .and. n_domain(1)==+1) remove = .True.
            if (any(mod(i_n, 56) == (/  9,10,11,12,  25,26,27,28,41,42,45,46,  49,50,53,54 /)) .and. .not.params%periodic_BC(2) .and. n_domain(2)==-1) remove = .True.
            if (any(mod(i_n, 56) == (/ 13,14,15,16,  29,30,31,32,43,44,47,48,  51,52,55,56 /)) .and. .not.params%periodic_BC(2) .and. n_domain(2)==+1) remove = .True.
            if (any(mod(i_n, 56) == (/ 17,18,19,20,  33,34,35,36,41,42,43,44,  49,50,51,52 /)) .and. .not.params%periodic_BC(3) .and. n_domain(3)==-1) remove = .True.
            if (any(mod(i_n, 56) == (/ 21,22,23,24,  37,38,39,40,45,46,47,48,  53,54,55,56 /)) .and. .not.params%periodic_BC(3) .and. n_domain(3)==+1) remove = .True.  

            if (remove) then
                ! deactivate this neighborhood relation
                hvy_neighbor(hvy_id, i_n) = -1
            endif
        enddo
    enddo
end subroutine
