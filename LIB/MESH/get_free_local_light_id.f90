!> \brief Locally (on the calling processor only) look for a free block and return the ID. If full, error.
!
!> \details Unlike get_free_light_id, this routine returns a free slot on a given proc only.
! ********************************************************************************************

subroutine get_free_local_light_id( params, mpirank, lgt_free_id, ignore_error )

    use module_params

    implicit none

    type (type_params), intent(in)                :: params             !> user defined parameter structure
    integer(kind=ik), intent(in)                  :: mpirank            !> On what rank do you want a free block ID?
    integer(kind=ik), intent(out)                 :: lgt_free_id        !> return a free ID in the range of the proc "mpirank"
    !> if there are no more free light ids, the code aborts, unless use set ignore_error=.true
    !> in which case it will return the -1 ligt id
    logical, optional, intent(in)                 :: ignore_error

    integer(kind=ik) :: k, first_light_id, last_light_id, i, tree_ID             ! local variables
    logical :: valid, ign_err, active_lists2

    lgt_free_id = -1

    ! first light of the mpirank we are looking at. Note mpirank is zero based but
    ! list is one based, so add 1
    call proc_to_lgt_data_start_id( first_light_id, mpirank, params%number_blocks )
    last_light_id = first_light_id + (params%number_blocks - 1)

    ! loop over the range of light IDs belonging to proc "mpirank"
    do k = first_light_id, last_light_id
        ! check: if the block is not active, then we found a free block to return
        if ( lgt_block(k,1) == -1 ) then
            lgt_free_id = k
            exit
        end if
    end do

    ign_err = .false.
    if (present(ignore_error)) ign_err = ignore_error

    ! if desired, the code will not exit here but instead return an invalid light
    ! id (-1)
    if ( .not. ign_err ) then
        ! error catching: is there no more free blocks on the list?
        if (lgt_free_id == -1) then
            write(*,'("rank=",i5)') params%rank
            call abort(4458110, "ERROR: We try to fetch a light free block ID from the list but all blocks are used on this CPU")
        end if
    endif

end subroutine get_free_local_light_id
