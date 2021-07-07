!> \brief Locally (on the calling processor only) look for a free block and return the ID. If full, error.
!
!> \details Unlike get_free_light_id, this routine returns a free slot on a given proc only. Note
!! it can sometimes be useful to return a free block which also is not on the active list. This may
!! seem odd (since by definition you shouldn't find any inactive block on the active list) but in
!! routines that modify the mesh (especially deleting blocks), the active list gets outdated.
! ********************************************************************************************

subroutine get_free_local_light_id( params, mpirank, lgt_block, lgt_free_id, lgt_active, lgt_n, ignore_error )

    use module_params

    implicit none

    type (type_params), intent(in)                :: params             !> user defined parameter structure
    integer(kind=ik), intent(in)                  :: mpirank            !> On what rank do you want a free block ID?
    integer(kind=ik), intent(inout)               :: lgt_block(:, :)    !> light data array (global list, redundant on mpi)
    integer(kind=ik), intent(out)                 :: lgt_free_id        !> return a free ID in the range of the proc "mpirank"
    integer(kind=ik), optional, intent(inout)     :: lgt_active(:)      !> list of active blocks (light data)
    integer(kind=ik), optional, intent(inout)     :: lgt_n              !> number of active blocks (light data)
    !> if there are no more free light ids, the code aborts, unless use set ignore_error=.true
    !> in which case it will return the -1 ligt id
    logical, optional, intent(in)                 :: ignore_error

    integer(kind=ik) :: k, first_light_id, last_light_id, i             ! local variables
    logical :: valid, ign_err

    lgt_free_id = -1

    ! first light of the mpirank we are looking at. Note mpirank is zero based but
    ! list is one based, so add 1
    call proc_to_lgt_data_start_id( first_light_id, mpirank, params%number_blocks )
    last_light_id = first_light_id + (params%number_blocks - 1)

    if (present(lgt_n) .and. present(lgt_active)) then
        ! if an active list is given, extra care is taken to return a free ID
        ! WHICH IS NOT ON THAT LIST. This may seem odd, but is required if the grid
        ! changes in a routine - it may then happen that a new block slipps into
        ! the active list.

        ! loop over the range of light IDs belonging to proc "mpirank"
        do k = first_light_id, last_light_id
            ! check: if the block is not active, then we found a free block to return
            if ( lgt_block(k,1) == -1 ) then
                ! we found d block inactive, but check the active list
                valid = .true.
                ! check if point is on active list:
                do i = 1, lgt_n
                    if ( lgt_active(i) == k) then
                        valid = .false.
                        exit
                    endif
                enddo
                if (valid) then
                    lgt_free_id = k
                    exit
                endif
            end if
        end do

    else
        ! without active list in the call, just take any inactive block

        ! loop over the range of light IDs belonging to proc "mpirank"
        do k = first_light_id, last_light_id
            ! check: if the block is not active, then we found a free block to return
            if ( lgt_block(k,1) == -1 ) then
                lgt_free_id = k
                exit
            end if
        end do

    endif

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
