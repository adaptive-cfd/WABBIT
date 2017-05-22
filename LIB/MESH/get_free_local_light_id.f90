!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name get_free_local_light_id.f90
!> \version 0.5
!> \author engels
!
!> \brief Locally (on the calling processor only) look for a free block and return the ID. If full, error.
!
!> \details Unlike get_free_light_id, this routine returns a free slot on a given proc only. Note
!! it can sometimes be useful to return a free block which also is not on the active list. This may
!! seem odd (since by definition you shouldn't find any inactive block on the active list) but in
!! routines that modify the mesh (especially deleting blocks), the active list gets outdated.
!!
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4
!
! ********************************************************************************************

subroutine get_free_local_light_id( params, mpirank, lgt_block, lgt_free_id, lgt_active, lgt_n )

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> On what rank do you want a free block ID?
    integer(kind=ik), intent(in)        :: mpirank
    !> light data array (global list, redundant on mpi)
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> return a free ID in the range of the proc "mpirank"
    integer(kind=ik), intent(out)       :: lgt_free_id
    !> list of active blocks (light data)
    integer(kind=ik), optional, intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), optional, intent(inout)     :: lgt_n

    ! local variables
    integer(kind=ik) :: k, first_light_id, last_light_id, i
    logical :: valid

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

    ! error catching: is there no more free blocks on the list?
    if (lgt_free_id == -1) then
      write(*,*) "ERROR: 4458110: We try to fetch a light free block ID from the list but all blocks are used on this CPU"
      stop
    end if

end subroutine get_free_local_light_id
