! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: get_free_light_id.f90
! version: 0.4
! author: msr, engels
!
! return light free block id
!
! input:    - first column of light data array
!           - length of input vector
! output:   - id of first passive block
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
!
! ********************************************************************************************

subroutine get_free_light_id( id, block_list, N )

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    integer(kind=ik), intent(in)        :: N
    ! light data array, first column
    integer(kind=ik), intent(in)        :: block_list(N)

    ! free id
    integer(kind=ik), intent(out)       :: id

    ! loop variables
    integer(kind=ik)                    :: k

!---------------------------------------------------------------------------------------------
! variables initialization

    id = -1

!---------------------------------------------------------------------------------------------
! main body

    ! loop over list
    do k = 1, N
      ! check: if the block is not active, then we found a free block to return
      if ( block_list(k) == -1 ) then
        id = k
        exit
      end if
    end do

    ! error catching: is there no more free blocks on the list?
    if (id == -1) then
      write(*,*) "ERROR: We try to fetch a light free block ID from the list but all blocks are used."
      stop
    end if

end subroutine get_free_light_id
