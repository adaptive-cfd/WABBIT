!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name get_free_light_id.f90
!> \version 0.4
!> \author msr, engels
!
!> \brief return light free block id
!
!> \details
!! input:
!!           - first column of light data array
!!           - length of input vector
!!
!! output:
!!           - id of first passive block
!!
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4
!
! ********************************************************************************************

subroutine get_free_light_id( id, lgt_block, N )

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    integer(kind=ik), intent(in)        :: N
    !> light data array, first column only (TODO FIXME: explain why??? does that make any difference??) (Thomas, 13/03/2018)
    integer(kind=ik), intent(in)        :: lgt_block(N)

    !> free id
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
      if ( lgt_block(k) == -1 ) then
        id = k
        exit
      end if
    end do

    ! error catching: is there no more free blocks on the list?
    if (id == -1) then
      call abort(6363709, "ERROR: We try to fetch a light free block ID from the list but all blocks are used.")
    end if

end subroutine get_free_light_id
