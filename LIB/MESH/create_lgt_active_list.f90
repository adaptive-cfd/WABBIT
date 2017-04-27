!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name create_lgt_active_list.f90
!> \version 0.4
!> \author msr
!
!> \brief create a list of all active blocks (light data)
!> \details dim 1: light data id
!! \n
!! input:    
!!           - light data
!!
!! output:
!!           - list of active blocks, light data id
!!           - number of active blocks
!!
!> = log ======================================================================================
!! \n
!! 23/11/16 - create subroutine
! ********************************************************************************************

subroutine create_lgt_active_list( lgt_block, lgt_active, lgt_n )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)

    !> list of active blocks (light data)
    integer(kind=ik), intent(out)       :: lgt_active(:)

    !> number of active blocks (light data)
    integer(kind=ik), intent(out)       :: lgt_n

    ! loop variables
    integer(kind=ik)                    :: k

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset active list
    lgt_active = -1

    ! list index
    lgt_n = 0

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all light data
    do k = 1, size(lgt_block, 1)

        ! block is active
        if ( lgt_block(k, 1) /= -1 ) then
            lgt_active( lgt_n + 1 ) = k
            lgt_n                   = lgt_n + 1
        end if

    end do

end subroutine create_lgt_active_list
