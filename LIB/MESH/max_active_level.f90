! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: max_active_level.f90
! version: 0.5
! author: engels
!
! return coarsest level in active block list
!
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
! ********************************************************************************************

function max_active_level( lgt_block, lgt_active, lgt_n )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    ! list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    ! return value
    integer(kind=ik)                    :: max_active_level
    ! loop variables
    integer(kind=ik)                    :: k, Jmax, max_treelevel

    max_treelevel = size( lgt_block, 2) -2
!---------------------------------------------------------------------------------------------
! main body
    Jmax = 0

    ! loop over all active blocks
    do k = 1, lgt_n
      Jmax = max(Jmax, lgt_block( lgt_active(k), max_treelevel+1))
    end do

    max_active_level = Jmax
end function
