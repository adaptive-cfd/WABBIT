! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: min_active_level.f90
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

function min_active_level( lgt_block, lgt_active, lgt_n )

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
    integer(kind=ik)                    :: min_active_level
    ! loop variables
    integer(kind=ik)                    :: k, Jmin, max_treelevel

    max_treelevel = size( lgt_block, 2) -2
!---------------------------------------------------------------------------------------------
! main body
    Jmin = max_treelevel

    ! loop over all active blocks
    do k = 1, lgt_n
      Jmin = min(Jmin, lgt_block( lgt_active(k), max_treelevel+1))
    end do

    min_active_level = Jmin
end function
