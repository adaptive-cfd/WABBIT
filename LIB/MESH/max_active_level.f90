!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name max_active_level.f90
!> \version 0.5
!> \author engels
!
!> \brief return finest level in active block list
!
!> \details
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4
! ********************************************************************************************

function max_active_level( lgt_block, lgt_active, lgt_n )

    !---------------------------------------------------------------------------------------------
    ! modules

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(in), optional :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in), optional :: lgt_n
    ! return value
    integer(kind=ik)                    :: max_active_level
    ! loop variables
    integer(kind=ik)                    :: k, Jmax, max_treelevel

    max_treelevel = size( lgt_block, 2) - extra_lgt_fields
    !---------------------------------------------------------------------------------------------
    ! main body
    Jmax = 0

    if (present(lgt_active) .and. present(lgt_n)) then
        ! call with active lists (to be preferred, much faster)
        ! loop over all active blocks
        do k = 1, lgt_n
            Jmax = max(Jmax, lgt_block( lgt_active(k), max_treelevel+IDX_MESH_LVL))
        end do

    else
        ! call without active lists
        do k = 1, size(lgt_block, 1)
            Jmax = max(Jmax, lgt_block( k, max_treelevel+IDX_MESH_LVL))
        end do
    end if

    max_active_level = Jmax
end function
