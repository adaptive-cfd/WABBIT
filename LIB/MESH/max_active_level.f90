!> \brief return finest level in active block list
! ********************************************************************************************

function max_active_level( lgt_block, lgt_active, lgt_n )

    implicit none

    integer(kind=ik), intent(in)        :: lgt_block(:, :)            !> light data array
    integer(kind=ik), intent(in), optional :: lgt_active(:)           !> list of active blocks (light data)
    integer(kind=ik), intent(in), optional :: lgt_n                   !> number of active blocks (light data)
    integer(kind=ik)                    :: max_active_level           ! return value
    integer(kind=ik)                    :: k, Jmax, max_treelevel     ! loop variables

    max_treelevel = size( lgt_block, 2) - EXTRA_LGT_FIELDS
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
