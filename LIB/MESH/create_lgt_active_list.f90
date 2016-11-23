! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: create_lgt_active_list.f90
! version: 0.4
! author: msr
!
! create a list of all active blocks (light data)
! dim 1: light data id
!
! input:    - light data
! output:   - list of active blocks, light data id
!           - number of active blocks
!
! = log ======================================================================================
!
! 23/11/16 - create subroutine
! ********************************************************************************************

subroutine create_lgt_active_list( block_list, lgt_active, n_lgt )

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! light data array
    integer(kind=ik), intent(in)        :: block_list(:, :)

    ! list of active blocks (light data)
    integer(kind=ik), intent(out)       :: lgt_active(:)

    ! number of active blocks (light data)
    integer(kind=ik), intent(out)       :: n_lgt

    ! loop variables
    integer(kind=ik)                    :: k

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset active list
    lgt_active = -1

    ! list index
    n_lgt = 0

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all light data
    do k = 1, size(block_list, 1)

        ! block is active
        if ( block_list(k, 1) /= -1 ) then
            lgt_active( n_lgt + 1 ) = k
            n_lgt = n_lgt + 1
        end if

    end do

end subroutine create_lgt_active_list
