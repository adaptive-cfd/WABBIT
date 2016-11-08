! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: respect_min_max_treelevel.f90
! version: 0.4
! author: msr
!
! unset refinement status in respect of min/max treelevel
!
! input:    - light data
!           - min/max treelevel
! output:   - light data arrays
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
! ********************************************************************************************

subroutine respect_min_max_treelevel( block_list, max_level, min_level)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! light data array
    integer(kind=ik), intent(inout)     :: block_list(:, :)
    ! heavy data array - block data
    integer(kind=ik), intent(in)        :: max_level, min_level

    ! loop variables
    integer(kind=ik)                    :: k, N

!---------------------------------------------------------------------------------------------
! variables initialization

    N = size(block_list, 1)

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all blocks
    do k = 1, N

        if ( block_list(k, 1) /= -1 ) then

            if ( ( block_list(k, max_level+2 ) ==  1) .and. ( block_list(k, max_level+1 ) >= max_level ) ) then
                ! can not refine
                block_list(k, max_level+2 ) = 0
            end if

            if ( ( block_list(k, max_level+2 ) == -1) .and. ( block_list(k, max_level+1 ) <= min_level ) ) then
                ! can not coarsen
                block_list(k, max_level+2 ) = 0
            end if

        end if

    end do

end subroutine respect_min_max_treelevel
