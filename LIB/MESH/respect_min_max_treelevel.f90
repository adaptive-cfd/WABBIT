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

subroutine respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)

    ! list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n

    ! treelevel restrictions
    integer(kind=ik)                    :: max_level, min_level

    ! loop variables
    integer(kind=ik)                    :: k

!---------------------------------------------------------------------------------------------
! variables initialization

    max_level = params%max_treelevel
    min_level = params%min_treelevel

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all active blocks
    do k = 1, lgt_n


        if ( ( lgt_block( lgt_active(k), max_level+2 ) ==  1) .and. ( lgt_block( lgt_active(k), max_level+1 ) >= max_level ) ) then
            ! can not refine
            lgt_block( lgt_active(k), max_level+2 ) = 0
        end if

        if ( ( lgt_block( lgt_active(k), max_level+2 ) == -1) .and. ( lgt_block( lgt_active(k), max_level+1 ) <= min_level ) ) then
            ! can not coarsen
            lgt_block( lgt_active(k), max_level+2 ) = 0
        end if

    end do

end subroutine respect_min_max_treelevel
