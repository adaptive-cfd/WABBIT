! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: does_block_exist.f90
! version: 0.4
! author: msr, engels
!
! return true, if block exist
!
! input:    - treecode
!           - light data array
!           - max treelvel
! output:   - .true. if block with treecode exists
!
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
! ********************************************************************************************

subroutine does_block_exist(treecode, lgt_block, max_treelevel, exists, light_id, lgt_active, lgt_n)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! max treelevel
    integer(kind=ik), intent(in)        :: max_treelevel
    ! block treecode
    integer(kind=ik), intent(in)        :: treecode(max_treelevel)

    ! light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)

    ! list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n

    ! .true. if block with treecode exists
    logical, intent(out)                :: exists
    ! light data id
    integer(kind=ik), intent(out)       :: light_id

    ! function to compare treecodes
    logical                             :: array_compare

    ! loop variables
    integer(kind=ik)                    :: k, N

!---------------------------------------------------------------------------------------------
! variables initialization

    N        = size( lgt_block, 1)
    exists   = .false.
    light_id = -1

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all active blocks
    do k = 1, lgt_n

        if ( array_compare( lgt_block( lgt_active(k) , 1:max_treelevel), treecode, max_treelevel) ) then
            exists   = .true.
            light_id = lgt_active(k)
        end if

    end do

end subroutine does_block_exist
