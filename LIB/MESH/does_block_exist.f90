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

subroutine does_block_exist(treecode, block_list, max_treelevel, exists, light_id)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! max treelevel
    integer(kind=ik), intent(in)        :: max_treelevel
    ! block treecode
    integer(kind=ik), intent(in)        :: treecode(max_treelevel)

    ! light data array
    integer(kind=ik), intent(in)        :: block_list(:, :)

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

    N        = size( block_list, 1)
    exists   = .false.
    light_id = -1

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all blocks
    do k = 1, N

        ! block is active
        if ( block_list(k, 1) /= -1 ) then

            if ( array_compare( block_list(k, 1:max_treelevel), treecode, max_treelevel) ) then
                exists   = .true.
                light_id = k
            end if

        end if

    end do

end subroutine does_block_exist
