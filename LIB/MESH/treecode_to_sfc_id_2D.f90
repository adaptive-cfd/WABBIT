!> \brief convert given treecode to position in sfc
!! input:    - treecode \n
!! output:   - position in sfc curve \n
! ********************************************************************************************

subroutine treecode_to_sfc_id_2D(sfc_id, treecode, n)

    implicit none
    integer(kind=ik), intent(out)       :: sfc_id
    integer(kind=ik), intent(in)        :: n              !> treecode size
    integer(kind=ik), intent(in)        :: treecode(n)
    integer(kind=ik)                    :: k              ! loop variable

    sfc_id = 0

    ! calculate sfc id
    ! change -1 in treecode to 0
    ! so every treecode is on highest level
    do k = 1, n
        if ( treecode(n-k+1) == -1 ) then
            ! nothing to do, treecode element is "0"
        else
            ! non zero element
            sfc_id = sfc_id + 4**(k-1) * treecode(n-k+1)
        end if
    end do

end subroutine treecode_to_sfc_id_2D
