!> \brief calculate treecode size, count elements which are not -1
!! input:    - treecode
!!           - length of treecode input vector
!! output:   - real treecode size (level of treecode)
! ********************************************************************************************

integer function treecode_size(treecode, N)

    use module_params     ! global parameters

    implicit none

    integer(kind=ik), intent(in)    :: N              !> length of treecode vector
    integer(kind=ik), intent(in)    :: treecode(N)    !> treecode vector
    integer(kind=ik)                :: i              ! loop variables

    treecode_size = 0

    do i = 1, N
        if ( treecode(i) /= -1 ) treecode_size = treecode_size + 1
    end do

end function treecode_size
