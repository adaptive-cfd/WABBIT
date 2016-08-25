! ********************************
! 2D AMR prototype
! --------------------------------
!
! calculate treecode size, count
! elements which are not -1
!
! name: treecode_size.f90
! date: 17.08.2016
! author: msr
! version: 0.1
!
! ********************************

integer function treecode_size(treecode)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), dimension(10)     :: treecode
    integer                             :: i

    treecode_size = 0

    do i = 1, 10
        if ( treecode(i) /= -1 ) treecode_size = treecode_size + 1
    end do

end function treecode_size
