!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name array_compare.f90
!> \version 0.4
!> \author msr
!
!> \brief compare 2 arrays
!
!> \details
!! input:    
!!           - 2 array
!!           - array size N
!!
!! output:   
!!           - .true. if arrays are equal
!!
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4
! ********************************************************************************************

logical function array_compare(array1, array2, N)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> array size
    integer(kind=ik), intent(in)        :: N
    !> arrays
    integer(kind=ik), intent(in)        :: array1(N), array2(N)

    ! loop variable
    integer(kind=ik)                    :: i

!---------------------------------------------------------------------------------------------
! variables initialization

    array_compare = .true.

!---------------------------------------------------------------------------------------------
! main body

    do i = 1, N
        if ( array1(i) /= array2(i) ) array_compare = .false.
    enddo

end function array_compare
