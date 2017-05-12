!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name treecode_size.f90
!> \version 0.4
!> \author msr
!
!> \brief calculate treecode size, count elements which are not -1
!
!> 
!! input:    
!!           - treecode
!!           - length of treecode input vector
!!
!! output:   
!!           - real treecode size (level of treecode)
!!
!! = log ======================================================================================
!! \n
!! 07/11/16 - switch to v0.4
! ********************************************************************************************

integer function treecode_size(treecode, N)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> length of treecode vector
    integer(kind=ik), intent(in)    :: N

    !> treecode vector
    integer(kind=ik), intent(in)    :: treecode(N)

    ! loop variables
    integer(kind=ik)                :: i

!---------------------------------------------------------------------------------------------
! variables initialization

    treecode_size = 0

!---------------------------------------------------------------------------------------------
! main body

    do i = 1, N
        if ( treecode(i) /= -1 ) treecode_size = treecode_size + 1
    end do

end function treecode_size
