!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name int_to_binary.f90
!> \version 0.4
!> \author msr
!
!> \brief convert a integer i to binary b \n
!! binary return as vector with length N
!
!> \details
!! input:    
!!           - integer to convert
!!           - length of output vector
!!
!! output:   
!!           - "binary" vector
!!
!! = log ======================================================================================
!! \n
!! 07/11/16 - switch to v0.4
! ********************************************************************************************

subroutine int_to_binary(i, N, b)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> integer to convert into binary
    integer(kind=ik), intent(in)    :: i

    !> length of binary output vector
    integer(kind=ik), intent(in)    :: N

    !> output vector
    integer(kind=ik), intent(out)   :: b(N)

    ! loop variables
    integer(kind=ik)                :: j, k

!---------------------------------------------------------------------------------------------
! variables initialization

    j = 1
    b = 0
    k = i

!---------------------------------------------------------------------------------------------
! main body

    do while (k > 0)
        b(j) = mod(k, 2)
        k = int(k/2)
        j = j + 1
    end do

end subroutine int_to_binary
