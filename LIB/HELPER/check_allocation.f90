! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: check_allocation.f90
! version: 0.5
! author: msr
!
! convert a integer i to binary b
! binary return as vector with length N
!
! input:    - allocation error
! output:   - stop program if error variable /= 0
!
! = log ======================================================================================
!
! 25/01/17 - create
! ********************************************************************************************

subroutine check_allocation(allocate_error)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! allocation error variable
    integer(kind=ik), intent(in) :: allocate_error

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    else
        ! nothing to do
    end if

end subroutine check_allocation
