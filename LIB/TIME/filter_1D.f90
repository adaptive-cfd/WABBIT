!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name filter_1D.f90
!> \version 0.5
!> \author msr
!
!> \brief 1D Filter subroutine
!
!> \details
!! input:    - filter stencil, data array, position of value to filter \n
!! output:   - filtered data
!! \n
!! = log ======================================================================================
!! \n
!! 27/03/17 - create
!
! ********************************************************************************************

subroutine filter_1D(phi, a, pos)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> datafield
    real(kind=rk), intent(inout)        :: phi(:)
    !> filter coefficients
    real(kind=rk), intent(in)           :: a(:)
    !> filter position
    integer(kind=ik), intent(in)        :: pos

    ! loop variable
    integer(kind=ik)                    :: k

    ! old values
    real(kind=rk)                       :: phi_old(size(phi,1))

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    phi_old = phi

!---------------------------------------------------------------------------------------------
! main body

    ! check filter stencil
    if ( size(phi) /= size(a) ) then
        write(*,'(80("_"))')
        print*, phi
        print*, a
        write(*,*) "ERROR: filter stencil has wrong size"
        stop
    end if

    ! filter data
    do k = 1, size(a)
        phi(pos) = phi(pos) + a(k)*phi_old(k)
    end do

end subroutine filter_1D
