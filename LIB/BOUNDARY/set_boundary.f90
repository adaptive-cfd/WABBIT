! ********************************
! 2D AMR prototype
! --------------------------------
!
! set boundary conditions
!
! name: set_boundary.f90
! date: 05.09.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine set_boundary()

    use module_params
    use module_blocks

    implicit none

    ! check boundary condition
    if (params%boundary == "periodic") then
         ! periodic condition
         call periodic_2D()
    else
         ! to do
    end if

end subroutine set_boundary
