! ********************************
! 2D AMR prototype
! --------------------------------
! 
! eigth order first derivative 
! nonperiodic
!
! name: D18j.f90
! date: 04.08.2016
! author: msr
! version: 0.1
! 
! ********************************
subroutine D18j(D, n, h)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)			    :: n
    real(kind=rk), intent(in)			        :: h
    real(kind=rk), dimension(n,n), intent(out)	:: D(n,n)
    real(kind=rk), dimension(9)			        :: stencil(9)

    stencil(1) = 1.0_rk/280.0_rk
    stencil(2) = -4.0_rk/105.0_rk
    stencil(3) = 1.0_rk/5.0_rk
    stencil(4) = -4.0_rk/5.0_rk
    stencil(5) = 0.0_rk
    stencil(6) = 4.0_rk/5.0_rk
    stencil(7) = -1.0_rk/5.0_rk
    stencil(8) = 4.0_rk/105.0_rk
    stencil(9) = -1.0_rk/280.0_rk

    call Dnonper(D, n, h, stencil, 9)

end subroutine D18j
