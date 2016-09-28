! ********************************
! 2D AMR prototype
! --------------------------------
! 
! sixth order second derivative
! periodic
!
! name: D26p.f90
! date: 04.08.2016
! author: msr
! version: 0.1
! 
! ********************************
subroutine D26p(D, n, h)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)		            :: n
    real(kind=rk), intent(in)		                :: h
    real(kind=rk), dimension(n,n), intent(out)	    :: D(n,n)
    real(kind=rk), dimension(7)		                :: stencil(7)

    stencil(1) = 1.0_rk/90.0_rk
    stencil(2) = -3.0_rk/20.0_rk
    stencil(3) = 3.0_rk/2.0_rk
    stencil(4) = -49.0_rk/18.0_rk
    stencil(5) = 3.0_rk/2.0_rk
    stencil(6) = -3.0_rk/20.0_rk
    stencil(7) = 1.0_rk/90.0_rk

    call Dper(D, n, h, stencil, 7)

end subroutine D26p
