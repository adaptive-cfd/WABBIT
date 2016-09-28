! ********************************
! 2D AMR prototype
! --------------------------------
!
! second order second derivative
! periodic
!
! name: D22p.f90
! date: 04.08.2016
! author: msr
! version: 0.1
!
! ********************************
subroutine D22p(D, n, h)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)                    :: n
    real(kind=rk), intent(in)                       :: h
    real(kind=rk), dimension(n,n), intent(out)      :: D(n,n)
    real(kind=rk), dimension(3)                     :: stencil(3)

    stencil(1) = 1.0_rk
    stencil(2) = -2.0_rk
    stencil(3) = 1.0_rk

    call Dper(D, n, h, stencil, 3)

end subroutine D22p
