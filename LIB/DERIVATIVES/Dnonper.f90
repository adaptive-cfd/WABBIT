! ********************************
! 2D AMR prototype
! --------------------------------
! 
! generalized function to create 
! non-periodic differentiation 
! matrices
!
! name: Dnonper.f90
! date: 04.08.2016
! author: msr
! version: 0.1
! 
! ********************************

subroutine Dnonper(D, n, h, stencil, m)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)			    :: n, m
    real(kind=rk), intent(in)			        :: h
    real(kind=rk), dimension(n,n), intent(out)	:: D(n,n)
    real(kind=rk), dimension(m), intent(in)	    :: stencil(m)

    integer(kind=ik)				            :: nrows, i
    real(kind=rk), dimension(m)			        :: ds

    call Dper(D, n, h, stencil, m)
 
    nrows = (m-1)/2

    ds = 0.0_rk

    do i = 1, nrows
        D(i,:) 			    = 0.0_rk
        call giveCertainOrder(ds, (0-(i-1)) , (m-(i-1)-1) , 1, h)
        D(i,1:m) 			= ds / h
        D(n-(i-1),:) 		= 0.0_rk
        call fliplr(ds, m)
        D(n-(i-1),n-m+1:n) 	= -ds / h
    end do

end subroutine Dnonper
