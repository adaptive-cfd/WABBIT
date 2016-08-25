! ********************************
! 2D AMR prototype
! --------------------------------
! 
! create derivative matrices with
! given stencil
!
! name: Dper.f90
! date: 04.08.2016
! author: msr
! version: 0.1
! 
! ********************************

subroutine Dper(D, n, h, stencil, m)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)			    :: n, m
    real(kind=rk), intent(in)			        :: h
    real(kind=rk), dimension(n,n), intent(out)	:: D(n,n)
    real(kind=rk), dimension(m), intent(in)	    :: stencil(m)

    integer(kind=ik)				            :: rows, i, pos1
    integer(kind=ik), dimension(m)		        :: pos
    real(kind=rk)				                :: m_2
   
    if (mod(m,2) /= 1 ) print*, "stencil have an even number of entries"

    m_2 = m / 2
    do i = 1, m
        pos(i) = -floor(m_2) + i - 1
    end do
  
    D(:,:) = 0.0_rk

    do rows = 1, n
        do i = 1, m
            pos1 = pos(i) + rows
            if (pos1 < 1) pos1 = pos1 + n
            if (pos1 > n) pos1 = pos1 - n
            D(rows,pos1) = stencil(i)
        end do
    end do

    D = D / h

end subroutine Dper
