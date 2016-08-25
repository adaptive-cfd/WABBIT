! ********************************
! 2D AMR prototype
! --------------------------------
! 
! giveCertainOrder script
!
! name: giveCertainOrder.f90
! date: 04.08.2016
! author: msr
! version: 0.1
! 
! ********************************

subroutine giveCertainOrder(ds, nMinus, nPlus, ableitung, h)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)						                :: nMinus, nPlus, ableitung
    real(kind=rk), intent(in)						                    :: h
    real(kind=rk), dimension(nPlus - nMinus + 1), intent(out)		    :: ds(nPlus - nMinus + 1)

    integer(kind=ik)							                        :: degrees, i, j, k, n, fac_n, inverse_error
    integer(kind=ik), dimension(nPlus - nMinus + 1)			            :: ipiv
    real(kind=rk), dimension(nPlus - nMinus + 1)				        :: js, Gs_row
    real(kind=rk), dimension(nPlus - nMinus + 1, nPlus - nMinus + 1)	:: Gs, work

    degrees = nPlus - nMinus + 1

    j = 1
    do i = nMinus, nPlus
        js(j) = i
        j = j + 1
    end do

    Gs(:,:) = 0.0_rk

    n = 0
    fac_n = 0

    do i = 1, degrees
        n = i-1
        fac_n = 0
        call factorial(n, fac_n)
        Gs_row = 0.0_rk
        do k = 1, degrees
            Gs(i,k) = h**(n) * js(k)**(n) / real(fac_n,8)
        end do
    end do

    inverse_error	= 0
    ipiv 		= 0
    ! LU factorization
    call dgetrf(degrees, degrees, Gs, degrees, ipiv, inverse_error)
    if (inverse_error /= 0) then
        stop 'Matrix factorization failed!'
    end if

    ! matrix inversion
    call dgetri(degrees, Gs, degrees, ipiv, work, degrees, inverse_error)
    if (inverse_error /= 0) then
        stop 'Matrix inversion failed!'
    end if
 
    ds = Gs(:,ableitung+1)

end subroutine giveCertainOrder
