! ********************************
! 2D AMR prototype
! --------------------------------
!
! refine the block by one level
!
! name: prediction_2D.f90
! date: 17.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine prediction_2D(coarse, fine, Bs)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)                                    :: Bs
    real(kind=rk), dimension(Bs,Bs), intent(out)                    :: fine
    real(kind=rk), dimension((Bs-1)/2+1,(Bs-1)/2+1), intent(in)     :: coarse

    integer(kind=ik)                                                :: i, j

    fine = 0.0_rk

    fine(1:Bs:2, 1:Bs:2) = coarse(:,:)

    ! second order interpolation
    ! y direction
    do i = 2, Bs, 2
        do j = 1, Bs, 2
            fine(i,j) = ( fine(i-1, j) + fine(i+1, j) ) / 2.0_rk
        end do
    end do

    ! x direction
    do i = 1, Bs
        do j = 2, Bs, 2
            fine(i,j) = ( fine(i, j-1) + fine(i, j+1) ) / 2.0_rk
        end do
    end do

end subroutine prediction_2D
