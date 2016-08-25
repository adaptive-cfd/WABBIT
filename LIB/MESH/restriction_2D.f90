! ********************************
! 2D AMR prototype
! --------------------------------
!
! coarsens the block by one level
!
! name: restriction_2D.f90
! date: 17.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine restriction_2D(fine, coarse, Bs)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)                                    :: Bs
    real(kind=rk), dimension(Bs,Bs), intent(in)                     :: fine
    real(kind=rk), dimension((Bs+1)/2, (Bs+1)/2), intent(out)       :: coarse

    coarse = 0.0_rk

    coarse(:, :) = fine(1:Bs:2, 1:Bs:2)

end subroutine restriction_2D
