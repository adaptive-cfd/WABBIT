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
! COMPILED IN THE MODULE "MODULE_INTERPOLATION.F90"
! ********************************

subroutine restriction_2D(fine, coarse)

    use module_params
    use module_blocks

    implicit none

    real(kind=rk), dimension(1:,1:), intent(in) :: fine
    real(kind=rk), dimension(1:,1:), intent(out) :: coarse
    integer(kind=ik) :: nfine, ncoarse

    ncoarse = size(coarse,1)
    nfine = size(fine,1)

    if ( 2*ncoarse-1 /= nfine ) then
      write(*,*) "restriction_2D: arrays wrongly sized.."
      stop
    endif

    coarse = 0.0_rk
    coarse(:, :) = fine(1:nfine:2,1:nfine:2)

end subroutine restriction_2D
