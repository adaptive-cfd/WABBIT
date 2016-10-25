! ********************************
! WABBIT
! --------------------------------
!
! refinement and coarsening subroutines
!
! name: module_interpolation.f90
! date: 30.09.2016
! author: msr
! version: 0.2
!
! ********************************

module module_interpolation
  use module_params
  use module_blocks

  implicit none

contains

    ! coarsen the block by one level
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

    ! refine the block by one level
    subroutine prediction_2D(coarse, fine)

        use module_params
        use module_blocks

        implicit none

        real(kind=rk), dimension(1:,1:), intent(out) :: fine
        real(kind=rk), dimension(1:,1:), intent(in)  :: coarse

        integer(kind=ik) :: i, j
        integer(kind=ik) :: ncoarse, nfine
        integer(kind=ik) :: icoarse, ifine

        ncoarse = size(coarse, 1)
        nfine   = size(fine, 1)

        if ( 2*ncoarse-1 /= nfine ) then
          write(*,*) "prediction_2D: arrays wrongly sized.."
          stop
        endif

        ! inititalize fine field (actually not necessary)
        fine = 9.0e9_rk
        ! fill matching points: the coarse and fine grid share a lot of points (as the
        ! fine grid results from insertion of one point between each coarse point)
        fine(1:nfine:2, 1:nfine:2) = coarse(:,:)


        if ( params%order_predictor == "multiresolution_2nd" ) then
            !-----------------------------------------------------------------------
            ! second order interpolation
            !-----------------------------------------------------------------------
            ! y direction
            do i = 2, nfine, 2
              do j = 1, nfine, 2
                fine(i,j) = ( fine(i-1, j) + fine(i+1, j) ) / 2.0_rk
              end do
            end do

            ! x direction
            do i = 1, nfine
              do j = 2, nfine, 2
                fine(i,j) = ( fine(i, j-1) + fine(i, j+1) ) / 2.0_rk
              end do
            end do
        elseif ( params%order_predictor == "multiresolution_4th"  ) then
            !-----------------------------------------------------------------------
            ! fourth order interpolation
            !-----------------------------------------------------------------------
             ! ist easier to use a 1D prediction operator and apply it just to rows
             ! and columns...

             ! along x
             do icoarse = 1, ncoarse ! we travel along the coarse grid
                 ifine = 2*(icoarse-1)+1
                 call prediction1D( coarse(:,icoarse), fine(:,ifine) )
             end do
             ! along y
             do icoarse = 1, ncoarse ! we travel along the coarse grid
                ifine = 2*(icoarse-1)+1
                call prediction1D( coarse(icoarse,:), fine(ifine,:) )
             end do
             ! between (this is also in x-direction, but we could also do it in y-dir)
             do i = 2, nfine, 2
                call prediction1D( fine(1:nfine:2,i ), fine( :, i) )
             end do
        else
             ! error case
             write(*,*) "prediction_2D: wrong method.."
             stop
        endif

    end subroutine prediction_2D


    subroutine prediction1D(coarse, fine)

        use module_params
        use module_blocks

        implicit none

        real(kind=rk), dimension(1:), intent(out) :: fine
        real(kind=rk), dimension(1:), intent(in) :: coarse

        integer(kind=ik) :: k, nfine, ncoarse
        real(kind=rk) :: a, b

        ncoarse = size(coarse,1)
        nfine = size(fine,1)

        if ( 2*ncoarse-1 /= nfine ) then
          write(*,*) "prediction1d: arrays wrongly sized.."
          stop
        endif

        ! this is the multiresolution predition operator.
        ! it pushes a signal from a coarser level to the next higher by
        ! interpolation

        fine(1:nfine:2) = coarse(:)

        ! fourth order:
        a = 9.0_rk/16.0_rk
        b =-1.0_rk/16.0_rk

        fine(2)     = (5.0_rk/16.0_rk)*coarse(1)+(15.0_rk/16.0_rk)*coarse(2)-(5.0_rk/16.0_rk)*coarse(3)+(1.0_rk/16.0_rk)*coarse(4)
        fine(nfine-1) = (1.0_rk/16.0_rk)*coarse(ncoarse-3) -(5.0_rk/16.0_rk)*coarse(ncoarse-2) +(15.0_rk/16.0_rk)*coarse(ncoarse-1) +(5.0_rk/16.0_rk)*coarse(ncoarse)

        do k = 2, ncoarse-2
            fine(2*k) = a*coarse(k)+a*coarse(k+1)+b*coarse(k-1)+b*coarse(k+2)
        end do

    end subroutine prediction1D

end module
