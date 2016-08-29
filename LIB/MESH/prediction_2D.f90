! ********************************
! 2D AMR prototype
! --------------------------------
!
! refine the block by one level
!
! name: prediction_2D.f90
! date: 17.08.2016
! author: msr, engels
! version: 0.1
!
! ********************************

subroutine prediction_2D(coarse, fine, Bs)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)                                    :: Bs
    real(kind=rk), dimension(1:Bs,1:Bs), intent(out)                    :: fine
    real(kind=rk), dimension(1:(Bs-1)/2+1,1:(Bs-1)/2+1), intent(in)     :: coarse

    integer(kind=ik) :: i, j
    integer(kind=ik) :: order, nc
    integer(kind=ik) :: icoarse, ifine

    ! temporarily choose the order here:
    order = 2
    ! we need the number of points in the coarse array (as fortran lacks Matlab's handy
    ! 'end' keyword)
    nc = (Bs-1)/2+1
    ! inititalize fine field with zeros (actually not necessary)
    fine = 9.0e9_rk
    ! fill matching points: the coarse and fine grid share a lot of points (as the
    ! fine grid results from insertion of one point between each coarse point)
    fine(1:Bs:2, 1:Bs:2) = coarse(:,:)


    if ( order == 2 ) then
        !-----------------------------------------------------------------------
        ! second order interpolation
        !-----------------------------------------------------------------------
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
    elseif ( order == 4 ) then
        !-----------------------------------------------------------------------
        ! fourth order interpolation
        !-----------------------------------------------------------------------
         ! ist easier to use a 1D prediction operator and apply it just to rows
         ! and columns...

         ! along x
         do icoarse = 1, nc ! we travel along the coarse grid
             ifine = 2*(icoarse-1)+1
             call prediction1D( coarse(:,icoarse), fine(:,ifine), bs )
         end do
         ! along y
         do icoarse = 1, nc ! we travel along the coarse grid
            ifine = 2*(icoarse-1)+1
            call prediction1D( coarse(icoarse,:), fine(ifine,:), bs )
         end do
         ! between (this is also in x-direction, but we could also do it in y-dir)
         do i = 2, bs, 2
            call prediction1D( fine(1:Bs:2,i ), fine( :, i) , bs )
         end do
    endif

end subroutine prediction_2D


subroutine prediction1D(coarse, fine, Bs)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)                       :: Bs
    real(kind=rk), dimension(1:Bs), intent(out)        :: fine
    real(kind=rk), dimension(1:(Bs-1)/2+1), intent(in) :: coarse

    integer(kind=ik) :: i, j, N, k
    integer(kind=ik) :: order, nc
    real(kind=rk) :: a,b


      ! this is the multiresolution predition operator.
      ! it pushes a signal from a coarser level to the next higher by
      ! interpolation

      N = size(coarse)
      ! note the size of the finer grid is not twice the coarse one, but one less!

      !fine=zeros(1,2*N-1); %

      fine(1:Bs:2) = coarse

      ! fourth order:
      a = 9.0_rk/16.0_rk
      b =-1.0_rk/16.0_rk

      fine(2)     = (5.0_rk/16.0_rk)*coarse(1)+(15.0_rk/16.0_rk)*coarse(2)-(5.0_rk/16.0_rk)*coarse(3)+(1.0_rk/16.0_rk)*coarse(4)
      fine(Bs-1) = (1.0_rk/16.0_rk)*coarse(N-3) -(5.0_rk/16.0_rk)*coarse(N-2) +(15.0_rk/16)*coarse(N-1) +(5.0_rk/16.0_rk)*coarse(N)

      do k = 2, N-2
          fine(2*k) = a*coarse(k)+a*coarse(k+1)+b*coarse(k-1)+b*coarse(k+2)
      end do

end subroutine
