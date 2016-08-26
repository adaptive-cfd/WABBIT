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
    real(kind=rk), dimension(Bs,Bs), intent(out)                    :: fine
    real(kind=rk), dimension((Bs-1)/2+1,(Bs-1)/2+1), intent(in)     :: coarse

    integer(kind=ik)                                                :: i, j
    integer(kind=ik) :: order, nc

    ! temporarily choose the order here:
    order = 4
    ! we need the number of points in the coarse array (as fortran lacks Matlab's handy
    ! 'end' keyword)
    nc = (Bs-1)/2+1
    ! inititalize fine field with zeros (actually not necessary)
    fine = 0.0_rk
    ! fill matching points: the coarse and fine grid share a lot of points (as the
    ! fine grid results from insertion of one point between each coarse point)
    fine(1:Bs:2, 1:Bs:2) = coarse(:,:)


    if ( order == 2) then
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
        fine(1:Bs:2,2)    = (5.0_rk/16.0_rk)*coarse(:,1)   +(15.0_rk/16.0_rk)*coarse(:,2)    -(5.0_rk/16.0_rk)* coarse(:,3)    +(1.0_rk/16.0_rk)*coarse(:,4)
        fine(1:Bs:2,Bs-1) = (1.0_rk/16.0_rk)*coarse(:,nc-3)-(5.0_rk/16.0_rk) *coarse(:,nc-2) +(15.0_rk/16.0_rk)*coarse(:,nc-1) +(5.0_rk/16.0_rk)*coarse(:,nc)

        fine(2   ,1:bs:2) = (5.0_rk/16.0_rk)*coarse(1,:)   +(15.0_rk/16.0_rk)*coarse(2,:)    -(5.0_rk/16.0_rk)* coarse(3,:)    +(1.0_rk/16.0_rk)*coarse(4,:)
        fine(bs-1,1:bs:2) = (1.0_rk/16.0_rk)*coarse(nc-3,:)-(5.0_rk/16.0_rk) *coarse(nc-2,:) +(15.0_rk/16.0_rk)*coarse(nc-1,:) +(5.0_rk/16.0_rk)*coarse(nc,:)

        ! along y-direction
        fine(4:bs-3:2,:) = (-1.0_rk/16.0_rk)*fine(1:bs-6:2,:) + (9.0_rk/16.0_rk)*fine(3:bs-4:2,:) + (9.0_rk/16.0_rk)*fine(5:bs-2:2,:) + (-1.0_rk/16.0_rk)*fine(7:bs:2,:);
        ! along x-direction
        fine(:,4:bs-3:2) = (-1.0_rk/16.0_rk)*fine(:,1:bs-6:2) + (9.0_rk/16.0_rk)*fine(:,3:bs-4:2) + (9.0_rk/16.0_rk)*fine(:,5:bs-2:2) + (-1.0_rk/16.0_rk)*fine(:,7:bs:2);
    endif

end subroutine prediction_2D
