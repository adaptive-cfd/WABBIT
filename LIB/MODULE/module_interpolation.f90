!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_interpolation.f90
!> \version 0.5
!> \author msr, engels
!
!> \brief refinement and coarsening subroutines
!
!>
!! = log ======================================================================================
!! \n
!! 04/11/16 - switch to v0.4 \n
!! 03/02/17 - create subroutines for 3D
!
! ********************************************************************************************

module module_interpolation

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    PRIVATE
    PUBLIC  :: restriction_2D,restriction_3D,prediction_2D,prediction_3D
!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    ! coarsen the block by one level
    subroutine restriction_2D(fine, coarse)

        implicit none

        real(kind=rk), dimension(1:,1:), intent(in) :: fine
        real(kind=rk), dimension(1:,1:), intent(out) :: coarse
        integer(kind=ik) :: nfine, ncoarse

        ncoarse = size(coarse,1)
        nfine = size(fine,1)

        if ( 2*ncoarse-1 /= nfine ) then
          call abort(888191,"ERROR: restriction_2D: arrays wrongly sized..")
        endif

        ! reseting not necessary
        !coarse = 0.0_rk
        coarse(:, :) = fine(1:nfine:2,1:nfine:2)

    end subroutine restriction_2D

    ! coarsen the block by one level
    subroutine restriction_3D(fine, coarse)

        implicit none

        real(kind=rk), dimension(1:,1:,1:), intent(in)  :: fine
        real(kind=rk), dimension(1:,1:,1:), intent(out) :: coarse
        integer(kind=ik) :: nfine, ncoarse

        ncoarse = size(coarse,1)
        nfine = size(fine,1)

        if ( 2*ncoarse-1 /= nfine ) then
          call abort(888192,"ERROR: restriction_3D: arrays wrongly sized..")
        endif

        coarse(:, :, :) = fine(1:nfine:2,1:nfine:2,1:nfine:2)

    end subroutine restriction_3D

    ! refine the block by one level
    subroutine prediction_2D(coarse, fine, order_predictor)

        implicit none

        real(kind=rk), dimension(1:,1:), intent(out) :: fine
        real(kind=rk), dimension(1:,1:), intent(in)  :: coarse

        character(len=80), intent(in)                :: order_predictor

        integer(kind=ik) :: i, j, l
        integer(kind=ik) :: ncoarse, nfine
        integer(kind=ik) :: icoarse, ifine
        ! interpolation coefficients
        ! a: one sided, b: central
        real(kind=rk) :: a(4), b(2)

        ncoarse = size(coarse, 1)
        nfine   = size(fine, 1)

        if ( 2*ncoarse-1 /= nfine ) then
          call abort(888193,"ERROR: prediction_2D: arrays wrongly sized..")
        endif

        ! fill matching points: the coarse and fine grid share a lot of points (as the
        ! fine grid results from insertion of one point between each coarse point)
        fine(1:nfine:2, 1:nfine:2) = coarse(:,:)


        if ( order_predictor == "multiresolution_2nd" ) then
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

        elseif ( order_predictor == "multiresolution_4th"  ) then
            !-----------------------------------------------------------------------
            ! fourth order interpolation
            !-----------------------------------------------------------------------
            ! init coefficients
            a(1) =  5.0_rk/16.0_rk
            a(2) =  15.0_rk/16.0_rk
            a(3) =  -5.0_rk/16.0_rk
            a(4) =  1.0_rk/16.0_rk

            b(1) =  9.0_rk/16.0_rk
            b(2) =  -1.0_rk/16.0_rk
            ! ist easier to use a 1D prediction operator and apply it just to rows
            ! and columns...
            ! along x
            do icoarse = 1, ncoarse
                ! we travel along the coarse grid
                ifine  = 2*(icoarse-1)+1
                fine( 2, ifine ) = a(1)*coarse(1,icoarse) + a(2)*coarse(2,icoarse) &
                                 + a(3)*coarse(3,icoarse) + a(4)*coarse(4,icoarse)

                fine( nfine-1, ifine ) = a(4)*coarse(ncoarse-3,icoarse) + a(3)*coarse(ncoarse-2,icoarse) &
                                       + a(2)*coarse(ncoarse-1,icoarse) + a(1)*coarse(ncoarse  ,icoarse)

                do l = 2, ncoarse-2
                    fine( 2*l, ifine ) = b(1)*coarse(l  ,icoarse) + b(1)*coarse(l+1,icoarse) &
                                       + b(2)*coarse(l-1,icoarse) + b(2)*coarse(l+2,icoarse)
                end do
            end do

            ! along y
            do icoarse = 1, ncoarse
                ! we travel along the coarse grid
                ifine  = 2*(icoarse-1)+1
                fine( ifine, 2  ) = a(1)*coarse(icoarse,1) + a(2)*coarse(icoarse,2) &
                                  + a(3)*coarse(icoarse,3) + a(4)*coarse(icoarse,4)

                fine( ifine, nfine-1  ) = a(4)*coarse(icoarse,ncoarse-3) + a(3)*coarse(icoarse,ncoarse-2) &
                                        + a(2)*coarse(icoarse,ncoarse-1) + a(1)*coarse(icoarse,ncoarse  )

                do l = 2, ncoarse-2
                    fine( ifine, 2*l ) = b(1)*coarse(icoarse,l  ) + b(1)*coarse(icoarse,l+1) &
                                       + b(2)*coarse(icoarse,l-1) + b(2)*coarse(icoarse,l+2)
                end do
            end do

            ! diagonal (this is also in x-direction, but we could also do it in y-dir)
            do ifine = 2, nfine, 2
                fine( 2, ifine  ) = a(1)*fine(1, ifine) + a(2)*fine(3, ifine) &
                                  + a(3)*fine(5, ifine) + a(4)*fine(7, ifine)

                fine( nfine-1, ifine ) = a(4)*fine(nfine-6, ifine) + a(3)*fine(nfine-4, ifine) &
                                       + a(2)*fine(nfine-2, ifine) + a(1)*fine(nfine, ifine  )

                do l = 4, nfine-3, 2
                    fine( l, ifine ) = b(1)*fine(l-1, ifine) + b(1)*fine(l+1, ifine) &
                                     + b(2)*fine(l-3, ifine) + b(2)*fine(l+3, ifine)
                end do
            end do
        else
             ! error case
             call abort(888194,"ERROR: prediction_2D: wrong method..")
        endif

    end subroutine prediction_2D

    ! refine the block by one level
    subroutine prediction_3D(coarse, fine, order_predictor)

        implicit none

        real(kind=rk), dimension(1:,1:,1:), intent(out) :: fine
        real(kind=rk), dimension(1:,1:,1:), intent(in)  :: coarse
        character(len=80), intent(in)                :: order_predictor

        integer(kind=ik) :: i, j, k, l
        integer(kind=ik) :: ncoarse, nfine
        integer(kind=ik) :: icoarse, ifine
        ! interpolation coefficients
        ! a: one sided, b: central
        real(kind=rk) :: a(4), b(2)

        ncoarse = size(coarse, 1)
        nfine   = size(fine, 1)

        if ( 2*ncoarse-1 /= nfine ) then
          call abort(888195,"ERROR: prediction_3D: arrays wrongly sized..")
        endif

        ! fill matching points: the coarse and fine grid share a lot of points (as the
        ! fine grid results from insertion of one point between each coarse point)
        fine(1:nfine:2, 1:nfine:2, 1:nfine:2) = coarse(:,:,:)

        if ( order_predictor == "multiresolution_2nd" ) then
            !-------------------------------------------------------------------
            ! second order interpolation
            !-------------------------------------------------------------------
            ! y direction
            do k = 1, nfine, 2
                do j = 1, nfine, 2
                    do i = 2, nfine, 2
                        fine(i,j,k) = ( fine(i-1, j, k) + fine(i+1, j, k) ) / 2.0_rk
                    end do
                end do
            end do

            ! x direction
            do k = 1, nfine, 2
                do j = 2, nfine, 2
                    do i = 1, nfine
                        fine(i,j,k) = ( fine(i, j-1, k) + fine(i, j+1, k) ) / 2.0_rk
                    end do
                end do
            end do

            ! z direction
            do k = 2, nfine, 2
                do j = 1, nfine
                    do i = 1, nfine
                        fine(i,j,k) = ( fine(i, j, k-1) + fine(i, j, k+1) ) / 2.0_rk
                    end do
                end do
            end do

        elseif ( order_predictor == "multiresolution_4th"  ) then
            !-----------------------------------------------------------------------
            ! fourth order interpolation
            !-----------------------------------------------------------------------
            ! init coefficients
            a(1) =  5.0_rk/16.0_rk
            a(2) =  15.0_rk/16.0_rk
            a(3) =  -5.0_rk/16.0_rk
            a(4) =  1.0_rk/16.0_rk

            b(1) =  9.0_rk/16.0_rk
            b(2) =  -1.0_rk/16.0_rk

            ! ist easier to use a 1D prediction operator and apply it just to rows
            ! and columns...
            ! loop over all z-level
            do k = 1, ncoarse
                ! along x
                do icoarse = 1, ncoarse
                    ! we travel along the coarse grid
                    ifine  = 2*(icoarse-1)+1

                    !call prediction1D( coarse(:,icoarse,k), fine(:,ifine,2*(k-1)+1) )

                    fine( 2, ifine, 2*(k-1)+1  ) = a(1)*coarse(1,icoarse,k) + a(2)*coarse(2,icoarse,k) &
                                                 + a(3)*coarse(3,icoarse,k) + a(4)*coarse(4,icoarse,k)

                    fine( nfine-1, ifine, 2*(k-1)+1  ) = a(4)*coarse(ncoarse-3,icoarse,k) + a(3)*coarse(ncoarse-2,icoarse,k) &
                                                       + a(2)*coarse(ncoarse-1,icoarse,k) + a(1)*coarse(ncoarse  ,icoarse,k)

                    do l = 2, ncoarse-2
                        fine( 2*l, ifine, 2*(k-1)+1 ) = b(1)*coarse(l  ,icoarse,k) + b(1)*coarse(l+1,icoarse,k) &
                                                      + b(2)*coarse(l-1,icoarse,k) + b(2)*coarse(l+2,icoarse,k)
                    end do
                end do

                ! along y
                do icoarse = 1, ncoarse
                    ! we travel along the coarse grid
                    ifine  = 2*(icoarse-1)+1

                    !call prediction1D( coarse(icoarse,:,k), fine(ifine,:,2*(k-1)+1) )

                    fine( ifine, 2, 2*(k-1)+1  )        = a(1)*coarse(icoarse,1,k) + a(2)*coarse(icoarse,2,k) &
                                                        + a(3)*coarse(icoarse,3,k) + a(4)*coarse(icoarse,4,k)

                    fine( ifine, nfine-1, 2*(k-1)+1  )  = a(4)*coarse(icoarse,ncoarse-3,k) + a(3)*coarse(icoarse,ncoarse-2,k) &
                                                        + a(2)*coarse(icoarse,ncoarse-1,k) + a(1)*coarse(icoarse,ncoarse  ,k)

                    do l = 2, ncoarse-2
                        fine( ifine, 2*l, 2*(k-1)+1 )   = b(1)*coarse(icoarse,l  ,k) + b(1)*coarse(icoarse,l+1,k) &
                                                        + b(2)*coarse(icoarse,l-1,k) + b(2)*coarse(icoarse,l+2,k)
                    end do
                end do

                ! between (this is also in x-direction, but we could also do it in y-dir)
                do ifine = 2, nfine, 2
                    !call prediction1D( fine(1:nfine:2,ifine,2*(k-1)+1 ), fine( :, ifine, 2*(k-1)+1) )

                    fine( 2, ifine, 2*(k-1)+1  ) = a(1)*fine(1, ifine, 2*(k-1)+1) + a(2)*fine(3, ifine, 2*(k-1)+1) &
                                                 + a(3)*fine(5, ifine, 2*(k-1)+1) + a(4)*fine(7, ifine, 2*(k-1)+1)

                    fine( nfine-1, ifine, 2*(k-1)+1 ) = a(4)*fine(nfine-6, ifine, 2*(k-1)+1) + a(3)*fine(nfine-4, ifine, 2*(k-1)+1) &
                                                      + a(2)*fine(nfine-2, ifine, 2*(k-1)+1) + a(1)*fine(nfine, ifine, 2*(k-1)+1  )

                    do l = 4, nfine-3, 2
                        fine( l, ifine, 2*(k-1)+1 ) = b(1)*fine(l-1, ifine, 2*(k-1)+1) + b(1)*fine(l+1, ifine, 2*(k-1)+1) &
                                                    + b(2)*fine(l-3, ifine, 2*(k-1)+1) + b(2)*fine(l+3, ifine, 2*(k-1)+1)
                    end do
                end do
            end do

            ! interpolate z
            do i = 1, nfine
                do j = 1, nfine
                    !call prediction1D( fine( i, j, 1:nfine:2 ), fine( i, j, : ) )

                    fine( i, j, 2  ) = a(1)*fine(i,j,1) + a(2)*fine(i,j,3) &
                                     + a(3)*fine(i,j,5) + a(4)*fine(i,j,7)

                    fine( i, j, nfine-1 ) = a(4)*fine(i,j,nfine-6) + a(3)*fine(i,j,nfine-4) &
                                          + a(2)*fine(i,j,nfine-2) + a(1)*fine(i,j,nfine  )

                    do l = 4, nfine-3, 2
                        fine( i, j, l ) = b(1)*fine(i,j,l-1) + b(1)*fine(i,j,l+1) &
                                        + b(2)*fine(i,j,l-3) + b(2)*fine(i,j,l+3)
                    end do
                end do
            end do

        else
            ! error case
            call abort(888196,"ERROR: prediction_2D: wrong method..")
        endif

    end subroutine prediction_3D


    subroutine prediction1D(coarse, fine)

        implicit none

        real(kind=rk), dimension(1:), intent(out) :: fine
        real(kind=rk), dimension(1:), intent(in) :: coarse

        integer(kind=ik) :: k, nfine, ncoarse
        real(kind=rk) :: a, b

        ncoarse = size(coarse,1)
        nfine = size(fine,1)

        if ( 2*ncoarse-1 /= nfine ) then
          call abort(888197,"ERROR: prediction1d: arrays wrongly sized..")
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
