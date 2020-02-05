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
    PUBLIC  :: restriction_2D,restriction_3D,prediction_2D,prediction_3D, restriction_prefilter_2D, restriction_prefilter_3D
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
        integer(kind=ik), dimension(2) :: nfine, ncoarse

        ncoarse(1) = size(coarse,1)
        nfine(1) = size(fine,1)
        ncoarse(2) = size(coarse,2)
        nfine(2) = size(fine,2)

        if ( 2*ncoarse(1)-1 /= nfine(1) .or. 2*ncoarse(2)-1 /= nfine(2)) then
            call abort(888191,"ERROR: restriction_2D: arrays wrongly sized..")
        endif

        coarse(:, :) = fine(1:nfine(1):2,1:nfine(2):2)

    end subroutine restriction_2D

    ! coarsen the block by one level
    subroutine restriction_prefilter_2D(fine, fine_filtered2, wavelet)
        implicit none

        real(kind=rk), dimension(1:,1:), intent(in) :: fine
        real(kind=rk), dimension(1:,1:), intent(out) :: fine_filtered2
        character(len=80), intent(in) :: wavelet
        real(kind=rk), dimension(1:size(fine,1),1:size(fine,2)) :: fine_filtered
        integer(kind=ik) :: ix, iy, shift, a, b, nx, ny
        real(kind=rk), allocatable, save :: HD(:)

        nx = size(fine,1)
        ny = size(fine,2)

        ! initialize filter according to wavelet
        if (.not. allocated(HD)) then
            select case(wavelet)
            case("CDF4,4")
                allocate( HD(-6:6) )

                ! H TILDE filter
                HD = (/ -2.0d0**(-9.d0), 0.0d0,  9.0d0*2.0d0**(-8.d0), -2.0d0**(-5.d0),  -63.0d0*2.0d0**(-9.d0),  9.0d0*2.0d0**(-5.d0), &
                87.0d0*2.0d0**(-7.d0), &
                9.0d0*2.0d0**(-5.d0), -63.0d0*2.0d0**(-9.d0), -2.0d0**(-5.d0), 9.0d0*2.0d0**(-8.d0), 0.0d0, -2.0d0**(-9.d0)/) ! TILDE

            case default
                call abort(0309192, "unkown biorothonal wavelet specified. Set course for adventure!")

            end select
        endif

        a = lbound(HD, dim=1)
        b = ubound(HD, dim=1)

        ! note: ghost node layer is returned as original values.

        ! copy (to fill ghost nodes with original data)
        fine_filtered(:,:) = fine
        fine_filtered(-a+1:nx-b,-a+1:ny-b) = 0.0_rk

        ! apply the filter (first in x-direction)
        do ix = -a+1, nx-b
            do shift = a, b
                fine_filtered(ix, -a+1:ny-b) = fine_filtered(ix, -a+1:ny-b) + fine(ix+shift, -a+1:ny-b)*HD(shift)
            enddo
        enddo

        ! copy (to fill ghost nodes with original data)
        fine_filtered2(:,:) = fine
        fine_filtered2(-a+1:nx-b,-a+1:ny-b) = 0.0_rk

        ! then in y-direction
        do iy = -a+1, ny-b
            do shift = a, b
                fine_filtered2(-a+1:nx-b, iy) = fine_filtered2(-a+1:nx-b, iy) + fine_filtered(-a+1:nx-b, iy+shift)*HD(shift)
            enddo
        enddo
    end subroutine


    subroutine restriction_prefilter_3D(fine, fine_filtered2, wavelet)
        implicit none

        real(kind=rk), dimension(1:,1:,1:), intent(in) :: fine
        real(kind=rk), dimension(1:,1:,1:), intent(out) :: fine_filtered2
        character(len=80), intent(in) :: wavelet
        real(kind=rk), allocatable, save :: fine_filtered(:,:,:)
        real(kind=rk), allocatable, save :: HD(:)
        integer(kind=ik) :: ix, iy, iz, shift, a, b, nx, ny, nz

        nx = size(fine,1)
        ny = size(fine,2)
        nz = size(fine,3)

        if (.not.allocated(fine_filtered)) allocate(fine_filtered(1:nx,1:ny,1:nz))

        ! initialize filter according to wavelet
        if (.not. allocated(HD)) then
            select case(wavelet)
            case("CDF4,4")
                allocate( HD(-6:6) )

                ! H TILDE filter ( from Sweldens 1996 paper)
                HD = (/ -2.0d0**(-9.d0), 0.0d0,  9.0d0*2.0d0**(-8.d0), -2.0d0**(-5.d0),  -63.0d0*2.0d0**(-9.d0),  9.0d0*2.0d0**(-5.d0), &
                87.0d0*2.0d0**(-7.d0), &
                9.0d0*2.0d0**(-5.d0), -63.0d0*2.0d0**(-9.d0), -2.0d0**(-5.d0), 9.0d0*2.0d0**(-8.d0), 0.0d0, -2.0d0**(-9.d0)/) ! TILDE

                ! attention. sweldens gives also the coefficients for CDF40, and there he does not have 1/16, but 1/32.
                ! his coefficients are thus divided by two. therefore, as we copy (g and h_tilde) from this paper
                ! and mix it with the 1/16 we had before, we need to multiply by TWO here.
                HD  = HD*2.0d0

            case default
                call abort(0309192, "unkown biorothonal wavelet specified. Set course for adventure!")

            end select
        endif

        a = lbound(HD, dim=1)
        b = ubound(HD, dim=1)

        ! note: ghost node layer is returned as original values.

        ! copy (to fill ghost nodes with original data)
        fine_filtered2(:, :, :) = fine
        fine_filtered2(-a+1:nx-b, -a+1:ny-b, -a+1:nz-b) = 0.0_rk

        ! apply the filter (first in x-direction)
        do ix = -a+1, nx-b
            do shift = a, b
                fine_filtered2(ix, -a+1:ny-b, -a+1:nz-b) = fine_filtered2(ix, -a+1:ny-b, -a+1:nz-b) + fine(ix+shift, -a+1:ny-b, -a+1:nz-b)*HD(shift)
            enddo
        enddo

        ! copy (to fill ghost nodes with original data)
        fine_filtered(:, :, :) = fine
        fine_filtered(-a+1:nx-b, -a+1:ny-b, -a+1:nz-b) = 0.0_rk

        ! then in y-direction
        do iy = -a+1, ny-b
            do shift = a, b
                fine_filtered(-a+1:nx-b, iy, -a+1:nz-b) = fine_filtered(-a+1:nx-b, iy, -a+1:nz-b) + fine_filtered2(-a+1:nx-b, iy+shift, -a+1:nz-b)*HD(shift)
            enddo
        enddo

        ! copy (to fill ghost nodes with original data)
        fine_filtered2(:, :, :) = fine
        fine_filtered2(-a+1:nx-b, -a+1:ny-b, -a+1:nz-b) = 0.0_rk

        ! then in z-direction
        do iz = -a+1, nz-b
            do shift = a, b
                fine_filtered2(-a+1:nx-b, -a+1:ny-b, iz) = fine_filtered2(-a+1:nx-b, -a+1:ny-b, iz) + fine_filtered(-a+1:nx-b, -a+1:ny-b, iz+shift)*HD(shift)
            enddo
        enddo
    end subroutine



    ! periodic index
    ! uses one-based indexing
    integer function perindex(i,n)
        implicit none
        integer, intent(in) :: i, n
        perindex = i
        do while (perindex<1 .or. perindex>n)
            if (perindex<1) then
                perindex = perindex+n
            elseif (perindex>n) then
                perindex = perindex-n
            endif
        end do

    end function


    ! coarsen the block by one level
    subroutine restriction_3D(fine, coarse)

        implicit none

        real(kind=rk), dimension(1:,1:,1:), intent(in)  :: fine
        real(kind=rk), dimension(1:,1:,1:), intent(out) :: coarse
        integer(kind=ik), dimension(3) :: nfine, ncoarse

        ncoarse(1) = size(coarse,1)
        nfine(1) = size(fine,1)
        ncoarse(2) = size(coarse,2)
        nfine(2) = size(fine,2)
        ncoarse(3) = size(coarse,3)
        nfine(3) = size(fine,3)

        if ( 2*ncoarse(1)-1 /= nfine(1) .or. 2*ncoarse(2)-1 /= nfine(2) .or. 2*ncoarse(3)-1 /= nfine(3) ) then
            call abort(888192,"ERROR: restriction_3D: arrays wrongly sized..")
        endif

        coarse(:, :, :) = fine(1:nfine(1):2,1:nfine(2):2,1:nfine(3):2)

    end subroutine restriction_3D

    ! refine the block by one level
    subroutine prediction_2D(coarse, fine, order_predictor)
        implicit none

        real(kind=rk), dimension(1:,1:), intent(inout) :: fine
        real(kind=rk), dimension(1:,1:), intent(inout) :: coarse
        character(len=*), intent(in)                :: order_predictor

        integer(kind=ik) :: i, j, l
        integer(kind=ik) :: nxcoarse, nxfine
        integer(kind=ik) :: nycoarse, nyfine
        integer(kind=ik) :: ixfine, iyfine
        ! interpolation coefficients
        ! a: one sided, b: central
        real(kind=rk) :: a(4), b(2)

        nxcoarse = size(coarse, 1)
        nxfine   = size(fine  , 1)
        nycoarse = size(coarse, 2)
        nyfine   = size(fine  , 2)

        if ( (2*nxcoarse-1 /= nxfine) .or. (2*nycoarse-1 /= nyfine)) then
            write(*,*) "coarse:", nxcoarse, nycoarse, "fine:", nxfine, nyfine
            call abort(888193,"ERROR: prediction_2D: arrays wrongly sized..")
        endif

        if ( ((nxfine<7) .or. (nyfine<7)).and.(order_predictor=="multiresolution_4th") ) then
            write(*,*) "coarse:", nxcoarse, nycoarse, "fine:", nxfine, nyfine
            call abort(888193,"ERROR: prediction_2D: not enough points for 4th order one-sided interp.")
        endif

        ! fill matching points: the coarse and fine grid share a lot of points (as the
        ! fine grid results from insertion of one point between each coarse point)
        fine(1:nxfine:2, 1:nyfine:2) = coarse(:,:)


        if ( order_predictor == "multiresolution_2nd" ) then
            !-------------------------------------------------------------------
            ! second order interpolation
            !-------------------------------------------------------------------
            ! y direction
            do i = 2, nxfine, 2
                do j = 1, nyfine, 2
                    fine(i,j) = ( fine(i-1, j) + fine(i+1, j) ) * 0.5_rk
                end do
            end do

            ! x direction
            do i = 1, nxfine
                do j = 2, nyfine, 2
                    fine(i,j) = ( fine(i, j-1) + fine(i, j+1) ) * 0.5_rk
                end do
            end do

        elseif ( order_predictor == "multiresolution_4th"  ) then
            !-------------------------------------------------------------------
            ! fourth order interpolation
            !-------------------------------------------------------------------
            ! one-side coeffs:
            a = (/ 5.0_rk/16.0_rk, 15.0_rk/16.0_rk, -5.0_rk/16.0_rk, 1.0_rk/16.0_rk /)
            ! centered coeffs (symmetric)
            b = (/ 9.0_rk/16.0_rk, -1.0_rk/16.0_rk /)

            ! step (a)
            ! first columns (x: const y: variable )
            ! these points require one-sided interpolation.
            fine( 2, 1:nyfine:2 ) = a(1)*fine( 1, 1:nyfine:2 ) &
            + a(2)*fine( 3, 1:nyfine:2 ) &
            + a(3)*fine( 5, 1:nyfine:2 ) &
            + a(4)*fine( 7, 1:nyfine:2 )

            ! last columns (same as above)
            fine( nxfine-1, 1:nyfine:2 ) = a(4)*fine( nxfine-6, 1:nyfine:2 ) &
            + a(3)*fine( nxfine-4, 1:nyfine:2 ) &
            + a(2)*fine( nxfine-2, 1:nyfine:2 ) &
            + a(1)*fine( nxfine,   1:nyfine:2 )

            ! interpolate regular columns
            do ixfine =  4, nxfine-3, 2
                fine( ixfine, 1:nyfine:2 ) = b(2)*fine( ixfine-3, 1:nyfine:2 ) &
                + b(1)*fine( ixfine-1, 1:nyfine:2 ) &
                + b(1)*fine( ixfine+1, 1:nyfine:2 ) &
                + b(2)*fine( ixfine+3, 1:nyfine:2 )
            enddo


            ! At this point, we have only every 2nd complete row missing
            ! so from now on, no step size 2 anymore

            ! first row
            ! these points requie one-sided interpolation.
            fine( 1:nxfine, 2 ) = a(1)*fine( 1:nxfine, 1 ) &
            + a(2)*fine( 1:nxfine, 3 ) &
            + a(3)*fine( 1:nxfine, 5 ) &
            + a(4)*fine( 1:nxfine, 7 )
            ! last row (same as above)
            fine( 1:nxfine, nyfine-1 ) = a(4)*fine( 1:nxfine, nyfine-6) &
            + a(3)*fine( 1:nxfine, nyfine-4 ) &
            + a(2)*fine( 1:nxfine, nyfine-2 ) &
            + a(1)*fine( 1:nxfine, nyfine )
            ! remaining interior rows
            do iyfine =  4, nyfine-3, 2
                fine( 1:nxfine, iyfine ) = b(2)*fine( 1:nxfine, iyfine-3 ) &
                + b(1)*fine( 1:nxfine, iyfine-1 ) &
                + b(1)*fine( 1:nxfine, iyfine+1 ) &
                + b(2)*fine( 1:nxfine, iyfine+3 )
            enddo

        else
            ! error case
            call abort(888194,"ERROR: prediction_2D: wrong method..")
        endif

    end subroutine prediction_2D


    ! refine the block by one level
    subroutine prediction_3D(coarse, fine, order_predictor)
        implicit none

        real(kind=rk), dimension(1:,1:,1:), intent(inout) :: fine
        real(kind=rk), dimension(1:,1:,1:), intent(inout) :: coarse
        character(len=*), intent(in) :: order_predictor

        integer(kind=ik) :: i, j, k
        integer(kind=ik) :: nxcoarse, nxfine, nycoarse, nyfine, nzcoarse, nzfine
        integer(kind=ik) :: ixfine, iyfine, izfine
        ! interpolation coefficients
        ! a: one sided, b: central
        real(kind=rk) :: a(4), b(2)

        nxcoarse = size(coarse, 1)
        nxfine   = size(fine  , 1)
        nycoarse = size(coarse, 2)
        nyfine   = size(fine  , 2)
        nzcoarse = size(coarse, 3)
        nzfine   = size(fine  , 3)

        if ( 2*nxcoarse-1 /= nxfine .or. 2*nycoarse-1 /= nyfine .or. 2*nzcoarse-1 /= nzfine ) then
            call abort(888195,"ERROR: prediction_3D: arrays wrongly sized..")
        endif

        ! fill matching points: the coarse and fine grid share a lot of points (as the
        ! fine grid results from insertion of one point between each coarse point)
        fine(1:nxfine:2, 1:nyfine:2, 1:nzfine:2) = coarse(:, :, :)

        if ( order_predictor == "multiresolution_2nd" ) then
            !-------------------------------------------------------------------
            ! second order interpolation
            !-------------------------------------------------------------------
            ! y direction
            do k = 1, nzfine, 2
                do j = 1, nyfine, 2
                    do i = 2, nxfine, 2
                        fine(i,j,k) = ( fine(i-1, j, k) + fine(i+1, j, k) ) * 0.5_rk
                    end do
                end do
            end do

            ! x direction
            do k = 1, nzfine, 2
                do j = 2, nyfine, 2
                    do i = 1, nxfine
                        fine(i,j,k) = ( fine(i, j-1, k) + fine(i, j+1, k) ) * 0.5_rk
                    end do
                end do
            end do

            ! z direction
            do k = 2, nzfine, 2
                do j = 1, nyfine
                    do i = 1, nxfine
                        fine(i,j,k) = ( fine(i, j, k-1) + fine(i, j, k+1) ) * 0.5_rk
                    end do
                end do
            end do

        elseif ( order_predictor == "multiresolution_4th"  ) then
            !-----------------------------------------------------------------------
            ! fourth order interpolation
            !-----------------------------------------------------------------------
            ! one-side coeffs:
            a = (/ 5.0_rk/16.0_rk, 15.0_rk/16.0_rk, -5.0_rk/16.0_rk, 1.0_rk/16.0_rk /)
            ! centered coeffs (symmetric)
            b = (/ 9.0_rk/16.0_rk, -1.0_rk/16.0_rk /)

            do izfine = 1, nzfine, 2
                ! --> in the planes, execute the 2d code
                ! step (a)
                ! first columns (x: const y: variable )
                ! these points requie one-sided interpolation.
                fine( 2, 1:nyfine:2, izfine ) = a(1)*fine( 1, 1:nyfine:2, izfine ) &
                + a(2)*fine( 3, 1:nyfine:2, izfine ) &
                + a(3)*fine( 5, 1:nyfine:2, izfine ) &
                + a(4)*fine( 7, 1:nyfine:2, izfine )

                ! last columns (same as above)
                fine( nxfine-1, 1:nyfine:2, izfine ) = a(4)*fine( nxfine-6, 1:nyfine:2, izfine ) &
                + a(3)*fine( nxfine-4, 1:nyfine:2, izfine ) &
                + a(2)*fine( nxfine-2, 1:nyfine:2, izfine ) &
                + a(1)*fine( nxfine,   1:nyfine:2, izfine )

                ! interpolate regular columns
                do ixfine =  4, nxfine-3, 2
                    fine( ixfine, 1:nyfine:2, izfine ) = b(2)*fine( ixfine-3, 1:nyfine:2, izfine ) &
                    + b(1)*fine( ixfine-1, 1:nyfine:2, izfine ) &
                    + b(1)*fine( ixfine+1, 1:nyfine:2, izfine ) &
                    + b(2)*fine( ixfine+3, 1:nyfine:2, izfine )
                enddo


                ! At this point, we have only every 2nd complete row missing
                ! so from now on, no step size 2 anymore

                ! first row
                ! these points requie one-sided interpolation.
                fine( 1:nxfine, 2, izfine ) = a(1)*fine( 1:nxfine, 1, izfine ) &
                + a(2)*fine( 1:nxfine, 3, izfine ) &
                + a(3)*fine( 1:nxfine, 5, izfine ) &
                + a(4)*fine( 1:nxfine, 7, izfine )
                ! last row (same as above)
                fine( 1:nxfine, nyfine-1, izfine ) = a(4)*fine( 1:nxfine, nyfine-6, izfine ) &
                + a(3)*fine( 1:nxfine, nyfine-4, izfine ) &
                + a(2)*fine( 1:nxfine, nyfine-2, izfine ) &
                + a(1)*fine( 1:nxfine, nyfine, izfine )
                ! remaining interior rows
                do iyfine =  4, nyfine-3, 2
                    fine( 1:nxfine, iyfine, izfine ) = b(2)*fine( 1:nxfine, iyfine-3, izfine ) &
                    + b(1)*fine( 1:nxfine, iyfine-1, izfine ) &
                    + b(1)*fine( 1:nxfine, iyfine+1, izfine ) &
                    + b(2)*fine( 1:nxfine, iyfine+3, izfine )
                enddo
            enddo

            ! interpolate the z-direction (completely missing planes, no step 2)
            ! first plane
            fine( :, :, 2 ) = a(1)*fine( :, :, 1 ) &
            + a(2)*fine( :, :, 3 ) &
            + a(3)*fine( :, :, 5 ) &
            + a(4)*fine( :, :, 7 )
            ! last plane
            fine( :, :, nzfine-1 ) = a(4)*fine( :, :, nzfine-6 ) &
            + a(3)*fine( :, :, nzfine-4 ) &
            + a(2)*fine( :, :, nzfine-2 ) &
            + a(1)*fine( :, :, nzfine )
            ! remaining planes
            do izfine =  4, nzfine-3, 2
                fine( :, :, izfine ) = b(2)*fine( :, :, izfine-3 ) &
                + b(1)*fine( :, :, izfine-1 ) &
                + b(1)*fine( :, :, izfine+1 ) &
                + b(2)*fine( :, :, izfine+3 )
            enddo


        else
            ! error case
            call abort(888196,"ERROR: prediction_2D: wrong method..")
        endif

    end subroutine prediction_3D


    ! subroutine prediction1D(coarse, fine)
    !
    !     implicit none
    !
    !     real(kind=rk), dimension(1:), intent(out) :: fine
    !     real(kind=rk), dimension(1:), intent(in) :: coarse
    !
    !     integer(kind=ik) :: k, nfine, ncoarse
    !     real(kind=rk) :: a, b
    !
    !     ncoarse = size(coarse,1)
    !     nfine = size(fine,1)
    !
    !     if ( 2*ncoarse-1 /= nfine ) then
    !       call abort(888197,"ERROR: prediction1d: arrays wrongly sized..")
    !     endif
    !
    !     ! this is the multiresolution predition operator.
    !     ! it pushes a signal from a coarser level to the next higher by
    !     ! interpolation
    !
    !     fine(1:nfine:2) = coarse(:)
    !
    !     ! fourth order:
    !     a = 9.0_rk/16.0_rk
    !     b =-1.0_rk/16.0_rk
    !
    !     fine(2)     = (5.0_rk/16.0_rk)*coarse(1)+(15.0_rk/16.0_rk)*coarse(2)-(5.0_rk/16.0_rk)*coarse(3)+(1.0_rk/16.0_rk)*coarse(4)
    !     fine(nfine-1) = (1.0_rk/16.0_rk)*coarse(ncoarse-3) -(5.0_rk/16.0_rk)*coarse(ncoarse-2) +(15.0_rk/16.0_rk)*coarse(ncoarse-1) +(5.0_rk/16.0_rk)*coarse(ncoarse)
    !
    !     do k = 2, ncoarse-2
    !         fine(2*k) = a*coarse(k)+a*coarse(k+1)+b*coarse(k-1)+b*coarse(k+2)
    !     end do
    !
    ! end subroutine prediction1D

end module
