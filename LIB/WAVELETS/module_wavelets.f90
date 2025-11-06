module module_wavelets
    use module_params
    use module_treelib
    use module_operators, only : setup_FD1_left_stencil, setup_FD2_stencil

    implicit none

    ! Everything is public by default
    PUBLIC
    SAVE

    ! Private interpolation stencils with added zeros
    private :: Stencil_Int_2nd, Stencil_Int_4th, Stencil_Int_6th, Stencil_Int_8th, Stencil_Int_10th, Stencil_Int_12th
    real(kind=rk), parameter :: Stencil_Int_2nd(-1:1) = (/ 1.0_rk, 2.0_rk, 1.0_rk /) / 2.0_rk
    real(kind=rk), parameter :: Stencil_Int_4th(-3:3) = (/ -1.0_rk, 0.0_rk, 9.0_rk, 16.0_rk, 9.0_rk, 0.0_rk, -1.0_rk /) / 2.0_rk**4
    real(kind=rk), parameter :: Stencil_Int_6th(-5:5) = (/ 3.0_rk, 0.0_rk, -25.0_rk, 0.0_rk, 150.0_rk, 2.0_rk**8, &
        150.0_rk, 0.0_rk, -25.0_rk, 0.0_rk, 3.0_rk /) / 2.0_rk**8
    real(kind=rk), parameter :: Stencil_Int_8th(-7:7) = (/ -5.0_rk, 0.0_rk, 49.0_rk, 0.0_rk, -245.0_rk, 0.0_rk, 1225.0_rk, &
        2.0_rk**11, 1225.0_rk, 0.0_rk, -245.0_rk, 0.0_rk, 49.0_rk, 0.0_rk, -5.0_rk /) / 2.0_rk**11
    real(kind=rk), parameter :: Stencil_Int_10th(-9:9) = (/ 35.0_rk, 0.0_rk, -405.0_rk, 0.0_rk, 2268.0_rk, 0.0_rk, -8820.0_rk, 0.0_rk, 39690.0_rk, &
        2.0_rk**16, 39690.0_rk, 0.0_rk, -8820.0_rk, 0.0_rk, 2268.0_rk, 0.0_rk, -405.0_rk, 0.0_rk, 35.0_rk /) / 2.0_rk**16
    real(kind=rk), parameter :: Stencil_Int_12th(-11:11) = (/ -63.0_rk, 0.0_rk, 847.0_rk, 0.0_rk, -5445.0_rk, 0.0_rk, 22869.0_rk, 0.0_rk, -76230.0_rk, 0.0_rk, 320166.0_rk, &
        2.0_rk**19, 320166.0_rk, 0.0_rk, -76230.0_rk, 0.0_rk, 22869.0_rk, 0.0_rk, -5445.0_rk, 0.0_rk, 847.0_rk, 0.0_rk, -63.0_rk /) / 2.0_rk**19

contains

#include "conversion_routines.f90"
#include "wavelet_decomposition_reconstruction.f90"


    ! this function is only used in merge_blocks, we could remove it
    subroutine restriction_prefilter_vct(params, u, u_filtered)
        implicit none
        type(type_params), intent(in) :: params
        real(kind=rk), dimension(:,:,:,:), intent(in) :: u
        real(kind=rk), dimension(:,:,:,:), intent(out) :: u_filtered

        if (.not. allocated(params%HD)) call abort(71717172, "wavelet not setup")

        call blockFilterXYZ_vct(params, u, u_filtered, params%HD, lbound(params%HD, dim=1), ubound(params%HD, dim=1))
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
    


    !> \brief Refine the block by one level using the prediciton operator
    ! NOTE:
    ! ======
    ! As of 07 Aug 2023, this routine does no longer include once-sided interpolation stencils
    ! because they must not be used anyways. One-sided interpolation corresponds to different
    ! wavelet functions near the boundaries. This functionality was used a long time ago with
    ! the "lazy wavelets" CDF20 and CDF40, but we have always seen to the code not using the one-sided
    ! stencils. Points that cannot be interpolated (where onse-sided interpolation would be used) are returned
    ! zero.
    ! o: matching points (checkerboard copy), x: zeros returned, i=interpolated properly
    !
    ! Input:                     Output (4th order):
    ! o   o   o   o   o          o x o x o x o x o
    !                            x x x x x x x x x
    ! o   o   o   o   o          o x o i o i o x o
    !                            x x i i i i i x x
    ! o   o   o   o   o          o x o i o i o x o
    !                            x x i i i i i x x
    ! o   o   o   o   o          o x o i o i o x o
    !                            x x x x x x x x x
    ! o   o   o   o   o          o x o x o x o x o
    ! Why do we have this routine?
    ! We could, after copying (checkerboard), also simply apply the wavelet%HR filter
    ! to all points. This yields the same result as this special routine. However,
    ! note how even for copied (odd indices) points, a number of multiplications and subsequent summation
    ! is required (although we multiply by zeros). For even points (which are indeed interpolated)
    ! we require as many multiplications as the filter HR, which contains zeros.
    ! This routine is thus more efficient (it skips odd points entirely and does not multiply by zero).
    ! As it is called for every ghost nodes patch (not just to upsample entire blocks), it is performance-critical.
    !
    ! inplace usage for a block:
    ! all prediction(hvy_block(1:params%Bs(1)+2*params%g:2,1:params%Bs(2)+2*params%g:2,:,ic,hvy_id), &
    !                hvy_block(1:params%Bs(1)+2*params%g-1,1:params%Bs(2)+2*params%g-1,:,ic,hvy_id), params%order_predictor)
    subroutine prediction(coarse, fine, order_predictor)
        implicit none

        real(kind=rk), dimension(:,:,:), intent(inout) :: fine
        real(kind=rk), dimension(:,:,:), intent(inout) :: coarse
        character(len=*), intent(in) :: order_predictor

        integer(kind=ik) :: i, j, k
        integer(kind=ik) :: nxcoarse, nxfine, nycoarse, nyfine, nzcoarse, nzfine
        integer(kind=ik) :: ixfine, iyfine, izfine, N, shift, shift_fine
        ! interpolation coefficients
        real(kind=rk), allocatable :: c(:)

        nxcoarse = size(coarse, 1)
        nxfine   = size(fine  , 1)
        nycoarse = size(coarse, 2)
        nyfine   = size(fine  , 2)
        nzcoarse = size(coarse, 3)
        nzfine   = size(fine  , 3)

        ! setup interpolation coefficients
        ! All orders from Taylor expansion
        ! Order 1-8 verified with Deriaz. JCP2023
        select case(order_predictor)
        case ("multiresolution_2nd")
            allocate(c(1:2))
            c = (/ 0.5_rk, 0.5_rk /)

        case ("multiresolution_4th")
            allocate(c(1:4))
            c = (/ -1.0_rk, 9.0_rk, 9.0_rk, -1.0_rk /) / 16.0_rk

        case ("multiresolution_6th")
            allocate(c(1:6))
            c = (/ 3.0_rk, -25.0_rk, 150.0_rk, 150.0_rk, -25.0_rk, 3.0_rk /) / 256.0_rk

        case ("multiresolution_8th")
            allocate(c(1:8))
            c = (/ -5.0_rk, 49.0_rk, -245.0_rk, 1225.0_rk, 1225.0_rk, -245.0_rk, 49.0_rk, -5.0_rk /) / 2048.0_rk

        case ("multiresolution_10th")
            allocate(c(1:10))
            c = (/ 35.0_rk, -405.0_rk, 2268.0_rk, -8820.0_rk, 39690.0_rk, 39690.0_rk, -8820.0_rk, 2268.0_rk, -405.0_rk, 35.0_rk /) / 65536.0_rk
        
        case ("multiresolution_12th")
            allocate(c(1:12))
            c = (/ -63.0_rk, 847.0_rk, -5445.0_rk, 22869.0_rk, -76230.0_rk, 320166.0_rk, 320166.0_rk, -76230.0_rk, 22869.0_rk, -5445.0_rk, 847.0_rk, -63.0_rk /) / 524288.0_rk

        case default
            call abort(23070811,"Error: unkown order_predictor="//trim(adjustl(order_predictor)))

        end select
        N = size(c,1)/2  ! which value can be interpolated first? 1-based, 4th and 6th skip 1 or 2 points

#ifdef DEV
        if ( 2*nxcoarse-1 /= nxfine .or. 2*nycoarse-1 /= nyfine .or. 2*nzcoarse-1 /= nzfine ) then
            call abort(888195,"ERROR: prediction - arrays wrongly sized..")
        endif
#endif

        ! prepare grid, set all values to zero or only those that will be overwritte (for inplace usage)
        ! as only edge points would be set to 0 that should not be used and the rest is overwritten, we do not need to do it
        ! but sometimes we init them as NaN and NaN is not fun so let's set them anyways
        fine(2:nxfine:2, :, :) = 0.0_rk
        fine(:, 2:nyfine:2, :) = 0.0_rk
        fine(:, :, 2:nzfine:2) = 0.0_rk
        ! fine = 0.0_rk

        ! fill matching points: the coarse and fine grid share a lot of points (as the
        ! fine grid results from insertion of one point between each coarse point).
        ! Sometimes called checkerboard copying
        fine(1:nxfine:2, 1:nyfine:2, 1:nzfine:2) = coarse(:, :, :)

        ! matrix operation version
        ! x-interpolation
        do shift = 1, size(c,1)
            shift_fine = -size(c,1)+(2*shift-1)
            fine(2*N:nxfine-(2*N-1):2,1:nyfine:2,1:nzfine:2) = fine(2*N:nxfine-(2*N-1):2,1:nyfine:2,1:nzfine:2) + c(shift)*fine(2*N+shift_fine:nxfine-(2*N-1)+shift_fine:2,1:nyfine:2,1:nzfine:2)
        end do
        ! y-interpolation
        do shift = 1, size(c,1)
            shift_fine = -size(c,1)+(2*shift-1)
            fine(:,2*N:nyfine-(2*N-1):2,1:nzfine:2) = fine(:,2*N:nyfine-(2*N-1):2,1:nzfine:2) + c(shift)*fine(:,2*N+shift_fine:nyfine-(2*N-1)+shift_fine:2,1:nzfine:2)
        end do
        ! z-interpolation
        do shift = 1, size(c,1)
            shift_fine = -size(c,1)+(2*shift-1)
            fine(:,:,2*N:nzfine-(2*N-1):2) = fine(:,:,2*N:nzfine-(2*N-1):2) + c(shift)*fine(:,:,2*N+shift_fine:nzfine-(2*N-1)+shift_fine:2)
        end do

        ! do izfine = 1, nzfine, 2
        !     ! in the z=const planes, we execute the 2D code.
        !     do ixfine= 2*N, nxfine-(2*N-1), 2
        !         do iyfine =  1, nyfine, 2
        !             ! note in this implementation, interp coeffs run
        !             ! c(1:2), c(1:4), c(1:6) for 2nd, 4th, 6th order respectively
        !             do shift = 1, size(c,1)
        !                 shift_fine = -size(c,1)+(2*shift-1)
        !                 fine(ixfine,iyfine,izfine) = fine(ixfine,iyfine,izfine) + c(shift)*fine(ixfine+shift_fine,iyfine,izfine)
        !             end do
        !         end do
        !     end do

        !     do ixfine = 1, nxfine, 1
        !         do iyfine =  2*N, nyfine-(2*N-1), 2
        !             ! note in this implementation, interp coeffs run
        !             ! c(1:2), c(1:4), c(1:6) for 2nd, 4th, 6th order respectively
        !             do shift = 1, size(c,1)
        !                 shift_fine = -size(c,1)+(2*shift-1)
        !                 fine(ixfine,iyfine,izfine) = fine(ixfine,iyfine,izfine) + c(shift)*fine(ixfine, iyfine+shift_fine,izfine)
        !             end do
        !         end do
        !     end do
        ! enddo

        ! ! finally, only 1D interpolation along z is missing.
        ! do izfine =  2*N, nzfine-(2*N-1), 2
        !     ! note in this implementation, interp coeffs run
        !     ! c(1:2), c(1:4), c(1:6) for 2nd, 4th, 6th order respectively
        !     do shift = 1, size(c,1)
        !         shift_fine = -size(c,1)+(2*shift-1)
        !         fine(:,:,izfine) = fine(:,:,izfine) + c(shift)*fine(:,:,izfine+shift_fine)
        !     end do
        ! enddo

    end subroutine


    ! Filters a block, all internal nodes, assuming g >= support of filter (crashes otherwise)
    ! Bs=6, g=2, HD=[-1:1]
    ! g g g g g g g g g g          g g g g g g g g g g          g g g g g g g g g g
    ! g g g g g g g g g g          g g p p p p p p g g          g g p g p g g g g g
    ! g g i i i i i i g g          g g f f f f f f g g          g g p i p i p i g g
    ! g g i i i i i i g g          g g f f f f f f g g          g g f i f i f i g g
    ! g g i i i i i i g g          g g f f f f f f g g          g g p i p i p i g g
    ! g g i i i i i i g g          g g f f f f f f g g          g g f i f i f i g g
    ! g g i i i i i i g g          g g f f f f f f g g          g g p i p i p i g g
    ! g g i i i i i i g g          g g f f f f f f g g          g g f i f i f i g g
    ! g g g g g g g g g g          g g p p p p p p g g          g g p g p g g g g g
    ! g g g g g g g g g g          g g g g g g g g g g          g g g g g g g g g g
    ! Fig1: input                    Fig2: out filtered           Fig3: filtered with restriction
    ! g=ghost i=internal f=filtered p=partially filtered
    ! For filtered: After the first dimension filtering we can ignore the edge-values (left and right side), this might not seem much for a full block
    !       but for ghost patches we sometimes have very small patches with large support.
    !       Additionally, if filter < g then some values at the edges need not be touched (top and bottom row)
    ! For filtered with restriction: Every dimension we only filter the restricted points (so NOT every second and at ghost point where filter)
    !       Afterwards, the next filtering considers for this dimension only those values so we have computing costs of (1/2 + 1/4 + 1/8) / 3 instead of 1
    ! As edges can be partially filtered or not they should not be used!
    subroutine blockFilterXYZ_vct( params, u, u_filtered, coefs_filter, fl_l, fl_r, do_restriction)
        implicit none
        type (type_params), intent(in) :: params
        real(kind=rk), dimension(:,:,:,:), intent(in) :: u
        real(kind=rk), dimension(:,:,:,:), intent(inout) :: u_filtered
        integer(kind=ik), intent(in) :: fl_l, fl_r  !< Filter length left and right
        real(kind=rk), intent(in) :: coefs_filter(fl_l:fl_r)
        logical, intent(in), optional :: do_restriction  !< when we do restriction anyways we can skip every second point

        integer(kind=ik) :: ix, iy, iz, ic, nx, ny, nz, nc, shift, ifs(1:3), ife(1:3), ils(1:3), ile(1:3)
        real(kind=rk), allocatable, save :: u_tmp(:,:,:)
        integer(kind=ik) :: s  ! sample_num, is either 1 or 2

        ! if the filter is just 1, then we copy and we're done.
        ! Yes, we use such stupid filters. They are in the CDFX0 wavelets (X=2,4,6,8,10,12)
        if (fl_l==0 .and. fl_r==0 .and. abs(coefs_filter(0)-1.0_rk)<=1.0e-10_rk) then
            u_filtered = u
            return
        endif


        nx = size(u, 1)
        ny = size(u, 2)
        nz = size(u, 3)
        nc = size(u, 4)

        ! compute loop boundaries - we have two different loop boundaries that we use
        ! 1) Filter has not yet been applied, we need g+1+fl_l and Bs+g+fl_r
        ! 2) Filter has already been applied for this dimension, we need only g+1 and Bs+g
        ifs = params%g+1+fl_l
        ife = params%Bs+params%g+fl_r
        ils = params%g+1
        ile = params%Bs+params%g
        if (params%dim == 2) then
            ifs(3) = 1
            ife(3) = 1
            ils(3) = 1
            ile(3) = 1
        endif

        if ((abs(fl_l) > params%g).or.(fl_r>params%g)) then
            write(*,'(A, i0, A, i0, A, i0)') "Filter size: ", fl_l, " / ", fl_r, ", but g= ", params%g
            call abort(202302208, "For applying the filter, not enough ghost nodes")
        endif

        ! stride, used if we do restriction
        s = 1
        if (present(do_restriction)) then
            if (do_restriction) s = 2
        endif

        ! it would be tempting to parallelize over ic here, but the temporary work array can become quite large.
        ! It is no benefit to parallelize anyway, as in memory they are very far apart.
        if (.not. allocated(u_tmp)) allocate( u_tmp(1:nx,1:ny,1:nz) )

        do ic = 1, nc
            
            u_tmp(:,:,:) = u(:,:,:,ic)
            u_filtered(:,:,:,ic) = u_tmp(:,:,:)

            ! filter results are consecutively added so we need to set to 0 at the beginning
            u_filtered(ils(1):ile(1):s, ifs(2):ife(2), ifs(3):ife(3), ic) = 0.0_rk
            do ix = ils(1), ile(1), s
                do shift = fl_l, fl_r
                    u_filtered(ix, ifs(2):ife(2), ifs(3):ife(3), ic) = u_filtered(ix, ifs(2):ife(2), ifs(3):ife(3), ic) &
                    + u_tmp(ix+shift, ifs(2):ife(2), ifs(3):ife(3))*coefs_filter(shift)
                enddo
            enddo

            u_tmp(:,:,:) = u_filtered(:,:,:,ic)
            ! filter results are consecutively added so we need to set to 0 at the beginning
            u_filtered(ils(1):ile(1):s, ils(2):ile(2):s, ifs(3):ife(3), ic) = 0.0_rk
            do iy = ils(2), ile(2), s
                do shift = fl_l, fl_r
                    u_filtered(ils(1):ile(1):s, iy, ifs(3):ife(3), ic) = u_filtered(ils(1):ile(1):s, iy, ifs(3):ife(3), ic) &
                    + u_tmp(ils(1):ile(1):s, iy+shift, ifs(3):ife(3))*coefs_filter(shift)
                enddo
            enddo

            if (params%dim == 3) then
                u_tmp(:,:,:) = u_filtered(:,:,:,ic)
                u_filtered(ils(1):ile(1):s, ils(2):ile(2):s, ils(3):ile(3):s, ic) = 0.0_rk
                do iz = ils(3), ile(3), s
                    do shift = fl_l, fl_r
                        u_filtered(ils(1):ile(1):s, ils(2):ile(2):s, iz, ic) = u_filtered(ils(1):ile(1):s, ils(2):ile(2):s, iz, ic) &
                        + u_tmp(ils(1):ile(1):s, ils(2):ile(2):s, iz+shift)*coefs_filter(shift)
                    enddo
                enddo
            endif
        enddo

    end subroutine

    ! filter everywhere where possible: do not start at the interior points, but
    ! start at the support of the filter
    ! Bs=6, g=2, HD=[-1:1]
    ! g g g g g g g g g g          g p p p p p p p p g          g p g p g p g p g g
    ! g g g g g g g g g g          g f f f f f f f f g          g p g p g p g p g g
    ! g g i i i i i i g g          g f f f f f f f f g          g f i f i f i f g g
    ! g g i i i i i i g g          g f f f f f f f f g          g p i p i p i p g g
    ! g g i i i i i i g g          g f f f f f f f f g          g f i f i f i f g g
    ! g g i i i i i i g g          g f f f f f f f f g          g p i p i p i p g g
    ! g g i i i i i i g g          g f f f f f f f f g          g f i f i f i f g g
    ! g g i i i i i i g g          g f f f f f f f f g          g p i p i p i p g g
    ! g g g g g g g g g g          g f f f f f f f f g          g f g f g f g f g g
    ! g g g g g g g g g g          g p p p p p p p p g          g p p p p p p p g g
    ! Fig1: input                    Fig2: out filtered           Fig3: filtered with restriction
    ! g= ghost i=internal f=filtered p=partially filtered
    ! For filtered: After the first dimension filtering we can ignore the edge-values (left and right side), this might not seem much for a full block
    !       but for ghost patches we sometimes have very small patches with large support
    ! For filtered with restriction: Every dimension we only filter the restricted points (so NOT every second and at edges where filter cannot be applied)
    !       Afterwards, the next filtering considers for this dimension only those values so we have computing costs of (1/2 + 1/4 + 1/8) / 3 instead of 1
    ! As edges can be partially filtered or not they should not be used!
    subroutine blockFilterXYZ_wherePossible_vct( params, u, u_filtered, coefs_filter, fl_l, fl_r, do_restriction)
        implicit none
        type (type_params), intent(in) :: params
        real(kind=rk), dimension(:,:,:,:), intent(in) :: u
        real(kind=rk), dimension(:,:,:,:), intent(inout) :: u_filtered
        integer(kind=ik), intent(in) :: fl_l, fl_r  !< Filter length left and right
        real(kind=rk), intent(in) :: coefs_filter(fl_l:fl_r)
        logical, intent(in), optional :: do_restriction  !< when we do restriction anyways we can skip every second point

        integer(kind=ik) :: ix, iy, iz, nx, ny, nz, nc, shift, g, Bs(1:3)
        real(kind=rk), allocatable, save :: u_tmp(:,:,:,:)
        integer(kind=ik) :: s  ! sample_num, is either 1 or 2

        ! if the filter is just 1, then we copy and we're done.
        ! Yes, we use such stupid filters. They are in the CDFX0 wavelets (X=2,4)
        if (fl_l==0 .and. fl_r==0 .and. abs(coefs_filter(0)-1.0_rk)<=1.0e-10_rk) then
            u_filtered = u
            return
        endif

        s = 1
        if (present(do_restriction)) then
            if (do_restriction) s = 2
        endif


        nx = size(u, 1)
        ny = size(u, 2)
        nz = size(u, 3)
        nc = size(u, 4)
        g  = params%g
        Bs = params%Bs

        if ((abs(fl_l) > g).or.(fl_r>g)) then
            write(*,'(A, i0, A, i0, A, i0)') "Filter size: ", fl_l, " / ", fl_r, ", but g= ", g
            call abort(202302209, "For applying the filter, not enough ghost nodes")
        endif

        if (.not. allocated(u_tmp)) allocate( u_tmp(1:nx,1:ny,1:nz,1:nc) )
        u_tmp = u

        u_filtered = u_tmp  ! copy everything
        u_filtered(-fl_l+1:nx-fl_r:s, :, :, :) = 0.0_rk  ! these values will be filtered and overwritten
        do ix = -fl_l+1, nx-fl_r, s
            do shift = fl_l, fl_r
                u_filtered(ix, :, :, :) = u_filtered(ix, :, :, :) + u_tmp(ix+shift, :, :, :)*coefs_filter(shift)
            enddo
        enddo

        u_tmp = u_filtered
        ! u_filtered = u_tmp  ! not necessary
        u_filtered(-fl_l+1:nx-fl_r:s, -fl_l+1:ny-fl_r:s, :, :) = 0.0_rk  ! these values will be filtered and overwritten again
        do iy = -fl_l+1, ny-fl_r, s
            do shift = fl_l, fl_r
                ! after previous filtering we can ignore some x-values which were not effected to decrease the computational load
                u_filtered(-fl_l+1:nx-fl_r:s, iy, :, :) = u_filtered(-fl_l+1:nx-fl_r:s, iy, :, :) &
                   + u_tmp(-fl_l+1:nx-fl_r:s, iy+shift, :, :)*coefs_filter(shift)
            enddo
        enddo


        if (nz /= 1) then

            u_tmp = u_filtered
            ! u_filtered = u_tmp  ! not necessary
            u_filtered(-fl_l+1:nx-fl_r:s, -fl_l+1:ny-fl_r:s, -fl_l+1:nz-fl_r:s, :) = 0.0_rk  ! these values will be filtered and overwritten again
            do iz = -fl_l+1, nz-fl_r, s
                do shift = fl_l, fl_r
                    ! after previous filtering we can ignore some x- and y-values which were not effected to decrease the computational load
                    u_filtered(-fl_l+1:nx-fl_r:s, -fl_l+1:ny-fl_r:s, iz, :) = u_filtered(-fl_l+1:nx-fl_r:s, -fl_l+1:ny-fl_r:s, iz, :) &
                    + u_tmp(-fl_l+1:nx-fl_r:s, -fl_l+1:ny-fl_r:s, iz+shift, :)*coefs_filter(shift)
                enddo
            enddo\

        endif

    end subroutine


    ! applies one of params% HD GD HR GR filters in each direction
    ! ToDo: for HD and GD filter, we might benefit from considering restriction, for HR from skipping points
    subroutine blockFilterCustom_vct( params, u, u_filtered, filter_x, filter_y, filter_z )
        implicit none
        type (type_params), intent(in) :: params
        real(kind=rk), dimension(:,:,:,:), intent(in) :: u
        real(kind=rk), dimension(:,:,:,:), intent(inout) :: u_filtered
        character(len=*), intent(in) :: filter_x, filter_y, filter_z
        integer(kind=ik) :: a, b
        real(kind=rk), allocatable :: coefs_filter(:)
        integer(kind=ik) :: ix, iy, iz, nx, ny, nz, nc, shift, g, Bs(1:3)
        real(kind=rk), allocatable, save :: u_tmp(:,:,:,:)

        nx = size(u, 1)
        ny = size(u, 2)
        nz = size(u, 3)
        nc = size(u, 4)
        g  = params%g
        Bs = params%Bs

        !-----------------------------------------------------------------------
        ! Filter x
        !-----------------------------------------------------------------------
        if (.not. allocated(u_tmp)) allocate( u_tmp(1:nx, 1:ny, 1:nz, 1:nc) )
        u_tmp = u

        select case(filter_x)
        case("HD")
            if (.not. allocated(params%HD)) call abort(202302203, "Wavelet setup not done yet?! HD")
            a = lbound(params%HD, dim=1)
            b = ubound(params%HD, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%HD

        case("GD")
            if (.not. allocated(params%GD)) call abort(202302203, "Wavelet setup not done yet?! GD")
            a = lbound(params%GD, dim=1)
            b = ubound(params%GD, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%GD

        case("HR")
            if (.not. allocated(params%HR)) call abort(202302203, "Wavelet setup not done yet?! HR")
            a = lbound(params%HR, dim=1)
            b = ubound(params%HR, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%HR

        case("GR")
            if (.not. allocated(params%GR)) call abort(202302203, "Wavelet setup not done yet?! GR")
            a = lbound(params%GR, dim=1)
            b = ubound(params%GR, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%GR

        case default
            call abort(202302201, "Unknown wavelet filter (x-dir), must be one of HD GD HR GR")

        end select

        if ((abs(a) > g).or.(b>g)) then
            write(*,*) a, b, "but g=", g
            call abort(202302211, "For applying the filter, not enough ghost nodes")
        endif

        u_filtered = u_tmp
        u_filtered(g+1:Bs(1)+g, :, :, :) = 0.0_rk
        do ix = g+1, Bs(1)+g
            do shift = a, b
                u_filtered(ix, :, :, :) = u_filtered(ix, :, :, :) + u_tmp(ix+shift, :, :, :)*coefs_filter(shift)
            enddo
        enddo

        !-----------------------------------------------------------------------
        ! Filter y
        !-----------------------------------------------------------------------
        select case(filter_y)
        case("HD")
            a = lbound(params%HD, dim=1)
            b = ubound(params%HD, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%HD

        case("GD")
            a = lbound(params%GD, dim=1)
            b = ubound(params%GD, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%GD

        case("HR")
            a = lbound(params%HR, dim=1)
            b = ubound(params%HR, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%HR

        case("GR")
            a = lbound(params%GR, dim=1)
            b = ubound(params%GR, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%GR

        case default
            call abort(202302201, "Unknown wavelet filter (y-dir), must be one of HD GD HR GR")

        end select

        if ((abs(a) > g).or.(b>g)) then
            write(*,*) a, b, "but g=", g
            call abort(202302211, "For applying the filter, not enough ghost nodes")
        endif

        u_tmp = u_filtered
        u_filtered = u_tmp
        u_filtered(:, g+1:Bs(2)+g, :, :) = 0.0_rk
        do iy = g+1, Bs(2)+g
            do shift = a, b
                u_filtered(:, iy, :, :) = u_filtered(:, iy, :, :) + u_tmp(:, iy+shift, :, :)*coefs_filter(shift)
            enddo
        enddo

        !-----------------------------------------------------------------------
        ! Filter z
        !-----------------------------------------------------------------------
        if (nz == 1) return

        select case(filter_z)
        case("HD")
            a = lbound(params%HD, dim=1)
            b = ubound(params%HD, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%HD

        case("GD")
            a = lbound(params%GD, dim=1)
            b = ubound(params%GD, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%GD

        case("HR")
            a = lbound(params%HR, dim=1)
            b = ubound(params%HR, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%HR

        case("GR")
            a = lbound(params%GR, dim=1)
            b = ubound(params%GR, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%GR

        case ("--")
            return

        case default
            call abort(202302201, "Unknown wavelet filter (z-dir), must be one of HD GD HR GR")

        end select

        if ((abs(a) > g).or.(b>g)) then
            write(*,*) a, b, "but g=", g
            call abort(202302211, "For applying the filter, not enough ghost nodes")
        endif

        u_tmp = u_filtered
        u_filtered = u_tmp
        u_filtered(:, :, g+1:Bs(3)+g, :) = 0.0_rk
        do iz = g+1, Bs(3)+g
            do shift = a, b
                u_filtered(:, :, iz, :) = u_filtered(:, :, iz, :) + u_tmp(:, :, iz+shift, :)*coefs_filter(shift)
            enddo
        enddo

    end subroutine

    !-------------------------------------------------------------------------------
    ! applies one of params% HD GD HR GR filters in one direction
    ! ToDo: for HD and GD filter, we might benefit from considering restriction, for HR from skipping points
    subroutine blockFilterCustom1_vct( params, u, u_filtered, filter, direction )
        implicit none
        type (type_params), intent(in) :: params
        real(kind=rk), dimension(:,:,:,:), intent(in) :: u
        real(kind=rk), dimension(:,:,:,:), intent(inout) :: u_filtered
        character(len=*), intent(in) :: filter, direction
        integer(kind=ik) :: a, b
        real(kind=rk), allocatable :: coefs_filter(:)
        integer(kind=ik) :: ix, iy, iz, nx, ny, nz, nc, shift, g, Bs(1:3)
        real(kind=rk), allocatable, save :: u_tmp(:,:,:,:)

        nx = size(u, 1)
        ny = size(u, 2)
        nz = size(u, 3)
        nc = size(u, 4)
        g  = params%g
        Bs = params%Bs

        if (.not. allocated(u_tmp)) allocate( u_tmp(1:nx, 1:ny, 1:nz, 1:nc) )
        u_tmp = u

        select case(filter)
        case("HD")
            if (.not. allocated(params%HD)) call abort(202302203, "Wavelet setup not done yet?! HD")
            a = lbound(params%HD, dim=1)
            b = ubound(params%HD, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%HD

        case("GD")
            if (.not. allocated(params%GD)) call abort(202302203, "Wavelet setup not done yet?! GD")
            a = lbound(params%GD, dim=1)
            b = ubound(params%GD, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%GD

        case("HR")
            if (.not. allocated(params%HR)) call abort(202302203, "Wavelet setup not done yet?! HR")
            a = lbound(params%HR, dim=1)
            b = ubound(params%HR, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%HR

        case("GR")
            if (.not. allocated(params%GR)) call abort(202302203, "Wavelet setup not done yet?! GR")
            a = lbound(params%GR, dim=1)
            b = ubound(params%GR, dim=1)
            if (allocated(coefs_filter)) deallocate(coefs_filter)
            allocate( coefs_filter(a:b) )
            coefs_filter = params%GR

        case default
            call abort(202302201, "Unknown wavelet filter, must be one of HD GD HR GR")

        end select

        if ((abs(a) > g).or.(b>g)) then
            write(*,*) a, b, "but g=", g
            call abort(202302211, "For applying the filter, not enough ghost nodes")
        endif

        select case(direction)
        case("x")
            u_filtered = u_tmp
            u_filtered(g+1:Bs(1)+g, :, :, :) = 0.0_rk
            do ix = g+1, Bs(1)+g
                do shift = a, b
                    u_filtered(ix, :, :, :) = u_filtered(ix, :, :, :) + u_tmp(ix+shift, :, :, :)*coefs_filter(shift)
                enddo
            enddo

        case("y")
            u_filtered = u_tmp
            u_filtered(:, g+1:Bs(2)+g, :, :) = 0.0_rk
            do iy = g+1, Bs(2)+g
                do shift = a, b
                    u_filtered(:, iy, :, :) = u_filtered(:, iy, :, :) + u_tmp(:, iy+shift, :, :)*coefs_filter(shift)
                enddo
            enddo

        case("z")
            u_filtered = u_tmp
            u_filtered(:, :, g+1:Bs(3)+g, :) = 0.0_rk
            do iz = g+1, Bs(3)+g
                do shift = a, b
                    u_filtered(:, :, iz, :) = u_filtered(:, :, iz, :) + u_tmp(:, :, iz+shift, :)*coefs_filter(shift)
                enddo
            enddo

        case default
            call abort(202302201, "Unknown direction [xyz]. You should go to bed more early, its better.")

        end select


    end subroutine

    !-------------------------------------------------------------------------------

    ! computes a one-level wavelet decomposition of a block.
    ! The computation of coefficients is possible on the entire block, but
    ! without sync'ing the reconstruction is not possible on the entire data.
    ! We have to sync wavelet-transformed blocks.
    ! Data are stored in Spaghetti-order (not Mallat-Order)
    subroutine waveletDecomposition_block(params, u)
        implicit none
        type (type_params), intent(in) :: params
        real(kind=rk), dimension(:,:,:,:), intent(inout) :: u

        real(kind=rk), allocatable, dimension(:,:,:,:), save :: sc, wc, test, ucopy
        integer(kind=ik) :: nx, ny, nz, nc, g, Bs(1:3), ii
        ! integer(kind=ik) :: ag, bg, ah, bh, ix, iy, iz, ic, shift
        ! real(kind=rk) :: ug, uh

        call WaveDecomposition_dim1( params, u )
    end subroutine

    !-------------------------------------------------------------------------------

    subroutine waveletReconstruction_block(params, u, SC_reconstruct, WC_reconstruct)
        implicit none
        type (type_params), intent(in) :: params
        real(kind=rk), dimension(:,:,:,:), intent(inout) :: u
        logical, optional, intent(in) :: SC_reconstruct, WC_reconstruct  !< En- or Disable one of the reconstructions, defaults to true

        call WaveReconstruction_dim1( params, u, SC_reconstruct, WC_reconstruct)
    end subroutine

    !-------------------------------------------------------------------------------

    ! ensures that 3/4 of the numbers are zero - required for reconstruction
    ! note when copying Spaghetti to Mallat, this is automatically done, but
    ! when manipulating coefficients, it may happen that we set nonzero values
    ! where a zero should be.
    subroutine setRequiredZerosWCSC_block(params, u)
        implicit none
        type (type_params), intent(in) :: params
        real(kind=rk), dimension(:,:,:,:), intent(inout) :: u

        integer(kind=ik) :: nx, ny

        nx = size(u, 1)
        ny = size(u, 2)

        if (modulo(params%g, 2) == 0) then
            ! even ghost nodes: 2nd coefficient will be set to zero
            !
            ! x 0 x 0 x 0
            ! 0 0 0 0 0 0
            ! x 0 x 0 x 0
            !
            u(2:nx:2, 2:ny:2, :, :) = 0.0_rk ! diagonal
            u(2:nx:2, 1:ny:2, :, :) = 0.0_rk ! to the right (+x)
            u(1:nx:2, 2:ny:2, :, :) = 0.0_rk ! to the top (+y)
        else
            ! odd ghost nodes: 1st coefficient will be set to zero
            !
            ! 0 0 0 0 0
            ! 0 x 0 x 0
            ! 0 0 0 0 0
            ! 0 x 0 x 0
            !
            u(1:nx:2, 1:ny:2, :, :) = 0.0_rk ! diagonal
            u(1:nx:2, 2:ny:2, :, :) = 0.0_rk ! +y
            u(2:nx:2, 1:ny:2, :, :) = 0.0_rk ! +x
        endif


    end subroutine



    ! manipulate wavelet coefficients in a neighborhood direction: applied
    ! if we find a coarser neighbor in this direction (coarse extension)
    subroutine coarseExtensionManipulateWC_block(params, wc, neighborhood, Nwcl_optional, Nwcr_optional, set_garbage)
        implicit none

        type (type_params), intent(in) :: params
        real(kind=rk), dimension(:,:,:,:), intent(inout) :: wc   !< input/output in WD spaghetti format
        integer(kind=ik), intent(in) :: neighborhood                 !< Which neighborhood to apply manipulation
        integer(kind=ik), intent(in), optional :: Nwcl_optional, Nwcr_optional
        !> If set to true, we set a really large value into the patch so that the simulation can crash on purpose
        !> This is used to make clear that those values are not correct and should not be used
        logical, intent(in), optional :: set_garbage

        integer(kind=ik) :: Nwcl, Nwcr, i_neighborhood
        integer(kind=ik) :: nx, ny, nz, nc, g, Bs(1:3), io(1:3), d, idx(2,3), i_dim, i_set
        logical :: setGarbage
        real(kind = rk) :: setNumber

        nx = size(wc, 1)
        ny = size(wc, 2)
        nz = size(wc, 3)
        nc = size(wc, 4)
        g = params%g
        Bs = params%bs
        Nwcl    = params%Nwcl
        Nwcr    = params%Nwcr
        d       = 2_ik ** params%dim

        ! sometimes we just need to delete the WC in the ghost nodes layer
        ! in which case we set Nwcl=Nwcr=g
        if (present(Nwcl_optional)) Nwcl = Nwcl_optional
        if (present(Nwcr_optional)) Nwcr = Nwcr_optional

        ! sometimes we are not allowed to access the WC (from finer level neighbors) so lets ensure the simulation crashes in that case
        setGarbage = .false.
        if (present(set_garbage)) setGarbage = set_garbage

        do i_set = 1,2
            idx(:, :) = 1

            ! 1: set inside of patch
            if (i_set == 1) then
                call get_indices_of_modify_patch(params%g, params%dim, neighborhood, idx, (/ nx, ny, nz/), (/Nwcl, Nwcl, Nwcl/), (/Nwcr, Nwcr, Nwcr/), &
                    g_p=(/ g, g, g/), g_m=(/ g, g, g/), lvl_diff=+1)
            ! 2: set ghost patch
            else
                call get_indices_of_ghost_patch(params%Bs, params%g, params%dim, neighborhood, idx, params%g, params%g, lvl_diff=+1)
            endif
            if (params%dim == 2) idx(:,3) = 1  ! make sure that in 2D we leave the third dimension undisturbed

            ! we need to know if the first point is a SC or WC for the patch and adapt the indices accordingly
            !     1 2 3 4 5 6 7 8 9 A B C
            !     G G G S W S W S W G G G
            !                 I I I
            ! 1-C - index numbering in hex format, G - ghost point, S - SC, W - WC, I - point of patch to be checked
            ! for g=odd, the SC are on even numbers; for g=even, the SC are on odd numbers
            ! depending on BS this can change as well
            io = 0
            io(1:params%dim) = modulo(g + idx(1, 1:params%dim)+1, 2)

            ! set really low number for ghost patch if we shoudln't access it
            setNumber = 0.0
            if (setGarbage .and. i_set == 2) setNumber = -9e200_rk

            ! set values, we have to skip the SC, so we delete first WX, then WY, WXY, then WZ, WXZ, WYZ, WXYZ
            wc(idx(1,1)+1-io(1):idx(2,1):2, idx(1,2)  +io(2):idx(2,2):2, idx(1,3)  +io(3):idx(2,3):2, 1:nc) = setNumber
            wc(idx(1,1)        :idx(2,1)  , idx(1,2)+1-io(2):idx(2,2):2, idx(1,3)  +io(3):idx(2,3):2, 1:nc) = setNumber
            if (params%dim == 3) then
                wc(idx(1,1)    :idx(2,1)  , idx(1,2)        :idx(2,2)  , idx(1,3)+1-io(3):idx(2,3):2, 1:nc) = setNumber
            endif

            ! wc(idx(1,1)+1-io(1):idx(2,1):2, idx(1,2)  +io(2):idx(2,2):2, idx(1,3)+io(3):idx(2,3):2, 1:nc) = setNumber  ! WX
            ! wc(idx(1,1)  +io(1):idx(2,1):2, idx(1,2)+1-io(2):idx(2,2):2, idx(1,3)+io(3):idx(2,3):2, 1:nc) = setNumber  ! WY
            ! wc(idx(1,1)+1-io(1):idx(2,1):2, idx(1,2)+1-io(2):idx(2,2):2, idx(1,3)+io(3):idx(2,3):2, 1:nc) = setNumber  ! WXY

            ! if (params%dim == 3) then
            !     wc(idx(1,1)  +io(1):idx(2,1):2, idx(1,2)  +io(2):idx(2,2):2, idx(1,3)+1-io(3):idx(2,3):2, 1:nc) = setNumber  ! WZ
            !     wc(idx(1,1)+1-io(1):idx(2,1):2, idx(1,2)  +io(2):idx(2,2):2, idx(1,3)+1-io(3):idx(2,3):2, 1:nc) = setNumber  ! WXZ
            !     wc(idx(1,1)  +io(1):idx(2,1):2, idx(1,2)+1-io(2):idx(2,2):2, idx(1,3)+1-io(3):idx(2,3):2, 1:nc) = setNumber  ! WYZ
            !     wc(idx(1,1)+1-io(1):idx(2,1):2, idx(1,2)+1-io(2):idx(2,2):2, idx(1,3)+1-io(3):idx(2,3):2, 1:nc) = setNumber  ! WXYZ
            ! endif
        enddo

        ! wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc, 2:d) = setNumber
    end subroutine



    subroutine coarseExtensionManipulateSC_block(params, wc, u_copy, neighborhood, ijk)
        implicit none

        type (type_params), intent(in) :: params
        real(kind=rk), dimension(:,:,:,:), intent(inout) :: wc   !< input/output in WD spaghetti format
        real(kind=rk), dimension(:,:,:,:), intent(in) :: u_copy  !< original values from which we copy
        integer(kind=ik), intent(in)   :: neighborhood               !< Which neighborhood to apply manipulation
        integer(kind=ik), intent(in), optional   :: ijk(2,3)         !< ijk of patch that we only care about, if not given it is the full block

        integer(kind=ik) :: nx, ny, nz, nc, g, Bs(1:3), i_set, i_neighborhood, io(1:3)
        integer(kind=ik) :: Nscl, Nscr, Nreconl, Nreconr, idx(2,3)
        logical          :: skip_copy  ! sometimes this is called but actually nothing will be changed, in this case we skip it

        nx = size(wc, 1)
        ny = size(wc, 2)
        nz = size(wc, 3)
        nc = size(wc, 4)
        g = params%g
        Bs = params%bs
        Nscl    = params%Nscl
        Nscr    = params%Nscr
        Nreconl = params%Nreconl
        Nreconr = params%Nreconr
        skip_copy = .false.

        ! only set inside patch, as the other should be updated with syncing
        do i_set = 1,1
            idx(:, :) = 1

            ! 1: set inside of patch
            if (i_set == 1) then
                call get_indices_of_modify_patch(params%g, params%dim, neighborhood, idx, (/ nx, ny, nz/), (/Nscl, Nscl, Nscl/), (/Nscr, Nscr, Nscr/), &
                    g_p=(/ g, g, g/), g_m=(/ g, g, g/), lvl_diff=+1)
            ! 2: set ghost patch
            else
                call get_indices_of_ghost_patch(params%Bs, params%g, params%dim, neighborhood, idx, params%g, params%g, lvl_diff=+1)
            endif

            ! we need to know if the first point is a SC or WC for the patch and adapt the indices accordingly
            !     1 2 3 4 5 6 7 8 9 A B C
            !     G G G S W S W S W G G G
            !                 I I I
            ! 1-C - index numbering in hex format, G - ghost point, S - SC, W - WC, I - point of patch to be checked
            ! for g=odd, the SC are on even numbers; for g=even, the SC are on odd numbers
            ! depending on BS this can change as well
            io = 0
            io(1:params%dim) = modulo(g + idx(1, 1:params%dim) + 1, 2)

            ! sometimes we do not want to look at the full block, then we can skip the copy where the indices do not matter
            if (present(ijk)) then
                ! check in each dimension if the indices lie outside the affected range
                ! copy patch has to end before beginning of affected region or begin before end of affected region
                if (idx(2,1) < ijk(1,1) .or. idx(1,1) > ijk(2,1)) skip_copy = .true.
                if (idx(2,2) < ijk(1,2) .or. idx(1,2) > ijk(2,2)) skip_copy = .true.
                if (idx(2,3) < ijk(1,3) .or. idx(1,3) > ijk(2,3) .and. params%dim == 3) skip_copy = .true.
            endif

            if (.not. skip_copy) then
                wc(idx(1,1)+io(1):idx(2,1):2, idx(1,2)+io(2):idx(2,2):2, idx(1,3)+io(3):idx(2,3):2, 1:nc) = &
            u_copy(idx(1,1)+io(1):idx(2,1):2, idx(1,2)+io(2):idx(2,2):2, idx(1,3)+io(3):idx(2,3):2, 1:nc)
            endif

        enddo

    end subroutine



    subroutine setup_wavelet(params, g_wavelet, g_RHS, verbose)
        implicit none
        type (type_params), intent(inout) :: params
        ! if called with g_wavelet, we return the number of ghost nodes
        ! required for the wavelet filters (used in postprocessing routines
        ! to decide which G is used.)
        integer(kind=ik), intent(out), optional :: g_wavelet, g_RHS
        logical, intent(in), optional :: verbose
        logical :: verbose1
        integer(kind=ik) :: i, j, g_min, a, block_min, diff_L, diff_R

        integer(kind=ik) :: CDFX, CDFY, FD1_size, FD2_size, FD_max_size
        real(kind=rk), allocatable :: h_ntilde(:)

        character(len=80) :: debug_file_name
        integer :: idx(2,3)

        if (allocated(params%GR)) deallocate(params%HD)
        if (allocated(params%GD)) deallocate(params%GD)
        if (allocated(params%HR)) deallocate(params%HR)
        if (allocated(params%GR)) deallocate(params%GR)
        if (allocated(params%MGR)) deallocate(params%MGR)

        verbose1 = .true.
        if (present(verbose)) verbose1 = verbose

        ! for non-lifted wavelets: (Donoho wavelets)
        ! Not really required.....doesn't hurt.
        params%Nscl = 0
        params%Nscr = 0
        params%Nwcl = 0
        params%Nwcr = 0
        params%Nreconl = 0
        params%Nreconr = 0

        ! put header first so output is easier to understand
        if (params%rank==0 .and. verbose1) then
            write(*, '("  ╭─╮      ╭─╮                           ╭─╮         ╭───╮           ╭─╮        ")')
            write(*, '("──╯ │ ╭────╯ │ ╭───   Wavelet-setup   ───╯ │ ╭───────╯   │   ╭───────╯ │ ╭──────")')
            write(*, '("    ╰─╯      ╰─╯                           ╰─╯           ╰───╯         ╰─╯      ")')
            write(*,'(2A)') "The wavelet is ", trim(adjustl(params%wavelet))
        endif

        ! the wavelet filter banks:
        ! HD - low pass decomposition filter, H_TILDE
        ! GD - high pass decomposition filter, G_TILDE - for CDF it is always HR with different sign for every second off-center value
        ! HR - low pass reconstruction filter, H
        ! GR - high pass reconstruction filter, G - for CDF it is always HD with different sign for every second off-center value

        select case(params%wavelet)
        case ("coiflet12")
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! coiflet - orthogonal, almost interpolating
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            ! Not fully tested -- do not use in production runs. (see below)
            ! Issue #1: the function refine_block uses the direct interpolation (polynomial), and not zero padding in wavelet space.
            !           consequently, the coarsen(refine(u)) unit test will fail.
            ! Issue #2: Coarse extension is to be clarified with coiflet -> copying of SC near interface may make less sense than for CDF ?
            allocate(params%HD(-4:7))
            allocate(params%GD(-6:5))  ! maybe needs to be shifted after reshifting wavelet filters to be symmetric

            allocate(params%HR(-7:4))
            allocate(params%GR(-5:6))  ! maybe needs to be shifted after reshifting wavelet filters to be symmetric

            ! copied from flusi coiflet (output)
            params%HD=(/1.638733646318000E-02, -4.146493678197000E-02, -6.737255472230000E-02, 3.861100668230900E-01, 8.127236354496100E-01, 4.170051844237800E-01,-7.648859907826000E-02,-5.943441864647000E-02, 2.368017194688000E-02, 5.611434819370000E-03,-1.823208870910000E-03,-7.205494453700000E-04/)
            params%GD=(/7.205494453700000E-04, -1.823208870910000E-03, -5.611434819370000E-03, 2.368017194688000E-02, 5.943441864647000E-02,-7.648859907826000E-02,-4.170051844237800E-01, 8.127236354496100E-01,-3.861100668230900E-01,-6.737255472230000E-02, 4.146493678197000E-02, 1.638733646318000E-02/)
            params%HR=(/-7.205494453700000E-04, -1.823208870910000E-03,  5.611434819370000E-03, 2.368017194688000E-02,-5.943441864647000E-02,-7.648859907826000E-02, 4.170051844237800E-01, 8.127236354496100E-01, 3.861100668230900E-01,-6.737255472230000E-02,-4.146493678197000E-02, 1.638733646318000E-02/)
            params%GR=(/1.638733646318000E-02,  4.146493678197000E-02, -6.737255472230000E-02,-3.861100668230900E-01, 8.127236354496100E-01,-4.170051844237800E-01,-7.648859907826000E-02, 5.943441864647000E-02, 2.368017194688000E-02,-5.611434819370000E-03,-1.823208870910000E-03, 7.205494453700000E-04/)

            ! As the coiflets are *not* strictly interpolating, only "almost-interpolating", the combination with 
            ! a fourth order interpolation in the ghost nodes is questionable. The "predictor" is used in the ghost nodes interpolation.
            params%order_predictor = "multiresolution_4th"
            CDFX = 4  ! This wavelet is quasi-interpolating of order 4

            do i = 1, 100
                write(*,*) "Warning: COIFLET12 is not yet fully tested and implemented, do not use for production runs!"
            enddo

        case default
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Cohen-Daubechies-Feauveau Wavelets (CDF) -- interpolating, biorthogonal
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            ! check if we have a CDF wavelet
            if (params%wavelet(1:3) /= "CDF") then
                call abort( 3006221, "Unkown bi-orthogonal wavelet specified. Set course for adventure! params%wavelet="//trim(adjustl(params%wavelet)) )
            endif

            ! Now read out the order
            read(params%wavelet(4:4), *) CDFX
            ! Thomas will probably hate me for paving the way for 10th and 12th order wavelets, but well...
            if (CDFX == 1) then
                read(params%wavelet(4:5), *) CDFX
            endif
            if (all(CDFX /= (/2,4,6,8,10,12/))) then
                call abort( 3006221, "Unkown bi-orthogonal wavelet specified. Set course for adventure! params%wavelet="//trim(adjustl(params%wavelet)) )
            endif
            if (CDFX < 10) then
                read(params%wavelet(5:5), *) CDFY
            else
                read(params%wavelet(6:6), *) CDFY
            endif
            if (CDFY == 1) then
                if (CDFX < 10) then
                    read(params%wavelet(5:6), *) CDFY
                else
                    read(params%wavelet(6:7), *) CDFY
                endif
            endif
            if (all(CDFY /= (/0,2,4,6,8,10,12/))) then
                call abort( 3006221, "Unkown bi-orthogonal wavelet specified. Set course for adventure! params%wavelet="//trim(adjustl(params%wavelet)) )
            endif
            ! Yes, Thomas, We could implement a CDF1212 with this, even though I doubt the filter convolution will be well representable with double precision
            ! But well, if someone wants to try it out, why not...

            if (params%rank==0) then
                write(*,'(A, I2, A, I2, A)') "Selected wavelet: CDF with interpolation order ", CDFX, " and ", CDFY, " vanishing moments."
            endif

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! The X in CDFXY specifies the H TILDE (HR) filter
            ! This will always be the same, taken from the interpolation stencils
            if (CDFX == 2) then
                allocate( params%HR(-1:1) )
                params%HR(:) = Stencil_Int_2nd(:)
            elseif (CDFX == 4) then
                allocate( params%HR(-3:3) )
                params%HR(:) = Stencil_Int_4th(:)
            elseif (CDFX == 6) then
                allocate( params%HR(-5:5) )
                params%HR(:) = Stencil_Int_6th(:)
            elseif (CDFX == 8) then
                allocate( params%HR(-7:7) )
                params%HR(:) = Stencil_Int_8th(:)
            elseif (CDFX == 10) then
                allocate( params%HR(-9:9) )
                params%HR(:) = Stencil_Int_10th(:)
            elseif (CDFX == 12) then
                allocate( params%HR(-11:11) )
                params%HR(:) = Stencil_Int_12th(:)
            else
                call abort( 3006221, "Unkown bi-orthogonal wavelet specified. Set course for adventure! params%wavelet="//trim(adjustl(params%wavelet)) )
            endif

            ! The Y in CDFXY specifies the H (HD) filter
            ! For unlifted wavelets, this is 1,
            ! For lifted wavelets, it is a convolution of H TILDE with X and a H TILDE with Y, in case X <= Y
            if (CDFY == 0) then
                ! unlifted wavelets, H is just 1 or simple downpass
                allocate( params%HD(0:0) )
                params%HD = (/1.0_rk/)
                ! multigrid restriction - second order central average for lowpass filtering
                allocate( params%MGR(-1:1) )
                params%MGR = (/1.0_rk/4.0_rk, 1.0_rk/2.0_rk, 1.0_rk/4.0_rk/)
            elseif (CDFX >= CDFY) then
                ! lifted wavelets
                ! For lifted wavelets, in case X <= Y, it is a convolution of H TILDE with X and a H TILDE with Y
                ! For Y > X, we only have hardcoded stencils for CDF24, CDF26, CDF28, CDF46

                ! At first, we will initialize the second H TILDE filter and set the center coefficient to 0
                allocate(h_ntilde(-(CDFY-1):CDFY-1))
                if (CDFY == 2) then
                    h_ntilde = Stencil_Int_2nd(:)
                elseif (CDFY == 4) then
                    h_ntilde = Stencil_Int_4th(:)
                elseif (CDFY == 6) then
                    h_ntilde = Stencil_Int_6th(:)
                elseif (CDFY == 8) then
                    h_ntilde = Stencil_Int_8th(:)
                elseif (CDFY == 10) then
                    h_ntilde = Stencil_Int_10th(:)
                elseif (CDFY == 12) then
                    h_ntilde = Stencil_Int_12th(:)
                else
                    call abort( 3006221, "Unkown bi-orthogonal wavelet specified. Set course for adventure! params%wavelet="//trim(adjustl(params%wavelet)) )
                endif
                h_ntilde(0) = 0.0_rk

                ! now initialize H filter filter (which is HD)
                allocate( params%HD( lbound(params%HR, dim=1) + lbound(h_ntilde, dim=1) : ubound(params%HR, dim=1) + ubound(h_ntilde, dim=1) ) )

                ! convolve both H TILDE filters to get H and add delta function
                ! h_k = \delta_{k,0} + \frac12\sum_m -1^{m}\cdot \tilde{h}_{N}(m) \cdot \tilde{h}_{\tilde{N}}(m-k)
                do i = lbound(params%HD, dim=1), ubound(params%HD, dim=1)
                    ! Initialize with delta function, this is only 1 on center and 0 anywhere else
                    if (i == 0) then
                        params%HD(i) = 1.0_rk
                    else
                        params%HD(i) = 0.0_rk
                    endif

                    ! now add the convolved part with second H TILDE
                    do j = lbound(params%HR, dim=1), ubound(params%HR, dim=1)
                        ! see if this contribution is inside bounds of h_ntilde
                        if (i-j < lbound(h_ntilde, dim=1) .or. i-j > ubound(h_ntilde, dim=1)) cycle

                        ! add contribution, we have to alternate the sign and multiply by 1/2
                        params%HD(i) = params%HD(i) + (-1.0_rk)**j * params%HR(j) * h_ntilde(i - j) / 2.0_rk
                    end do
                end do
            ! Now follows some hardcoded stencils for Y > X
            ! These are lifted wavelets, but I have no clue how they are constructed
            elseif (params%wavelet == "CDF24") then
                allocate( params%HD(-4:4) )
                ! from Daubechies - Ten lectures on wavelets, Table 8.2
                params%HD = (/3.0_rk, -6.0_rk, -16.0_rk, 38.0_rk, 90.0_rk, 38.0_rk, -16.0_rk, -6.0_rk, 3.0_rk/) / 128.0_rk
            elseif (params%wavelet == "CDF26") then
                allocate( params%HD(-6:6) )
                ! from Daubechies - Ten lectures on wavelets, Table 8.2
                params%HD = (/-5.0_rk, 10.0_rk, 34.0_rk, -78.0_rk, -123.0_rk, 324.0_rk, 700.0_rk, 324.0_rk, -123.0_rk, -78.0_rk, 34.0_rk, 10.0_rk, -5.0_rk/) / 1024.0_rk
            elseif (params%wavelet == "CDF28") then
                allocate( params%HD(-8:8) )
                ! from Daubechies - Ten lectures on wavelets, Table 8.2
                params%HD = (/35.0_rk, -70.0_rk, -300.0_rk, 670.0_rk, 1228.0_rk, -3126.0_rk, -3796.0_rk, 10718.0_rk, 22050.0_rk, &
                    10718.0_rk, -3796.0_rk, -3126.0_rk, 1228.0_rk, 670.0_rk, -300.0_rk, -70.0_rk, 35.0_rk/) / 32768.0_rk
            elseif (params%wavelet == "CDF46") then
                ! Sweldens paper, "The Lifting Scheme: A Custom-Design
                ! Construction of Biorthogonal Wavelets" table 2 for N_tilde=6
                allocate( params%HD(-8:8) )
                params%HD = (/ 9.0_rk*2.0_rk**(-14.0_rk),          0.0_rk, &
                            -35.0_rk*2.0_rk**(-12.0_rk),   9.0_rk*2.0_rk**(-10.0_rk), &
                            189.0_rk*2.0_rk**(-12.0_rk), -59.0_rk*2.0_rk**(-10.0_rk), &
                            -477.0_rk*2.0_rk**(-12.0_rk), 153.0_rk*2.0_rk**( -9.0_rk), &
                            5379.0_rk*2.0_rk**(-13.0_rk), 153.0_rk*2.0_rk**( -9.0_rk), &
                            -477.0_rk*2.0_rk**(-12.0_rk), -59.0_rk*2.0_rk**(-10.0_rk), &
                            189.0_rk*2.0_rk**(-12.0_rk),   9.0_rk*2.0_rk**(-10.0_rk), &
                            -35.0_rk*2.0_rk**(-12.0_rk),          0.0_rk, &
                            9.0_rk*2.0_rk**(-14.0_rk) /)  
            else 
                call abort( 3006221, "For lifted CDF wavelets, first number must be smaller or equal to second number or you use CDF24, CDF26, CDF28 or CDF46. Set course for adventure! params%wavelet="//trim(adjustl(params%wavelet)) )
            endif
                
            ! for multigrid restriction, we use the same filter as HD if not an unlifted wavelet
            if (CDFY /= 0) then
                allocate( params%MGR( lbound(params%HD, dim=1) : ubound(params%HD, dim=1) ) )
                params%MGR(:) = params%HD(:)
            endif

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! G TILDE filter - HR filter with different sign for every second off-center value
            allocate( params%GD( lbound(params%HR, dim=1):ubound(params%HR, dim=1)) )
            do i = lbound(params%GD, dim=1), ubound(params%GD, dim=1)
                params%GD(i) = (-1.0_rk)**(i) * params%HR(i)
            enddo
                        
            ! G filter - HD filter with different sign for every second off-center value
            ! if (params%wavelet(5:5) == "0") then  ! for unlifted wavelets GR filter is set larger than HD filter
            !     allocate( params%GR(-1:1) )
            !     params%GR = (/ 0.0_rk, 1.0_rk, 0.0_rk /)
            ! else
                allocate( params%GR( lbound(params%HD, dim=1):ubound(params%HD, dim=1)) )
                do i = lbound(params%GR, dim=1), ubound(params%GR, dim=1)
                    params%GR(i) = (-1.0_rk)**(i) * params%HD(i)
                enddo
            ! endif

            ! Unlifted or lifted - every CDFX0er wavelet is considered unlifted and the rest lifted
            params%isLiftedWavelet = CDFY /= 0

            ! order predictor is decided by X in CDFXY (it is coinciding with the HR filter actually)
            if (CDFX == 2) then
                params%order_predictor = "multiresolution_2nd"
            elseif (CDFX == 4) then
                params%order_predictor = "multiresolution_4th"
            elseif (CDFX == 6) then
                params%order_predictor = "multiresolution_6th"
            elseif (CDFX == 8) then
                params%order_predictor = "multiresolution_8th"
            elseif (CDFX == 10) then
                params%order_predictor = "multiresolution_10th"
            elseif (CDFX == 12) then
                params%order_predictor = "multiresolution_12th"
            endif
        end select

        ! retrieve FD operator stencil sizes as needed for some setups, use h_ntilde as dummy variable and i,j as temporary integers
        if (params%order_discretization /= "not-initialized") then
            call setup_FD1_left_stencil(params%order_discretization, h_ntilde, i,j)
            FD1_size = max(i,j)
            call setup_FD2_stencil(params%order_discretization, h_ntilde, i,j)
            FD2_size = max(i,j)
            FD_max_size = max(FD1_size, FD2_size)
        else
            FD_max_size = 0
        endif

        ! minimum number of ghost nodes required for this wavelet
        ! value is determined by largest filter bank, let's write all out for clarity
        g_min = 0
        g_min = max(g_min, abs(lbound(params%HD, dim=1)))
        g_min = max(g_min, abs(ubound(params%HD, dim=1)))
        g_min = max(g_min, abs(lbound(params%GD, dim=1)))
        g_min = max(g_min, abs(ubound(params%GD, dim=1)))
        g_min = max(g_min, abs(lbound(params%HR, dim=1)))
        g_min = max(g_min, abs(ubound(params%HR, dim=1)))
        g_min = max(g_min, abs(lbound(params%GR, dim=1)))
        g_min = max(g_min, abs(ubound(params%GR, dim=1)))

        ! g_RHS is usually decided by FD order
        ! g_RHS is also dependent on X in CDFXY of the wavelet due to how we synch each stage:
        ! the third stage (prediction) needs points from the boundary to correctly interpolate values
        if (present(g_RHS)) then
            g_RHS = CDFX / 2
            if (g_RHS < FD_max_size) then
                g_RHS = FD_max_size
            endif
        endif

        ! compute coarse extension parameters. (see inkscape drawing 1d-ghostnodes-coarseextension-V2.svg)
        !
        ! Nscl, Nscr:
        ! determined with (h_tilde filter) not reaching into ghost nodes.
        ! Remaining points are copied (orange arrow)
        !
        ! Pink zero wc: reconstruction of last ghost node with (g filter)
        ! must use all zero wc (YOU CAN SKIP THIS!)
        !
        ! Red zero WC: wc that is computed from copied SC is zero
        ! (g_tilde filter)
        !
        ! Reconstruction: any modified zero or copy (g filter & h filter)
        ! as the g filter is wider, it determines the length (but check!)
        ! -> it is always the g filter because it is wider AND we set more WC
        ! to zero
        !
        params%Nscl    = max(abs(lbound(params%HD, dim=1)) - 1, 0)
        params%Nwcl    = params%Nscl + abs(lbound(params%GD, dim=1))
        params%Nreconl = params%Nwcl + abs(lbound(params%GR, dim=1))

        params%Nscr    = ubound(params%HD, dim=1)
        params%Nwcr    = params%Nscr + ubound(params%GD, dim=1)
        params%Nreconr = params%Nwcr + ubound(params%GR, dim=1)

        ! ! for unlifted wavelets no SC are copied, however some WC do have to be wiped in case we need to reconstruct (CVS and image denoising)
        ! ! This is, because at a coarse-fine interface the filter to compute the WC stretches over the border, creating the same dilemma as for the CE
        ! ! normally, if no WC are altered we can simply recopy the values that we had
        ! if (.not. params%isLiftedWavelet) then
        !     params%Nwcl = params%Nwcl + 2
        !     params%Nwcr = params%Nwcr + 2
        !     params%Nreconl = params%Nreconl + 2
        !     params%Nreconr = params%Nreconr + 2
        ! endif

        ! for wavelets with regularity higher than the wavelets NWC needs to be increased, the reasoning for me is yet unclear
        ! this was investigated and found using the invertibility test
        if (params%isLiftedWavelet) then
            a = 0
            if (params%wavelet(4:5) == "24" .or. params%wavelet(4:5) == "46") a = 2
            if (params%wavelet(4:5) == "26") a = 4
            if (params%wavelet(4:5) == "28") a = 6
            params%Nwcl = params%Nwcl + a
            params%Nwcr = params%Nwcr + a
            params%Nreconl = params%Nreconl + a
            params%Nreconr = params%Nreconr + a
        endif

        !--------------------------------------------------------------------------------------------------------
        ! If NWC of CE is smaller than the size of the FD-stencils, this can create strong divergence peaks.
        ! In order to reduce this, we set the minimum Nwc to the size of the FD stencils.
        ! This mainly affects the unlifted wavelets (or you pair CDF22 with FD >= 6th order, you weirdo)
        i = 2*FD_max_size

        diff_L         = max(i - params%Nwcl, 0)
        params%Nwcl    = params%Nwcl + diff_L
        params%Nreconl = params%Nreconl + diff_L
        
        diff_R         = max(i - params%Nwcr, 0)
        params%Nwcr    = params%Nwcr + diff_R
        params%Nreconr = params%Nreconr + diff_R

        if (params%rank==0) then
            write(*,'(A)') "Wavelet setup is adjusted to discretization: "//trim(adjustl(params%order_discretization))
        endif
        
        
        if ((params%useCoarseExtension .or. params%useSecurityZone) .and. params%rank==0 ) then
           write(*, '(A,i3,1x,i3," L/R")') "Increased Nwc to consider FD-stencil size by ", diff_L, diff_R
        endif

        ! significant refinement without coarse extension can cause trouble, let's give a warning to the user (that no-one will probably read ever)
        if (params%refinement_indicator == 'significant' .and. .not. params%useCoarseExtension) then
           write(*, '(A)') 'WARNING: Significant refinement are prone to grid instabilities of our discrete operators. You should use the coarse extension in order to filter coarse-fine grid interfaces!'
        endif
        !--------------------------------------------------------------------------------------------------------

        ! we are debugging the patches for Coarse Extension SC Copy and WC Zero to a file so that we can check that they are correct
#ifdef DEV
        open(16,file="CE_bounds.dat",status='replace')
        write(16,'(3(A, i0), 2(A))') "% dim=", params%dim, ", Bs=", params%Bs(1), ", g=", params%g, ", CDF=", params%wavelet
        write(16,'(A, A13, 6(A15))') "% ", "Neighborhood" , "idx_1", "idx_2", "idy_1", "idy_2", "idz_1", "idz_2"
        do i = 1, 56*3
            call get_indices_of_modify_patch(params%g, params%dim, i, idx, (/ params%Bs(1)+2*params%g, params%Bs(2)+2*params%g, merge(1, params%Bs(3)+2*params%g, params%dim==2)/), (/params%Nscl, params%Nscl, params%Nscl/), (/params%Nscr, params%Nscr, params%Nscr/), &
                g_p=(/ params%g, params%g, params%g/), g_m=(/ params%g, params%g, params%g/), lvl_diff=+1)
            write(16,'(7(i15))') i, idx(:, :)
            call get_indices_of_modify_patch(params%g, params%dim, i, idx, (/ params%Bs(1)+2*params%g, params%Bs(2)+2*params%g, merge(1, params%Bs(3)+2*params%g, params%dim==2)/), (/params%Nwcl, params%Nwcl, params%Nwcl/), (/params%Nwcr, params%Nwcr, params%Nwcr/), &
                g_p=(/ params%g, params%g, params%g/), g_m=(/ params%g, params%g, params%g/), lvl_diff=+1)
            write(16,'(7(i15))') i, idx(:, :)
        enddo
        close(16)
#endif
        
        !---------------------------------------------------------------------------------------------------------


        if (present(g_wavelet)) then
            ! if we return the minimal value of ghosts, we assume that you are going to use it
            ! instead of ignoring it -> we can omit checking 
            g_wavelet = g_min
        else
            if (params%g < g_min) then
                write(*,'(A, A, i3, A, i3)') trim(adjustl(params%wavelet)), " g=", params%g, " < g_min=", g_min
                call abort(8888881, "The selected wavelet requires more ghost nodes.")
            endif
        endif

        if (params%rank==0 .and. verbose1) then
            if (params%useCoarseExtension) then
                write(*,'(A55, i4, i4)') "During coarse extension, we will copy SC (L,R):", params%Nscl, params%Nscr
                write(*,'(A55, i4, i4)') "During coarse extension, we will delete WC (L,R):", params%Nwcl, params%Nwcr
            endif

            ! For the leaf-first loop, we need 3*h as minimum blocksize, as we have an upwards dependency for the leaf-decomposition
            ! So, for lower BS we do level-wise loop (which performs worse) and if the BS is high enough, we do the more optimized leaf-first level-wise loop
            block_min = 0
            if (params%isLiftedWavelet .and. maxval(params%Bs(:)) /= 0) then
                block_min = 3* max(abs(lbound(params%HD, dim=1)), abs(ubound(params%HD, dim=1)))
                if (any(params%Bs(:params%dim) < block_min)) then
                    write(*, '(A, i3, A, i3, A)') 'Bs=', minval(params%Bs(1:params%dim)), " < 3*h=", block_min,", not using optimized wavelet decomposition algorithm"
                else
                    write(*, '(A, i3, A, i3, A)') 'Bs=', minval(params%Bs(1:params%dim)), " >= 3*h=", block_min,", using optimized wavelet decomposition algorithm"
                endif
            endif
            write(*,'(2A)') "The predictor is: ", trim(adjustl(params%order_predictor))
            write(*,'(A,"[",i3,":",i2,"]=",14(es12.4,1x))') "HD", lbound(params%HD, dim=1), ubound(params%HD, dim=1), params%HD
            write(*,'(A,"[",i3,":",i2,"]=",14(es12.4,1x))') "GD", lbound(params%GD, dim=1), ubound(params%GD, dim=1), params%GD
            write(*,'(A,"[",i3,":",i2,"]=",14(es12.4,1x))') "HR", lbound(params%HR, dim=1), ubound(params%HR, dim=1), params%HR
            write(*,'(A,"[",i3,":",i2,"]=",14(es12.4,1x))') "GR", lbound(params%GR, dim=1), ubound(params%GR, dim=1), params%GR
            write(*, '(20("╭─╮ "))')
            write(*, '(20("╯ ╰─"))')
            ! this code has been used to plot our wavelets in PYTHON.
            a = maxval( (/&
            abs(lbound(params%HD, dim=1)), &
            abs(lbound(params%GD, dim=1)), &
            abs(lbound(params%HR, dim=1)), &
            abs(lbound(params%GR, dim=1)), &
            abs(ubound(params%HD, dim=1)), &
            abs(ubound(params%GD, dim=1)), &
            abs(ubound(params%HR, dim=1)), &
            abs(ubound(params%GR, dim=1)) &
            /) )

            write(*,'(A)', advance='no') "HD=["
            do i = -a, +a
                if ((i < lbound(params%HD, dim=1) ).or.(i > ubound(params%HD, dim=1) )) then
                    write(*,'(es15.8,", ")', advance='no') 0.0_rk
                else
                    write(*,'(es15.8,", ")', advance='no') params%HD(i)
                endif
            enddo
            write(*,'(A)') "]"

            write(*,'(A)', advance='no') "GD=["
            do i = -a, +a
                if ((i < lbound(params%GD, dim=1) ).or.(i > ubound(params%GD, dim=1) )) then
                    write(*,'(es15.8,", ")', advance='no') 0.0_rk
                else
                    write(*,'(es15.8,", ")', advance='no') params%GD(i)
                endif
            enddo
            write(*,'(A)') "]"

            write(*,'(A)', advance='no') "HR=["
            do i = -a, +a
                if ((i < lbound(params%HR, dim=1) ).or.(i > ubound(params%HR, dim=1) )) then
                    write(*,'(es15.8,", ")', advance='no') 0.0_rk
                else
                    write(*,'(es15.8,", ")', advance='no') params%HR(i)
                endif
            enddo
            write(*,'(A)') "]"

            write(*,'(A)', advance='no') "GR=["
            do i = -a, +a
                if ((i < lbound(params%GR, dim=1) ).or.(i > ubound(params%GR, dim=1) )) then
                    write(*,'(es15.8,", ")', advance='no') 0.0_rk
                else
                    write(*,'(es15.8,", ")', advance='no') params%GR(i)
                endif
            enddo
            write(*,'(A)') "]"

            write(*,'("sum(",A2,")=",es15.8, A, es10.3)') "HD", sum(params%HD), "  , Diff to 1.0=", sum(params%HD)-1.0_rk
            write(*,'("sum(",A2,")=",es15.8)')            "GD", sum(params%GD)
            write(*,'("sum(",A2,")=",es15.8, A, es10.3)') "HR", sum(params%HR), "  , Diff to 2.0=", sum(params%HR)-2.0_rk
            write(*,'("sum(",A2,")=",es15.8)')            "GR", sum(params%GR)

            write(*, '(20("╭─╮ "))')
            write(*, '(20("╯ ╰─"))')
        endif
    end subroutine



    ! 1D filter, pass skip_g=params%g to skip ghost points, pass do_restriction to filter only every second point
    subroutine filter1dim(params, u, u_filtered, filter, fl_l, fl_r, skip_g, sampling)
        implicit none
        type (type_params), intent(in) :: params                        !< good ol' params
        real(kind=rk), dimension(:), intent(in) :: u                   !< the signal to be filtered
        real(kind=rk), dimension(:), intent(inout) :: u_filtered       !< the resulting filtered signal
        integer, intent(in) :: fl_l, fl_r                               !< filter bound indices left and right
        real(kind=rk), dimension(fl_l:fl_r), intent(in) :: filter       !< the actual filter
        integer, intent(in) :: skip_g(1:2)                              !< 0 to filter everything, (/params%g,/params%g) to skip ghost points
        !> sampling rate, 1 for normal mode, 2 for restriction, no optional parameter as this is inside critical point-loop and I want to avoid if-clauses
        integer, intent(in) :: sampling                                        

        integer(kind=ik) :: N, i, j

        N = size(u)

        ! ! apply filter f to periodic signal u, i.e. convolute with filter

        ! classical pointwise
        ! u_filtered = 0.0_rk
        ! do i =  1+skip_g, N-skip_g, sampling
        !     do j = fl_l, fl_r
        !         u_filtered(i) = u_filtered(i) + u(i+j) * filter(j)
        !     enddo
        ! enddo

        ! vectorized filter with sum
        do i =  1+skip_g(1), N-skip_g(2), sampling
            u_filtered(i) = sum(u(i+fl_l:i+fl_r) * filter(:))
        enddo

        ! ! vectorize filter application
        ! u_filtered = 0.0_rk
        ! do j = fl_l, fl_r
        !     u_filtered(1+skip_g: N-skip_g: sampling) = u_filtered(1+skip_g: N-skip_g: sampling) + u(1+skip_g+j: N-skip_g+j: sampling) * filter(j)
        ! enddo
    end subroutine


    ! multi-dimensional FWT using 1D non-periodic filtering routines
    ! Inplace
    subroutine WaveDecomposition_dim1( params, u_wc )
        implicit none
        type (type_params), intent(in) :: params
        ! Input: data Output: FWT(DATA) in Spagghetti-ordering
        real(kind=rk), dimension(:, :, :, :), intent(inout) :: u_wc
        real(kind=rk), dimension(:), allocatable, save :: buffer1, buffer2
        integer(kind=ik) :: ix, iy, iz, ic, g, Bs(1:3), nx, ny, nz, nc, maxn
        nx = size(u_wc, 1)
        ny = size(u_wc, 2)
        nz = size(u_wc, 3)
        nc = size(u_wc, 4)
        g  = params%g
        Bs = params%Bs

        if (.not. allocated(params%HD)) call abort(1717229, "Wavelet setup not called?!")
        if (.not. allocated(params%GD)) call abort(1717231, "Wavelet setup not called?!")

        maxn = maxval((/ nx, ny, nz /))
        if (allocated(buffer1)) then
            if (size(buffer1, dim=1)<maxn) deallocate(buffer1)
        endif
        if (.not.allocated(buffer1)) allocate(buffer1(1:maxn))
        if (allocated(buffer2)) then
            if (size(buffer2, dim=1)<=maxn) deallocate(buffer2)
        endif
        if (.not.allocated(buffer2)) allocate(buffer2(1:maxn))

        ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        do ic = 1, nc
            do iy = 1, ny
                do iz = 1, nz
                    ! low-pass filter (scaling function)
                    call filter1dim(params, u_wc(:,iy,iz,ic), buffer1(1:nx), params%HD, lbound(params%HD,dim=1), ubound(params%HD,dim=1), skip_g=(/params%g,params%g/), sampling=2)

                    ! high-pass filter (these guys are the details)
                    call filter1dim(params, u_wc(:,iy,iz,ic), buffer2(1:nx), params%GD, lbound(params%GD,dim=1), ubound(params%GD,dim=1), skip_g=(/params%g+1,params%g/), sampling=2)

                    ! decimation by 2, sort into array in spaghetti form SC WC
                    u_wc((g+1):(Bs(1)+g):2,iy,iz,ic) = buffer1( (g+1):(Bs(1)+g):2 )
                    u_wc((g+2):(Bs(1)+g):2,iy,iz,ic) = buffer2( (g+2):(Bs(1)+g):2 )
                enddo
            enddo
        enddo

        ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        do ic = 1, nc
            do ix = (g+1), (Bs(1)+g)
                do iz = 1, nz
                    ! low-pass filter (scaling function)
                    call filter1dim(params, u_wc(ix,:,iz,ic), buffer1(1:ny), params%HD, lbound(params%HD,dim=1), ubound(params%HD,dim=1), skip_g=(/params%g,params%g/), sampling=2)

                    ! high-pass filter (these guys are the details)
                    call filter1dim(params, u_wc(ix,:,iz,ic), buffer2(1:ny), params%GD, lbound(params%GD,dim=1), ubound(params%GD,dim=1), skip_g=(/params%g+1,params%g/), sampling=2)

                    ! decimation by 2, sort into array in spaghetti form SC WC
                    u_wc(ix,(g+1):(Bs(2)+g):2,iz,ic) = buffer1( (g+1):(Bs(2)+g):2 )
                    u_wc(ix,(g+2):(Bs(2)+g):2,iz,ic) = buffer2( (g+2):(Bs(2)+g):2 )
                enddo
            enddo
        enddo

        ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (params%dim == 3) then
            do ic = 1, nc
                do ix = (g+1), (Bs(1)+g)
                    do iy = (g+1), (Bs(2)+g)
                        ! low-pass filter (scaling function)
                        call filter1dim(params, u_wc(ix,iy,:,ic), buffer1(1:nz), params%HD, lbound(params%HD,dim=1), ubound(params%HD,dim=1), skip_g=(/params%g,params%g/), sampling=2)

                        ! high-pass filter (these guys are the details)
                        call filter1dim(params, u_wc(ix,iy,:,ic), buffer2(1:nz), params%GD, lbound(params%GD,dim=1), ubound(params%GD,dim=1), skip_g=(/params%g+1,params%g/), sampling=2)

                        ! decimation by 2, sort into array in spaghetti form SC WC
                        u_wc(ix,iy,(g+1):(Bs(3)+g):2,ic) = buffer1( (g+1):(Bs(3)+g):2 )
                        u_wc(ix,iy,(g+2):(Bs(3)+g):2,ic) = buffer2( (g+2):(Bs(3)+g):2 )
                    enddo
                enddo
            enddo ! loop over components
        endif
    end subroutine


    !-----------------------------------------------------------------------------
    ! Reconstruction from low- and high pass filtered coefficients.
    ! Data is first upsampled, then filtered with the reconstruction filters.
    ! Note reconstruction filters are reverse of decomposition filters.
    !
    ! Input: hh gh hh gh hh gh hh gh
    !        hg gg hg gg hg gg hg gg
    !        hh gh hh gh hh gh hh gh
    !        hg gg hg gg hg gg hg gg
    !        hh gh hh gh hh gh hh gh
    !        hg gg hg gg hg gg hg gg
    ! Note: in input in spaghetti ordering is synced
    !-----------------------------------------------------------------------------
    subroutine WaveReconstruction_dim1( params, u_wc, SC_reconstruct, WC_reconstruct)
        implicit none
        type (type_params), intent(in) :: params
        !> Input is in spaghetti ordering (synchronized, ie with ghost nodes)
        real(kind=rk), dimension(:, :, :, :), intent(inout) :: u_wc
        real(kind=rk), dimension(:), allocatable, save :: buffer1, buffer2, buffer3
        integer(kind=ik) :: ix, iy, iz, ic, g, Bs(1:3), nx, ny, nz, nc, maxn
        integer(kind=ik) :: io  ! if g is odd the SCs start from the second point

        logical, optional, intent(in) :: SC_reconstruct, WC_reconstruct  !< En- or Disable one of the reconstructions, defaults to true
        logical :: SC_rec, WC_rec 
        SC_rec = .true.
        WC_rec = .true.
        if (present(SC_reconstruct)) SC_rec = SC_reconstruct
        if (present(WC_reconstruct)) WC_rec = WC_reconstruct

        nx = size(u_wc, 1)
        ny = size(u_wc, 2)
        nz = size(u_wc, 3)
        nc = size(u_wc, 4)
        g  = params%g
        Bs = params%Bs

        ! we need to shift if g is even as SCs start from second point
        io = 0
        if (modulo(g,2)==1) io = 1

        if (.not. allocated(params%HD)) call abort(1717229, "Wavelet setup not called?!")
        if (.not. allocated(params%GD)) call abort(1717231, "Wavelet setup not called?!")

        maxn = maxval((/ nx, ny, nz /))
        if (allocated(buffer1)) then
            if (size(buffer1, dim=1)<maxn) deallocate(buffer1)
        endif
        if (.not.allocated(buffer1)) allocate(buffer1(1:maxn))
        if (allocated(buffer2)) then
            if (size(buffer2, dim=1)<=maxn) deallocate(buffer2)
        endif
        if (.not.allocated(buffer2)) allocate(buffer2(1:maxn))
        if (allocated(buffer3)) then
            if (size(buffer3, dim=1)<=maxn) deallocate(buffer3)
        endif
        if (.not.allocated(buffer3)) allocate(buffer3(1:maxn+2))

        ! if SC should not be reconstructed, we still apply all filters (to have cross-effects) but wipe the SC
        if (.not. SC_rec) then
            if (nz==1) then
                u_wc(1+io:Bs(1)+2*g:2, 1+io:Bs(2)+2*g:2, 1, 1:nc) = 0.0_rk
            else
                u_wc(1+io:Bs(1)+2*g:2, 1+io:Bs(2)+2*g:2, 1+io:Bs(3)+2*g:2, 1:nc) = 0.0_rk
            endif
        endif

        ! If only the SC should be reconstructed, then GR is skipped all-together

        ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (WC_rec) then
            do ic = 1, nc; do iy = 1, ny; do iz = 1, nz
                ! fill upsampling buffer for low-pass filter: every second point
                ! apply low-pass filter to upsampled signal
                buffer3 = 0.0_rk
                buffer3(1+io:nx:2) = u_wc(1+io:Bs(1)+2*g:2, iy, iz, ic) ! SC
                call filter1dim(params, buffer3(1:nx), buffer1(1:nx), params%HR, lbound(params%HR,dim=1), ubound(params%HR,dim=1), skip_g=(/params%g,params%g/), sampling=1)

                ! fill upsampling buffer for high-pass filter: every second point
                ! GR locations shifted by one point
                buffer3 = 0.0_rk
                buffer3(2-io:nx:2) = u_wc(2-io:Bs(1)+2*g:2, iy, iz, ic) ! WC
                call filter1dim(params, buffer3(1:nx), buffer2(1:nx), params%GR, lbound(params%GR,dim=1), ubound(params%GR,dim=1), skip_g=(/params%g,params%g/), sampling=1)

                u_wc(:, iy, iz, ic) = buffer1(1:nx) + buffer2(1:nx)
            enddo; enddo; enddo
        elseif (SC_rec) then
            do ic = 1, nc; do iy = 1, ny; do iz = 1, nz
                ! fill upsampling buffer for low-pass filter: every second point
                ! apply low-pass filter to upsampled signal
                buffer3 = 0.0_rk
                buffer3(1+io:nx:2) = u_wc(1+io:Bs(1)+2*g:2, iy, iz, ic) ! SC
                call filter1dim(params, buffer3(1:nx), u_wc(1:nx, iy, iz, ic), params%HR, lbound(params%HR,dim=1), ubound(params%HR,dim=1), skip_g=(/params%g,params%g/), sampling=1)
            enddo; enddo; enddo
        endif

        ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! ignore ghost points for dimensions that were already treated
        if (WC_rec) then
            do ic = 1, nc; do ix = (g+1), (Bs(1)+g); do iz = 1, nz
                ! fill upsampling buffer for low-pass filter: every second point
                buffer3 = 0.0_rk
                buffer3(1+io:nx:2) = u_wc(ix, 1+io:Bs(2)+2*g:2, iz, ic) ! SC
                ! buffer3(1:ny:2) = u_wc(ix, 1:n(2), iz, ic)
                call filter1dim(params, buffer3(1:ny), buffer1(1:ny), params%HR, lbound(params%HR,dim=1), ubound(params%HR,dim=1), skip_g=(/params%g,params%g/), sampling=1)

                ! fill upsampling buffer for high-pass filter: every second point
                ! shifted by one point to account for lbound(GR) possibly being larger than g
                buffer3 = 0.0_rk
                buffer3(2-io:nx:2) = u_wc(ix, 2-io:Bs(2)+2*g:2, iz, ic) ! WC
                ! buffer3(1:ny:2) = u_wc(ix, n(2)+1:2*n(2), iz, ic)
                call filter1dim(params, buffer3(1:ny), buffer2(1:ny), params%GR, lbound(params%GR,dim=1), ubound(params%GR,dim=1), skip_g=(/params%g,params%g/), sampling=1)

                u_wc(ix, :, iz, ic) = buffer1(1:ny) + buffer2(1:ny)
            enddo; enddo; enddo
        elseif (SC_rec) then
            do ic = 1, nc; do ix = (g+1), (Bs(1)+g); do iz = 1, nz
                ! fill upsampling buffer for low-pass filter: every second point
                buffer3 = 0.0_rk
                buffer3(1+io:nx:2) = u_wc(ix, 1+io:Bs(2)+2*g:2, iz, ic) ! SC
                ! buffer3(1:ny:2) = u_wc(ix, 1:n(2), iz, ic)
                call filter1dim(params, buffer3(1:ny), u_wc(ix, 1:ny, iz, ic), params%HR, lbound(params%HR,dim=1), ubound(params%HR,dim=1), skip_g=(/params%g,params%g/), sampling=1)
            enddo; enddo; enddo
        endif

        !!!!!!!!!!!!!!!!!
        if (nz==1) return
        !!!!!!!!!!!!!!!!!

        ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! ignore ghost points for dimensions that were already treated
        if (WC_rec) then
            do ic = 1, nc; do ix = (g+1), (Bs(1)+g); do iy = (g+1), (Bs(2)+g)
                ! fill upsampling buffer for low-pass filter: every second point
                buffer3 = 0.0_rk
                buffer3(1+io:nx:2) = u_wc(ix, iy, 1+io:Bs(3)+2*g:2, ic) ! SC
                call filter1dim(params, buffer3(1:nz), buffer1(1:nz), params%HR, lbound(params%HR,dim=1), ubound(params%HR,dim=1), skip_g=(/params%g,params%g/), sampling=1)

                ! fill upsampling buffer for high-pass filter: every second point
                ! shifted by one point to account for lbound(GR) possibly being larger than g
                buffer3 = 0.0_rk
                buffer3(2-io:nx:2) = u_wc(ix, iy, 2-io:Bs(3)+2*g:2, ic) ! WC
                call filter1dim(params, buffer3(1:nz), buffer2(1:nz), params%GR, lbound(params%GR,dim=1), ubound(params%GR,dim=1), skip_g=(/params%g,params%g/), sampling=1)

                u_wc(ix, iy, :, ic) = buffer1(1:nz) + buffer2(1:nz)
            enddo; enddo; enddo
        elseif (SC_rec) then
            do ic = 1, nc; do ix = (g+1), (Bs(1)+g); do iy = (g+1), (Bs(2)+g)
                ! fill upsampling buffer for low-pass filter: every second point
                buffer3 = 0.0_rk
                buffer3(1+io:nx:2) = u_wc(ix, iy, 1+io:Bs(3)+2*g:2, ic) ! SC
                call filter1dim(params, buffer3(1:nz), u_wc(ix, iy, 1:nz, ic), params%HR, lbound(params%HR,dim=1), ubound(params%HR,dim=1), skip_g=(/params%g,params%g/), sampling=1)
            enddo; enddo; enddo
        endif

    end subroutine


    !> \brief Renorm the wavelet coefficients of a spaghetti-decomposed block from L_infty to L_infty, L_1, L_2 or dotH_1 norm
    subroutine wavelet_renorm_block(params, val_block, val_renormed, level_block, level_ref, indices, verbose_check)
        implicit none

        !> user defined parameter structure
        type (type_params), intent(in)         :: params
        !> heavy data for one block (hence 4D) expected in spaghetti-decomposed ordering
        real(kind=rk), intent(inout)           :: val_block(:, :, :, :)
        !> renormed data with SC set to 0
        real(kind=rk), intent(inout)           :: val_renormed(:, :, :, :)
        !> If we use L2 or H1 normalization, the threshold eps is level-dependent, hence
        !! we pass the level of this block to this routine
        integer(kind=ik), intent(in)           :: level_block
        !> If we use L2 or H1 normalization, the threshold eps is level-dependent, hence
        !! we pass a reference level to this routine, it defaults to params%Jmax
        integer(kind=ik), intent(in), optional :: level_ref
        !> Indices of patch if not the whole interior block should be thresholded, used for securityZone
        integer(kind=ik), intent(in), optional :: indices(1:2, 1:3)
        logical, intent(in), optional          :: verbose_check  !< No matter the value, if this is present we debug

        integer(kind=ik)                       :: idx(2,3), nc, g, Jref, i_dim, Bs(1:3), io(1:3)

        nc     = size(val_block, 4)
        Bs     = params%Bs
        g      = params%g
        Jref   = 0
        if (present(level_ref)) Jref = level_ref

        ! set the indices we want to threshold
        idx(:, :) = 1
        if (present(indices)) then
            idx(:, :) = indices(:, :)
        else  ! full interior block
            idx(1, 1) = g+1
            idx(2, 1) = Bs(1)+g
            idx(1, 2) = g+1
            idx(2, 2) = Bs(2)+g
            if (params%dim == 3) then
                idx(1, 3) = g+1
                idx(2, 3) = Bs(3)+g
            endif
        endif

        ! we need to know if the first point is a SC or WC for the patch and adapt the indices accordingly
        !     1 2 3 4 5 6 7 8 9 A B C
        !     G G G S W S W S W G G G
        !                 I I I
        ! 1-C - index numbering in hex format, G - ghost point, S - SC, W - WC, I - point of patch to be checked
        ! for g=odd, the SC are on even numbers; for g=even, the SC are on odd numbers
        ! depending on BS this can change as well
        io = 0
        io(1:params%dim) = modulo(g + idx(1, 1:params%dim) + 1, 2)

        ! copy only part we need
        val_renormed(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) = &
           val_block(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc)

        ! set sc to zero to more easily compute the maxval, use offset if first index is not a SC
        val_renormed(idx(1,1)+io(1):idx(2,1):2, idx(1,2)+io(2):idx(2,2):2, idx(1,3)+io(3):idx(2,3):2, 1:nc) = 0.0_rk

        ! We renorm by multiplying all values by level shifts+1 (SC) and then on the level change all SC-factors to WC-factors (apply WC-factor^2)
        ! -1 from the fact, that wavelet decomposed values on a block of level J will be decomposed values on level J-1
        select case(params%eps_norm)
        case ("Linfty")
            ! Wavelets are directly infinity normed, so we do nothing
        case ("L1")
            ! Wavelets get factor 2^-j per SC and 2^j per WC
            ! apply level shift
            val_renormed(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) = \
            val_renormed(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) * 2.0_rk**(dble((Jref-level_block-1)*params%dim))
            ! change SC to WC factors on this level - factor^2
            val_renormed(idx(1,1)+io(1):idx(2,1):2, idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) = \
            val_renormed(idx(1,1)+io(1):idx(2,1):2, idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) / 4.0_rk
            val_renormed(idx(1,1):idx(2,1), idx(1,2)+io(2):idx(2,2):2, idx(1,3):idx(2,3), 1:nc) = \
            val_renormed(idx(1,1):idx(2,1), idx(1,2)+io(2):idx(2,2):2, idx(1,3):idx(2,3), 1:nc) / 4.0_rk
            if (params%dim == 3) then
                val_renormed(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3)+io(3):idx(2,3):2, 1:nc) = \
                val_renormed(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3)+io(3):idx(2,3):2, 1:nc) / 4.0_rk
            endif
        case ("L2")
            ! Wavelets get factor 2^-j/2 per SC and 2^j/2 per WC
            ! apply level shift
            val_renormed(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) = \
            val_renormed(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) * 2.0_rk**(dble((Jref-level_block-1)*params%dim)/2.0_rk)
            ! change SC to WC factors on this level - factor^2
            val_renormed(idx(1,1)+io(1):idx(2,1):2, idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) = \
            val_renormed(idx(1,1)+io(1):idx(2,1):2, idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) / 2.0_rk
            val_renormed(idx(1,1):idx(2,1), idx(1,2)+io(2):idx(2,2):2, idx(1,3):idx(2,3), 1:nc) = \
            val_renormed(idx(1,1):idx(2,1), idx(1,2)+io(2):idx(2,2):2, idx(1,3):idx(2,3), 1:nc) / 2.0_rk
            if (params%dim == 3) then
                val_renormed(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3)+io(3):idx(2,3):2, 1:nc) = \
                val_renormed(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3)+io(3):idx(2,3):2, 1:nc) / 2.0_rk
            endif
        case ("H1")
            ! Wavelets get factor 2^(2-d)j/(2d) per SC and 2^(d-2)j/(2d) per WC, for d=2 this is equivalent to Linfty, for d=3 to 2^-j/6 (between Linfty and L2)
            ! JB ToDo - this needs to be checked, for 2D it is Linfty norm
            ! JB - for 3D, WC and SC seem not to be inverse of each other, I think it is 2^(d-2)j/(d) per WC and WC = SC**2
            ! JB - and then another mystery to solve, but the factor between each seems to be 2^(-2/3)
            if (params%dim == 3) then
                ! apply level shift
                val_renormed(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) = \
                val_renormed(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) * 2.0_rk**(dble((Jref-level_block)*(2.0_rk-params%dim))/2.0_rk)
                ! change SC to WC factors on this level - factor^2
                val_renormed(idx(1,1)+io(1):idx(2,1):2, idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) = \
                val_renormed(idx(1,1)+io(1):idx(2,1):2, idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc) * 2.0_rk**(2.0_rk*(params%dim-2.0_rk)/3.0_rk)
                val_renormed(idx(1,1):idx(2,1), idx(1,2)+io(2):idx(2,2):2, idx(1,3):idx(2,3), 1:nc) = \
                val_renormed(idx(1,1):idx(2,1), idx(1,2)+io(2):idx(2,2):2, idx(1,3):idx(2,3), 1:nc) * 2.0_rk**(2.0_rk*(params%dim-2.0_rk)/3.0_rk)
                val_renormed(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3)+io(3):idx(2,3):2, 1:nc) = \
                val_renormed(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3)+io(3):idx(2,3):2, 1:nc) * 2.0_rk**(2.0_rk*(params%dim-2.0_rk)/3.0_rk)
            endif
        case default
            call abort(241024, "ERROR:Unknown wavelet normalization!")
        end select
    end subroutine




end module
