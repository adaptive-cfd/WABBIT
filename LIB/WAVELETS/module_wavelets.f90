module module_wavelets
    use module_params
    use module_treelib

    implicit none

contains

#include "conversion_routines.f90"

    ! coarsen the block by one level
    subroutine restriction_2D(fine, coarse)
        implicit none

        real(kind=rk), dimension(1:,1:), intent(in) :: fine
        real(kind=rk), dimension(1:,1:), intent(out) :: coarse
        integer(kind=ik), dimension(2) :: nfine, ncoarse

        ncoarse(1) = size(coarse,1)
        ncoarse(2) = size(coarse,2)
        nfine(1)   = size(fine,1)
        nfine(2)   = size(fine,2)

#ifdef DEV
        if ( 2*ncoarse(1)-1 /= nfine(1) .or. 2*ncoarse(2)-1 /= nfine(2)) then
            write(*,*) shape(coarse), ":", shape(fine)
            call abort(888191,"ERROR: restriction_2D: arrays wrongly sized..")
        endif
#endif

        coarse(:, :) = fine(1:nfine(1):2, 1:nfine(2):2)

    end subroutine

    ! ! Version without vectors
    ! subroutine restriction_prefilter(params, u, u_filtered)
    !     implicit none
    !     type(type_params), intent(in) :: params
    !     real(kind=rk), dimension(1:,1:,1:), intent(in) :: u
    !     real(kind=rk), dimension(1:,1:,1:), intent(out) :: u_filtered

    !     if (.not. allocated(params%HD)) call abort(71717172, "wavelet not setup")

    !     call blockFilterXYZ(params, u, u_filtered, params%HD, lbound(params%HD, dim=1), ubound(params%HD, dim=1))
    ! end subroutine


    subroutine restriction_prefilter_vct(params, u, u_filtered)
        implicit none
        type(type_params), intent(in) :: params
        real(kind=rk), dimension(1:,1:,1:,1:), intent(in) :: u
        real(kind=rk), dimension(1:,1:,1:,1:), intent(out) :: u_filtered

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


    ! coarsen the block by one level
    subroutine restriction_3D(fine, coarse)
        implicit none

        real(kind=rk), dimension(1:,1:,1:), intent(in)  :: fine
        real(kind=rk), dimension(1:,1:,1:), intent(out) :: coarse
        integer(kind=ik), dimension(3) :: nfine, ncoarse

        ncoarse(1) = size(coarse,1)
        ncoarse(2) = size(coarse,2)
        ncoarse(3) = size(coarse,3)
        nfine(1) = size(fine,1)
        nfine(2) = size(fine,2)
        nfine(3) = size(fine,3)

#ifdef DEV
        if ( 2*ncoarse(1)-1 /= nfine(1) .or. 2*ncoarse(2)-1 /= nfine(2) .or. 2*ncoarse(3)-1 /= nfine(3) ) then
            call abort(888192,"ERROR: restriction_3D: arrays wrongly sized..")
        endif
#endif

        coarse(:, :, :) = fine(1:nfine(1):2,1:nfine(2):2,1:nfine(3):2)
    end subroutine
    


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
    subroutine prediction(coarse, fine, order_predictor)
        implicit none

        real(kind=rk), dimension(1:,1:,1:), intent(inout) :: fine
        real(kind=rk), dimension(1:,1:,1:), intent(inout) :: coarse
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

        case default
            call abort(23070811,"Error: unkown order_predictor="//trim(adjustl(order_predictor)))

        end select
        N = size(c,1)/2  ! which value can be interpolated first? 1-based, 4th and 6th skip 1 or 2 points

#ifdef DEV
        if ( 2*nxcoarse-1 /= nxfine .or. 2*nycoarse-1 /= nyfine .or. 2*nzcoarse-1 /= nzfine ) then
            call abort(888195,"ERROR: prediction - arrays wrongly sized..")
        endif
#endif

        ! fill matching points: the coarse and fine grid share a lot of points (as the
        ! fine grid results from insertion of one point between each coarse point).
        ! Sometimes called checkerboard copying
        fine = 0.0_rk
        fine(1:nxfine:2, 1:nyfine:2, 1:nzfine:2) = coarse(:, :, :)


        do izfine = 1, nzfine, 2
            ! in the z=const planes, we execute the 2D code.
            do ixfine= 2*N, nxfine-(2*N-1), 2
                do iyfine =  1, nyfine, 2
                    ! note in this implementation, interp coeffs run
                    ! c(1:2), c(1:4), c(1:6) for 2nd, 4th, 6th order respectively
                    do shift = 1, size(c,1)
                        shift_fine = -size(c,1)+(2*shift-1)
                        fine(ixfine,iyfine,izfine) = fine(ixfine,iyfine,izfine) + c(shift)*fine(ixfine+shift_fine,iyfine,izfine)
                    end do
                end do
            end do

            do ixfine = 1, nxfine, 1
                do iyfine =  2*N, nyfine-(2*N-1), 2
                    ! note in this implementation, interp coeffs run
                    ! c(1:2), c(1:4), c(1:6) for 2nd, 4th, 6th order respectively
                    do shift = 1, size(c,1)
                        shift_fine = -size(c,1)+(2*shift-1)
                        fine(ixfine,iyfine,izfine) = fine(ixfine,iyfine,izfine) + c(shift)*fine(ixfine, iyfine+shift_fine,izfine)
                    end do
                end do
            end do
        enddo

        ! finally, only 1D interpolation along z is missing.
        do izfine =  2*N, nzfine-(2*N-1), 2
            ! note in this implementation, interp coeffs run
            ! c(1:2), c(1:4), c(1:6) for 2nd, 4th, 6th order respectively
            do shift = 1, size(c,1)
                shift_fine = -size(c,1)+(2*shift-1)
                fine(:,:,izfine) = fine(:,:,izfine) + c(shift)*fine(:,:,izfine+shift_fine)
            end do
        enddo

    end subroutine


    ! ! Please note applying a filter requires also manipulating the ghost nodes
    ! ! Version without vectors
    ! subroutine blockFilterXYZ( params, u, u_filtered, coefs_filter, a, b )
    !     implicit none
    !     type (type_params), intent(in) :: params
    !     real(kind=rk), dimension(1:,1:,1:), intent(in) :: u
    !     real(kind=rk), dimension(1:,1:,1:), intent(inout) :: u_filtered
    !     integer(kind=ik) :: a, b
    !     real(kind=rk), intent(in) :: coefs_filter(a:b)
    !     integer(kind=ik) :: ix, iy, iz, nx, ny, nz, shift, g, Bs(1:3)
    !     real(kind=rk), allocatable, save :: u_tmp(:,:,:)

    !     ! if the filter is just 1, then we copy and we're done.
    !     ! Yes, we use such stupid filters. They are in the CDFX0 wavelets (X=2,4)
    !     if (a==0 .and. b==0 .and. abs(coefs_filter(0)-1.0_rk)<=1.0e-10_rk) then
    !         u_filtered = u
    !         return
    !     endif

    !     nx = size(u, 1)
    !     ny = size(u, 2)
    !     nz = size(u, 3)
    !     g  = params%g
    !     Bs = params%Bs

    !     if ((abs(a) > g).or.(b>g)) then
    !         write(*,*) a, b, "but g=", g
    !         call abort(202302209, "For applying the filter, not enough ghost nodes")
    !     endif

    !     if (.not. allocated(u_tmp)) allocate( u_tmp(1:nx,1:ny,1:nz) )
    !     u_tmp = u

    !     u_filtered = u_tmp
    !     u_filtered(g+1:Bs(1)+g, :, :) = 0.0_rk
    !     do ix = g+1, Bs(1)+g
    !         ! do ix = -a+1, nx-b
    !         do shift = a, b
    !             u_filtered(ix, :, :) = u_filtered(ix, :, :) + u_tmp(ix+shift, :, :)*coefs_filter(shift)
    !         enddo
    !     enddo

    !     u_tmp = u_filtered
    !     u_filtered = u_tmp
    !     u_filtered(:, g+1:Bs(2)+g, :) = 0.0_rk
    !     do iy = g+1, Bs(2)+g
    !         ! do iy = -a+1, ny-b
    !         do shift = a, b
    !             u_filtered(:, iy, :) = u_filtered(:, iy, :) + u_tmp(:, iy+shift, :)*coefs_filter(shift)
    !         enddo
    !     enddo

    !     if (nz == 1) return

    !     u_tmp = u_filtered
    !     u_filtered = u_tmp
    !     u_filtered(:, :, g+1:Bs(3)+g) = 0.0_rk
    !     do iz = g+1, Bs(3)+g
    !         do shift = a, b
    !             u_filtered(:, :, iz) = u_filtered(:, :, iz) + u_tmp(:, :, iz+shift)*coefs_filter(shift)
    !         enddo
    !     enddo

    ! end subroutine


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
        real(kind=rk), dimension(1:,1:,1:,1:), intent(in) :: u
        real(kind=rk), dimension(1:,1:,1:,1:), intent(inout) :: u_filtered
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


        nx = size(u, 1)
        ny = size(u, 2)
        nz = size(u, 3)
        nc = size(u, 4)
        g  = params%g
        Bs = params%Bs

        if ((abs(fl_l) > g).or.(fl_r>g)) then
            write(*,*) fl_l, fl_r, "but g=", g
            call abort(202302209, "For applying the filter, not enough ghost nodes")
        endif

        s = 1
        if (present(do_restriction)) then
            if (do_restriction) s = 2
        endif

        if (.not. allocated(u_tmp)) allocate( u_tmp(1:nx,1:ny,1:nz,1:nc) )
        u_tmp = u

        u_filtered = u_tmp
        if (nz == 1) then
            ! filter results are consecutively added so we need to set to 0 at the beginning
            u_filtered(g+1:Bs(1)+g:s, g+1+fl_l:Bs(2)+g+fl_r, :, :) = 0.0_rk
            do ix = g+1, Bs(1)+g, s
                do shift = fl_l, fl_r
                    ! ignore the edge parts where we do not need a filter value from for further dimensions: g+1+fl_l:Bs+g+fl_r
                    u_filtered(ix, g+1+fl_l:Bs(2)+g+fl_r, :, :) = u_filtered(ix, g+1+fl_l:Bs(2)+g+fl_r, :, :) &
                    + u_tmp(ix+shift, g+1+fl_l:Bs(2)+g+fl_r, :, :)*coefs_filter(shift)
                enddo
            enddo
        else
            ! filter results are consecutively added so we need to set to 0 at the beginning
            u_filtered(g+1:Bs(1)+g:s, g+1+fl_l:Bs(2)+g+fl_r, g+1+fl_l:Bs(3)+g+fl_r, :) = 0.0_rk
            do ix = g+1, Bs(1)+g, s
                do shift = fl_l, fl_r
                    ! ignore the edge parts where we do not need a filter value from for further dimensions: g+1+fl_l:Bs+g+fl_r
                    u_filtered(ix, g+1+fl_l:Bs(2)+g+fl_r, g+1+fl_l:Bs(3)+g+fl_r, :) = u_filtered(ix, g+1+fl_l:Bs(2)+g+fl_r, g+1+fl_l:Bs(3)+g+fl_r, :) &
                    + u_tmp(ix+shift, g+1+fl_l:Bs(2)+g+fl_r, g+1+fl_l:Bs(3)+g+fl_r, :)*coefs_filter(shift)
                enddo
            enddo
        endif

        u_tmp = u_filtered
        if (nz == 1) then
            ! filter results are consecutively added so we need to set to 0 at the beginning
            u_filtered(g+1:Bs(1)+g:s, g+1:Bs(2)+g:s, :, :) = 0.0_rk
            do iy = g+1, Bs(2)+g, s
                do shift = fl_l, fl_r
                    ! ignore the edge parts where we do not need a filter value from for further dimensions: g+1+fl_l:Bs+g+fl_r
                    u_filtered(g+1:Bs(1)+g:s, iy, :, :) = u_filtered(g+1:Bs(1)+g:s, iy, :, :) &
                       + u_tmp(g+1:Bs(1)+g:s, iy+shift, :, :)*coefs_filter(shift)
                enddo
            enddo
        else
            ! filter results are consecutively added so we need to set to 0 at the beginning
            u_filtered(g+1:Bs(1)+g:s, g+1:Bs(2)+g:s, g+1+fl_l:Bs(3)+g+fl_r, :) = 0.0_rk
            do iy = g+1, Bs(2)+g, s
                do shift = fl_l, fl_r
                    ! ignore the edge parts where we do not need a filter value from for further dimensions: g+1+fl_l:Bs+g+fl_r
                    u_filtered(g+1:Bs(1)+g:s, iy, g+1+fl_l:Bs(3)+g+fl_r, :) = u_filtered(g+1:Bs(1)+g:s, iy, g+1+fl_l:Bs(3)+g+fl_r, :) &
                       + u_tmp(g+1:Bs(1)+g:s, iy+shift, g+1+fl_l:Bs(3)+g+fl_r, :)*coefs_filter(shift)
                enddo
            enddo
        endif


        if (nz == 1) return

        u_tmp = u_filtered
        ! filter results are consecutively added so we need to set to 0 at the beginning
        u_filtered(g+1:Bs(1)+g:s, g+1:Bs(2)+g:s, g+1:Bs(3)+g:s, :) = 0.0_rk
        do iz = g+1, Bs(3)+g, s
            do shift = fl_l, fl_r
                u_filtered(g+1:Bs(1)+g:s, g+1:Bs(2)+g:s, iz, :) = u_filtered(g+1:Bs(1)+g:s, g+1:Bs(2)+g:s, iz, :) &
                   + u_tmp(g+1:Bs(1)+g:s, g+1:Bs(2)+g:s, iz+shift, :)*coefs_filter(shift)
            enddo
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
        real(kind=rk), dimension(1:,1:,1:,1:), intent(in) :: u
        real(kind=rk), dimension(1:,1:,1:,1:), intent(inout) :: u_filtered
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
            write(*,*) fl_l, fl_r, "but g=", g
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


        if (nz == 1) return

        u_tmp = u_filtered
        ! u_filtered = u_tmp  ! not necessary
        u_filtered(-fl_l+1:nx-fl_r:s, -fl_l+1:ny-fl_r:s, -fl_l+1:nz-fl_r:s, :) = 0.0_rk  ! these values will be filtered and overwritten again
        do iz = -fl_l+1, nz-fl_r, s
            do shift = fl_l, fl_r
                ! after previous filtering we can ignore some x- and y-values which were not effected to decrease the computational load
                u_filtered(-fl_l+1:nx-fl_r:s, -fl_l+1:ny-fl_r:s, iz, :) = u_filtered(-fl_l+1:nx-fl_r:s, -fl_l+1:ny-fl_r:s, iz, :) &
                   + u_tmp(-fl_l+1:nx-fl_r:s, -fl_l+1:ny-fl_r:s, iz+shift, :)*coefs_filter(shift)
            enddo
        enddo

    end subroutine


    ! applies one of params% HD GD HR GR filters in each direction
    ! ToDo: for HD and GD filter, we might benefit from considering restriction, for HR from skipping points
    subroutine blockFilterCustom_vct( params, u, u_filtered, filter_x, filter_y, filter_z )
        implicit none
        type (type_params), intent(in) :: params
        real(kind=rk), dimension(1:,1:,1:,1:), intent(in) :: u
        real(kind=rk), dimension(1:,1:,1:,1:), intent(inout) :: u_filtered
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
            call abort(202302201, "Unknown wavelet filter, must be one of HD GD HR GR")

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
            call abort(202302201, "Unknown wavelet filter, must be one of HD GD HR GR")

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
            call abort(202302201, "Unknown wavelet filter, must be one of HD GD HR GR")

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
        real(kind=rk), dimension(1:,1:,1:,1:), intent(in) :: u
        real(kind=rk), dimension(1:,1:,1:,1:), intent(inout) :: u_filtered
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
        real(kind=rk), dimension(1:,1:,1:,1:), intent(inout) :: u

        real(kind=rk), allocatable, dimension(:,:,:,:), save :: sc, wc, test, ucopy
        integer(kind=ik) :: nx, ny, nz, nc, g, Bs(1:3), ii
        ! integer(kind=ik) :: ag, bg, ah, bh, ix, iy, iz, ic, shift
        ! real(kind=rk) :: ug, uh

        call WaveDecomposition_dim1( params, u )
    end subroutine


    ! computes a one-level wavelet decomposition of a block.
    ! The computation of coefficients is possible on the entire block, but
    ! without sync'ing the reconstruction is not possible on the entire data.
    ! We have to sync wavelet-transformed blocks.
    ! Data are stored in Spaghetti-order (not Mallat-Order)
    subroutine waveletDecomposition_block_old(params, u)
        implicit none
        type (type_params), intent(in) :: params
        real(kind=rk), dimension(1:,1:,1:,1:), intent(inout) :: u

        real(kind=rk), allocatable, dimension(:,:,:,:), save :: sc, wc, test, ucopy
        integer(kind=ik) :: nx, ny, nz, nc, g, Bs(1:3), ii
        ! integer(kind=ik) :: ag, bg, ah, bh, ix, iy, iz, ic, shift
        ! real(kind=rk) :: ug, uh

        nx = size(u, 1)
        ny = size(u, 2)
        nz = size(u, 3)
        nc = size(u, 4)
        g  = params%g
        Bs = params%Bs

#ifdef DEV
        if (nz /= 1) call abort(7223839, "currently 2D only")
        if (modulo(Bs(1),2)/=0) call abort(7223139, "only even Bs is possible with biorthogonal wavelets")
        if (modulo(Bs(2),2)/=0) call abort(7223139, "only even Bs is possible with biorthogonal wavelets")
#endif

        if (allocated(sc)) then
            if ((size(sc,1)/=nx).or.(size(sc,2)/=ny).or.(size(sc,3)/=nz).or.(size(sc,4)/=nc)) deallocate(sc)
        endif
        if (allocated(wc)) then
            if ((size(wc,1)/=nx).or.(size(wc,2)/=ny).or.(size(wc,3)/=nz).or.(size(wc,4)/=nc)) deallocate(wc)
        endif

        if (.not. allocated(sc)) allocate( sc(1:nx, 1:ny, 1:nz, 1:nc) )
        if (.not. allocated(wc)) allocate( wc(1:nx, 1:ny, 1:nz, 1:nc) )
        ! if (.not. allocated(u_wc)) allocate( u_wc(1:nx, 1:ny, 1:nz, 1:nc) )

        ! alternative algorithm (true Mallat ordering, but with ghost nodes)
        call blockFilterCustom1_vct( params, u, sc, "HD", "x" )
        call blockFilterCustom1_vct( params, u, wc, "GD", "x" )

        u( 1:Bs(1)/2, :, :, :)       = sc( (g+1):(Bs(1)+g):2, :, :, :)
        u( Bs(1)/2+1:Bs(1), :, :, :) = wc( (g+1):(Bs(1)+g):2, :, :, :)

        call blockFilterCustom1_vct( params, u, sc, "HD", "y" )
        call blockFilterCustom1_vct( params, u, wc, "GD", "y" )

        u(:, 1:Bs(2)/2, :, :) =  sc(:, (g+1):(Bs(2)+g):2, :, :)
        u(:, Bs(2)/2+1:Bs(2), :, :) = wc(:, (g+1):(Bs(2)+g):2, :, :)

        ! Note at this point U contains SC/WC in "true Mallat ordering", but note
        ! that data includes ghost nodes.

        ! copy to Spaghetti ordering
        sc=0.0_rk
        sc( (g+1):(Bs(1)+g):2, (g+1):(Bs(1)+g):2, :, :) = u(1:Bs(1)/2, 1:Bs(2)/2, :, :)
        sc( (g+2):(Bs(1)+g):2, (g+1):(Bs(1)+g):2, :, :) = u(1:Bs(1)/2, Bs(2)/2+1:Bs(2), :, :)
        sc( (g+1):(Bs(1)+g):2, (g+2):(Bs(1)+g):2, :, :) = u(Bs(1)/2+1:Bs(1), 1:Bs(2)/2, :, :)
        sc( (g+2):(Bs(1)+g):2, (g+2):(Bs(1)+g):2, :, :) = u(Bs(1)/2+1:Bs(1), Bs(2)/2+1:Bs(2), :, :)

        ! copy to the back spaghetti-ordered coefficients to the block
        u = sc
    end subroutine

    !-------------------------------------------------------------------------------

    subroutine waveletReconstruction_block(params, u)
        implicit none
        type (type_params), intent(in) :: params
        real(kind=rk), dimension(1:,1:,1:,1:), intent(inout) :: u

        call WaveReconstruction_dim1( params, u )
    end subroutine

    !-------------------------------------------------------------------------------

!     subroutine waveletReconstruction_block_old(params, u)
!         implicit none
!         type (type_params), intent(in) :: params
!         real(kind=rk), dimension(1:,1:,1:,1:), intent(inout) :: u
!
!         real(kind=rk), allocatable, dimension(:,:,:,:,:), save :: wc
!         integer(kind=ik) :: nx, ny, nz, nc, g, Bs(1:3)
!
!         nx = size(u, 1)
!         ny = size(u, 2)
!         nz = size(u, 3)
!         nc = size(u, 4)
!         g  = params%g
!         Bs = params%Bs
!
!         if (allocated(sc)) then
!             if ((size(sc,1)/=nx).or.(size(sc,2)/=ny).or.(size(sc,3)/=nz).or.(size(sc,4)/=nc)) deallocate(sc)
!         endif
!         if (allocated(wcx)) then
!             if ((size(wcx,1)/=nx).or.(size(wcx,2)/=ny).or.(size(wcx,3)/=nz).or.(size(wcx,4)/=nc)) deallocate(wcx)
!         endif
!         if (allocated(wcy)) then
!             if ((size(wcy,1)/=nx).or.(size(wcy,2)/=ny).or.(size(wcy,3)/=nz).or.(size(wcy,4)/=nc)) deallocate(wcy)
!         endif
!         if (allocated(wcxy)) then
!             if ((size(wcxy,1)/=nx).or.(size(wcxy,2)/=ny).or.(size(wcxy,3)/=nz).or.(size(wcxy,4)/=nc)) deallocate(wcxy)
!         endif
!
!         if (.not. allocated(sc  )) allocate(   sc(1:nx, 1:ny, 1:nz, 1:nc) )
!         if (.not. allocated(wcx )) allocate(  wcx(1:nx, 1:ny, 1:nz, 1:nc) )
!         if (.not. allocated(wcy )) allocate(  wcy(1:nx, 1:ny, 1:nz, 1:nc) )
!         if (.not. allocated(wcxy)) allocate( wcxy(1:nx, 1:ny, 1:nz, 1:nc) )
!
! #ifdef DEV
!         if (nz /= 1) call abort(7223839, "currently 2D only")
!         if (modulo(Bs(1), 2) /= 0) call abort(99111,"This code requires Bs even")
!         if (modulo(Bs(2), 2) /= 0) call abort(99111,"This code requires Bs even")
! #endif
!
!         call spaghetti2inflatedMallat_block(params, u, sc, wcx, wcy, wcxy)
!
!         call blockFilterCustom_vct( params, sc  , sc  , "HR", "HR", "--" ) ! inplace should work, only a copy statement from the input
!         call blockFilterCustom_vct( params, wcx , wcx , "HR", "GR", "--" ) ! inplace should work, only a copy statement from the input
!         call blockFilterCustom_vct( params, wcy , wcy , "GR", "HR", "--" ) ! inplace should work, only a copy statement from the input
!         call blockFilterCustom_vct( params, wcxy, wcxy, "GR", "GR", "--" ) ! inplace should work, only a copy statement from the input
!
!         u = sc + wcx + wcy + wcxy
!     end subroutine

    !-------------------------------------------------------------------------------

    ! ensures that 3/4 of the numbers are zero - required for reconstruction
    ! note when copying Spaghetti to Mallat, this is automatically done, but
    ! when manipulating coefficients, it may happen that we set nonzero values
    ! where a zero should be.
    subroutine setRequiredZerosWCSC_block(params, u)
        implicit none
        type (type_params), intent(in) :: params
        real(kind=rk), dimension(1:,1:,1:,1:), intent(inout) :: u

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
        real(kind=rk), dimension(1:,1:,1:,1:), intent(inout) :: wc   !< input/output in WD spaghetti format
        integer(kind=ik), intent(in) :: neighborhood                 !< Which neighborhood to apply manipulation
        integer(kind=ik), intent(in), optional :: Nwcl_optional, Nwcr_optional
        !> If set to true, we set a really large value into the patch so that the simulation can crash on purpose
        !> This is used to make clear that those values are not correct and should not be used
        logical, intent(in), optional :: set_garbage

        integer(kind=ik) :: Nwcl, Nwcr, i_neighborhood
        integer(kind=ik) :: nx, ny, nz, nc, g, Bs(1:3), d, idx(2,3), i_dim, i_set
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

            ! we need to know if the first point is a SC or WC for the patch we check and skip it if it is a WC
            !     1 2 3 4 5 6 7 8 9 A B C
            !     G G G S W S W S W G G G
            !                 I I I
            ! 1-C - index numbering in hex format, G - ghost point, S - SC, W - WC, I - point of patch to be checked
            ! Patch I is checked, but we need to know that index 7 has a WC and should be skipped
            ! this is for parity with inflatedMallat version where SC and WC are situated on the SC indices of the spaghetti format
            ! for g=odd, the SC are on even numbers; for g=even, the SC are on odd numbers
            idx(1, 1:params%dim) = idx(1, 1:params%dim) + modulo(g + idx(1, 1:params%dim) + 1, 2)
            ! also, when the patch size to be copied is odd, the last point is only partially included but its WC have to be considered
            ! this gives problem if the last point is a SC so we need to handle this special case
            do i_dim = 1, params%dim
                if (idx(2, i_dim) /= size(wc, i_dim)) then
                    idx(2, i_dim) = idx(2, i_dim) + modulo(g + idx(2, i_dim), 2)
                endif
            enddo

            ! set really low number for ghost patch if we shoudln't access it
            setNumber = 0.0
            if (setGarbage .and. i_set == 2) setNumber = -9e200_rk

            ! set values, we have to skip the SC
            wc(idx(1,1)+1:idx(2,1):2, idx(1,2)  :idx(2,2):2, idx(1,3)  :idx(2,3):2, 1:nc) = setNumber
            wc(idx(1,1)  :idx(2,1)  , idx(1,2)+1:idx(2,2):2, idx(1,3)  :idx(2,3):2, 1:nc) = setNumber
            if (params%dim == 3) then
                wc(idx(1,1)  :idx(2,1)  , idx(1,2)  :idx(2,2)  , idx(1,3)+1:idx(2,3):2, 1:nc) = setNumber
            endif
        enddo

        ! wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc, 2:d) = setNumber
    end subroutine



    subroutine coarseExtensionManipulateSC_block(params, wc, u_copy, neighborhood, ijk)
        implicit none

        type (type_params), intent(in) :: params
        real(kind=rk), dimension(1:,1:,1:,1:), intent(inout) :: wc   !< input/output in WD spaghetti format
        real(kind=rk), dimension(1:,1:,1:,1:), intent(in) :: u_copy  !< original values from which we copy
        integer(kind=ik), intent(in)   :: neighborhood               !< Which neighborhood to apply manipulation
        integer(kind=ik), intent(in), optional   :: ijk(2,3)         !< ijk of patch that we only care about, if not given it is the full block

        integer(kind=ik) :: nx, ny, nz, nc, g, Bs(1:3), i_set, i_neighborhood
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

            ! we need to know if the first point is a SC or WC for the patch we check and skip it if it is a WC
            !     1 2 3 4 5 6 7 8 9 A B C
            !     G G G S W S W S W G G G
            !                 I I I
            ! 1-C - index numbering in hex format, G - ghost point, S - SC, W - WC, I - point of patch to be checked
            ! Patch I is checked, but we need to know that index 7 has a WC and should be skipped
            ! this is for parity with inflatedMallat version where SC and WC are situated on the SC indices of the spaghetti format
            ! for g=odd, the SC are on even numbers; for g=even, the SC are on odd numbers
            idx(1, 1:params%dim) = idx(1, 1:params%dim) + modulo(g + idx(1, 1:params%dim) + 1, 2)
            ! also, when the patch size to be copied is odd, the last point is only partially included but its WC have to be considered
            ! However, here we treat SC so we can ignore it completely

            ! sometimes we do not want to look at the full block, then we can skip the copy where the indices do not matter
            if (present(ijk)) then
                ! check in each dimension if the indices lie outside the affected range
                ! copy patch has to end before beginning of affected region or begin before end of affected region
                if (idx(2,1) < ijk(1,1) .or. idx(1,1) > ijk(2,1)) skip_copy = .true.
                if (idx(2,2) < ijk(1,2) .or. idx(1,2) > ijk(2,2)) skip_copy = .true.
                if (idx(2,3) < ijk(1,3) .or. idx(1,3) > ijk(2,3) .and. params%dim == 3) skip_copy = .true.
            endif

            if (.not. skip_copy) then
                wc(idx(1,1):idx(2,1):2, idx(1,2):idx(2,2):2, idx(1,3):idx(2,3):2, 1:nc) = &
                    u_copy(idx(1,1):idx(2,1):2, idx(1,2):idx(2,2):2, idx(1,3):idx(2,3):2, 1:nc)
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
        integer(kind=ik) :: i, g_min, a, block_min

        if (allocated(params%GR)) deallocate(params%HD)
        if (allocated(params%GD)) deallocate(params%GD)
        if (allocated(params%HR)) deallocate(params%HR)
        if (allocated(params%GR)) deallocate(params%GR)

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

        ! the wavelet filter banks:
        ! HD - low pass decomposition filter, H_TILDE
        ! GD - high pass decomposition filter, G_TILDE - for CDF it is always HR with different sign for every second off-center value
        ! HR - low pass reconstruction filter, H
        ! GR - high pass reconstruction filter, G - for CDF it is always HD with different sign for every second off-center value

        ! check if we have a CDF wavelet
        if (params%wavelet(1:3) /= "CDF") then
            call abort( 3006221, "Unkown bi-orthogonal wavelet specified. Set course for adventure! params%wavelet="//trim(adjustl(params%wavelet)) )
        endif
        
        ! The HR filter is always defined by the first number of the CDFXY filter
        ! The HD filter changes depending on both numbers
        if (params%wavelet(4:4) == "2") then
            ! H TILDE filter
            allocate( params%HR(-1:1) )
            params%HR = (/ 0.5_rk, 1.0_rk, 0.5_rk  /)

            ! H filter
            if (params%wavelet(5:5) .eq. "0") then
                allocate( params%HD(0:0) )
                params%HD = (/1.0_rk/)
            elseif (params%wavelet(5:5) == "2") then
                allocate( params%HD(-2:2) )
                params%HD = (/-1.0_rk, +2.0_rk, +6.0_rk, +2.0_rk, -1.0_rk/) / 8.0_rk
            elseif (params%wavelet(5:5) == "4") then
                allocate( params%HD(-4:4) )
                ! from Daubechies - Ten lectures on wavelets, Table 8.2
                params%HD = (/3.0_rk, -6.0_rk, -16.0_rk, 38.0_rk, 90.0_rk, 38.0_rk, -16.0_rk, -6.0_rk, 3.0_rk/) / 128.0_rk
            elseif (params%wavelet(5:5) == "6") then
                allocate( params%HD(-6:6) )
                ! from Daubechies - Ten lectures on wavelets, Table 8.2
                params%HD = (/-5.0_rk, 10.0_rk, 34.0_rk, -78.0_rk, -123.0_rk, 324.0_rk, 700.0_rk, 324.0_rk, -123.0_rk, -78.0_rk, 34.0_rk, 10.0_rk, -5.0_rk/) / 1024.0_rk
            elseif (params%wavelet(5:5) == "8") then
                allocate( params%HD(-8:8) )
                ! from Daubechies - Ten lectures on wavelets, Table 8.2
                params%HD = (/35.0_rk, -70.0_rk, -300.0_rk, 670.0_rk, 1228.0_rk, -3126.0_rk, -3796.0_rk, 10718.0_rk, 22050.0_rk, &
                    10718.0_rk, -3796.0_rk, -3126.0_rk, 1228.0_rk, 670.0_rk, -300.0_rk, -70.0_rk, 35.0_rk/) / 32768.0_rk
            else
                call abort( 3006221, "Unkown bi-orthogonal wavelet specified. Set course for adventure! params%wavelet="//trim(adjustl(params%wavelet)) )
            endif

        elseif (params%wavelet(4:4) == "4") then
            ! H TILDE filter
            allocate( params%HR(-3:3) )
            params%HR = (/ -1.0_rk/16.0_rk, 0.0_rk, 9.0_rk/16.0_rk, 1.0_rk, 9.0_rk/16.0_rk, 0.0_rk, -1.0_rk/16.0_rk  /)

            ! H filter
            if (params%wavelet(5:5) == "0") then
                allocate( params%HD(0:0) )
                params%HD = (/1.0_rk/)
            elseif (params%wavelet(5:5) == "2") then
                ! Sweldens paper, "The Lifting Scheme: A Custom-Design
                ! Construction of Biorthogonal Wavelets" table 3 for N_tilde=2
                allocate( params%HD(-4:4) )
                params%HD = (/ 2.0_rk**(-6.0_rk), 0.0_rk, &
                              -2.0_rk**(-3.0_rk), 2.0_rk**(-2.0_rk), &
                       23.0_rk*2.0_rk**(-5.0_rk), 2.0_rk**(-2.0_rk), &
                              -2.0_rk**(-3.0_rk), 0.0_rk, &
                               2.0_rk**(-6.0_rk) /)
            elseif (params%wavelet(5:5) == "4") then
                ! Sweldens paper, "The Lifting Scheme: A Custom-Design
                ! Construction of Biorthogonal Wavelets" table 2 for N_tilde=4
                allocate( params%HD(-6:6) )
                params%HD = (/ -2.0_rk**(-9.0_rk),         0.0_rk, &
                         9.0_rk*2.0_rk**(-8.0_rk),        -2.0_rk**(-5.0_rk), &
                       -63.0_rk*2.0_rk**(-9.0_rk),  9.0_rk*2.0_rk**(-5.0_rk), &
                        87.0_rk*2.0_rk**(-7.0_rk),  9.0_rk*2.0_rk**(-5.0_rk), &
                       -63.0_rk*2.0_rk**(-9.0_rk),        -2.0_rk**(-5.0_rk), &
                         9.0_rk*2.0_rk**(-8.0_rk),         0.0_rk, &
                               -2.0_rk**(-9.0_rk)/)
            elseif (params%wavelet(5:5) == "6") then
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
                call abort( 3006221, "Unkown bi-orthogonal wavelet specified. Set course for adventure! params%wavelet="//trim(adjustl(params%wavelet)) )      
            endif
        elseif (params%wavelet(4:4) == "6") then
            ! H TILDE filter  
            allocate( params%HR(-5:5) )
            params%HR = (/ 3.0_rk, 0.0_rk, -25.0_rk, 0.0_rk, 150.0_rk, 256.0_rk, 150.0_rk, 0.0_rk, -25.0_rk, 0.0_rk, 3.0_rk /) / 256.0_rk

            ! H filter
            if (params%wavelet(5:5) == "0") then
                allocate( params%HD(0:0) )
                params%HD = (/1.0_rk/)
            elseif (params%wavelet(5:5) == "2") then
                ! Sweldens paper, "The Lifting Scheme: A Custom-Design
                ! Construction of Biorthogonal Wavelets" table 3 for N_tilde=2
                allocate( params%HD(-6:6) )
                params%HD = (/ -3.0_rk*2.0_rk**(-10.0_rk), 0.0_rk, &
                               11.0_rk*2.0_rk**( -9.0_rk), 0.0_rk, &
                             -125.0_rk*2.0_rk**(-10.0_rk), 2.0_rk**(-2.0_rk), &
                              181.0_rk*2.0_rk**( -8.0_rk), 2.0_rk**(-2.0_rk), &
                             -125.0_rk*2.0_rk**(-10.0_rk), 0.0_rk, &
                               11.0_rk*2.0_rk**( -9.0_rk), 0.0_rk, &
                               -3.0_rk*2.0_rk**(-10.0_rk)/)
            elseif (params%wavelet(5:5) == "4") then
                ! Sweldens paper, "The Lifting Scheme: A Custom-Design
                ! Construction of Biorthogonal Wavelets" table 3 for N_tilde=4
                allocate( params%HD(-8:8) )
                params%HD = (/ 3.0_rk*2.0_rk**(-13.0_rk),        0.0_rk, &
                             -13.0_rk*2.0_rk**(-11.0_rk),        0.0_rk, &
                              87.0_rk*2.0_rk**(-11.0_rk),       -2.0_rk**(-5.0_rk), &
                            -243.0_rk*2.0_rk**(-11.0_rk), 9.0_rk*2.0_rk**(-5.0_rk), &
                            2721.0_rk*2.0_rk**(-12.0_rk), 9.0_rk*2.0_rk**(-5.0_rk), &
                            -243.0_rk*2.0_rk**(-11.0_rk),       -2.0_rk**(-5.0_rk), &
                              87.0_rk*2.0_rk**(-11.0_rk),        0.0_rk, &
                             -13.0_rk*2.0_rk**(-11.0_rk),        0.0_rk, &
                               3.0_rk*2.0_rk**(-13.0_rk) /)  
            elseif (params%wavelet(5:5) == "6") then
                ! Sweldens paper, "The Lifting Scheme: A Custom-Design
                ! Construction of Biorthogonal Wavelets" table 3 for N_tilde=6
                allocate( params%HD(-10:10) )
                params%HD = (/ -9.0_rk*2.0_rk**(-17.0_rk),          0.0_rk, &
                               75.0_rk*2.0_rk**(-16.0_rk),          0.0_rk, &
                            -1525.0_rk*2.0_rk**(-17.0_rk),   3.0_rk*2.0_rk**(-9.0_rk), &
                              825.0_rk*2.0_rk**(-14.0_rk), -25.0_rk*2.0_rk**(-9.0_rk), &
                            -7425.0_rk*2.0_rk**(-16.0_rk),  75.0_rk*2.0_rk**(-8.0_rk), &
                            21201.0_rk*2.0_rk**(-15.0_rk),  75.0_rk*2.0_rk**(-8.0_rk), &
                            -7425.0_rk*2.0_rk**(-16.0_rk), -25.0_rk*2.0_rk**(-9.0_rk), &
                              825.0_rk*2.0_rk**(-14.0_rk),   3.0_rk*2.0_rk**(-9.0_rk), &
                            -1525.0_rk*2.0_rk**(-17.0_rk),          0.0_rk, &
                               75.0_rk*2.0_rk**(-16.0_rk),          0.0_rk, &
                               -9.0_rk*2.0_rk**(-17.0_rk) /)
            else
                call abort( 3006221, "Unkown bi-orthogonal wavelet specified. Set course for adventure! params%wavelet="//trim(adjustl(params%wavelet)) )
            endif    
        else
            call abort( 3006221, "Unkown bi-orthogonal wavelet specified. Set course for adventure! params%wavelet="//trim(adjustl(params%wavelet)) )
        endif

        ! G TILDE filter - HR filter with different sign for every second off-center value
        allocate( params%GD( lbound(params%HR, dim=1)+1:ubound(params%HR, dim=1)+1) )
        do i = lbound(params%GD, dim=1), ubound(params%GD, dim=1)
            params%GD(i) = (-1.0_rk)**(i-1) * params%HR(i-1)
        enddo
                    
        ! G filter - HD filter with different sign for every second off-center value
        if (params%wavelet(5:5) == "0") then  ! for unlifted wavelets GR filter is set larger than HD filter
            allocate( params%GR(-2:0) )
            params%GR = (/ 0.0_rk, 1.0_rk, 0.0_rk /)
        else
            allocate( params%GR( lbound(params%HD, dim=1)-1:ubound(params%HD, dim=1)-1) )
            do i = lbound(params%GR, dim=1), ubound(params%GR, dim=1)
                params%GR(i) = (-1.0_rk)**(i+1) * params%HD(i+1)
            enddo
        endif

        ! Unlifted or lifted - every CDFX0er wavelet is considered unlifted and the rest lifted
        params%isLiftedWavelet = params%wavelet(5:5) /= "0"

        ! order predictor is decided by X in CDFXY (it is coinciding with the HR filter actually)
        if (params%wavelet(4:4) == "2") then
            params%order_predictor = "multiresolution_2nd"
        elseif (params%wavelet(4:4) == "4") then
            params%order_predictor = "multiresolution_4th"
        elseif (params%wavelet(4:4) == "6") then
            params%order_predictor = "multiresolution_6th"
        endif

        ! minimum number of ghost nodes required for this wavelet
        ! value is determined by largest filter bank
        g_min = 0
        g_min = max(g_min, abs(lbound(params%HD, dim=1)))
        g_min = max(g_min, abs(ubound(params%HD, dim=1)))
        g_min = max(g_min, abs(lbound(params%GD, dim=1)))
        g_min = max(g_min, abs(ubound(params%GD, dim=1)))
        g_min = max(g_min, abs(lbound(params%HR, dim=1)))
        g_min = max(g_min, abs(ubound(params%HR, dim=1)))
        g_min = max(g_min, abs(lbound(params%GR, dim=1)))
        g_min = max(g_min, abs(ubound(params%GR, dim=1)))

        ! g_RHS is decided by X in CDFXY
        if (present(g_RHS)) then
            if (params%wavelet(4:4) == "2") then
                g_RHS = 1
            elseif (params%wavelet(4:4) == "4") then
                g_RHS = 2
            elseif (params%wavelet(4:4) == "6") then
                g_RHS = 3
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

        ! for unlifted wavelets no SC are copied, however some WC do have to be wiped in case we ned to reconstruct (CVS and image denoising)
        ! normally, if no WC are altered we can simply recopy the valeus that we had
        if (.not. params%isLiftedWavelet) then
            params%Nwcl = params%Nwcl - 2
            params%Nreconl = params%Nreconl - 2
            params%Nwcr = params%Nwcr - 2
            params%Nreconr = params%Nreconr - 2
        endif

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

        ! ! conditions for minimum blocksize for lifted wavelets arising from coarse extension:
        ! !    1. coarse extension reconstructs needs SC or WC from finer neighbors
        ! !    2. finer neighbors need values for their reconstruction that we altered with our reconstruction
        ! ! the Bs check is because sometimes this is called and Bs is not set yet, I have no idea when to check it then
        ! block_min = 0
        ! if (params%isLiftedWavelet .and. maxval(params%Bs(:)) /= 0 .and. .not. params%CVS) then
        !     ! first condition: we need ghost points in wavelet decomposed form from finer neighbours for our reconstruction which we cannot get
        !     ! this is critical and we currently cannot get those values
        !     block_min = max(block_min, params%Nreconr + max(abs(lbound(params%GR, dim=1))-2, abs(lbound(params%HR, dim=1))-2))
        !     block_min = max(block_min, params%Nreconl + max(ubound(params%GR, dim=1), ubound(params%HR, dim=1)))
        !     ! second condition: our finer neighbors need our reconstructed values to reconstruct itself it's values
        !     ! if that is the case we have an upwards dependency which is currently not handled
        !     block_min = max(block_min, params%Nreconr + (ubound(params%HR, dim=1)+1)/2)
        !     block_min = max(block_min, params%Nreconl + abs(lbound(params%HR, dim=1))/2)

        !     ! block needs to be larger as constraints
        !     block_min = block_min+1

        !     if (any(params%Bs(:params%dim) < block_min)) then
        !         write(*,'(A, A, 3(i3), A, i3)') trim(adjustl(params%wavelet)), " Bs=", params%Bs(:), " < block_min=", block_min
        !         call abort(8888881, "The selected wavelet requires larger blocksizes to do the correct coarse extension.")
        !     endif
        ! endif

        if (params%rank==0 .and. verbose1) then
            write(*, '("  ╭─╮      ╭─╮                           ╭─╮         ╭───╮           ╭─╮        ")')
            write(*, '("──╯ │ ╭────╯ │ ╭───   Wavelet-setup   ───╯ │ ╭───────╯   │   ╭───────╯ │ ╭──────")')
            write(*, '("    ╰─╯      ╰─╯                           ╰─╯           ╰───╯         ╰─╯      ")')
            write(*,'(2A)') "The wavelet is ", trim(adjustl(params%wavelet))
            write(*,'(A55, i4, i4)') "During coarse extension, we will copy SC (L,R):", params%Nscl, params%Nscr
            write(*,'(A55, i4, i4)') "During coarse extension, we will delete WC (L,R):", params%Nwcl, params%Nwcr
            write(*,'(A55, i4, i4)') "During coarse extension, we will reconstruct u (L,R):", params%Nreconl, params%Nreconr
            ! if (block_min /= 0 .and. .not. params%CVS) write(*,'(A55, i4)') "From coarse extension we have a minimum blocksize of:", block_min
            write(*,'(2A)') "The predictor is: ", trim(adjustl(params%order_predictor))
            write(*,'(A,"[",i2,":",i1,"]=",14(es12.4,1x))') "HD", lbound(params%HD, dim=1), ubound(params%HD, dim=1), params%HD
            write(*,'(A,"[",i2,":",i1,"]=",14(es12.4,1x))') "GD", lbound(params%GD, dim=1), ubound(params%GD, dim=1), params%GD
            write(*,'(A,"[",i2,":",i1,"]=",14(es12.4,1x))') "HR", lbound(params%HR, dim=1), ubound(params%HR, dim=1), params%HR
            write(*,'(A,"[",i2,":",i1,"]=",14(es12.4,1x))') "GR", lbound(params%GR, dim=1), ubound(params%GR, dim=1), params%GR
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
            write(*, '(20("╭─╮ "))')
            write(*, '(20("╯ ╰─"))')
        endif
    end subroutine



    ! 1D filter, pass skip_g=params%g to skip ghost points, pass do_restriction to filter only every second point
    subroutine filter1dim(params, u, u_filtered, filter, fl_l, fl_r, skip_g, sampling)
        implicit none
        type (type_params), intent(in) :: params                        !< good ol' params
        real(kind=rk), dimension(1:), intent(in) :: u                   !< the signal to be filtered
        real(kind=rk), dimension(1:), intent(inout) :: u_filtered       !< the resulting filtered signal
        integer, intent(in) :: fl_l, fl_r                               !< filter bound indices left and right
        real(kind=rk), dimension(fl_l:fl_r), intent(in) :: filter       !< the actual filter
        integer, intent(in) :: skip_g                                   !< 0 to filter everything, params%g to skip ghost points
        !> sampling rate, 1 for normal mode, 2 for restriction, no optional parameter as this is inside critical point-loop and I want to avoid if-clauses
        integer, intent(in) :: sampling                                        

        integer(kind=ik) :: N, i, j

        N = size(u)
        u_filtered = 0.0_rk

        ! apply filter f to periodic signal u, i.e. convolute with filter
        do i =  1+skip_g, N-skip_g, sampling
            do j = fl_l, fl_r
                u_filtered(i) = u_filtered(i) + u(i+j) * filter(j)
            enddo
        enddo
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
                    call filter1dim(params, u_wc(:,iy,iz,ic), buffer1(1:nx), params%HD, lbound(params%HD,dim=1), ubound(params%HD,dim=1), skip_g=params%g, sampling=2)

                    ! high-pass filter (these guys are the details)
                    call filter1dim(params, u_wc(:,iy,iz,ic), buffer2(1:nx), params%GD, lbound(params%GD,dim=1), ubound(params%GD,dim=1), skip_g=params%g, sampling=2)

                    ! decimation by 2, sort into array in spaghetti form SC WC
                    u_wc((g+1):(Bs(1)+g):2,iy,iz,ic) = buffer1( (g+1):(Bs(1)+g):2 )
                    u_wc((g+2):(Bs(1)+g):2,iy,iz,ic) = buffer2( (g+1):(Bs(1)+g):2 )
                enddo
            enddo
        enddo

        ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        do ic = 1, nc
            do ix = (g+1), (Bs(1)+g)
                do iz = 1, nz
                    ! low-pass filter (scaling function)
                    call filter1dim(params, u_wc(ix,:,iz,ic), buffer1(1:ny), params%HD, lbound(params%HD,dim=1), ubound(params%HD,dim=1), skip_g=params%g, sampling=2)

                    ! high-pass filter (these guys are the details)
                    call filter1dim(params, u_wc(ix,:,iz,ic), buffer2(1:ny), params%GD, lbound(params%GD,dim=1), ubound(params%GD,dim=1), skip_g=params%g, sampling=2)

                    ! decimation by 2, sort into array in spaghetti form SC WC
                    u_wc(ix,(g+1):(Bs(2)+g):2,iz,ic) = buffer1( (g+1):(Bs(2)+g):2 )
                    u_wc(ix,(g+2):(Bs(2)+g):2,iz,ic) = buffer2( (g+1):(Bs(2)+g):2 )
                enddo
            enddo
        enddo

        ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (params%dim == 3) then
            do ic = 1, nc
                do ix = (g+1), (Bs(1)+g)
                    do iy = (g+1), (Bs(2)+g)
                        ! low-pass filter (scaling function)
                        call filter1dim(params, u_wc(ix,iy,:,ic), buffer1(1:nz), params%HD, lbound(params%HD,dim=1), ubound(params%HD,dim=1), skip_g=params%g, sampling=2)

                        ! high-pass filter (these guys are the details)
                        call filter1dim(params, u_wc(ix,iy,:,ic), buffer2(1:nz), params%GD, lbound(params%GD,dim=1), ubound(params%GD,dim=1), skip_g=params%g, sampling=2)

                        ! decimation by 2, sort into array in spaghetti form SC WC
                        u_wc(ix,iy,(g+1):(Bs(3)+g):2,ic) = buffer1( (g+1):(Bs(3)+g):2 )
                        u_wc(ix,iy,(g+2):(Bs(3)+g):2,ic) = buffer2( (g+1):(Bs(3)+g):2 )
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
    subroutine WaveReconstruction_dim1( params, u_wc )
        implicit none
        type (type_params), intent(in) :: params
        !> Input is in spaghetti ordering (synchronized, ie with ghost nodes)
        real(kind=rk), dimension(:, :, :, :), intent(inout) :: u_wc
        real(kind=rk), dimension(:), allocatable, save :: buffer1, buffer2, buffer3
        integer(kind=ik) :: ix, iy, iz, ic, g, Bs(1:3), nx, ny, nz, nc, maxn
        integer(kind=ik) :: io  ! if g is odd the SCs start from the second point
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
        if (.not.allocated(buffer3)) allocate(buffer3(1:maxn))

        ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        do ic = 1, nc
            do iy = 1, ny
                do iz = 1, nz
                    ! fill upsampling buffer for low-pass filter: every second point
                    ! apply low-pass filter to upsampled signal
                    buffer3 = 0.0_rk
                    buffer3(1+io:nx-io:2) = u_wc(1+io:Bs(1)+2*g-io:2, iy, iz, ic) ! SC
                    call filter1dim(params, buffer3(1:nx), buffer1(1:nx), params%HR, lbound(params%HR,dim=1), ubound(params%HR,dim=1), skip_g=params%g, sampling=1)

                    ! fill upsampling buffer for high-pass filter: every second point
                    buffer3 = 0.0_rk
                    buffer3(1+io:nx-io:2) = u_wc(2+io:Bs(1)+2*g-io:2, iy, iz, ic) ! WC
                    call filter1dim(params, buffer3(1:nx), buffer2(1:nx), params%GR, lbound(params%GR,dim=1), ubound(params%GR,dim=1), skip_g=params%g, sampling=1)

                    u_wc(:, iy, iz, ic) = buffer1(1:nx) + buffer2(1:nx)
                enddo
            enddo
        enddo

        ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        do ic = 1, nc
            do ix = (g+1), (Bs(1)+g)  ! ignore ghost points
                do iz = 1, nz
                    ! fill upsampling buffer for low-pass filter: every second point
                    buffer3 = 0.0_rk
                    buffer3(1+io:nx-io:2) = u_wc(ix, 1+io:Bs(2)+2*g-io:2, iz, ic) ! SC
                    ! buffer3(1:ny:2) = u_wc(ix, 1:n(2), iz, ic)
                    call filter1dim(params, buffer3(1:ny), buffer1(1:ny), params%HR, lbound(params%HR,dim=1), ubound(params%HR,dim=1), skip_g=params%g, sampling=1)

                    ! fill upsampling buffer for high-pass filter: every second point
                    buffer3 = 0.0_rk
                    buffer3(1+io:nx-io:2) = u_wc(ix, 2+io:Bs(2)+2*g-io:2, iz, ic) ! WC
                    ! buffer3(1:ny:2) = u_wc(ix, n(2)+1:2*n(2), iz, ic)
                    call filter1dim(params, buffer3(1:ny), buffer2(1:ny), params%GR, lbound(params%GR,dim=1), ubound(params%GR,dim=1), skip_g=params%g, sampling=1)

                    u_wc(ix, :, iz, ic) = buffer1(1:ny) + buffer2(1:ny)
                enddo
            enddo
        enddo

        !!!!!!!!!!!!!!!!!
        if (nz==1) return
        !!!!!!!!!!!!!!!!!

        ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        do ic = 1, nc
            do ix = (g+1), (Bs(1)+g)  ! ignore ghost points
                do iy = (g+1), (Bs(2)+g)  ! ignore ghost points
                    ! fill upsampling buffer for low-pass filter: every second point
                    buffer3 = 0.0_rk
                    buffer3(1+io:nx-io:2) = u_wc(ix, iy, 1+io:Bs(3)+2*g-io:2, ic) ! SC
                    call filter1dim(params, buffer3(1:nz), buffer1(1:nz), params%HR, lbound(params%HR,dim=1), ubound(params%HR,dim=1), skip_g=params%g, sampling=1)

                    ! fill upsampling buffer for high-pass filter: every second point
                    buffer3 = 0.0_rk
                    buffer3(1+io:nx-io:2) = u_wc(ix, iy, 2+io:Bs(3)+2*g-io:2, ic) ! WC
                    call filter1dim(params, buffer3(1:nz), buffer2(1:nz), params%GR, lbound(params%GR,dim=1), ubound(params%GR,dim=1), skip_g=params%g, sampling=1)

                    u_wc(ix, iy, :, ic) = buffer1(1:nz) + buffer2(1:nz)
                enddo
            enddo
        enddo

    end subroutine






end module
