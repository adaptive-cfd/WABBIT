!-------------------------------------------------------------------------------
!> 
!! Periodic, level-0 wavelet decomposition for one sub-level.
!!
!! This routine applies one decomposition pass along x/y(/z) on the interior block
!! only, assuming level-0 periodicity. Instead of reading ghost nodes, each filter
!! tap index is mapped back into the interior domain via periodic projection.
!!
!! The selected sub-level controls the filter dilation:
!!   stride = 2^(sub_level)
!! so higher sub-levels reuse the same filter coefficients at larger spacing.
!!
!! For sub_level>0, only the parent SC lattice is decomposed (consecutive decomposition).
!! Existing coarser WC coefficients are left untouched.
!!
!! Example (1D, Bs=8): decomposition from L-1 to L-2 (sub_level=2)
!!   indices:  1  2  3  4  5  6  7  8
!!   before :  S1 W0 W1 W0 S1 W0 W1 W0
!!   after  :  S2 W0 W1 W0 W2 W0 W1 W0
!!   (only parent-S1 positions 1 and 5 are decomposed; W1 and W0 are untouched)
!! Here, S is SC, W is WC, and the subscript is the level (0 for L0, 1 for L-1, etc.)
!-------------------------------------------------------------------------------
subroutine sub0_wavelet_decompose_sublevel(params, u, u_d, sub_level)
    use module_globals
    use module_params

    implicit none

    type(type_params), intent(in) :: params
    real(kind=rk), intent(in)     :: u(:, :, :, :)
    real(kind=rk), intent(inout)  :: u_d(:, :, :, :)
    integer(kind=ik), intent(in)  :: sub_level

    integer(kind=ik) :: nx, ny, nz, nc
    integer(kind=ik) :: ix, iy, iz, ic
    integer(kind=ik) :: ix0, ix1, iy0, iy1, iz0, iz1
    real(kind=rk), allocatable :: line_in(:), line_out(:)

    ! basic contract checks
    if (sub_level < 1) call abort(260001_ik, "sub0_wavelet_decompose_sublevel: sub_level must be >= 1")
    if (.not. allocated(params%HD)) call abort(260002_ik, "sub0_wavelet_decompose_sublevel: setup_wavelet missing HD")
    if (.not. allocated(params%GD)) call abort(260003_ik, "sub0_wavelet_decompose_sublevel: setup_wavelet missing GD")

    nx = size(u, 1)
    ny = size(u, 2)
    nz = size(u, 3)
    nc = size(u, 4)

    if (size(u_d, 1) /= nx .or. size(u_d, 2) /= ny .or. size(u_d, 3) /= nz .or. size(u_d, 4) /= nc) then
        call abort(260004_ik, "sub0_wavelet_decompose_sublevel: u and u_d must have same shape")
    endif

    ! interior index ranges (no ghost update here)
    ix0 = params%g + 1
    ix1 = params%g + params%Bs(1)
    iy0 = params%g + 1
    iy1 = params%g + params%Bs(2)
    iz0 = 1
    iz1 = 1
    if (params%dim == 3) then
        iz0 = params%g + 1
        iz1 = params%g + params%Bs(3)
    endif

    ! start from input and overwrite interior progressively dimension by dimension
    u_d = u

    allocate(line_in(1:max(nx, ny, nz)), line_out(1:max(nx, ny, nz)))

    ! X-direction pass
    do ic = 1, nc
        do iz = iz0, iz1
            do iy = iy0, iy1
                line_in(1:nx) = u_d(1:nx, iy, iz, ic)
                call apply_periodic_decompose_line(line_in, line_out, ix0, ix1, sub_level, params%HD, params%GD)
                u_d(ix0:ix1, iy, iz, ic) = line_out(ix0:ix1)
            enddo
        enddo
    enddo

    ! Y-direction pass
    do ic = 1, nc
        do iz = iz0, iz1
            do ix = ix0, ix1
                line_in(1:ny) = u_d(ix, 1:ny, iz, ic)
                call apply_periodic_decompose_line(line_in, line_out, iy0, iy1, sub_level, params%HD, params%GD)
                u_d(ix, iy0:iy1, iz, ic) = line_out(iy0:iy1)
            enddo
        enddo
    enddo

    ! Z-direction pass (3D only)
    if (params%dim == 3) then
        do ic = 1, nc
            do iy = iy0, iy1
                do ix = ix0, ix1
                    line_in(1:nz) = u_d(ix, iy, 1:nz, ic)
                    call apply_periodic_decompose_line(line_in, line_out, iz0, iz1, sub_level, params%HD, params%GD)
                    u_d(ix, iy, iz0:iz1, ic) = line_out(iz0:iz1)
                enddo
            enddo
        enddo
    endif

    deallocate(line_in, line_out)

contains

    !-----------------------------------------------------------------------
    ! Apply one periodic 1D decomposition pass on a line segment [i0:i1].
    ! Only parent-SC positions are updated; parent-WC positions remain untouched.
    !-----------------------------------------------------------------------
    subroutine apply_periodic_decompose_line(line_i, line_o, i0, i1, level_local, HD, GD)
        real(kind=rk), intent(in)    :: line_i(:)
        real(kind=rk), intent(inout) :: line_o(:)
        integer(kind=ik), intent(in) :: i0, i1, level_local
        real(kind=rk), intent(in)    :: HD(:), GD(:)

        integer(kind=ik) :: i, j, idx, stride, parent_step, child_step
        real(kind=rk) :: acc

        ! dilated-filter spacing for this sub-level
        stride = 2_ik**(level_local)
        parent_step = stride
        child_step = 2_ik * parent_step
        line_o = line_i

        do i = i0, i1
            ! only parent-SC lattice is decomposed at this stage
            if (modulo(i-i0, parent_step) /= 0) cycle

            acc = 0.0_rk
            if (modulo(i-i0, child_step) == 0) then
                do j = lbound(HD, 1), ubound(HD, 1)
                    ! periodic projection instead of ghost-node access
                    idx = periodic_project_index(i + stride * j, i0, i1)
                    acc = acc + line_i(idx) * HD(j)
                enddo
            else
                do j = lbound(GD, 1), ubound(GD, 1)
                    ! periodic projection instead of ghost-node access
                    idx = periodic_project_index(i + stride * j, i0, i1)
                    acc = acc + line_i(idx) * GD(j)
                enddo
            endif
            line_o(i) = acc
        enddo
    end subroutine

    !-----------------------------------------------------------------------
    ! Map any integer index back to the periodic interior interval [i0:i1].
    !-----------------------------------------------------------------------
    pure integer(kind=ik) function periodic_project_index(i, i0, i1)
        integer(kind=ik), intent(in) :: i, i0, i1
        integer(kind=ik) :: n

        n = i1 - i0 + 1
        periodic_project_index = i0 + modulo(i - i0, n)
    end function

end subroutine


!-------------------------------------------------------------------------------
!>
!! Periodic, level-0 wavelet reconstruction for one sub-level.
!!
!! Reconstruction is performed dimension-wise on the interior block with periodic
!! projection of filter accesses. For each target location, HR contributions are
!! taken from SC-parity source points and GR contributions from WC-parity points.
!!
!! The selected sub-level controls filter dilation via stride = 2^(sub_level).
!!
!! Example (1D, Bs=8): reconstruction from L-2 to L-1 (sub_level=2)
!!   indices:  1  2  3  4  5  6  7  8
!!   before :  S2 W0 W1 W0 W2 W0 W1 W0
!!   after  :  S1 W0 W1 W0 S1 W0 W1 W0
!!   (only parent-S1 positions 1 and 5 are reconstructed; W1 and W0)
!! Here, S is SC, W is WC, and the subscript is the level (0 for L0, 1 for L-1, etc.)
!-------------------------------------------------------------------------------
subroutine sub0_wavelet_reconstruct_sublevel(params, u_d, u_r, sub_level)
    use module_globals
    use module_params

    implicit none

    type(type_params), intent(in) :: params
    real(kind=rk), intent(in)     :: u_d(:, :, :, :)
    real(kind=rk), intent(inout)  :: u_r(:, :, :, :)
    integer(kind=ik), intent(in)  :: sub_level

    integer(kind=ik) :: nx, ny, nz, nc
    integer(kind=ik) :: ix, iy, iz, ic
    integer(kind=ik) :: ix0, ix1, iy0, iy1, iz0, iz1
    real(kind=rk), allocatable :: line_in(:), line_out(:)

    ! basic contract checks
    if (sub_level < 1) call abort(260011_ik, "sub0_wavelet_reconstruct_sublevel: sub_level must be >= 1")
    if (.not. allocated(params%HR)) call abort(260012_ik, "sub0_wavelet_reconstruct_sublevel: setup_wavelet missing HR")
    if (.not. allocated(params%GR)) call abort(260013_ik, "sub0_wavelet_reconstruct_sublevel: setup_wavelet missing GR")

    nx = size(u_d, 1)
    ny = size(u_d, 2)
    nz = size(u_d, 3)
    nc = size(u_d, 4)

    if (size(u_r, 1) /= nx .or. size(u_r, 2) /= ny .or. size(u_r, 3) /= nz .or. size(u_r, 4) /= nc) then
        call abort(260014_ik, "sub0_wavelet_reconstruct_sublevel: u_d and u_r must have same shape")
    endif

    ! interior index ranges (no ghost update here)
    ix0 = params%g + 1
    ix1 = params%g + params%Bs(1)
    iy0 = params%g + 1
    iy1 = params%g + params%Bs(2)
    iz0 = 1
    iz1 = 1
    if (params%dim == 3) then
        iz0 = params%g + 1
        iz1 = params%g + params%Bs(3)
    endif

    ! start from decomposed input and overwrite interior progressively
    u_r = u_d

    allocate(line_in(1:max(nx, ny, nz)), line_out(1:max(nx, ny, nz)))

    ! X-direction pass
    do ic = 1, nc
        do iz = iz0, iz1
            do iy = iy0, iy1
                line_in(1:nx) = u_r(1:nx, iy, iz, ic)
                call apply_periodic_reconstruct_line(line_in, line_out, ix0, ix1, sub_level, params%HR, params%GR)
                u_r(ix0:ix1, iy, iz, ic) = line_out(ix0:ix1)
            enddo
        enddo
    enddo

    ! Y-direction pass
    do ic = 1, nc
        do iz = iz0, iz1
            do ix = ix0, ix1
                line_in(1:ny) = u_r(ix, 1:ny, iz, ic)
                call apply_periodic_reconstruct_line(line_in, line_out, iy0, iy1, sub_level, params%HR, params%GR)
                u_r(ix, iy0:iy1, iz, ic) = line_out(iy0:iy1)
            enddo
        enddo
    enddo

    ! Z-direction pass (3D only)
    if (params%dim == 3) then
        do ic = 1, nc
            do iy = iy0, iy1
                do ix = ix0, ix1
                    line_in(1:nz) = u_r(ix, iy, 1:nz, ic)
                    call apply_periodic_reconstruct_line(line_in, line_out, iz0, iz1, sub_level, params%HR, params%GR)
                    u_r(ix, iy, iz0:iz1, ic) = line_out(iz0:iz1)
                enddo
            enddo
        enddo
    endif

    deallocate(line_in, line_out)

contains

    !-----------------------------------------------------------------------
    ! Apply one periodic 1D reconstruction pass on [i0:i1].
    ! Only parent-SC lattice points are reconstructed; parent-WC stays untouched.
    !-----------------------------------------------------------------------
    subroutine apply_periodic_reconstruct_line(line_i, line_o, i0, i1, level_local, HR, GR)
        real(kind=rk), intent(in)    :: line_i(:)
        real(kind=rk), intent(inout) :: line_o(:)
        integer(kind=ik), intent(in) :: i0, i1, level_local
        real(kind=rk), intent(in)    :: HR(:), GR(:)

        integer(kind=ik) :: i, j, idx, stride, parent_step, child_step
        real(kind=rk) :: acc

        ! dilated-filter spacing for this sub-level
        stride = 2_ik**(level_local)
        parent_step = stride
        child_step = 2_ik * parent_step
        line_o = line_i

        do i = i0, i1
            ! only reconstruct parent-SC lattice points
            if (modulo(i-i0, parent_step) /= 0) cycle

            acc = 0.0_rk

            do j = lbound(HR, 1), ubound(HR, 1)
                ! periodic projection instead of ghost-node access
                idx = periodic_project_index(i + stride * j, i0, i1)
                if (modulo(idx-i0, child_step) == 0) then
                    acc = acc + line_i(idx) * HR(j)
                endif
            enddo

            do j = lbound(GR, 1), ubound(GR, 1)
                ! periodic projection instead of ghost-node access
                idx = periodic_project_index(i + stride * j, i0, i1)
                if (modulo(idx-i0, child_step) == parent_step) then
                    acc = acc + line_i(idx) * GR(j)
                endif
            enddo

            line_o(i) = acc
        enddo
    end subroutine

    !-----------------------------------------------------------------------
    ! Map any integer index back to the periodic interior interval [i0:i1].
    !-----------------------------------------------------------------------
    pure integer(kind=ik) function periodic_project_index(i, i0, i1)
        integer(kind=ik), intent(in) :: i, i0, i1
        integer(kind=ik) :: n

        n = i1 - i0 + 1
        periodic_project_index = i0 + modulo(i - i0, n)
    end function

end subroutine


!-------------------------------------------------------------------------------
!>
!! Apply periodic decomposition passes on level 0 for a sub-level range.
!!
!! Pass k uses sub-level k for k=n_in...n_out (with automatic direction).
!! `u_io` is updated in-place, `u_work` is a required scratch array.
!-------------------------------------------------------------------------------
subroutine sub0_wavelet_decompose_n(params, u_io, n_in, n_out, u_work)
    use module_globals
    use module_params

    implicit none

    type(type_params), intent(in) :: params
    real(kind=rk), intent(inout)  :: u_io(:, :, :, :)
    integer(kind=ik), intent(in)  :: n_in, n_out
    real(kind=rk), intent(inout)  :: u_work(:, :, :, :)

    integer(kind=ik) :: k

    if (n_in < 0 .or. n_out <= n_in) call abort(260021_ik, "sub0_wavelet_decompose_n: n_in must be >= 0 and n_out must be > n_in")

    do k = n_in+1, n_out, 1
        call sub0_wavelet_decompose_sublevel(params, u_io, u_work, k)
        u_io = u_work
    enddo
end subroutine


!-------------------------------------------------------------------------------
!>
!! Apply periodic reconstruction passes on level 0 for a sub-level range.
!!
!! Typical inverse call is n_in=sub_max, n_out=1.
!! `u_io` is updated in-place, `u_work` is a required scratch array.
!-------------------------------------------------------------------------------
subroutine sub0_wavelet_reconstruct_n(params, u_io, n_in, n_out, u_work)
    use module_globals
    use module_params

    implicit none

    type(type_params), intent(in) :: params
    real(kind=rk), intent(inout)  :: u_io(:, :, :, :)
    integer(kind=ik), intent(in)  :: n_in, n_out
    real(kind=rk), intent(inout)  :: u_work(:, :, :, :)

    integer(kind=ik) :: k, dk

    if (n_out < 0 .or. n_out >= n_in) call abort(260031_ik, "sub0_wavelet_reconstruct_n: n_out must be >= 0 and n_out must be < n_in")

    do k = n_in, n_out+1, -1
        call sub0_wavelet_reconstruct_sublevel(params, u_io, u_work, k)
        u_io = u_work
    enddo
end subroutine


!-------------------------------------------------------------------------------
!>
!! Delete coefficients on a selected sub-level and type ("SC" or "WC").
!!
!! `level_local = 0` refers to L0 (SC0/WC0), `1` refers to L-1, etc.
!!
!! In d dimensions, level-local coefficient lattice points satisfy that each
!! coordinate belongs to {0, 2^level_local} modulo 2^(level_local+1) relative
!! to the interior start index. The all-zero combination is SC, all other
!! combinations on this lattice are WC components.
!!
!! Example (1D, Bs=8): delete WC on L-1 -> level_local=1, coeff_kind="WC"
!!   indices:  1  2  3  4  5  6  7  8
!!   before :  S1 W0 W1 W0 S1 W0 W1 W0
!!   after  :  S1 W0  . W0 S1 W0  . W0
!!
!! Example in L-2 context (same delete operation):
!!   indices:  1  2  3  4  5  6  7  8
!!   before :  S2 W0 W1 W0 W2 W0 W1 W0
!!   after  :  S2 W0  . W0 W2 W0  . W0
!!   (deleting WC on L-1 removes index 3 and 7 entries only)
!! Here, S is SC, W is WC, and the subscript is the level (0 for L0, 1 for L-1, etc.)
!-------------------------------------------------------------------------------
subroutine sub0_delete_coefficients(params, u_io, level_local, coeff_kind)
    use module_globals
    use module_params

    implicit none

    type(type_params), intent(in) :: params
    real(kind=rk), intent(inout)  :: u_io(:, :, :, :)
    integer(kind=ik), intent(in)  :: level_local
    character(len=*), intent(in)  :: coeff_kind

    integer(kind=ik) :: ix, iy, iz, ic
    integer(kind=ik) :: ix0, ix1, iy0, iy1, iz0, iz1
    integer(kind=ik) :: step_parent, step_child
    integer(kind=ik) :: px, py, pz
    logical :: on_level_lattice, is_sc

    if (level_local < 0) call abort(260041_ik, "sub0_delete_coefficients: level_local must be >= 0")
    if (coeff_kind /= "SC" .and. coeff_kind /= "WC") then
        call abort(260042_ik, "sub0_delete_coefficients: coeff_kind must be SC or WC")
    endif

    ix0 = params%g + 1
    ix1 = params%g + params%Bs(1)
    iy0 = params%g + 1
    iy1 = params%g + params%Bs(2)
    iz0 = 1
    iz1 = 1
    if (params%dim == 3) then
        iz0 = params%g + 1
        iz1 = params%g + params%Bs(3)
    endif

    step_parent = 2_ik**(level_local)
    step_child = 2_ik * step_parent

    do ic = 1, size(u_io, 4)
        do iz = iz0, iz1
            do iy = iy0, iy1
                do ix = ix0, ix1
                    px = modulo(ix-ix0, step_child)
                    py = modulo(iy-iy0, step_child)
                    pz = modulo(iz-iz0, step_child)

                    on_level_lattice = (px == 0_ik .or. px == step_parent) .and. (py == 0_ik .or. py == step_parent)
                    if (params%dim == 3) then
                        on_level_lattice = on_level_lattice .and. (pz == 0_ik .or. pz == step_parent)
                    endif

                    if (.not. on_level_lattice) cycle

                    is_sc = (px == 0_ik .and. py == 0_ik)
                    if (params%dim == 3) is_sc = is_sc .and. (pz == 0_ik)

                    if (coeff_kind == "SC" .and. is_sc) then
                        u_io(ix, iy, iz, ic) = 0.0_rk
                    elseif (coeff_kind == "WC" .and. .not. is_sc) then
                        u_io(ix, iy, iz, ic) = 0.0_rk
                    endif
                enddo
            enddo
        enddo
    enddo
end subroutine