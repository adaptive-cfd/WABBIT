! Probe sampling and output (wide ASCII file probes.t)
!
! This file is included into module_mesh.

subroutine init_probes_file(params, overwrite)
    implicit none
    type(type_params), intent(inout) :: params
    logical, intent(in) :: overwrite

    logical, save :: initialized = .false.
    integer, save :: iu = -1
    integer(kind=ik), save :: digits = 1_ik, digits_lines = 1_ik, digits_line_points = 1_ik

    logical :: enabled, exists
    integer(kind=ik) :: required_g
    integer(kind=ik) :: p, il, ip, iv
    integer(kind=ik) :: mpierr
    logical :: all_probes_ok
    character(len=cshort) :: column_name, column_format

    enabled = ((params%n_probes > 0 .or. params%n_probe_lines > 0) .and. params%N_probe_variables > 0) .and. &
        (params%nsave_probes /= 99999999_ik .or. abs(params%tsave_probes-9999999.9_rk) > 1e-12_rk)

    if (.not. enabled) return
    if (initialized) return

    if (params%physics_type == 'navier_stokes') then
        call abort(260520, 'ERROR: Probes currently not implemented for navier_stokes')
    endif

    if (params%probe_interpolation_order < 0 .or. params%probe_interpolation_order > 2) then
        call abort(260520, 'ERROR: probe_interpolation_order must be 0, 1, or 2')
    endif

    p = params%probe_interpolation_order
    select case (p)
    case (0)
        required_g = 0
    case (1)
        required_g = 1
    case (2)
        required_g = 3
    end select
    if (params%g < required_g) then
        call abort(260520, 'ERROR: Not enough ghost nodes for requested probe_interpolation_order')
    endif

    ! validate probe locations against domain [0, domain_size]
    all_probes_ok = .true.
    if (any(params%probe_x < 0.0_rk) .or. any(params%probe_x > params%domain_size(1))) all_probes_ok = .false.
    if (any(params%probe_y < 0.0_rk) .or. any(params%probe_y > params%domain_size(2))) all_probes_ok = .false.
    if (params%dim == 3) then
        if (any(params%probe_z < 0.0_rk) .or. any(params%probe_z > params%domain_size(3))) all_probes_ok = .false.
    endif

    if (.not. all_probes_ok) then
        call abort(260520, 'ERROR: At least one probe point is outside domain [0, domain_size].')
    endif

    ! number of digits for probe id
    digits = 1
    if (params%n_probes >= 10) then
        digits = int(log10(real(params%n_probes, rk)), kind=ik) + 1_ik
    endif
    digits_lines = 1
    if (params%n_probe_lines >= 10) digits_lines = int(log10(real(params%n_probe_lines, rk)), kind=ik) + 1_ik
    digits_line_points = 1
    if (any(params%probe_line_npoints(1:params%n_probe_lines) >= 10)) then
        digits_line_points = int(log10(real(maxval(params%probe_line_npoints(1:params%n_probe_lines)), rk)), kind=ik) + 1_ik
    endif

    if (params%rank == 0) then
        inquire(file='probes.t', exist=exists)
        if (overwrite .or. .not. exists) then
            open(newunit=iu, file='probes.t', status='replace', action='write')

            ! write time
            write(iu, '(A15)', advance='no') '%          time'
            ! write entries for every variable as "probeID:varname", right aligned in 15 characters
            do ip = 1, params%n_probes
                do iv = 1, params%N_probe_variables
                    write(column_format, '(A,I0,A)') '(i0.', digits, ',A,A)'
                    write(column_name, column_format) ip, ':', trim(adjustl(params%probe_variables(iv)))
                    write(iu, '(1x, A15)', advance='no') column_name(1:min(len_trim(column_name), 15))
                enddo
            enddo
            ! write entries for every line probe point as "lineID:probeID:varname", right aligned in 15 characters
            do il = 1, params%n_probe_lines
                do ip = 1, params%probe_line_npoints(il)
                    do iv = 1, params%N_probe_variables
                        write(column_format, '(A,I0,A,I0,A)') '(i0.', digits_lines, ',A,i0.', digits_line_points, ',A,A)'
                        write(column_name, column_format) il, ':', ip, ':', trim(adjustl(params%probe_variables(iv)))
                        write(iu, '(1x, A15)', advance='no') column_name(1:min(len_trim(column_name), 15))
                    enddo
                enddo
            enddo
            ! end line
            write(iu,*)

            close(iu)
            iu = -1
        else
            ! nothing to do here; we will append lazily when first writing data
        endif
    endif

    initialized = .true.

end subroutine init_probes_file


subroutine finalize_probes_file(params)
    implicit none
    type(type_params), intent(inout) :: params

    logical, save :: initialized = .false.
    integer, save :: iu = -1

    ! This is a no-op placeholder: we keep the file open for performance.
    ! The OS will close it at program termination.
    ! If you want explicit close, make init_probes_file expose its unit.

end subroutine finalize_probes_file


logical function it_is_time_to_probe(time, iteration, params)
    implicit none
    real(kind=rk), intent(in) :: time
    integer(kind=ik), intent(in) :: iteration
    type(type_params), intent(in) :: params

    real(kind=rk), parameter :: tol = 1e-12_rk
    logical :: enabled

    enabled = ((params%n_probes > 0 .or. params%n_probe_lines > 0) .and. params%N_probe_variables > 0) .and. &
        (params%nsave_probes /= 99999999_ik .or. abs(params%tsave_probes-9999999.9_rk) > 1e-12_rk)

    it_is_time_to_probe = .false.
    if (.not. enabled) return
    if (time + tol < params%probe_start_time) return

    if (params%nsave_probes /= 99999999_ik) then
        if (modulo(iteration, params%nsave_probes) == 0) then
            it_is_time_to_probe = .true.
            return
        endif
    endif

    if (abs(params%tsave_probes-9999999.9_rk) > 1e-12_rk) then
        if (abs(mod(time-params%probe_start_time, params%tsave_probes)) < tol .or. &
            abs(mod(time-params%probe_start_time, params%tsave_probes)-params%tsave_probes) < tol) then
            it_is_time_to_probe = .true.
            return
        endif
    endif

end function it_is_time_to_probe


subroutine probes_wrapper(time, params, hvy_block, hvy_tmp, hvy_mask, tree_ID)
    implicit none

    real(kind=rk), intent(in) :: time
    type(type_params), intent(inout) :: params
    real(kind=rk), intent(inout) :: hvy_block(:, :, :, :, :)
    real(kind=rk), intent(inout) :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), intent(inout) :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), intent(in) :: tree_ID

    integer(kind=ik) :: k, hvy_id, lgt_id
    real(kind=rk) :: x0(3), dx(3)
    integer(kind=2) :: n_domain(3)

    integer, save :: iu_local = -1  !< unit for writing probes on this rank; only rank 0 will have it, but we keep track of it to avoid reopening the file multiple times
    logical, save :: have_unit = .false.  !< whether this rank has a unit open for writing probes; only rank 0 will have it, but we keep track of it to avoid reopening the file multiple times
    integer(kind=ik), save :: write_counter = 0_ik  !< counts how many lines have been written since last flush, to control flushing frequency
    integer(kind=ik) :: il, ip, iv  !< loop indices for probe lines, points and variables
    real(kind=rk), allocatable, save :: vals(:)  !< array with probe values, one-dimensionalized
    real(kind=rk) :: xq(3)  !< query point coordinates
    real(kind=rk) :: t0
    integer(kind=ik) :: mpierr
    character(len=cshort) :: write_format

    if (.not. allocated(vals)) allocate(vals(1:params%n_probes * params%N_probe_variables + sum(params%probe_line_npoints(1:params%n_probe_lines)) * params%N_probe_variables))

    if (params%physics_type /= 'ACM-new') then
        call abort(2505206, 'Probes currently implemented for ACM-new only')
    endif

    ! we need to be sure that variables are synced - RHS sync is enough for derivative quantities
    t0 = MPI_Wtime()
    call sync_ghosts_RHS_tree(params, hvy_block, tree_ID)
    call toc("probes_wrapper (sync)", 95, MPI_Wtime()-t0)

    ! mask might be needed for probes
    t0 = MPI_Wtime()
    call createMask_tree(params, time, hvy_mask, hvy_tmp)
    ! extra preparatory step. Some things have to be prepared on grid-level and not on block-level. This is especially important when deriving any quantities by solving the poisson-equation
    if (params%physics_type == "NSPP") then
        call pressure_from_velocity(params, time, hvy_block, hvy_tmp, hvy_mask, tree_ID_flow)
    endif
    call toc("probes_wrapper (prepare variables on grid-level)", 96, MPI_Wtime()-t0)

    t0 = MPI_Wtime()
    n_domain = 0
    ! compute requested probe variables on the grid into hvy_tmp
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin(params, lgt_id, x0, dx)
        if (.not. all(params%periodic_BC)) call get_adjacent_boundary_surface_normal(params, lgt_id, n_domain)

        call PREPARE_SAVE_DATA_meta(params%physics_type, time, hvy_block(:,:,:,:,hvy_id), &
            params%g, x0, dx, hvy_tmp(:,:,:,:,hvy_id), &
            hvy_mask(:,:,:,:,merge(1, hvy_id, size(hvy_mask,5) == 1)), n_domain, names_override=params%probe_variables)
    enddo
    call toc("probes_wrapper (prepare variables on block)", 97, MPI_Wtime()-t0)

    ! We only need ghost values for the compact interpolation stencils used by p=1 and p=2.
    ! p=0 samples the lower grid neighbor directly so it's always on the grid
    if (params%probe_interpolation_order > 0) then
        t0 = MPI_Wtime()
        call sync_ghosts_RHS_tree(params, hvy_tmp(:,:,:,1:params%N_probe_variables,:), tree_ID, g_minus=max(merge(1,3, params%probe_interpolation_order == 1), params%g_RHS), g_plus=max(merge(1,3, params%probe_interpolation_order == 1), params%g_RHS))
        call toc("probes_wrapper (sync)", 95, MPI_Wtime()-t0)
    endif

    ! all processors set all values to -Inf, then when we do max reduction, any value that is written will dominate and be populated for rank0 to write it
    vals = -huge(1.0_rk)

    ! loop over all probe points
    t0 = MPI_Wtime()
    do ip = 1, params%n_probes
        xq = 0.0_rk
        xq(1) = params%probe_x(ip)
        xq(2) = params%probe_y(ip)
        xq(3) = params%probe_z(ip)

        ! find a local block that contains the query point
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            call get_block_spacing_origin(params, lgt_id, x0, dx)

            if (point_in_block(params, xq, x0, dx)) then
                do iv = 1, params%N_probe_variables
                    vals((ip-1)*params%N_probe_variables + iv) = interpolate_probe_tensor(params, hvy_tmp(:,:,:,iv,hvy_id), xq, x0, dx, params%probe_interpolation_order)
                enddo
                exit
            endif
        enddo
    enddo

    ! loop over all line probe points
    do il = 1, params%n_probe_lines
        do ip = 1, params%probe_line_npoints(il)
            xq = 0.0_rk
            xq(1) = params%probe_line_x1(il) + real(ip-1, rk) * (params%probe_line_x2(il) - params%probe_line_x1(il)) / real(params%probe_line_npoints(il)-1, rk)
            xq(2) = params%probe_line_y1(il) + real(ip-1, rk) * (params%probe_line_y2(il) - params%probe_line_y1(il)) / real(params%probe_line_npoints(il)-1, rk)
            xq(3) = params%probe_line_z1(il) + real(ip-1, rk) * (params%probe_line_z2(il) - params%probe_line_z1(il)) / real(params%probe_line_npoints(il)-1, rk)

            ! find a local block that contains the query point
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                call get_block_spacing_origin(params, lgt_id, x0, dx)

                if (point_in_block(params, xq, x0, dx)) then
                    do iv = 1, params%N_probe_variables
                        vals((params%n_probes + sum(params%probe_line_npoints(1:il-1)) + ip-1)*params%N_probe_variables + iv) = &
                            interpolate_probe_tensor(params, hvy_tmp(:,:,:,iv,hvy_id), xq, x0, dx, params%probe_interpolation_order)
                    enddo
                    exit
                endif
            enddo
        enddo
    enddo

    call MPI_Allreduce(MPI_IN_PLACE, vals, size(vals), MPI_REAL8, MPI_MAX, WABBIT_COMM, mpierr)
    call toc("probes_wrapper (probe interpolation)", 98, MPI_Wtime()-t0)

    if (params%rank == 0) then
        if (.not. have_unit) then
            ! We rely on the fact that init_probes_file opened the file with a unit on rank 0.
            ! Unfortunately newunit is not globally accessible here; reopen in append mode.
            open(newunit=iu_local, file='probes.t', status='old', position='append', action='write')
            have_unit = .true.
        endif

        write(write_format, '(A,I0,A)') '(ES15.8, ', params%n_probes * params%N_probe_variables + sum(params%probe_line_npoints(1:params%n_probe_lines)) * params%N_probe_variables, '(",", ES15.8))'
        write(iu_local, write_format) time, vals(1:params%n_probes * params%N_probe_variables + sum(params%probe_line_npoints(1:params%n_probe_lines)) * params%N_probe_variables)

        write_counter = write_counter + 1
        if (modulo(write_counter, flush_frequency) == 0) call flush(iu_local)
    endif

contains

    logical function point_in_block(params_local, xq_local, x0_local, dx_local)
        implicit none
        type(type_params), intent(in) :: params_local
        real(kind=rk), intent(in) :: xq_local(3)
        real(kind=rk), intent(in) :: x0_local(3)
        real(kind=rk), intent(in) :: dx_local(3)

        real(kind=rk) :: xmin, xmax
        integer(kind=ik) :: d

        point_in_block = .true.
        do d = 1, params_local%dim
            xmin = x0_local(d)
            xmax = x0_local(d) + dx_local(d) * real(params_local%Bs(d), rk)
            if (xq_local(d) < xmin - 1e-14_rk .or. xq_local(d) > xmax + 1e-14_rk) then
                point_in_block = .false.
                return
            endif
        enddo
    end function point_in_block


    real(kind=rk) function interpolate_probe_tensor(params_local, f, xq_local, x0_local, dx_local, p)
        implicit none
        type(type_params), intent(in) :: params_local   !< params
        real(kind=rk), intent(in) :: f(:,:,:)           !< block data
        real(kind=rk), intent(in) :: xq_local(3)        !< query point coordinates
        real(kind=rk), intent(in) :: x0_local(3)        !< block origin coordinates
        real(kind=rk), intent(in) :: dx_local(3)        !< block spacing
        integer(kind=ik), intent(in) :: p               !< probe interpolation order: 0=floor-neighbor, 1=linear kernel, 2=delta kernel

        integer(kind=ik) :: ix0, iy0, iz0
        integer(kind=ik) :: ix, iy, iz
        integer(kind=ik) :: support
        integer(kind=ik) :: g_local
        real(kind=rk) :: x_grid, y_grid, z_grid
        real(kind=rk) :: wx, wy, wz

        g_local = params_local%g

        ! compute center point index as lower nearest point
        ix0 = int(floor((xq_local(1) - x0_local(1)) / dx_local(1)), kind=ik) + g_local + 1
        iy0 = int(floor((xq_local(2) - x0_local(2)) / dx_local(2)), kind=ik) + g_local + 1
        if (params_local%dim == 3) then
            iz0 = int(floor((xq_local(3) - x0_local(3)) / dx_local(3)), kind=ik) + g_local + 1
        else
            iz0 = 1
        endif
        interpolate_probe_tensor = 0.0_rk

        select case (p)
        case (0)
            ! Floor-neighbor sampling: take the lower grid point in each active
            ! direction. This avoids any need for ghost synchronization.
            interpolate_probe_tensor = f(ix0, iy0, iz0)

        case (1)
            ! Compact linear interpolation. Each axis contributes a tent kernel
            ! with support radius 1 grid spacing, so only the immediate neighbors
            ! around the query point can contribute.
            support = 1

            do iz = iz0 - merge(0, support, params_local%dim==2), iz0 + merge(0, support, params_local%dim==2)
                z_grid = x0_local(3) + real(iz - (g_local + 1), rk) * dx_local(3)
                if (params_local%dim == 2) then
                    wz = 1.0_rk
                else
                    wz = linear_interpolation(xq_local(3) - z_grid, dx_local(3))
                endif
                do iy = iy0 - support, iy0 + support
                    y_grid = x0_local(2) + real(iy - (g_local + 1), rk) * dx_local(2)
                    wy = linear_interpolation(xq_local(2) - y_grid, dx_local(2))
                    do ix = ix0 - support, ix0 + support
                        x_grid = x0_local(1) + real(ix - (g_local + 1), rk) * dx_local(1)
                        wx = linear_interpolation(xq_local(1) - x_grid, dx_local(1))
                        interpolate_probe_tensor = interpolate_probe_tensor + wx * wy * wz * f(ix, iy, iz)
                    enddo
                enddo
            enddo

        case (2)
            ! Compact delta interpolation. This uses the smoother discrete delta
            ! kernel with support radius 3 grid spacings, so we sum over a wider
            ! tensor-product stencil than for linear interpolation.
            support = 3

            do iz = iz0 - merge(0, support, params_local%dim==2), iz0 + merge(0, support, params_local%dim==2)
                z_grid = x0_local(3) + real(iz - (g_local + 1), rk) * dx_local(3)
                if (params_local%dim == 2) then
                    wz = 1.0_rk
                else
                    wz = delta_interpolation(xq_local(3) - z_grid, dx_local(3))
                endif
                do iy = iy0 - support, iy0 + support
                    y_grid = x0_local(2) + real(iy - (g_local + 1), rk) * dx_local(2)
                    wy = delta_interpolation(xq_local(2) - y_grid, dx_local(2))
                    do ix = ix0 - support, ix0 + support
                        x_grid = x0_local(1) + real(ix - (g_local + 1), rk) * dx_local(1)
                        wx = delta_interpolation(xq_local(1) - x_grid, dx_local(1))
                        interpolate_probe_tensor = interpolate_probe_tensor + wx * wy * wz * f(ix, iy, iz)
                    enddo
                enddo
            enddo

        case default
            call abort(2505208, 'probe_interpolation_order must be 0, 1, or 2')
        end select

    end function interpolate_probe_tensor

end subroutine probes_wrapper
