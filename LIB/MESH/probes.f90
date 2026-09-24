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
    if (params%n_probes > 0) then
        if (any(params%probe_x < 0.0_rk) .or. any(params%probe_x > params%domain_size(1))) all_probes_ok = .false.
        if (any(params%probe_y < 0.0_rk) .or. any(params%probe_y > params%domain_size(2))) all_probes_ok = .false.
        if (params%dim == 3) then
            if (any(params%probe_z < 0.0_rk) .or. any(params%probe_z > params%domain_size(3))) all_probes_ok = .false.
        endif
    endif
    if (params%n_probe_lines > 0) then
        if (any(params%probe_line_x1 < 0.0_rk) .or. any(params%probe_line_x1 > params%domain_size(1))) all_probes_ok = .false.
        if (any(params%probe_line_x2 < 0.0_rk) .or. any(params%probe_line_x2 > params%domain_size(1))) all_probes_ok = .false.
        if (any(params%probe_line_y1 < 0.0_rk) .or. any(params%probe_line_y1 > params%domain_size(2))) all_probes_ok = .false.
        if (any(params%probe_line_y2 < 0.0_rk) .or. any(params%probe_line_y2 > params%domain_size(2))) all_probes_ok = .false.
        if (params%dim == 3) then
            if (any(params%probe_line_z1 < 0.0_rk) .or. any(params%probe_line_z1 > params%domain_size(3))) all_probes_ok = .false.
            if (any(params%probe_line_z2 < 0.0_rk) .or. any(params%probe_line_z2 > params%domain_size(3))) all_probes_ok = .false.
        endif
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
            write(iu, '(A)', advance='no') '% time'
            ! write entries for every variable as "probeID:varname", right aligned in 15 characters
            do ip = 1, params%n_probes
                do iv = 1, params%N_probe_variables
                    write(column_format, '(A,I0,A)') '(A, i0.', digits, ',A,A)'
                    write(column_name, column_format) 'Point', ip, ':', trim(adjustl(params%probe_variables(iv)))
                    write(iu, '(A,A)', advance='no') tfile_separator, trim(adjustl(column_name))
                enddo
            enddo
            ! write entries for every line probe point as "lineID:probeID:varname", right aligned in 15 characters
            do il = 1, params%n_probe_lines
                do ip = 1, params%probe_line_npoints(il)
                    do iv = 1, params%N_probe_variables
                        write(column_format, '(A,I0,A,I0,A)') '(A, i0.', digits_lines, ',A,i0.', digits_line_points, ',A,A)'
                        write(column_name, column_format) 'Line', il, ':Point', ip, ':', trim(adjustl(params%probe_variables(iv)))
                        write(iu, '(A,A)', advance='no') tfile_separator, trim(adjustl(column_name))
                    enddo
                enddo
            enddo
            ! end line
            write(iu,'(A)') ""

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

    integer(kind=ik) :: k, hvy_id, lgt_id, probe_var_0, probe_var_E, tmp_size
    real(kind=rk) :: x0(3), dx(3)
    integer(kind=2) :: n_domain(3)

    integer, save :: iu_local = -1  !< unit for writing probes on this rank; only rank 0 will have it, but we keep track of it to avoid reopening the file multiple times
    logical, save :: have_unit = .false.  !< whether this rank has a unit open for writing probes; only rank 0 will have it, but we keep track of it to avoid reopening the file multiple times
    integer(kind=ik), save :: write_counter = 0_ik  !< counts how many lines have been written since last flush, to control flushing frequency
    integer(kind=ik) :: il, ip, iv, pt, n_points_total, nvars_batch  !< loop indices for probe lines, points and variables
    real(kind=rk), allocatable, save :: vals(:)  !< array with probe values, one-dimensionalized
    real(kind=rk), allocatable, save :: xq_all(:,:)  !< query point coordinates for all probes and line points
    real(kind=rk), allocatable, save :: fq(:,:)  !< interpolated values, one column per query point, sized for the largest variable batch
    character(len=cshort) :: interpolationMethod
    real(kind=rk) :: t0
    character(len=cshort) :: write_format

    n_points_total = params%n_probes + sum(params%probe_line_npoints(1:params%n_probe_lines))

    if (.not. allocated(vals)) allocate(vals(1:n_points_total * params%N_probe_variables))
    if (.not. allocated(xq_all)) allocate(xq_all(1:3, 1:n_points_total))
    if (.not. allocated(fq)) allocate(fq(1:size(hvy_tmp, 4), 1:n_points_total))

    select case (params%probe_interpolation_order)
    case (0)
        interpolationMethod = "floor"
    case (1)
        interpolationMethod = "linear"
    case (2)
        interpolationMethod = "delta"
    end select

    ! assemble the query points once: first all single probes, then all line probes
    ! (this order matches the layout of "vals" and thus of the probes.t output file)
    xq_all = 0.0_rk
    do ip = 1, params%n_probes
        xq_all(1, ip) = params%probe_x(ip)
        xq_all(2, ip) = params%probe_y(ip)
        xq_all(3, ip) = params%probe_z(ip)
    enddo
    do il = 1, params%n_probe_lines
        do ip = 1, params%probe_line_npoints(il)
            pt = params%n_probes + sum(params%probe_line_npoints(1:il-1)) + ip
            xq_all(1, pt) = params%probe_line_x1(il) + real(ip-1, rk) * (params%probe_line_x2(il) - params%probe_line_x1(il)) / real(params%probe_line_npoints(il)-1, rk)
            xq_all(2, pt) = params%probe_line_y1(il) + real(ip-1, rk) * (params%probe_line_y2(il) - params%probe_line_y1(il)) / real(params%probe_line_npoints(il)-1, rk)
            xq_all(3, pt) = params%probe_line_z1(il) + real(ip-1, rk) * (params%probe_line_z2(il) - params%probe_line_z1(il)) / real(params%probe_line_npoints(il)-1, rk)
        enddo
    enddo

    if (params%physics_type /= 'ACM-new' .and. params%physics_type /= 'NSPP') then
        call abort(2505206, 'Probes currently implemented for ACM-new and NSPP only')
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

    ! it could be that we would like to save more probe variables than we have available slots in hvy_tmp
    ! in that case we do several loops
    tmp_size = size(hvy_tmp, 4)
    t0 = MPI_Wtime()
    n_domain = 0
    ! all processors initialize all values to -Inf, then when we do max reduction, any value that is written will dominate and be populated for rank0 to write it
    vals = -huge(1.0_rk)
    do probe_var_0 = 1, params%N_probe_variables, tmp_size
        ! this is the end index for this batch of probe variables
        probe_var_E = min(probe_var_0 + tmp_size - 1, params%N_probe_variables)

        ! compute requested probe variables on the grid into hvy_tmp
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            call get_block_spacing_origin(params, lgt_id, x0, dx)
            if (.not. all(params%periodic_BC)) call get_adjacent_boundary_surface_normal(params, lgt_id, n_domain)

            call PREPARE_SAVE_DATA_meta(params%physics_type, time, hvy_block(:,:,:,:,hvy_id), &
                params%g, x0, dx, hvy_tmp(:,:,:,:,hvy_id), &
                hvy_mask(:,:,:,:,merge(1, hvy_id, size(hvy_mask,5) == 1)), n_domain, names_override=params%probe_variables(probe_var_0:probe_var_E))
        enddo
        call toc("probes_wrapper (prepare variables on block)", 97, MPI_Wtime()-t0)

        ! We only need ghost values for the compact interpolation stencils used by p=1 and p=2.
        ! p=0 samples the lower grid neighbor directly so it's always on the grid
        if (params%probe_interpolation_order > 0) then
            t0 = MPI_Wtime()
            call sync_ghosts_RHS_tree(params, hvy_tmp(:,:,:,1:probe_var_E-probe_var_0+1,:), tree_ID, g_minus=max(merge(1,3, params%probe_interpolation_order == 1), params%g_RHS), g_plus=max(merge(1,3, params%probe_interpolation_order == 1), params%g_RHS))
            call toc("probes_wrapper (sync)", 95, MPI_Wtime()-t0)
        endif

        ! interpolate this batch of variables at all query points (probes + line points) at once.
        ! Ghost nodes have already been synced above according to probe_interpolation_order, so
        ! we can skip the (redundant) sync inside interpolatePointCloud_tree.
        t0 = MPI_Wtime()

        nvars_batch = probe_var_E - probe_var_0 + 1

        call interpolatePointCloud_tree(params, hvy_tmp(:,:,:,1:nvars_batch,:), tree_ID, xq_all, fq(1:nvars_batch,:), interpolationMethod, sync=.false.)

        do pt = 1, n_points_total
            do iv = probe_var_0, probe_var_E
                vals((pt-1)*params%N_probe_variables + iv) = fq(iv-probe_var_0+1, pt)
            enddo
        enddo
    enddo

    call toc("probes_wrapper (probe interpolation)", 98, MPI_Wtime()-t0)

    if (params%rank == 0) then
        if (.not. have_unit) then
            ! We rely on the fact that init_probes_file opened the file with a unit on rank 0.
            ! Unfortunately newunit is not globally accessible here; reopen in append mode.
            open(newunit=iu_local, file='probes.t', status='old', position='append', action='write')
            have_unit = .true.
        endif

        write(write_format, '(A,I0,A)') '(ES15.8, ', params%n_probes * params%N_probe_variables + sum(params%probe_line_npoints(1:params%n_probe_lines)) * params%N_probe_variables, '(";", ES15.8))'
        write(iu_local, write_format) time, vals(1:params%n_probes * params%N_probe_variables + sum(params%probe_line_npoints(1:params%n_probe_lines)) * params%N_probe_variables)

        write_counter = write_counter + 1
        if (modulo(write_counter, flush_frequency) == 0) call flush(iu_local)
    endif

end subroutine probes_wrapper
