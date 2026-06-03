!> \brief this is a prototype - it solves the pressure poisson equation
!-----------------------------------------------------------------------------------------------------


subroutine proto_pressure_multigrid(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData
    use module_fft
    use module_poisson

    implicit none

    interface
        subroutine compute_NonLinear(params, hvy_u, hvy_NL, order_discretization, treeID)
            use module_globals
            use module_params
            type(type_params), intent(in) :: params
            real(kind=rk), intent(inout) :: hvy_u(:,:,:,:,:)
            real(kind=rk), intent(out) :: hvy_NL(:,:,:,:,:)
            character(len=cshort), intent(in) :: order_discretization
            integer(kind=ik), intent(in) :: treeID
        end subroutine compute_NonLinear

        subroutine compute_divergence_tree(params, hvy_u, hvy_div, order_discretization, treeID)
            use module_globals
            use module_params
            type(type_params), intent(in) :: params
            real(kind=rk), intent(in) :: hvy_u(:,:,:,:,:)
            real(kind=rk), intent(out) :: hvy_div(:,:,:,:)
            character(len=cshort), intent(in) :: order_discretization
            integer(kind=ik), intent(in) :: treeID
        end subroutine compute_divergence_tree
    end interface

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: file_ux, file_uy, file_uz, file_p
    real(kind=rk)                      :: dt, t, time
    integer(kind=ik)                   :: k_block, lgt_ID, hvy_id, Bs(1:3)

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :), hvy_work(:, :, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, Jmin, ic, nc, i_cycle, it, tc_length, mpierr, g(1:3)

    character(len=cshort)              :: fname, cycle_type
    logical                            :: exist_p
    real(kind=rk)                      :: x0(1:3), dx(1:3), domain(1:3), norm(1:9), volume
    integer(kind=tsize)                :: treecode

    real(kind=rk)        :: t_block, t_loop, t_cycle

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (fname and other options)
    call get_command_argument(2, file_ux)

    ! does the user need help?
    if (file_ux=='--help' .or. file_ux=='--h' .or. file_ux=='-h') then
        if (params%rank==0) then
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Wabbit prototype: Multi-Grid Gauss-Seidel solver"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Read in a data field (2D or 3D) b and solves the system laplacian u = b"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') "./wabbit-post --proto-GS-multigrid ux.h5 uy.h5 uz.h5 [p.h5] WAVELET LAPLACE-ORDER IT_CYCLE IT_GS"
            write(*, '(A)') "    ux.h5 = Velocity in X-direction"
            write(*, '(A)') "    uy.h5 = Velocity in Y-direction"
            write(*, '(A)') "    uz.h5 = Velocity in Z-direction (optional, for 2D pass ux or uy again"
            write(*, '(A)') "    p.h5 = Pressure field for comparison (optional)"
            write(*, '(A)') "    WAVELET = CDFXY"
            write(*, '(A)') "    LAPLACE-ORDER = 2, 4, 6, 8 for the order of the Laplace operator"
            write(*, '(A)') "    IT_CYCLE = Number of cycles to perform"
            write(*, '(A)') "    IT_GS = Number of GS-sweeps for upwards iterations"
            write(*, '(A)') "-----------------------------------------------------------"
        end if
        return
    endif
    call check_file_exists(trim(file_ux))
    call get_command_argument(3, file_uy)
    call check_file_exists(trim(file_uy))
    call get_command_argument(4, file_uz)
    call check_file_exists(trim(file_uz))

    call get_command_argument(5, file_p)
    inquire( file=file_p, exist=exist_p )
    if (exist_p) then
        if (params%rank == 0) write(*, '(A, A, A)') "Found file ", trim(file_p), " for comparison."
    else
        if (params%rank == 0) write(*, '(A)') "No comparison will be made, as no file was provided."
    endif

    call get_command_argument(5 + merge(1,0, exist_p), params%wavelet)

    ! set laplace order
    call get_command_argument(6 + merge(1,0, exist_p), params%poisson_order)
    ! set number of cycles
    call get_command_argument(7 + merge(1,0, exist_p), cycle_type)
    read (cycle_type, *) params%poisson_cycle_it
    call get_command_argument(8 + merge(1,0, exist_p), cycle_type)
    read (cycle_type, *) params%poisson_GS_it
    params%poisson_Sync_it = 2

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file_ux, lgt_n(tree_ID), time, it, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    ! ! odd BS test - but that doesn't work so great for now
    ! BS(1:params%dim) = Bs(1:params%dim) + 1

    params%Jmax = tc_length
    params%Jmin = 0   ! yes, I want to go down to one block only!
    params%n_eqn = params%dim
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs
    allocate(params%butcher_tableau(1,1))

    params%n_eqn = 3 + merge(1,0, exist_p)
    nc = 1

    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component(1) = "0"
    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component = 1
    select case(params%poisson_order)
    case ("FD_4th_comp_2_2")
        params%order_discretization = "FD_4th_central"
    case ("FD_6th_comp_3_3")
        params%order_discretization = "FD_6th_central"
    case default
        params%order_discretization = params%poisson_order
    end select
    params%poisson_coarsest = "FFT"
    params%FFT_accuracy = "FD"  ! FD or spectral
    params%poisson_Sync_it = 2

    Bs = params%Bs

    call setup_wavelet(params, params%g)
    call setup_laplacian_stencils(params, params%g)

    t_block = MPI_Wtime()
    call fft_initialize(params)
    call toc( "fft initialize", 10101, MPI_Wtime()-t_block )

    ! we need to decompose values, so we need atleast 2**dim/2**(dim-1) blocks
    params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) * 2**params%dim/(2**params%dim-1) + 8 )

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_work=hvy_work, nrhs_slots1=1, neqn_hvy_tmp=8)

    ! read input data
    if (exist_p) then
        call readHDF5vct_tree( (/file_ux, file_uy, file_uz, file_p/), params, hvy_block, tree_ID)
    else
        call readHDF5vct_tree( (/file_ux, file_uy, file_uz/), params, hvy_block, tree_ID)
    endif

    ! maye we didn't sync
    call sync_ghosts_tree( params, hvy_block, tree_ID )

    ! balance the load
    call balanceLoad_tree(params, hvy_block, tree_ID)

    ! ------------------------------------------------------------
    ! Gauss-Seidel Multi-Grid - solve the system Ax = b
    !    When going to lower levels, we compute the residual r = b - Ax and solve on lower levels Ae = r
    !    This is repeated until we reach the lowest level, then we go back upwards, completing one V-cycle.
    ! ------------------------------------------------------------
    !    for now I assume equidistant grid, so I loop level-wise
    !    in order to avoid leaf-wise refinement status hacking for now
    ! ------------------------------------------------------------

    call init_t_file('multigrid-cycle.t', .true., (/'    residual L2', '    residual L1', 'residual Linfty', '           time'/))
    call init_t_file('multigrid-iteration.t', .true., (/'direction', 'iteration', '     time'/))
    if (exist_p) then
        call init_t_file('multigrid-compare.t', .true., (/'             iteration', '    diff 2 spectral L2', '    diff 2 spectral L1', 'diff 2 spectral Linfty'/))
    endif
    ! first entries are always ignored when reading in files, as this is t=0?
    call append_t_file('multigrid-cycle.t', (/0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk/))
    call append_t_file('multigrid-iteration.t', (/0.0_rk, 0.0_rk, 0.0_rk/))
    if (exist_p) call append_t_file('multigrid-compare.t', (/0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk/))

    call componentWiseNorm_tree(params, hvy_block(:,:,:,1:size(hvy_block,4),:), tree_ID, "Mean", norm(1:size(hvy_block,4)), threshold_state_vector=.false.)
    if (params%rank == 0) write(*, '(A, 3(1x,es10.3))') "--- Mean values of read ux, uy, (uz), p: ", norm(1:size(hvy_block,4))

    ! compute RHS = -div ( 1/2 ((grad u + u grad)u))
    ! in 2D: -1/2 * ((uu_dx + uv_dy + u*u_dx + v*u_dy)_dx +
    !                (vu_dx + vv_dy + u*v_dx + v*v_dy)_dy)
    ! in 3D: -1/2 * ((uu_dx + uv_dy + uw_dz + u*u_dx + v*u_dy + w*u_dz)_dx +
    !                (vu_dx + vv_dy + vw_dz + u*v_dx + v*v_dy + w*v_dz)_dy +
    !                (wu_dx + wv_dy + ww_dz + u*w_dx + v*w_dy + w*w_dz)_dz)
    t_block = MPI_Wtime()
    call compute_NonLinear(params, hvy_block(:,:,:,1:3,:), hvy_tmp(:,:,:,1:3,:), params%order_discretization, tree_ID)
    call sync_ghosts_tree(params, hvy_tmp(:,:,:,1:3,:), tree_ID)
    call toc( "Compute Non-Linear term", 9998, MPI_Wtime()-t_block )

    t_block = MPI_Wtime()
    call compute_divergence_tree(params, hvy_tmp(:,:,:,1:3,:), hvy_block(:,:,:,1,:), params%order_discretization, tree_ID)
    call sync_ghosts_tree(params, hvy_block(:,:,:,1:1,:), tree_ID)
    call toc( "Compute Divergence", 9999, MPI_Wtime()-t_block )

    call componentWiseNorm_tree(params, hvy_block(:,:,:,1:size(hvy_block,4),:), tree_ID, "Mean", norm(1:size(hvy_block,4)), threshold_state_vector=.false.)
    if (params%rank == 0) write(*, '(A, 3(1x,es10.3))') "--- Mean values of RHS and p: ", norm(1), norm(size(hvy_block,4))
    do ic = 1,size(hvy_block,4)
        do k_block = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k_block, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            hvy_block(:,:,:,ic,hvy_id) = hvy_block(:,:,:,ic,hvy_id) - norm(ic)
        enddo
    enddo

    ! For the Mehrstellenverfahren, we actually do solve the system Au = Bb
    ! So before doing anything, we apply the matrix B on b, this is a large tensorial matrix, but we only have to apply this operation once
    if (params%poisson_order == "FD_6th_mehrstellen") then
        t_block = MPI_Wtime()
        do k_block = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k_block, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            call get_block_spacing_origin( params, lgt_ID, x0, dx )
            ! apply the stencil to the RHS
            call GS_compute_Ax(params, hvy_block(:,:,:,1,hvy_id), hvy_tmp(:,:,:,1,hvy_id), dx, apply_B_RHS=.true.)
            hvy_block(:,:,:,1,hvy_id) = hvy_tmp(:,:,:,1,hvy_id)
        enddo
        call toc( "RHS Preparation Bb", 10000, MPI_Wtime()-t_block )

        t_block = MPI_Wtime()
        call sync_ghosts_tree(params, hvy_block(:,:,:,1:1,:), tree_ID)
        call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )
    endif

    ! prepare full tree grid, this will be populated along the way
    ! blocks are practically empty, but we fill them along the way so ref flag will be set to 0 to allow synching
    call init_full_tree(params, tree_ID, set_ref=0)

    ! init solution as zero
    do k_block = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        hvy_tmp(:,:,:,1,hvy_id) = 0.0_rk
    enddo

    ! compute several v- or f-cycles
    do i_cycle = 1, params%poisson_cycle_it

        ! ! usually we want spectral resolution for the coarsest level, for debugging we can disable this
        ! fft_order_discretization = "spectral"

        call multigrid_vcycle(params, hvy_tmp(:,:,:,1:1,:), hvy_block(:,:,:,1:1,:), hvy_tmp(:,:,:,2:size(hvy_tmp,4),:), tree_ID, residual_out=norm, verbose=.true.)

        ! laplacian is invariant to shifts of constant values
        ! our values are defined with zero mean for comparison
        ! as multigrid might accidently introduce a constant offset, we remove it
        call componentWiseNorm_tree(params, hvy_tmp(:,:,:,1:1,:), tree_ID, "Mean", norm(1:1), threshold_state_vector=.false.)
        do k_block = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k_block, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            hvy_tmp(:,:,:,1,hvy_id) = hvy_tmp(:,:,:,1,hvy_id) - norm(1)
        enddo
        if (params%rank == 0) write(*, '(A, es10.3, A)') "--- Mean value: ", norm(1), " ---"

        ! compare to spectral solutions with spectral or FD accuracy
        if (exist_p) then
            ! compute difference to spectral accuracy
            do k_block = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k_block, tree_ID)
            if (.not. block_is_leaf(params, hvy_id)) cycle
            hvy_work(:,:,:,1,hvy_id,1) = hvy_tmp(:,:,:,1,hvy_id) - hvy_block(:,:,:,4,hvy_id)
            enddo
            ! compute norms
            call componentWiseNorm_tree(params, hvy_work(:,:,:,1:1,:,1), tree_ID, "L2", norm(1:1), threshold_state_vector=.false.)
            if (params%rank == 0) write(*, '(A, es10.3, A)') "--- Diff spectral L2: ", norm(1), " ---"
            call componentWiseNorm_tree(params, hvy_work(:,:,:,1:1,:,1), tree_ID, "L1", norm(2:2), threshold_state_vector=.false.)
            if (params%rank == 0) write(*, '(A, es10.3, A)') "--- Diff spectral L1: ", norm(2), " ---"
            call componentWiseNorm_tree(params, hvy_work(:,:,:,1:1,:,1), tree_ID, "Linfty", norm(3:3), threshold_state_vector=.false.)
            if (params%rank == 0) write(*, '(A, es10.3, A)') "--- Diff spectral Linfty: ", norm(3), " ---"

            call append_t_file('multigrid-compare.t', (/dble(i_cycle), norm(1), norm(2), norm(3)/))
        endif

        ! for not saving the intermediate values, also change all save times to 0
        enddo

        ! delete all non-leaf blocks with daughters as we for now do not have any use for them
        call prune_fulltree2leafs(params, tree_ID)
        
        ! save file under new name
        write(fname, '(A, i4.4, A)') "p_00", 0, "00000.h5"
        call saveHDF5_tree(fname, dble(0), 0, 1, params, hvy_tmp, tree_ID )

        ! save file under new name
        write(fname, '(A, i4.4, A)') "res_00", 0, "00000.h5"
        call saveHDF5_tree(fname, dble(0), 0, 1, params, hvy_tmp(:,:,:,2:size(hvy_tmp,4),:), tree_ID )

        ! save file under new name
        write(fname, '(A, I4.4, A)') "RHS_00", 0, "00000.h5"
        call saveHDF5_tree(fname, dble(0), 0, 1, params, hvy_block, tree_ID )

        ! ! for saving intermediate value, also change all save times to i_cycle
        !     ! save spectral u at new time position
        !     if (exist_p) then
        !         write(fname, '(A, I4.4, A)') "u-spectral_00", i_cycle, "00000.h5"
        !         call saveHDF5_tree(fname, dble(i_cycle), i_cycle, 1, params, hvy_block(:,:,:,2,:), tree_ID )
        !     endif
        !     call init_full_tree(params, tree_ID, set_ref=0)
        ! enddo

    t_block = MPI_Wtime()
    call fft_destroy(params)
    call toc( "fft destroy", 10100, MPI_Wtime()-t_block )

    call summarize_profiling( WABBIT_COMM )

end subroutine proto_pressure_multigrid


subroutine proto_NSI_EE(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData
    use module_fft
    use module_poisson
    use module_ini_files_parser_mpi
    use module_time_step  ! statistics wrapper
    ! this module is the saving wrapper (e.g. save state vector or vorticity)
    ! it exists to disentangle module_forest and module_IO
    use module_saving

    implicit none

    interface
        subroutine compute_divergence_tree(params, hvy_u, hvy_div, order_discretization, treeID)
            use module_globals
            use module_params
            type(type_params), intent(in) :: params
            real(kind=rk), intent(in) :: hvy_u(:,:,:,:,:)
            real(kind=rk), intent(out) :: hvy_div(:,:,:,:)
            character(len=cshort), intent(in) :: order_discretization
            integer(kind=ik), intent(in) :: treeID
        end subroutine compute_divergence_tree

        subroutine compute_NSI_RHS(params, hvy_u, hvy_mask, hvy_RHS, order_discretization, treeID, nu, C_eta)
            use module_globals
            use module_params
            type(type_params), intent(in) :: params
            real(kind=rk), intent(in) :: hvy_u(:,:,:,:,:)
            real(kind=rk), intent(in) :: hvy_mask(:,:,:,:,:)
            real(kind=rk), intent(out) :: hvy_RHS(:,:,:,:,:)
            character(len=cshort), intent(in) :: order_discretization
            integer(kind=ik), intent(in) :: treeID
            real(kind=rk), intent(in), optional :: nu, C_eta
        end subroutine compute_NSI_RHS

        subroutine compute_projection(params, hvy_u, hvy_p, hvy_u_div0, order_discretization, treeID, dt)
            use module_globals
            use module_params
            type(type_params), intent(in) :: params
            real(kind=rk), intent(in) :: hvy_u(:,:,:,:,:)
            real(kind=rk), intent(in) :: hvy_p(:,:,:,:,:)
            real(kind=rk), intent(out) :: hvy_u_div0(:,:,:,:,:)
            character(len=cshort), intent(in) :: order_discretization
            integer(kind=ik), intent(in) :: treeID
            real(kind=rk), intent(in) :: dt
        end subroutine compute_projection
    end interface


    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: file_params
    real(kind=rk)                      :: dt, time
    integer(kind=ik)                   :: k_block, lgt_ID, hvy_id, Bs(1:3), g(1:3), iteration, j, l

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :), hvy_mask(:, :, :, :, :), hvy_work(:, :, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, Jmin, ic, nc, i_cycle, it, tc_length, mpierr

    character(len=cshort)              :: fname, order_disc_nonlinear, order_disc_pressure, order_laplacian
    real(kind=rk)                      :: x0(1:3), dx(1:3), domain(1:3), norm(1:6), volume
    integer(kind=tsize)                :: treecode
    logical                            :: it_is_time_to_save_data=.false., overwrite, error_OOM, is_equidistant
    integer(kind=ik)                   :: Nblocks_rhs, Nblocks, lgt_n_tmp, mpicode, Jmin1, Jmax1

    type(inifile) :: FILE

    ! hack variables needed from ACM
    real(kind=rk) :: nu, C_eta
    logical :: penalization

    real(kind=rk)        :: t_block, t_loop, t_cycle, t2, t3, t4

    ! ! this routine works only on one tree
    ! allocate( hvy_n(1), lgt_n(1) )

    ! init time loop
    time          = 0.0_rk
    dt            = 0.0_rk
    iteration     = 0

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (fname and other options)
    call get_command_argument(2, file_params)

    ! does the user need help?
    if (file_params=='--help' .or. file_params=='--h' .or. file_params=='-h') then
        if (params%rank==0) then
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Wabbit prototype: Multi-Grid Gauss-Seidel solver"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Read in a data field (2D or 3D) b and solves the system laplacian u = b"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') "./wabbit-post --proto-GS-multigrid PARAMS.ini"
            write(*, '(A)') "    PARAMS.ini = ini-file to read in, has to be set to ACM for now"
            write(*, '(A)') "-----------------------------------------------------------"
        end if
        return
    endif
    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, file_params )
    ! have the pysics module read their own parameters
    call init_physics_modules( params, file_params, params%N_mask_components  )

    order_laplacian = params%poisson_order
    order_disc_nonlinear = params%poisson_order
    order_disc_pressure = params%poisson_order
    params%order_discretization = params%poisson_order
    
    ! some special snowflakes
    if (params%poisson_order == "FD_6th_mehrstellen") then
        order_disc_nonlinear = "FD_6th_central"
        order_disc_pressure = "FD_6th_central"
        params%order_discretization = "FD_6th_central"
    ! elseif (params%poisson_order == "FD_4th_comp_1_3") then
    !     order_disc_nonlinear = "FD_4th_comp_1_3"
    !     order_disc_pressure = "FD_4th_comp_1_3"
    !     params%order_discretization = "FD_4th_comp_1_3"
    endif

    ! HACK - read in values that are only read by ACM module but that we need
    call read_ini_file_mpi(FILE, file_params, .true.)
    call read_param_mpi(FILE, 'ACM-new', 'nu', nu, 5.0e-5_rk )
    call read_param_mpi(FILE, 'VPM', 'C_eta', C_eta, 1.0_rk )
    call read_param_mpi(FILE, 'VPM', 'penalization', penalization, .false.)

    call setup_wavelet(params, params%g)
    ! order predictor is decided by X in CDFXY (it is coinciding with the HR filter actually) - but should be 2 higher for poisson solver
    if (params%wavelet(4:4) == "2") then
        params%order_predictor = "multiresolution_4th"
    elseif (params%wavelet(4:4) == "4") then
        params%order_predictor = "multiresolution_6th"
    elseif (params%wavelet(4:4) == "6") then
        params%order_predictor = "multiresolution_8th"
    elseif (params%wavelet(4:4) == "8") then
        params%order_predictor = "multiresolution_10th"
    elseif (params%wavelet(4:5) == "10" .or. params%wavelet(4:5) == "12") then
        params%order_predictor = "multiresolution_12th"
    endif
    call setup_laplacian_stencils(params, params%g)

    t_block = MPI_Wtime()
    call fft_initialize(params)
    call toc( "fft initialize", 10101, MPI_Wtime()-t_block )

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_mask=hvy_mask, hvy_work=hvy_work, neqn_hvy_tmp=8)

    ! On all blocks, set the initial condition (incl. synchronize ghosts)
    call setInitialCondition_tree( params, hvy_block, tree_ID_flow, params%adapt_inicond, time, iteration, hvy_mask, hvy_tmp, hvy_work=hvy_work)


    ! initialize t-files
    overwrite = .false.
    if (params%rank==0 .and. iteration==0) then
        overwrite = .true.
    endif

    ! the physics modules should initialize the ascii files they require, e.g. for
    ! saving kinetic energy over time, etc.
    call INITIALIZE_ASCII_FILES_meta( params%physics_type, time, overwrite )

    ! a few files have to be intialized by wabbit, because they are logfiles produced
    ! by wabbit independent of the physics modules.
    call init_t_file('dt.t', overwrite)
    call init_t_file('performance.t', overwrite)
    call init_t_file('eps_norm.t', overwrite)
    call init_t_file('krylov_err.t', overwrite)
    call init_t_file('balancing.t', overwrite)
    call init_t_file('block_xfer.t', overwrite)
    call init_t_file('thresholding.t', overwrite)

    if (params%rank==0) then
        call Initialize_runtime_control_file()
    endif

    ! shorten indices
    Bs   = params%Bs
    g(:) = params%g
    if (params%dim==2) g(3) = 0

    ! ------------------------------------------------------------
    ! Gauss-Seidel Multi-Grid - solve the system Ax = b
    !    When going to lower levels, we compute the residual r = b - Ax and solve on lower levels Ae = r
    !    This is repeated until we reach the lowest level, then we go back upwards, completing one V-cycle.
    ! ------------------------------------------------------------
    !    for now I assume equidistant grid, so I loop level-wise
    !    in order to avoid leaf-wise refinement status hacking for now
    ! ------------------------------------------------------------

    ! Maybe u is not completely divergence free in the discrete sense, so we do one projection in order to get a divergence free velocity field
    call sync_ghosts_tree( params, hvy_block(:,:,:,1:params%dim,:), tree_ID_flow )
    call compute_divergence_tree(params, hvy_block(:,:,:,1:params%dim,:), hvy_tmp(:,:,:,params%dim+1,:), order_disc_pressure, tree_ID_flow)
    call sync_ghosts_tree( params, hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), tree_ID_flow )
    call multigrid_solve(params, hvy_tmp(:,:,:,params%dim+2:params%dim+2,:), hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), hvy_tmp(:,:,:,params%dim+3:size(hvy_tmp,4),:), tree_ID_flow, init_0=.true., verbose=.false., hvy_full=hvy_tmp)
    call sync_ghosts_tree( params, hvy_tmp(:,:,:,params%dim+2:params%dim+2,:), tree_ID_flow )
    call compute_projection(params, hvy_block(:,:,:,1:params%dim,:), hvy_tmp(:,:,:,params%dim+2:params%dim+2,:), hvy_block(:,:,:,1:params%dim,:), order_disc_pressure, tree_ID_flow, 1.0_rk)

    ! --- following block is for debugging the helmholtz projection ----
    ! call sync_ghosts_tree( params, hvy_block(:,:,:,1:params%dim,:), tree_ID_flow )
    ! write(fname, '(A, i12.12, A)') "ux_", int(100*1.0e6),".h5"
    ! call saveHDF5_tree(fname, 100.0_rk, 0, 1, params, hvy_block, tree_ID )
    ! write(fname, '(A, i12.12, A)') "uy_", int(100*1.0e6),".h5"
    ! call saveHDF5_tree(fname, 100.0_rk, 0, 2, params, hvy_block, tree_ID )

    ! ! Maybe u is not completely divergence free in the discrete sense, so we do one projection in order to get a divergence free velocity field
    ! call sync_ghosts_tree( params, hvy_block(:,:,:,1:params%dim,:), tree_ID_flow )
    ! call compute_divergence_tree(params, hvy_block(:,:,:,1:params%dim,:), hvy_tmp(:,:,:,params%dim+1,:), order_disc_pressure, tree_ID_flow)
    ! call sync_ghosts_tree( params, hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), tree_ID_flow )
    ! call multigrid_solve(params, hvy_block(:,:,:,params%dim+1:params%dim+1,:), hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), hvy_tmp(:,:,:,params%dim+2:size(hvy_tmp,4),:), tree_ID_flow, init_0=.true., verbose=.false., hvy_full=hvy_tmp)
    ! call sync_ghosts_tree( params, hvy_block(:,:,:,params%dim+1:params%dim+1,:), tree_ID_flow )
    ! call compute_projection(params, hvy_block(:,:,:,1:params%dim,:), hvy_block(:,:,:,params%dim+1:params%dim+1,:), hvy_block(:,:,:,1:params%dim,:), order_disc_pressure, tree_ID_flow, 1.0_rk)

    ! call sync_ghosts_tree( params, hvy_block(:,:,:,1:params%dim,:), tree_ID_flow )
    ! write(fname, '(A, i12.12, A)') "u2x_", int(100*1.0e6),".h5"
    ! call saveHDF5_tree(fname, 100.0_rk, 0, 1, params, hvy_block, tree_ID )
    ! write(fname, '(A, i12.12, A)') "u2y_", int(100*1.0e6),".h5"
    ! call saveHDF5_tree(fname, 100.0_rk, 0, 2, params, hvy_block, tree_ID )
    ! -----------------------------------------------------------------

    ! we have to recompute p
    call sync_ghosts_tree( params, hvy_block(:,:,:,1:params%dim,:), tree_ID_flow )
    call compute_NSI_RHS(params, hvy_block(:,:,:,1:params%dim,:), hvy_mask, hvy_tmp(:,:,:,1:params%dim,:), order_disc_nonlinear, tree_ID_flow, nu, C_eta)
    call sync_ghosts_tree( params, hvy_tmp(:,:,:,1:params%dim,:), tree_ID_flow )
    call compute_divergence_tree(params, hvy_tmp(:,:,:,1:params%dim,:), hvy_tmp(:,:,:,params%dim+1,:), order_disc_pressure, tree_ID_flow)
    call sync_ghosts_tree( params, hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), tree_ID_flow )
    call multigrid_solve(params, hvy_tmp(:,:,:,params%dim+2:params%dim+2,:), hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), hvy_tmp(:,:,:,params%dim+3:size(hvy_tmp,4),:), tree_ID_flow, hvy_full=hvy_tmp)
    ! multigrid solve needs consistent usage of hvy_tmp, so we write data from hvy_tmp to hvy_block for pressure afterwards
    do k_block = 1, hvy_n(tree_ID_flow)
        hvy_id = hvy_active(k_block, tree_ID_flow)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        hvy_block(:,:,:,params%dim+1,hvy_id) = hvy_tmp(:,:,:,params%dim+2,hvy_id)
    enddo

    !*******************************************************************
    ! statistics
    !*******************************************************************
    t_block = MPI_wtime()
    ! we need to sync ghost nodes for some derived qtys, for sure
    call sync_ghosts_RHS_tree( params, hvy_block, tree_ID_flow )
    call statistics_wrapper(time, dt, params, hvy_block, hvy_tmp, hvy_mask, tree_ID_flow)
    call toc( "TOPLEVEL: statistics", 13, MPI_wtime()-t_block)
    !***********************************************************************
    ! Write fields to HDF5 file
    !***********************************************************************
    ! do not save any output before this time (so maybe revoke the previous decision)
    if (time+1e-12_rk>=params%write_time_first) then
        call save_data( iteration, time, params, hvy_block, hvy_tmp, hvy_mask, tree_ID_flow )
    endif

    !---------------------------------------------------------------------------
    ! main time loop
    !---------------------------------------------------------------------------
    if (params%rank==0) then
        write(*,*) ""
        write(*,'(10(" "), "╔", 48("═"), "╗")') 
        write(*,'(10(" "), A)') "║ On your marks, ready, starting main time loop! ║"
        write(*,'(10(" "), "╚", 48("═"), "╝")')
        write(*,*) ""
    endif

    ! ------------------------------------------------------------
    ! Explicit Euler timestepping
    !    we assume equidistant grid for now, so no REC loop
    do while ( time<params%time_max .and. iteration<params%nt)
        t2 = MPI_wtime()
        !***********************************************************************
        ! Refine - make sure the solution at t+dt can be resolved on the grid
        !***********************************************************************
        t4 = MPI_wtime()
        if ( params%adapt_tree ) then
            ! synchronization before refinement (because the interpolation takes place on the extended blocks
            ! including the ghost nodes)
            ! Note: at this point the grid is rather coarse (fewer blocks), and the sync step is rather cheap.
            ! Snych'ing becomes much more expensive one the grid is refined.
            call sync_ghosts_tree( params, hvy_block, tree_ID_flow )

            ! refine the mesh after refinement_indicator, usually "everywhere" or "significant". When using
            ! "significant", the refinement flags from the last call to adapt_tree call are reused. 
            ! This might not be given for the first iteration so we just skip this (as adapt_tree was not called yet)
            ! and use the "everywhere" indicator. Note: after resuming a run from backup, this works as well, because
            ! the refinement_flag is 0 and this results in "significant" refining in fact all blocks.
            if (params%refinement_indicator == "significant" .and. iteration == 0) then
                call refine_tree( params, hvy_block, "everywhere", tree_ID=tree_ID_flow, error_OOM=error_OOM, check_full_tree=.true.)
            else
                call refine_tree( params, hvy_block, params%refinement_indicator, tree_ID=tree_ID_flow, error_OOM=error_OOM, check_full_tree=.true.)
            endif
            ! if refine_tree runs out-of-memory (OOM), it does not actually do the refinement, but returns after realizing
            ! there won't be enough mem. We can thus jump to 17 to save a backup and terminate.
            if (error_OOM) goto 17
        endif
        call toc( "TOPLEVEL: refinement", 10, MPI_wtime()-t4)

        Nblocks_rhs = lgt_n(tree_ID_flow)
        Jmin1 = minActiveLevel_tree(tree_ID_flow)
        Jmax1 = maxActiveLevel_tree(tree_ID_flow)
        is_equidistant = Jmin1==Jmax1

        !***********************************************************************
        ! Timestep - Euler Explicit or RK2 / Heun's method
        !***********************************************************************
        t4 = MPI_wtime()

        ! RK_generic
        t3 = MPI_wtime()
        call sync_ghosts_tree( params, hvy_block(:,:,:,1:params%dim,:), tree_ID_flow )
        call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)

        ! caluclate timestep, here we use that of ACM and set the speed of sound to close to 0
        call calculate_time_step(params, time, iteration, hvy_block, dt, tree_ID)

        ! compute RHS for the first stage without pressure gradient
        t3 = MPI_wtime()
        call compute_NSI_RHS(params, hvy_block(:,:,:,1:params%dim,:), hvy_mask, hvy_work(:,:,:,1:params%dim,:,1), order_disc_nonlinear, tree_ID_flow, nu, C_eta )
        call toc( "timestep: compute_NSI_RHS", 22, MPI_wtime()-t3)

        ! compute RHS for the pressure Poisson equation, solve the Poisson equation and add pressure gradient
        t3 = MPI_wtime()
        call sync_ghosts_tree( params, hvy_work(:,:,:,1:params%dim,:,1), tree_ID_flow )
        call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
        t3 = MPI_wtime()
        call compute_divergence_tree(params, hvy_work(:,:,:,1:params%dim,:,1), hvy_tmp(:,:,:,1,:), order_disc_pressure, tree_ID_flow)
        call toc( "timestep: compute_divergence_tree", 23, MPI_wtime()-t3)
        t3 = MPI_wtime()
        ! if the grid doesn't change, then we can reuse the previous solution as initial guess
        ! if not, then we should start from zero, as blocks might have been moved or new ones activated, resulting in rubbish initial guess
        if (params%adapt_tree) then
            call multigrid_solve(params, hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), hvy_tmp(:,:,:,1:1,:), &
                                hvy_tmp(:,:,:,params%dim+2:size(hvy_tmp,4),:), tree_ID_flow, init_0=.true., &
                                verbose=.false., hvy_full=hvy_tmp)
        else
            call multigrid_solve(params, hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), hvy_tmp(:,:,:,1:1,:), &
                                hvy_tmp(:,:,:,params%dim+2:size(hvy_tmp,4),:), tree_ID_flow, init_0=.false., &
                                verbose=.false., hvy_full=hvy_tmp)
        endif
        call toc( "timestep: multigrid_solve", 24, MPI_wtime()-t3)
        t3 = MPI_wtime()
        call compute_projection(params, hvy_work(:,:,:,1:params%dim,:,1), hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), hvy_work(:,:,:,1:params%dim,:,1), order_disc_pressure, tree_ID_flow, 1.0_rk)
        call toc( "timestep: compute_projection", 25, MPI_wtime()-t3)

        ! compute k_1, k_2, .... (coefficients for final stage)
        do j = 2, size(params%butcher_tableau, 1) - 1
            ! build intermediate velocity as input for RHS of this stage
            do k_block = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k_block, tree_ID)
                hvy_tmp(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),1:params%dim,hvy_id) = hvy_block(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),1:params%dim,hvy_id)
            end do
            do l = 2, j
                ! check if coefficient is zero - if so, avoid loop over all components and active blocks
                if (abs(params%butcher_tableau(j,l)) < 1.0e-8_rk) cycle

                do k_block = 1, hvy_n(tree_ID)
                    hvy_id = hvy_active(k_block, tree_ID)
                    ! new input for computation of k-coefficients
                    ! k_j = RHS((t+dt*c_j, data_field(t) + sum(a_jl*k_l))
                    hvy_tmp(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),1:params%dim,hvy_id) = hvy_tmp(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),1:params%dim,hvy_id) &
                    + dt * params%butcher_tableau(j,l) * hvy_work(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),1:params%dim,hvy_id,l-1)
                end do
            end do
            t3 = MPI_wtime()
            call sync_ghosts_tree( params, hvy_tmp(:,:,:,1:params%dim,:), tree_ID_flow )
            call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
            ! compute RHS for the current stage without pressure gradient
            t3 = MPI_wtime()
            call compute_NSI_RHS(params, hvy_tmp(:,:,:,1:params%dim,:), hvy_mask, hvy_work(:,:,:,1:params%dim,:,j), order_disc_nonlinear, tree_ID_flow, nu, C_eta )
            ! incremental version, add old pressure gradient so that new pressure to be solved is only Dp - this reduces the amount of iterations for projection needed to converge, however it only really helps for adaptive grids
            if (.not. is_equidistant) then
                call compute_projection(params, hvy_work(:,:,:,1:params%dim,:,j), hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), hvy_work(:,:,:,1:params%dim,:,j), order_disc_pressure, tree_ID_flow, 1.0_rk)
            endif
            call toc( "timestep: compute_NSI_RHS", 22, MPI_wtime()-t3)
            ! compute RHS for the pressure Poisson equation, solve the Poisson equation and add pressure gradient
            t3 = MPI_wtime()
            call sync_ghosts_tree( params, hvy_work(:,:,:,1:params%dim,:,j), tree_ID_flow )
            call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
            t3 = MPI_wtime()
            call compute_divergence_tree(params, hvy_work(:,:,:,1:params%dim,:,j), hvy_tmp(:,:,:,1,:), order_disc_pressure, tree_ID_flow)
            call toc( "timestep: compute_divergence_tree", 23, MPI_wtime()-t3)
            t3 = MPI_wtime()
            ! interestingly, for adaptive grids we better initialize with zero, as it needs less iterations to converge
            call multigrid_solve(params, hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), hvy_tmp(:,:,:,1:1,:), hvy_tmp(:,:,:,params%dim+2:size(hvy_tmp,4),:), tree_ID_flow, init_0= .not. is_equidistant, verbose=.false., hvy_full=hvy_tmp)
            call toc( "timestep: multigrid_solve", 24, MPI_wtime()-t3)
            t3 = MPI_wtime()
            call compute_projection(params, hvy_work(:,:,:,1:params%dim,:,j), hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), hvy_work(:,:,:,1:params%dim,:,j), order_disc_pressure, tree_ID_flow, 1.0_rk)
            call toc( "timestep: compute_projection", 25, MPI_wtime()-t3)
        enddo

        ! final stage (actual final update of state vector)
        ! for the RK4 the final stage looks like this:
        ! data_field(t+dt) = data_field(t) + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4)
        do k_block = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k_block, tree_ID)
            do j = 2, size(params%butcher_tableau, 2)
                ! check if coefficient is zero - if so, avoid loop over all components and active blocks
                if ( abs(params%butcher_tableau(size(params%butcher_tableau, 1),j)) < 1.0e-8_rk) then
                    cycle
                endif

                ! ... dt*(b1*k1 + b2*k2+ ..)
                ! params%butcher_tableau(size(params%butcher_tableau,1)), since we want to access last line,
                ! e.g. b1 = butcher(last line,2)
                hvy_block(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),1:params%dim,hvy_id) = hvy_block(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),1:params%dim,hvy_id) &
                    + dt*params%butcher_tableau(size(params%butcher_tableau,1),j) * hvy_work(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),1:params%dim,hvy_id,j-1)
            end do
        end do
        t3 = MPI_wtime()
        call sync_ghosts_tree( params, hvy_block(:,:,:,1:params%dim,:), tree_ID_flow )
        call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)

        ! write(fname, '(A, i12.12, A)') "ux_", int(time*1.0e6),".h5"
        ! call saveHDF5_tree(fname, time, 0, 1, params, hvy_block, tree_ID )
        ! write(fname, '(A, i12.12, A)') "uy_", int(time*1.0e6),".h5"
        ! call saveHDF5_tree(fname, time, 0, 2, params, hvy_block, tree_ID )
        ! write(fname, '(A, i12.12, A)') "p_", int(time*1.0e6),".h5"
        ! call saveHDF5_tree(fname, time, 0, params%dim+1, params, hvy_block, tree_ID )

        ! increase time variables after all RHS substeps
        time = time + dt
        iteration = iteration + 1

        ! determine if it is time to save data
        it_is_time_to_save_data = .false.
        if ((params%write_method=='fixed_freq' .and. modulo(iteration, params%write_freq)==0) .or. &
            (params%write_method=='fixed_time' .and. abs(mod(time, params%write_time))<1.0e-12_rk)) then
            it_is_time_to_save_data = .true.
        endif
        ! do not save any output before this time (so maybe revoke the previous decision)
        if (time+1e-12_rk<params%write_time_first) then
            it_is_time_to_save_data = .false.
        endif
        ! it can rarely happen that not all proc arrive at the same time at the above condition, then some decide to
        ! save data and others do not. this is a rare but severe problem, to solve it, synchronize:
        call MPI_BCAST( it_is_time_to_save_data, 1, MPI_LOGICAL, 0, WABBIT_COMM, mpicode )

        ! let's do some regular projections - this can be changed by a parameter
        if (modulo(iteration, params%nprojection_NSI)==0) then
            ! Maybe u is not completely divergence free in the discrete sense, so we do one projection in order to get a divergence free velocity field
            t3 = MPI_wtime()
            call sync_ghosts_tree( params, hvy_block(:,:,:,1:params%dim,:), tree_ID_flow )
            call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
            t3 = MPI_wtime()
            call compute_divergence_tree(params, hvy_block(:,:,:,1:params%dim,:), hvy_tmp(:,:,:,params%dim+2,:), order_disc_pressure, tree_ID_flow)
            call toc( "timestep: compute_divergence_tree", 23, MPI_wtime()-t3)
            t3 = MPI_wtime()
            call sync_ghosts_tree( params, hvy_tmp(:,:,:,params%dim+2:params%dim+2,:), tree_ID_flow )
            call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
            t3 = MPI_wtime()
            call multigrid_solve(params, hvy_tmp(:,:,:,params%dim+3:params%dim+3,:), hvy_tmp(:,:,:,params%dim+2:params%dim+2,:), hvy_tmp(:,:,:,params%dim+4:size(hvy_tmp,4),:), tree_ID_flow, init_0=.true., verbose=.false., hvy_full=hvy_tmp)
            call toc( "timestep: multigrid_solve", 24, MPI_wtime()-t3)
            t3 = MPI_wtime()
            call sync_ghosts_tree( params, hvy_tmp(:,:,:,params%dim+3:params%dim+3,:), tree_ID_flow )
            call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
            t3 = MPI_wtime()
            call compute_projection(params, hvy_block(:,:,:,1:params%dim,:), hvy_tmp(:,:,:,params%dim+3:params%dim+3,:), hvy_block(:,:,:,1:params%dim,:), order_disc_pressure, tree_ID_flow, 1.0_rk)
            call toc( "timestep: compute_projection", 25, MPI_wtime()-t3)
        endif

        ! ! we reconstruct the pressure if we do statistics or save data to have it in respective accuracy
        ! ! if (it_is_time_to_save_data .or. (modulo(iteration, params%nsave_stats)==0).or.(abs(mod(time, params%tsave_stats))<1e-12_rk)) then
        ! ! pressure is not needed in statistics, so we skip it here
        ! if (it_is_time_to_save_data) then
        !     t3 = MPI_wtime()
        !     call sync_ghosts_tree( params, hvy_block(:,:,:,1:params%dim,:), tree_ID_flow )
        !     call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
        !     t3 = MPI_wtime()
        !     call compute_NSI_RHS(params, hvy_block(:,:,:,1:params%dim,:), hvy_mask, hvy_tmp(:,:,:,1:params%dim,:), order_disc_nonlinear, tree_ID_flow, nu, C_eta)
        !     call toc( "timestep: compute_NSI_RHS", 22, MPI_wtime()-t3)
        !     t3 = MPI_wtime()
        !     call sync_ghosts_tree( params, hvy_tmp(:,:,:,1:params%dim,:), tree_ID_flow )
        !     call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
        !     t3 = MPI_wtime()
        !     call compute_divergence_tree(params, hvy_tmp(:,:,:,1:params%dim,:), hvy_tmp(:,:,:,params%dim+2,:), order_disc_pressure, tree_ID_flow)
        !     call toc( "timestep: compute_divergence_tree", 23, MPI_wtime()-t3)
        !     t3 = MPI_wtime()
        !     call sync_ghosts_tree( params, hvy_tmp(:,:,:,params%dim+2:params%dim+2,:), tree_ID_flow )
        !     call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
        !     t3 = MPI_wtime()
        !     call multigrid_solve(params, hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), hvy_tmp(:,:,:,params%dim+2:params%dim+2,:), hvy_tmp(:,:,:,params%dim+3:size(hvy_tmp,4),:), tree_ID_flow, hvy_full=hvy_tmp)
        !     call toc( "timestep: multigrid_solve", 24, MPI_wtime()-t3)
        !     do k_block = 1, hvy_n(tree_ID_flow)
        !         hvy_id = hvy_active(k_block, tree_ID)
        !         hvy_block(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),params%dim+1,hvy_id) = hvy_tmp(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),4,hvy_id)
        !     enddo
        ! endif

        call toc( "TOPLEVEL: time stepper", 11, MPI_wtime()-t4)

        !*******************************************************************
        ! statistics
        !*******************************************************************
        if ( (modulo(iteration, params%nsave_stats)==0).or.(abs(mod(time, params%tsave_stats))<1e-12_rk) .or. abs(time - params%time_max) < 1e-12_rk ) then
            t4 = MPI_wtime()
            ! we need to sync ghost nodes for some derived qtys, for sure
            call sync_ghosts_RHS_tree( params, hvy_block, tree_ID_flow )

            call statistics_wrapper(time, dt, params, hvy_block, hvy_tmp, hvy_mask, tree_ID_flow)
            call toc( "TOPLEVEL: statistics", 13, MPI_wtime()-t4)
        endif

        !***********************************************************************
        ! Adapt mesh - (coarsening where possible)
        !***********************************************************************
        t4 = MPI_wtime()
        ! adapt the mesh
        if ( params%adapt_tree ) then
            ! some coarsening indicators require us to know the mask function (if
            ! it is considered as secondary criterion, e.g.). Creating the mask is a high-level
            ! routine that relies on forests and pruned trees, which are not available in the module_mesh.
            ! Hence the mask is created here.
            if (params%threshold_mask) then
                ! create mask function at current time
                call createMask_tree(params, time, hvy_mask, hvy_tmp)

                ! actual coarsening (including the mask function)
                call adapt_tree( time, params, hvy_block, tree_ID_flow, params%coarsening_indicator, hvy_tmp, &
                    hvy_mask=hvy_mask, hvy_work=hvy_work)
            else
                ! actual coarsening (no mask function is required)
                call adapt_tree( time, params, hvy_block, tree_ID_flow, params%coarsening_indicator, hvy_tmp, hvy_work=hvy_work)
            endif
        endif
        call toc( "TOPLEVEL: adapt mesh", 14, MPI_wtime()-t4)
        Nblocks = lgt_n(tree_ID_flow)

        !***********************************************************************
        ! Write fields to HDF5 file
        !***********************************************************************
        if (it_is_time_to_save_data) then
            ! JB comment: Tests show, that this will not affec the solution, but it breaks the significant refinement somehow, so it is disabled for now
            ! ------------------------------------------------------------
            ! ! before saving, we will do another projection to make sure the velocity field is as divergence free as possible
            ! ! this should get rid of divergence error from adaptive mesh changes
            ! ! Maybe u is not completely divergence free in the discrete sense, so we do one projection in order to get a divergence free velocity field
            ! t3 = MPI_wtime()
            ! call sync_ghosts_tree( params, hvy_block(:,:,:,1:params%dim,:), tree_ID_flow )
            ! call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
            ! t3 = MPI_wtime()
            ! call compute_divergence_tree(params, hvy_block(:,:,:,1:params%dim,:), hvy_tmp(:,:,:,params%dim+2,:), order_disc_pressure, tree_ID_flow)
            ! call toc( "timestep: compute_divergence_tree", 23, MPI_wtime()-t3)
            ! t3 = MPI_wtime()
            ! call sync_ghosts_tree( params, hvy_tmp(:,:,:,params%dim+2:params%dim+2,:), tree_ID_flow )
            ! call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
            ! t3 = MPI_wtime()
            ! call multigrid_solve(params, hvy_tmp(:,:,:,params%dim+3:params%dim+3,:), hvy_tmp(:,:,:,params%dim+2:params%dim+2,:), hvy_tmp(:,:,:,params%dim+4:size(hvy_tmp,4),:), tree_ID_flow, init_0=.true., verbose=.false., hvy_full=hvy_tmp)
            ! call toc( "timestep: multigrid_solve", 24, MPI_wtime()-t3)
            ! t3 = MPI_wtime()
            ! call sync_ghosts_tree( params, hvy_tmp(:,:,:,params%dim+3:params%dim+3,:), tree_ID_flow )
            ! call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
            ! t3 = MPI_wtime()
            ! call compute_projection(params, hvy_block(:,:,:,1:params%dim,:), hvy_tmp(:,:,:,params%dim+3:params%dim+3,:), hvy_block(:,:,:,1:params%dim,:), order_disc_pressure, tree_ID_flow, 1.0_rk)
            ! call toc( "timestep: compute_projection", 25, MPI_wtime()-t3)

            ! now at this position we compute the pressure field to have it in the saved data
            t3 = MPI_wtime()
            call sync_ghosts_tree( params, hvy_block(:,:,:,1:params%dim,:), tree_ID_flow )
            call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
            t3 = MPI_wtime()
            call compute_NSI_RHS(params, hvy_block(:,:,:,1:params%dim,:), hvy_mask, hvy_tmp(:,:,:,1:params%dim,:), order_disc_nonlinear, tree_ID_flow, nu, C_eta)
            call toc( "timestep: compute_NSI_RHS", 22, MPI_wtime()-t3)
            t3 = MPI_wtime()
            call sync_ghosts_tree( params, hvy_tmp(:,:,:,1:params%dim,:), tree_ID_flow )
            call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
            t3 = MPI_wtime()
            call compute_divergence_tree(params, hvy_tmp(:,:,:,1:params%dim,:), hvy_tmp(:,:,:,params%dim+2,:), order_disc_pressure, tree_ID_flow)
            call toc( "timestep: compute_divergence_tree", 23, MPI_wtime()-t3)
            t3 = MPI_wtime()
            call sync_ghosts_tree( params, hvy_tmp(:,:,:,params%dim+2:params%dim+2,:), tree_ID_flow )
            call toc( "timestep: sync_ghost_tree", 21, MPI_wtime()-t3)
            t3 = MPI_wtime()
            call multigrid_solve(params, hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), hvy_tmp(:,:,:,params%dim+2:params%dim+2,:), hvy_tmp(:,:,:,params%dim+3:size(hvy_tmp,4),:), tree_ID_flow, init_0=.true., verbose=.false., hvy_full=hvy_tmp)
            call toc( "timestep: multigrid_solve", 24, MPI_wtime()-t3)
            do k_block = 1, hvy_n(tree_ID_flow)
                hvy_id = hvy_active(k_block, tree_ID_flow)
                hvy_block(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),params%dim+1,hvy_id) = hvy_tmp(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),4,hvy_id)
            enddo


            ! NOTE new versions (>16/12/2017) call physics module routines call prepare_save_data. These
            ! routines create the fields to be stored in the work array hvy_tmp in the first 1:params%N_fields_saved
            ! slots. the state vector (hvy_block) is copied if desired.
            call save_data( iteration, time, params, hvy_block, hvy_tmp, hvy_mask, tree_ID_flow )
        endif

        ! output on screen
        t2 = MPI_wtime() - t2
        if (params%rank==0) then
            write(*, '("RUN: it=",i7)', advance='no') iteration
            write(*, '(" time=",f16.9, " t_wall=",es10.3)', advance='no') time, t2
            write(*, '(" Nb=(",i6,"/",i6,")")', advance='no') Nblocks_rhs, Nblocks
            write(*, '(" J=(",i2,":",i2,"/",i2,":",i2, ")")', advance='no') Jmin1, Jmax1, minActiveLevel_tree(tree_ID_flow), maxActiveLevel_tree(tree_ID_flow)
            ! the goofy factor (1.0+1.0/(2**params%dim)) is to take into account that if we have this many blocks at the RHS level,
            ! we still need a few more for the full tree transform (the virtual nodes which are not leafs) 1.25 in 2D and 1.125 in 3D
            write(*, '(" dt=",es8.1," mem=",i3,"%")', advance='no') dt, nint(((1.0+1.0/(2**params%dim))*dble(Nblocks_rhs+lgt_n(tree_ID_mask))/dble(size(lgt_block,1)))*100.0_rk)
            write(*, '("")')  ! line break

             ! prior to 11/04/2019, this file was called timesteps_info.t but it was missing some important
             ! information, so I renamed it when adding those (since post-scripts would no longer be compatible
             ! it made sense to me to change the name)
             call append_t_file( 'performance.t', (/time, dble(iteration), t2, dble(Nblocks_rhs), dble(Nblocks), &
             dble(minActiveLevel_tree(tree_ID_flow)), &
             dble(maxActiveLevel_tree(tree_ID_flow)), &
             dble(params%number_procs), dble(size(lgt_block,1)) /) )
        end if

    enddo

 17 t_block = MPI_Wtime()
    call fft_destroy(params)
    call toc( "fft destroy", 10100, MPI_Wtime()-t_block )

    call summarize_profiling( WABBIT_COMM )

end subroutine


!> \brief computes the non-linear term of the incompressible NS equations in skew symmetric form
subroutine compute_NonLinear(params, hvy_u, hvy_NL, order_discretization, treeID)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    type(type_params), intent(in) :: params
    real(kind=rk), intent(inout) :: hvy_u(:,:,:,:,:)
    real(kind=rk), intent(out) :: hvy_NL(:,:,:,:,:)
    character(len=cshort), intent(in) :: order_discretization
    integer(kind=ik), intent(in) :: treeID

    integer(kind=ik) :: k_block, hvy_id, lgt_ID, ix, iy, iz
    real(kind=rk) :: dx(1:3), x0(1:3), dx_inv, dy_inv, dz_inv
    real(kind=rk) :: u_dx, v_dx, w_dx, u_dy, v_dy, w_dy, u_dz, v_dz, w_dz
    real(kind=rk) :: uu_dx, uv_dy, uw_dz, vu_dx, vv_dy, vw_dz, wu_dx, wv_dy, ww_dz

    !> parameters for FD1_l, FD1_r, FD2 operator
    real(kind=rk), allocatable, dimension(:) :: FD1_l, FD1_r
    integer(kind=ik) :: FD1_ls, FD1_le, FD1_rs, FD1_re

    ! Setup finite difference stencils
    call setup_FD1_left_stencil(order_discretization, FD1_l, FD1_ls, FD1_le)
    call setup_FD1_right_stencil(order_discretization, FD1_r, FD1_rs, FD1_re)

    ! loop over all blocks
    do k_block = 1, hvy_n(treeID)
        hvy_id = hvy_active(k_block, treeID)
        call hvy2lgt(lgt_ID, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin(params, lgt_ID, x0, dx)
        dx_inv = 1.0_rk / dx(1)
        dy_inv = 1.0_rk / dx(2)
        if (params%dim == 3) then
            dz_inv = 1.0_rk / dx(3)
        else
            dz_inv = 0.0_rk
        endif

        ! compute non-linear term
        if (params%dim == 2) then
            ! 2D: NL_x = 1/2 * (uu_dx + uv_dy + u*u_dx + v*u_dy)
            !     NL_y = 1/2 * (vu_dx + vv_dy + u*v_dx + v*v_dy)
            do iy = params%g+1, params%Bs(2)+params%g
                do ix = params%g+1, params%Bs(1)+params%g
                    ! First derivatives
                    u_dx   = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,1,1,hvy_id)) * dx_inv
                    v_dx   = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,1,2,hvy_id)) * dx_inv

                    u_dy   = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,1,1,hvy_id)) * dy_inv
                    v_dy   = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,1,2,hvy_id)) * dy_inv

                    ! Non-linear energy terms
                    uu_dx = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix+FD1_rs:ix+FD1_re,iy,1,1,hvy_id) * hvy_u(ix+FD1_rs:ix+FD1_re,iy,1,1,hvy_id)) * dx_inv
                    uv_dy = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy+FD1_rs:iy+FD1_re,1,1,hvy_id) * hvy_u(ix,iy+FD1_rs:iy+FD1_re,1,2,hvy_id)) * dy_inv

                    vu_dx = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix+FD1_rs:ix+FD1_re,iy,1,2,hvy_id) * hvy_u(ix+FD1_rs:ix+FD1_re,iy,1,1,hvy_id)) * dx_inv
                    vv_dy = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy+FD1_rs:iy+FD1_re,1,2,hvy_id) * hvy_u(ix,iy+FD1_rs:iy+FD1_re,1,2,hvy_id)) * dy_inv

                    hvy_NL(ix,iy,1,1,hvy_id) = -0.5_rk*(uu_dx + uv_dy   + hvy_u(ix,iy,1,1,hvy_id)*u_dx + hvy_u(ix,iy,1,2,hvy_id)*u_dy )
                    hvy_NL(ix,iy,1,2,hvy_id) = -0.5_rk*(vu_dx + vv_dy   + hvy_u(ix,iy,1,1,hvy_id)*v_dx + hvy_u(ix,iy,1,2,hvy_id)*v_dy )
                enddo
            enddo
        elseif (params%dim == 3) then
            ! 3D: NL_x = 1/2 * (uu_dx + uv_dy + uw_dz + u*u_dx + v*u_dy + w*u_dz)
            !     NL_y = 1/2 * (vu_dx + vv_dy + vw_dz + u*v_dx + v*v_dy + w*v_dz)
            !     NL_z = 1/2 * (wu_dx + wv_dy + ww_dz + u*w_dx + v*w_dy + w*w_dz)
            do iz = params%g+1, params%Bs(3)+params%g
                do iy = params%g+1, params%Bs(2)+params%g
                    do ix = params%g+1, params%Bs(1)+params%g
                        ! First derivatives
                        u_dx = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix+FD1_rs:ix+FD1_re,iy,iz,1,hvy_id)) * dx_inv
                        v_dx = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix+FD1_rs:ix+FD1_re,iy,iz,2,hvy_id)) * dx_inv
                        w_dx = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix+FD1_rs:ix+FD1_re,iy,iz,3,hvy_id)) * dx_inv

                        u_dy = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy+FD1_rs:iy+FD1_re,iz,1,hvy_id)) * dy_inv
                        v_dy = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy+FD1_rs:iy+FD1_re,iz,2,hvy_id)) * dy_inv
                        w_dy = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy+FD1_rs:iy+FD1_re,iz,3,hvy_id)) * dy_inv

                        u_dz = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy,iz+FD1_rs:iz+FD1_re,1,hvy_id)) * dz_inv
                        v_dz = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy,iz+FD1_rs:iz+FD1_re,2,hvy_id)) * dz_inv
                        w_dz = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy,iz+FD1_rs:iz+FD1_re,3,hvy_id)) * dz_inv

                        ! Non-linear energy terms
                        uu_dx = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,iz,1,hvy_id) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,iz,1,hvy_id)) * dx_inv
                        uv_dy = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,iz,1,hvy_id) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,iz,2,hvy_id)) * dy_inv
                        uw_dz = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy,iz+FD1_ls:iz+FD1_le,1,hvy_id) * hvy_u(ix,iy,iz+FD1_ls:iz+FD1_le,3,hvy_id)) * dz_inv

                        vu_dx = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,iz,2,hvy_id) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,iz,1,hvy_id)) * dx_inv
                        vv_dy = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,iz,2,hvy_id) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,iz,2,hvy_id)) * dy_inv
                        vw_dz = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy,iz+FD1_ls:iz+FD1_le,2,hvy_id) * hvy_u(ix,iy,iz+FD1_ls:iz+FD1_le,3,hvy_id)) * dz_inv

                        wu_dx = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,iz,3,hvy_id) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,iz,1,hvy_id)) * dx_inv
                        wv_dy = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,iz,3,hvy_id) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,iz,2,hvy_id)) * dy_inv
                        ww_dz = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy,iz+FD1_ls:iz+FD1_le,3,hvy_id) * hvy_u(ix,iy,iz+FD1_ls:iz+FD1_le,3,hvy_id)) * dz_inv

                        hvy_NL(ix,iy,iz,1,hvy_id) = -0.5_rk * (uu_dx + uv_dy + uw_dz + &
                            hvy_u(ix,iy,iz,1,hvy_id)*u_dx + hvy_u(ix,iy,iz,2,hvy_id)*u_dy + hvy_u(ix,iy,iz,3,hvy_id)*u_dz)
                        hvy_NL(ix,iy,iz,2,hvy_id) = -0.5_rk * (vu_dx + vv_dy + vw_dz + &
                            hvy_u(ix,iy,iz,1,hvy_id)*v_dx + hvy_u(ix,iy,iz,2,hvy_id)*v_dy + hvy_u(ix,iy,iz,3,hvy_id)*v_dz)
                        hvy_NL(ix,iy,iz,3,hvy_id) = -0.5_rk * (wu_dx + wv_dy + ww_dz + &
                            hvy_u(ix,iy,iz,1,hvy_id)*w_dx + hvy_u(ix,iy,iz,2,hvy_id)*w_dy + hvy_u(ix,iy,iz,3,hvy_id)*w_dz)
                    enddo
                enddo
            enddo
        endif
    enddo

end subroutine compute_NonLinear


!> \brief computes the divergence of a vector field
subroutine compute_divergence_tree(params, hvy_u, hvy_div, order_discretization, treeID)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    type(type_params), intent(in) :: params
    real(kind=rk), intent(in) :: hvy_u(:,:,:,:,:)
    real(kind=rk), intent(out) :: hvy_div(:,:,:,:)
    character(len=cshort), intent(in) :: order_discretization
    integer(kind=ik), intent(in) :: treeID

    integer(kind=ik) :: k_block, hvy_id, lgt_ID
    real(kind=rk) :: dx(1:3), x0(1:3)

    ! loop over all blocks
    do k_block = 1, hvy_n(treeID)
        hvy_id = hvy_active(k_block, treeID)
        call hvy2lgt(lgt_ID, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin(params, lgt_ID, x0, dx)

        call compute_divergence(hvy_u(:,:,:,1:params%dim,hvy_id), dx, params%Bs, params%g, order_discretization, hvy_div(:,:,:,hvy_id))
    enddo
end subroutine compute_divergence_tree


subroutine compute_NSI_RHS(params, hvy_u, hvy_mask, hvy_RHS, order_discretization, treeID, nu, C_eta)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    type(type_params), intent(in) :: params
    real(kind=rk), intent(in) :: hvy_u(:,:,:,:,:)
    real(kind=rk), intent(in) :: hvy_mask(:,:,:,:,:)
    real(kind=rk), intent(out) :: hvy_RHS(:,:,:,:,:)
    character(len=cshort), intent(in) :: order_discretization
    integer(kind=ik), intent(in) :: treeID
    real(kind=rk), intent(in), optional :: nu, C_eta

    integer(kind=ik) :: k_block, hvy_id, lgt_ID, ix, iy, iz
    real(kind=rk) :: dx(1:3), x0(1:3), dx_inv, dy_inv, dz_inv, dx2_inv, dy2_inv, dz2_inv, C_eta_inv, nu_set
    real(kind=rk) :: u_dx, u_dy, u_dz, u_dxdx, u_dydy, u_dzdz, u_dxdy, u_dxdz, &
                     v_dx, v_dy, v_dz, v_dxdx, v_dydy, v_dzdz, v_dxdy, v_dydz, &
                     w_dx, w_dy, w_dz, w_dxdx, w_dydy, w_dzdz, w_dxdz, w_dydz, &
                     penalx, penaly, penalz, u, v, w, &
                     uu_dx, uv_dy, uw_dz, vu_dx, vv_dy, vw_dz, wu_dx, wv_dy, ww_dz
    
    !> parameters for FD1_l operator
    real(kind=rk), allocatable, dimension(:) :: FD1_l, FD1_r, FD2
    integer(kind=ik) :: FD1_ls, FD1_le, FD1_rs, FD1_re, FD2_s, FD2_e

    ! Determine the viscosity value to use
    if (present(nu)) then
        nu_set = nu
    else
        nu_set = 1.0e-5_rk
    endif

    ! Setup stencils using the unified interface from module_operators
    call setup_FD1_left_stencil(order_discretization, FD1_l, FD1_ls, FD1_le)
    call setup_FD1_right_stencil(order_discretization, FD1_r, FD1_rs, FD1_re)
    
    ! Performance optimization: skip FD2 setup if viscosity is zero
    if (abs(nu_set) < 1.0e-14_rk) then
        allocate(FD2(0:0))
        FD2(0) = 0.0_rk
        FD2_s = 0
        FD2_e = 0
    else
        call setup_FD2_stencil(order_discretization, FD2, FD2_s, FD2_e)
    endif
    
    ! loop over all blocks
    do k_block = 1, hvy_n(treeID)
        hvy_id = hvy_active(k_block, treeID)
        call hvy2lgt(lgt_ID, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin(params, lgt_ID, x0, dx)
        dx_inv = 1.0_rk / dx(1)
        dy_inv = 1.0_rk / dx(2)
        dx2_inv = 1.0_rk / (dx(1)**2)
        dy2_inv = 1.0_rk / (dx(2)**2)
        if (params%dim == 3) then
            dz_inv = 1.0_rk / dx(3)
            dz2_inv = 1.0_rk / (dx(3)**2)
        else
            dz_inv = 0.0_rk ! no z direction in 2D
            dz2_inv = 0.0_rk ! no z direction in 2D
        endif
        ! HACK ! This is part of module ACM and cannot be read normaly, so we read it in by hand
        ! C_eta_inv = 1.0_rk / params%C_eta
        ! nu = params%nu
        if (present(C_eta)) then
            C_eta_inv = 1.0_rk / C_eta
        else
            C_eta_inv = 1.0_rk
        endif

        ! compute RHS -1/2 (\nabla \cdot u + u \cdot \nabla) u + \nu \Delta u - mask/C_eta (u - u_s)
        if (params%dim == 2) then
            do iy = params%g+1, params%Bs(2)+params%g
                do ix = params%g+1, params%Bs(1)+params%g
                    ! Velocity
                    u = hvy_u(ix, iy, 1, 1, hvy_id)
                    v = hvy_u(ix, iy, 1, 2, hvy_id)

                    ! First derivatives
                    u_dx   = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix+FD1_rs:ix+FD1_re,iy,1,1,hvy_id)) * dx_inv
                    v_dx   = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix+FD1_rs:ix+FD1_re,iy,1,2,hvy_id)) * dx_inv

                    u_dy   = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy+FD1_rs:iy+FD1_re,1,1,hvy_id)) * dy_inv
                    v_dy   = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy+FD1_rs:iy+FD1_re,1,2,hvy_id)) * dy_inv

                    ! Non-linear energy terms
                    uu_dx = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,1,1,hvy_id) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,1,1,hvy_id)) * dx_inv
                    uv_dy = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,1,1,hvy_id) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,1,2,hvy_id)) * dy_inv

                    vu_dx = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,1,2,hvy_id) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,1,1,hvy_id)) * dx_inv
                    vv_dy = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,1,2,hvy_id) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,1,2,hvy_id)) * dy_inv

                    ! Second derivatives
                    u_dxdx = sum(FD2(FD2_s:FD2_e) * hvy_u(ix+FD2_s:ix+FD2_e,iy,1,1,hvy_id)) * dx2_inv
                    v_dxdx = sum(FD2(FD2_s:FD2_e) * hvy_u(ix+FD2_s:ix+FD2_e,iy,1,2,hvy_id)) * dx2_inv

                    u_dydy = sum(FD2(FD2_s:FD2_e) * hvy_u(ix,iy+FD2_s:iy+FD2_e,1,1,hvy_id)) * dy2_inv
                    v_dydy = sum(FD2(FD2_s:FD2_e) * hvy_u(ix,iy+FD2_s:iy+FD2_e,1,2,hvy_id)) * dy2_inv

                    ! Penalization terms
                    penalx = -hvy_mask(ix,iy,1,1,hvy_id) * C_eta_inv * (u - hvy_mask(ix,iy,1,2,hvy_id))
                    penaly = -hvy_mask(ix,iy,1,1,hvy_id) * C_eta_inv * (v - hvy_mask(ix,iy,1,3,hvy_id))

                    ! Fill RHS
                    hvy_RHS(ix,iy,1,1,hvy_id) = -0.5_rk*(uu_dx + uv_dy   + u*u_dx + v*u_dy ) + nu_set*(u_dxdx + u_dydy ) + penalx
                    hvy_RHS(ix,iy,1,2,hvy_id) = -0.5_rk*(vu_dx + vv_dy   + u*v_dx + v*v_dy ) + nu_set*(v_dxdx + v_dydy ) + penaly
                enddo
            enddo
        elseif (params%dim == 3) then
            do iz = params%g+1, params%Bs(3)+params%g
                do iy = params%g+1, params%Bs(2)+params%g
                    do ix = params%g+1, params%Bs(1)+params%g
                        ! Velocity
                        u = hvy_u(ix, iy, iz, 1, hvy_id)
                        v = hvy_u(ix, iy, iz, 2, hvy_id)
                        w = hvy_u(ix, iy, iz, 3, hvy_id)

                        ! First derivatives
                        u_dx = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix+FD1_rs:ix+FD1_re,iy,iz,1,hvy_id)) * dx_inv
                        v_dx = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix+FD1_rs:ix+FD1_re,iy,iz,2,hvy_id)) * dx_inv
                        w_dx = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix+FD1_rs:ix+FD1_re,iy,iz,3,hvy_id)) * dx_inv

                        u_dy = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy+FD1_rs:iy+FD1_re,iz,1,hvy_id)) * dy_inv
                        v_dy = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy+FD1_rs:iy+FD1_re,iz,2,hvy_id)) * dy_inv
                        w_dy = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy+FD1_rs:iy+FD1_re,iz,3,hvy_id)) * dy_inv

                        u_dz = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy,iz+FD1_rs:iz+FD1_re,1,hvy_id)) * dz_inv
                        v_dz = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy,iz+FD1_rs:iz+FD1_re,2,hvy_id)) * dz_inv
                        w_dz = sum(FD1_r(FD1_rs:FD1_re) * hvy_u(ix,iy,iz+FD1_rs:iz+FD1_re,3,hvy_id)) * dz_inv

                        ! Non-linear energy terms
                        uu_dx = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,iz,1,hvy_id) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,iz,1,hvy_id)) * dx_inv
                        uv_dy = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,iz,1,hvy_id) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,iz,2,hvy_id)) * dy_inv
                        uw_dz = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy,iz+FD1_ls:iz+FD1_le,1,hvy_id) * hvy_u(ix,iy,iz+FD1_ls:iz+FD1_le,3,hvy_id)) * dz_inv

                        vu_dx = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,iz,2,hvy_id) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,iz,1,hvy_id)) * dx_inv
                        vv_dy = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,iz,2,hvy_id) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,iz,2,hvy_id)) * dy_inv
                        vw_dz = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy,iz+FD1_ls:iz+FD1_le,2,hvy_id) * hvy_u(ix,iy,iz+FD1_ls:iz+FD1_le,3,hvy_id)) * dz_inv

                        wu_dx = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,iz,3,hvy_id) * hvy_u(ix+FD1_ls:ix+FD1_le,iy,iz,1,hvy_id)) * dx_inv
                        wv_dy = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,iz,3,hvy_id) * hvy_u(ix,iy+FD1_ls:iy+FD1_le,iz,2,hvy_id)) * dy_inv
                        ww_dz = sum(FD1_l(FD1_ls:FD1_le) * hvy_u(ix,iy,iz+FD1_ls:iz+FD1_le,3,hvy_id) * hvy_u(ix,iy,iz+FD1_ls:iz+FD1_le,3,hvy_id)) * dz_inv

                        ! Second derivatives
                        u_dxdx = sum(FD2(FD2_s:FD2_e) * hvy_u(ix+FD2_s:ix+FD2_e,iy,iz,1,hvy_id)) * dx2_inv
                        u_dydy = sum(FD2(FD2_s:FD2_e) * hvy_u(ix,iy+FD2_s:iy+FD2_e,iz,1,hvy_id)) * dy2_inv
                        u_dzdz = sum(FD2(FD2_s:FD2_e) * hvy_u(ix,iy,iz+FD2_s:iz+FD2_e,1,hvy_id)) * dz2_inv

                        v_dxdx = sum(FD2(FD2_s:FD2_e) * hvy_u(ix+FD2_s:ix+FD2_e,iy,iz,2,hvy_id)) * dx2_inv
                        v_dydy = sum(FD2(FD2_s:FD2_e) * hvy_u(ix,iy+FD2_s:iy+FD2_e,iz,2,hvy_id)) * dy2_inv
                        v_dzdz = sum(FD2(FD2_s:FD2_e) * hvy_u(ix,iy,iz+FD2_s:iz+FD2_e,2,hvy_id)) * dz2_inv

                        w_dxdx = sum(FD2(FD2_s:FD2_e) * hvy_u(ix+FD2_s:ix+FD2_e,iy,iz,3,hvy_id)) * dx2_inv
                        w_dydy = sum(FD2(FD2_s:FD2_e) * hvy_u(ix,iy+FD2_s:iy+FD2_e,iz,3,hvy_id)) * dy2_inv
                        w_dzdz = sum(FD2(FD2_s:FD2_e) * hvy_u(ix,iy,iz+FD2_s:iz+FD2_e,3,hvy_id)) * dz2_inv

                        ! Penalization terms
                        penalx = -hvy_mask(ix,iy,iz,1,hvy_id) * C_eta_inv * (u - hvy_mask(ix,iy,iz,2,hvy_id))
                        penaly = -hvy_mask(ix,iy,iz,1,hvy_id) * C_eta_inv * (v - hvy_mask(ix,iy,iz,3,hvy_id))
                        penalz = -hvy_mask(ix,iy,iz,1,hvy_id) * C_eta_inv * (w - hvy_mask(ix,iy,iz,4,hvy_id))

                        ! Fill RHS
                        hvy_RHS(ix,iy,iz,1,hvy_id) = -0.5_rk * (uu_dx + uv_dy + uw_dz + u*u_dx + v*u_dy + w*u_dz) &
                            + nu_set*(u_dxdx + u_dydy + u_dzdz) + penalx

                        hvy_RHS(ix,iy,iz,2,hvy_id) = -0.5_rk * (vu_dx + vv_dy + vw_dz + u*v_dx + v*v_dy + w*v_dz) &
                            + nu_set*(v_dxdx + v_dydy + v_dzdz) + penaly

                        hvy_RHS(ix,iy,iz,3,hvy_id) = -0.5_rk * (wu_dx + wv_dy + ww_dz + u*w_dx + v*w_dy + w*w_dz) &
                            + nu_set*(w_dxdx + w_dydy + w_dzdz) + penalz
                    enddo
                enddo
            enddo
        endif
    enddo

end subroutine compute_NSI_RHS


!> \brief computes the divergence of a vector field
subroutine compute_projection(params, hvy_u, hvy_p, hvy_u_div0, order_discretization, treeID, dt)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    type(type_params), intent(in) :: params
    real(kind=rk), intent(in) :: hvy_u(:,:,:,:,:)
    real(kind=rk), intent(in) :: hvy_p(:,:,:,:,:)
    real(kind=rk), intent(out) :: hvy_u_div0(:,:,:,:,:)
    character(len=cshort), intent(in) :: order_discretization
    integer(kind=ik), intent(in) :: treeID
    real(kind=rk), intent(in) :: dt

    integer(kind=ik) :: k_block, hvy_id, lgt_ID, ix, iy, iz
    real(kind=rk) :: dx(1:3), x0(1:3), dx_inv, dy_inv, dz_inv
    real(kind=rk) :: p_dx, p_dy, p_dz

    !> parameters for FD1_r operator
    real(kind=rk), allocatable, dimension(:) :: FD1_l, FD1_r
    integer(kind=ik) :: FD1_ls, FD1_le, FD1_rs, FD1_re

    ! Setup stencils using the unified interface from module_operators
    call setup_FD1_left_stencil(order_discretization, FD1_l, FD1_ls, FD1_le)
    call setup_FD1_right_stencil(order_discretization, FD1_r, FD1_rs, FD1_re)

    ! write(*,'(A,10(es10.3))') "Left stencil FD1: ", FD1_l(FD1_ls:FD1_le)
    ! write(*,'(A,10(es10.3))') "Right stencil FD1: ", FD1_r(FD1_rs:FD1_re)
    
    ! loop over all blocks
    do k_block = 1, hvy_n(treeID)
        hvy_id = hvy_active(k_block, treeID)
        call hvy2lgt(lgt_ID, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin(params, lgt_ID, x0, dx)
        dx_inv = 1.0_rk / dx(1)
        dy_inv = 1.0_rk / dx(2)
        if (params%dim == 3) then
            dz_inv = 1.0_rk / dx(3)
        else
            dz_inv = 0.0_rk ! no z direction in 2D
        endif

        ! compute projection u = u* - dt * \nabla p
        if (params%dim == 2) then
            do iy = params%g+1, params%Bs(2)+params%g
                do ix = params%g+1, params%Bs(1)+params%g
                    p_dx = sum(FD1_r(FD1_rs:FD1_re) * hvy_p(ix+FD1_rs:ix+FD1_re,iy,1,1,hvy_id)) * dx_inv
                    p_dy = sum(FD1_r(FD1_rs:FD1_re) * hvy_p(ix,iy+FD1_rs:iy+FD1_re,1,1,hvy_id)) * dy_inv

                    hvy_u_div0(ix,iy,1,1,hvy_id) = hvy_u(ix,iy,1,1,hvy_id) - dt*p_dx
                    hvy_u_div0(ix,iy,1,2,hvy_id) = hvy_u(ix,iy,1,2,hvy_id) - dt*p_dy
                enddo
            enddo
        elseif (params%dim == 3) then
            do iz = params%g+1, params%Bs(3)+params%g
                do iy = params%g+1, params%Bs(2)+params%g
                    do ix = params%g+1, params%Bs(1)+params%g
                        p_dx = sum(FD1_r(FD1_rs:FD1_re) * hvy_p(ix+FD1_rs:ix+FD1_re,iy,iz,1,hvy_id)) * dx_inv
                        p_dy = sum(FD1_r(FD1_rs:FD1_re) * hvy_p(ix,iy+FD1_rs:iy+FD1_re,iz,1,hvy_id)) * dy_inv
                        p_dz = sum(FD1_r(FD1_rs:FD1_re) * hvy_p(ix,iy,iz+FD1_rs:iz+FD1_re,1,hvy_id)) * dz_inv

                        hvy_u_div0(ix,iy,iz,1,hvy_id) = hvy_u(ix,iy,iz,1,hvy_id) - dt*p_dx
                        hvy_u_div0(ix,iy,iz,2,hvy_id) = hvy_u(ix,iy,iz,2,hvy_id) - dt*p_dy
                        hvy_u_div0(ix,iy,iz,3,hvy_id) = hvy_u(ix,iy,iz,3,hvy_id) - dt*p_dz
                    enddo
                enddo
            enddo
        endif
    enddo
end subroutine compute_projection