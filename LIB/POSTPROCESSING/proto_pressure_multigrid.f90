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
            character(len=*), intent(in) :: order_discretization
            integer, intent(in) :: treeID
        end subroutine compute_NonLinear

        subroutine compute_divergence(params, hvy_u, hvy_div, order_discretization, treeID)
            use module_globals
            use module_params
            type(type_params), intent(in) :: params
            real(kind=rk), intent(in) :: hvy_u(:,:,:,:,:)
            real(kind=rk), intent(out) :: hvy_div(:,:,:,:)
            character(len=*), intent(in) :: order_discretization
            integer, intent(in) :: treeID
        end subroutine compute_divergence
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
    real(kind=rk)                      :: x0(1:3), dx(1:3), domain(1:3), norm(1:6), volume
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
            write(*, '(A)') "    ux.h5 = Velocity in Y-direction"
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
    call get_command_argument(6 + merge(1,0, exist_p), params%laplacian_order)
    if (params%laplacian_order == "CFD_2nd") then
        params%laplacian_stencil_size = 1
    elseif (params%laplacian_order == "CFD_4th") then
        params%laplacian_stencil_size = 2
    elseif (params%laplacian_order == "CFD_6th") then
        params%laplacian_stencil_size = 3
    elseif (params%laplacian_order == "CFD_8th") then
        params%laplacian_stencil_size = 4
    elseif (params%laplacian_order == "MST_6th") then
        params%laplacian_stencil_size = 1
    else
        call abort(250612, "ERROR: laplace order not recognized. Use CFD_2nd, CFD_4th, CFD_6th, CFD_8th or MST_6th")
    endif
        ! set number of cycles
    call get_command_argument(7 + merge(1,0, exist_p), cycle_type)
    read (cycle_type, *) params%laplacian_cycle_it
    call get_command_argument(8 + merge(1,0, exist_p), cycle_type)
    read (cycle_type, *) params%laplacian_GS_it

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
    select case(params%laplacian_order)
    case ("CFD_2nd", "CFD_4th", "CFD_6th", "CFD_8th")
        params%order_discretization = "FD_" // params%laplacian_order(5:7) // "_central"
    case ("MST_6th")
        params%order_discretization = "FD_6th_central"
    end select
    params%laplacian_coarsest = "FFT"
    params%FFT_accuracy = "FD"  ! FD or spectral

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

    ! read input data, if comparison file is given, we read it in as well
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

    call componentWiseNorm_tree(params, hvy_block(:,:,:,1:size(hvy_block,4),:), tree_ID, "Mean", norm(1:size(hvy_block,4)))
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
    call compute_divergence(params, hvy_tmp(:,:,:,1:3,:), hvy_block(:,:,:,1,:), params%order_discretization, tree_ID)
    call sync_ghosts_tree(params, hvy_block(:,:,:,1:1,:), tree_ID)
    call toc( "Compute Divergence", 9999, MPI_Wtime()-t_block )

    call componentWiseNorm_tree(params, hvy_block(:,:,:,1:size(hvy_block,4),:), tree_ID, "Mean", norm(1:size(hvy_block,4)))
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
    if (params%laplacian_order == "MST_6th") then
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
    do i_cycle = 1, params%laplacian_cycle_it

        ! ! usually we want spectral resolution for the coarsest level, for debugging we can disable this
        ! fft_order_discretization = "spectral"

        call multigrid_vcycle(params, hvy_tmp(:,:,:,1:1,:), hvy_block(:,:,:,1:1,:), hvy_tmp(:,:,:,2:size(hvy_tmp,4),:), tree_ID)

        ! laplacian is invariant to shifts of constant values
        ! our values are defined with zero mean for comparison
        ! as multigrid might accidently introduce a constant offset, we remove it
        call componentWiseNorm_tree(params, hvy_tmp(:,:,:,1:1,:), tree_ID, "Mean", norm(1:1))
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
            call componentWiseNorm_tree(params, hvy_work(:,:,:,1:1,:,1), tree_ID, "L2", norm(1:1))
            if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Diff spectral L2: ", norm(1), " ---"
            call componentWiseNorm_tree(params, hvy_work(:,:,:,1:1,:,1), tree_ID, "L1", norm(2:2))
            if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Diff spectral L1: ", norm(2), " ---"
            call componentWiseNorm_tree(params, hvy_work(:,:,:,1:1,:,1), tree_ID, "Linfty", norm(3:3))
            if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Diff spectral Linfty: ", norm(3), " ---"

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

    integer(kind=ik) :: k_block, hvy_id, lgt_ID, ix, iy, iz, a, ia
    real(kind=rk) :: dx(1:3), x0(1:3), dx_inv, dy_inv, dz_inv
    real(kind=rk) :: u_dx, v_dx, w_dx, u_dy, v_dy, w_dy, u_dz, v_dz, w_dz
    real(kind=rk) :: uu_dx, uv_dy, uw_dz, vu_dx, vv_dy, vw_dz, wu_dx, wv_dy, ww_dz

    real(kind=rk), parameter :: a_TW4(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    real(kind=rk), parameter :: a_TWR4(-3:3) = (/-0.020843142770_rk, +0.166705904415_rk, -0.770882380518_rk, 0.0_rk, 0.770882380518_rk, -0.166705904415_rk, 0.020843142770_rk/)
    real(kind=rk), parameter :: a_FD2(-1:1) = (/ -0.5_rk, 0.0_rk, +0.5_rk /)
    real(kind=rk), parameter :: a_FD4(-2:2) = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)
    real(kind=rk), parameter :: a_FD6(-3:3) = (/-1.0_rk/60.0_rk, 3.0_rk/20.0_rk, -3.0_rk/4.0_rk, 0.0_rk, 3.0_rk/4.0_rk, -3.0_rk/20.0_rk, 1.0_rk/60.0_rk/)
    real(kind=rk), parameter :: a_FD8(-4:4) = (/1.0_rk/280.0_rk, -4.0_rk/105.0_rk, 1.0_rk/5.0_rk, -4.0/5.0_rk, 0.0_rk, 4.0_rk/5.0_rk, -1.0_rk/5.0_rk, 4.0_rk/105.0_rk, -1.0_rk/280.0_rk/)

    real(kind=rk), allocatable :: a_FD(:)

    select case(order_discretization)
        case("FD_2nd_central")
            a = 1
            allocate(a_FD(-a:a))
            a_FD(:) = a_FD2(:)
        case("FD_4th_central")
            a = 2
            allocate(a_FD(-a:a))
            a_FD(:) = a_FD4(:)
        case("FD_6th_central")
            a = 3
            allocate(a_FD(-a:a))
            a_FD(:) = a_FD6(:)
        case("FD_8th_central")
            a = 4
            allocate(a_FD(-a:a))
            a_FD(:) = a_FD8(:)
        case("TW_4th_central")
            a = 3
            allocate(a_FD(-a:a))
            a_FD(:) = a_TW4(:)
        case("TWR_4th_central")
            a = 3
            allocate(a_FD(-a:a))
            a_FD(:) = a_TWR4(:)
        case default
            call abort(250615, "ERROR: order of discretization not known: " // trim(order_discretization))
    end select

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
                    u_dx   = sum(a_FD(-a:a) * hvy_u(ix-a:ix+a,iy,1,1,hvy_id)) * dx_inv
                    v_dx   = sum(a_FD(-a:a) * hvy_u(ix-a:ix+a,iy,1,2,hvy_id)) * dx_inv

                    u_dy   = sum(a_FD(-a:a) * hvy_u(ix,iy-a:iy+a,1,1,hvy_id)) * dy_inv
                    v_dy   = sum(a_FD(-a:a) * hvy_u(ix,iy-a:iy+a,1,2,hvy_id)) * dy_inv

                    uu_dx = sum(a_FD(-a:a) * hvy_u(ix-a:ix+a,iy,1,1,hvy_id) * hvy_u(ix-a:ix+a,iy,1,1,hvy_id)) * dx_inv
                    uv_dy = sum(a_FD(-a:a) * hvy_u(ix,iy-a:iy+a,1,1,hvy_id) * hvy_u(ix,iy-a:iy+a,1,2,hvy_id)) * dy_inv

                    vu_dx = sum(a_FD(-a:a) * hvy_u(ix-a:ix+a,iy,1,2,hvy_id) * hvy_u(ix-a:ix+a,iy,1,1,hvy_id)) * dx_inv
                    vv_dy = sum(a_FD(-a:a) * hvy_u(ix,iy-a:iy+a,1,2,hvy_id) * hvy_u(ix,iy-a:iy+a,1,2,hvy_id)) * dy_inv

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
                        u_dx = sum(a_FD(-a:a) * hvy_u(ix-a:ix+a,iy,iz,1,hvy_id)) * dx_inv
                        v_dx = sum(a_FD(-a:a) * hvy_u(ix-a:ix+a,iy,iz,2,hvy_id)) * dx_inv
                        w_dx = sum(a_FD(-a:a) * hvy_u(ix-a:ix+a,iy,iz,3,hvy_id)) * dx_inv

                        u_dy = sum(a_FD(-a:a) * hvy_u(ix,iy-a:iy+a,iz,1,hvy_id)) * dy_inv
                        v_dy = sum(a_FD(-a:a) * hvy_u(ix,iy-a:iy+a,iz,2,hvy_id)) * dy_inv
                        w_dy = sum(a_FD(-a:a) * hvy_u(ix,iy-a:iy+a,iz,3,hvy_id)) * dy_inv

                        u_dz = sum(a_FD(-a:a) * hvy_u(ix,iy,iz-a:iz+a,1,hvy_id)) * dz_inv
                        v_dz = sum(a_FD(-a:a) * hvy_u(ix,iy,iz-a:iz+a,2,hvy_id)) * dz_inv
                        w_dz = sum(a_FD(-a:a) * hvy_u(ix,iy,iz-a:iz+a,3,hvy_id)) * dz_inv

                        uu_dx = sum(a_FD(-a:a) * hvy_u(ix-a:ix+a,iy,iz,1,hvy_id) * hvy_u(ix-a:ix+a,iy,iz,1,hvy_id)) * dx_inv
                        uv_dy = sum(a_FD(-a:a) * hvy_u(ix,iy-a:iy+a,iz,1,hvy_id) * hvy_u(ix,iy-a:iy+a,iz,2,hvy_id)) * dy_inv
                        uw_dz = sum(a_FD(-a:a) * hvy_u(ix,iy,iz-a:iz+a,1,hvy_id) * hvy_u(ix,iy,iz-a:iz+a,3,hvy_id)) * dz_inv

                        vu_dx = sum(a_FD(-a:a) * hvy_u(ix-a:ix+a,iy,iz,2,hvy_id) * hvy_u(ix-a:ix+a,iy,iz,1,hvy_id)) * dx_inv
                        vv_dy = sum(a_FD(-a:a) * hvy_u(ix,iy-a:iy+a,iz,2,hvy_id) * hvy_u(ix,iy-a:iy+a,iz,2,hvy_id)) * dy_inv
                        vw_dz = sum(a_FD(-a:a) * hvy_u(ix,iy,iz-a:iz+a,2,hvy_id) * hvy_u(ix,iy,iz-a:iz+a,3,hvy_id)) * dz_inv

                        wu_dx = sum(a_FD(-a:a) * hvy_u(ix-a:ix+a,iy,iz,3,hvy_id) * hvy_u(ix-a:ix+a,iy,iz,1,hvy_id)) * dx_inv
                        wv_dy = sum(a_FD(-a:a) * hvy_u(ix,iy-a:iy+a,iz,3,hvy_id) * hvy_u(ix,iy-a:iy+a,iz,2,hvy_id)) * dy_inv
                        ww_dz = sum(a_FD(-a:a) * hvy_u(ix,iy,iz-a:iz+a,3,hvy_id) * hvy_u(ix,iy,iz-a:iz+a,3,hvy_id)) * dz_inv

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
subroutine compute_divergence(params, hvy_u, hvy_div, order_discretization, treeID)
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

    integer(kind=ik) :: k_block, hvy_id, lgt_ID, ix, iy, iz, a
    real(kind=rk) :: dx(1:3), x0(1:3), dx_inv, dy_inv, dz_inv
    real(kind=rk) :: u_dx, v_dy, w_dz

    real(kind=rk), parameter :: a_TW4(-3:3) = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    real(kind=rk), parameter :: a_TWR4(-3:3) = (/-0.020843142770_rk, +0.166705904415_rk, -0.770882380518_rk, 0.0_rk, 0.770882380518_rk, -0.166705904415_rk, 0.020843142770_rk/)
    real(kind=rk), parameter :: a_FD2(-1:1) = (/ -0.5_rk, 0.0_rk, +0.5_rk /)
    real(kind=rk), parameter :: a_FD4(-2:2) = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)
    real(kind=rk), parameter :: a_FD6(-3:3) = (/-1.0_rk/60.0_rk, 3.0_rk/20.0_rk, -3.0_rk/4.0_rk, 0.0_rk, 3.0_rk/4.0_rk, -3.0_rk/20.0_rk, 1.0_rk/60.0_rk/)
    real(kind=rk), parameter :: a_FD8(-4:4) = (/1.0_rk/280.0_rk, -4.0_rk/105.0_rk, 1.0_rk/5.0_rk, -4.0/5.0_rk, 0.0_rk, 4.0_rk/5.0_rk, -1.0_rk/5.0_rk, 4.0_rk/105.0_rk, -1.0_rk/280.0_rk/)
    real(kind=rk), allocatable :: a_FD(:)
    select case(order_discretization)
        case("FD_2nd_central")
            a = 1
            allocate(a_FD(-a:a))
            a_FD(:) = a_FD2(:)
        case("FD_4th_central")
            a = 2
            allocate(a_FD(-a:a))
            a_FD(:) = a_FD4(:)
        case("FD_6th_central")
            a = 3
            allocate(a_FD(-a:a))
            a_FD(:) = a_FD6(:)
        case("FD_8th_central")
            a = 4
            allocate(a_FD(-a:a))
            a_FD(:) = a_FD8(:)
        case("TW_4th_central")
            a = 3
            allocate(a_FD(-a:a))
            a_FD(:) = a_TW4(:)
        case("TWR_4th_central")
            a = 3
            allocate(a_FD(-a:a))
            a_FD(:) = a_TWR4(:)
        case default
            call abort(250615, "ERROR: order of discretization not known: " // trim(order_discretization))
    end select
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

        ! compute divergence
        if (params%dim == 2) then
            ! 2D: div u = u_dx + v_dy
            do iy = params%g+1, params%Bs(2)+params%g
                do ix = params%g+1, params%Bs(1)+params%g
                    u_dx = sum(a_FD(-a:a) * hvy_u(ix-a:ix+a,iy,1,1,hvy_id)) * dx_inv
                    v_dy = sum(a_FD(-a:a) * hvy_u(ix,iy-a:iy+a,1,2,hvy_id)) * dy_inv

                    hvy_div(ix,iy,1,hvy_id) = u_dx + v_dy
                enddo
            enddo
        elseif (params%dim == 3) then
            ! 3D: div u = u_dx + v_dy + w_dz
            do iz = params%g+1, params%Bs(3)+params%g
                do iy = params%g+1, params%Bs(2)+params%g
                    do ix = params%g+1, params%Bs(1)+params%g
                        u_dx = sum(a_FD(-a:a) * hvy_u(ix-a:ix+a,iy,iz,1,hvy_id)) * dx_inv
                        v_dy = sum(a_FD(-a:a) * hvy_u(ix,iy-a:iy+a,iz,2,hvy_id)) * dy_inv
                        w_dz = sum(a_FD(-a:a) * hvy_u(ix,iy,iz-a:iz+a,3,hvy_id)) * dz_inv

                        hvy_div(ix,iy,iz,hvy_id) = u_dx + v_dy + w_dz
                    enddo
                enddo
            enddo
        endif
    enddo
end subroutine compute_divergence