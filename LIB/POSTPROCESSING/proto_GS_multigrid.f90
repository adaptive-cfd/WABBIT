!> \brief this is a prototype - it tries out the Gauss-Seidel Multigrid for equidistant grids
!-----------------------------------------------------------------------------------------------------


subroutine proto_GS_multigrid(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData
    use module_fft
    use module_poisson

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: file_b, file_u, file_uFD
    real(kind=rk)                      :: dt, t, time
    integer(kind=ik)                   :: k_block, lgt_ID, hvy_id, Bs(1:3)

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :), hvy_work(:, :, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, Jmin, ic, nc, a, i_cycle, it, tc_length, mpierr, g(1:3)

    character(len=cshort)              :: fname, cycle_type, fft_order_discretization
    logical                            :: exist_u, exist_uFD
    integer(kind=ik)                   :: laplace_order, it_cycle, it_down, it_up, it_coarsest
    real(kind=rk)                      :: x0(1:3), dx(1:3), domain(1:3), norm(1:6), volume, tol_cg
    integer(kind=tsize)                :: treecode

    real(kind=rk)        :: t_block, t_loop, t_cycle

    ! FD central differences for 2nd derivative of order 2,4,6,8 
    real(kind=rk), allocatable :: dx2_FD(:)

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (fname and other options)
    call get_command_argument(2, file_b)

    ! does the user need help?
    if (file_b=='--help' .or. file_b=='--h' .or. file_b=='-h') then
        if (params%rank==0) then
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Wabbit prototype: Multi-Grid Gauss-Seidel solver"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Read in a data field (2D or 3D) b and solves the system laplacian u = b"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') "./wabbit-post --proto-GS-multigrid b.h5 [u.h5] [u_FD.h5] WAVELET LAPLACE-ORDER CYCLE IT_CYCLE IT_DOWN IT_UP"
            write(*, '(A)') "    b.h5 = RHS field to be solved"
            write(*, '(A)') "    u.h5 = solution field for comparison (optional)"
            write(*, '(A)') "    u_FD.h5 = solution field with FD accuracy for comparison (optional)"
            write(*, '(A)') "    WAVELET = CDFXY"
            write(*, '(A)') "    LAPLACE-ORDER = 2, 4, 6, 8 for the order of the Laplace operator"
            write(*, '(A)') "    CYCLE = 'V' or 'F' for V-cycle or F-cycle"
            write(*, '(A)') "    IT_CYCLE = Number of cycles to perform"
            write(*, '(A)') "    IT_DOWN = Number of GS-sweeps for downwards iterations"
            write(*, '(A)') "    IT_UP = Number of GS-sweeps for upwards iterations"
            write(*, '(A)') "-----------------------------------------------------------"
        end if
        return
    endif
    call check_file_exists(trim(file_b))

    call get_command_argument(3, file_u)
    inquire( file=file_u, exist=exist_u )
    if (exist_u) then
        if (params%rank == 0) write(*, '(A, A, A)') "Found file ", trim(file_u), " for comparison."
        call get_command_argument(4, file_uFD)
        inquire( file=file_uFD, exist=exist_uFD )
        if (exist_uFD) then
            if (params%rank == 0) write(*, '(A, A, A)') "Found file ", trim(file_uFD), " for FD accuracy comparison."
        else
            if (params%rank == 0) write(*, '(A)') "No FD accuracy comparison will be made, as no file was provided."
        endif
    else
        if (params%rank == 0) write(*, '(A)') "No comparison will be made, as no file was provided."
        exist_uFD = .false.
    endif

    call get_command_argument(3 + merge(1,0, exist_u) + merge(1,0, exist_uFD), params%wavelet)

    ! set laplace order
    call get_command_argument(4 + merge(1,0, exist_u) + merge(1,0, exist_uFD), cycle_type)
    read (cycle_type, *) laplace_order
    a = laplace_order / 2
    allocate( dx2_FD(-a:a) )
    if (laplace_order == 2) then
        fft_order_discretization = "FD_2nd_central"
        dx2_FD = (/ 1.0_rk, -2.0_rk, 1.0_rk /)
    elseif (laplace_order == 4) then
        fft_order_discretization = "FD_4th_central"
        dx2_FD = (/-1.0_rk/12.0_rk,  4.0_rk/3.0_rk,  -5.0_rk/2.0_rk,  4.0_rk/3.0_rk,  -1.0_rk/12.0_rk /)
    elseif (laplace_order == 6) then
        fft_order_discretization = "FD_6th_central"
        dx2_FD = (/ 1.0_rk/90.0_rk, -3.0_rk/20.0_rk,  3.0_rk/2.0_rk,-49.0_rk/18.0_rk,  3.0_rk/2.0_rk, -3.0_rk/20.0_rk, 1.0_rk/90.0_rk/)
    elseif (laplace_order == 8) then
        fft_order_discretization = "FD_8th_central"
        dx2_FD = (/-1.0_rk/560.0_rk, 8.0_rk/315.0_rk,-1.0_rk/5.0_rk,  8.0_rk/5.0_rk,-205.0_rk/72.0_rk, 8.0_rk/5.0_rk, -1.0_rk/5.0_rk, 8.0_rk/315.0_rk, -1.0_rk/560.0_rk /)
    endif
        ! set number of cycles
    call get_command_argument(6 + merge(1,0, exist_u) + merge(1,0, exist_uFD), cycle_type)
    read (cycle_type, *) it_cycle
    call get_command_argument(7 + merge(1,0, exist_u) + merge(1,0, exist_uFD), cycle_type)
    read (cycle_type, *) it_down
    call get_command_argument(8 + merge(1,0, exist_u) + merge(1,0, exist_uFD), cycle_type)
    read (cycle_type, *) it_up
    ! set cycle type
    call get_command_argument(5 + merge(1,0, exist_u) + merge(1,0, exist_uFD), cycle_type)

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file_b, lgt_n(tree_ID), time, it, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    params%Jmax = tc_length
    params%Jmin = 0   ! yes, I want to go down to one block only!
    params%n_eqn = params%dim
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs
    allocate(params%butcher_tableau(1,1))

    ! set only one variable for now
    params%n_eqn = 1 + merge(1,0, exist_u) + merge(1,0, exist_uFD)
    nc = 1

    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component(1) = "0"
    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component = 1
    params%order_discretization = "FD_4th_central"

    Bs = params%Bs

    call setup_wavelet(params, params%g)

    t_block = MPI_Wtime()
    call fft_initialize(params)
    call toc( "fft initialize", 10101, MPI_Wtime()-t_block )

    ! we need to decompose values, so we need atleast 2**dim/2**(dim-1) blocks
    params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) * 2**params%dim/(2**params%dim-1) + 8 )

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_work=hvy_work, nrhs_slots1=1, neqn_hvy_tmp=5)

    ! read input data, if comparison file is given, we read it in as well
    if (exist_u) then
        if (exist_uFD) then
            call readHDF5vct_tree( (/file_b, file_u, file_uFD/), params, hvy_block, tree_ID)
        else
            call readHDF5vct_tree( (/file_b, file_u/), params, hvy_block, tree_ID)
        endif
    else
        call readHDF5vct_tree( (/file_b/), params, hvy_block, tree_ID)
    endif

    ! debug max-min
    ! do k = 1, hvy_n(tree_ID)
    !     hvy_id = hvy_active(k, tree_ID)
    !     call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
    !     write(*, '(A, i0, 2(A, es9.2))') "Block ", lgt_id, " : ", minval(hvy_block(:,:,:,:,hvy_id)), " - ", maxval(hvy_block(:,:,:,:,hvy_id))
    ! enddo

    ! maye we didn't sync
    call sync_ghosts_tree( params, hvy_block, tree_ID )

    ! balance the load
    call balanceLoad_tree(params, hvy_block, tree_ID)

    ! prepare full tree grid, this will be populated along the way
    ! blocks are practically empty, but we fill them along the way so ref flag will be set to 0 to allow synching
    call init_full_tree(params, tree_ID, set_ref=0)

    ! ------------------------------------------------------------
    ! Gauss-Seidel Multi-Grid - solve the system Ax = b
    !    When going to lower levels, we compute the residual r = b - Ax and solve on lower levels Ae = r
    !    This is repeated until we reach the lowest level, then we go back upwards, completing one V-cycle.
    ! ------------------------------------------------------------
    !    for now I assume equidistant grid, so I loop level-wise
    !    in order to avoid leaf-wise refinement status hacking for now
    ! ------------------------------------------------------------
    Jmin = 0
    it_coarsest = -1  ! JMin=0 : if set to -1, FFT is used, elsewise CG, for JMin>0 GS sweeps are always used
    tol_cg = 1e-4

    call init_t_file('multigrid-cycle.t', .true., (/'    residual L2', '    residual L1', 'residual Linfty', '           time'/))
    call init_t_file('multigrid-iteration.t', .true., (/'direction', 'iteration', '     time'/))
    if (exist_u) then
        if (exist_uFD) then
            call init_t_file('multigrid-compare.t', .true., (/'             iteration', '    diff 2 spectral L2', '    diff 2 spectral L1', 'diff 2 spectral Linfty', '          diff 2 FD L2', '          diff 2 FD L1', '      diff 2 FD Linfty'/))
        else
            call init_t_file('multigrid-compare.t', .true., (/'             iteration', '    diff 2 spectral L2', '    diff 2 spectral L1', 'diff 2 spectral Linfty'/))
        endif
    endif
    ! first entries are always ignored when reading in files, as this is t=0?
    call append_t_file('multigrid-cycle.t', (/0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk/))
    call append_t_file('multigrid-iteration.t', (/0.0_rk, 0.0_rk, 0.0_rk/))
    if (exist_uFD .and. exist_u) call append_t_file('multigrid-compare.t', (/0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk/))
    if (exist_u .and. .not. exist_uFD) call append_t_file('multigrid-compare.t', (/0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk/))

    ! init solution as zero
    do k_block = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        hvy_tmp(:,:,:,1:nc,hvy_id) = 0.0_rk
    enddo

    ! compute several v- or f-cycles
    do i_cycle = 1,it_cycle

        ! ! usually we want spectral resolution for the coarsest level, for debugging we can disable this
        ! fft_order_discretization = "spectral"

        if (trim(cycle_type) == "V") then
            call multigrid_vcycle(params, hvy_tmp(:,:,:,1:nc,:), hvy_block(:,:,:,1:nc,:), hvy_tmp(:,:,:,nc+1:size(hvy_tmp,4),:), tree_ID, a, dx2_FD, it_down, it_up, fft_order_discretization)
        elseif (trim(cycle_type) == "F") then
            call abort(250514, 'F-Cycle is discontinued')
        else
            call abort(250429, "ERROR: cycle type not recognized. Use 'V' or 'F'")
        endif

        ! laplacian is invariant to shifts of constant values
        ! our values are defined with zero mean for comparison
        ! as multigrid might accidently introduce a constant offset, we remove it
        call componentWiseNorm_tree(params, hvy_tmp(:,:,:,1:nc,:), tree_ID, "Mean", norm(1:nc))
        do ic = 1,nc
            do k_block = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k_block, tree_ID)
                call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
                hvy_tmp(:,:,:,ic,hvy_id) = hvy_tmp(:,:,:,ic,hvy_id) - norm(ic)
            enddo
        enddo
        if (params%rank == 0) write(*, '(A, es10.3, A)') "--- Mean value: ", norm(1), " ---"

        ! compare to spectral solutions with spectral or FD accuracy
        if (exist_u) then
            ! compute difference to spectral accuracy
            do k_block = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k_block, tree_ID)
                if (.not. block_is_leaf(params, hvy_id)) cycle
                hvy_work(:,:,:,1:nc,hvy_id,1) = hvy_tmp(:,:,:,1:nc,hvy_id) - hvy_block(:,:,:,nc+1:2*nc,hvy_id)
            enddo
            ! compute norms
            call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:,1), tree_ID, "L2", norm)
            if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Diff spectral L2: ", norm(1), " ---"
            call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:,1), tree_ID, "L1", norm(nc+1:2*nc))
            if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Diff spectral L1: ", norm(nc+1), " ---"
            call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:,1), tree_ID, "Linfty", norm(2*nc+1:3*nc))
            if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Diff spectral Linfty: ", norm(2*nc+1), " ---"

            if (exist_uFD) then
                ! compute difference to FD accuracy
                do k_block = 1, hvy_n(tree_ID)
                    hvy_id = hvy_active(k_block, tree_ID)
                    if (.not. block_is_leaf(params, hvy_id)) cycle
                    hvy_work(:,:,:,1:nc,hvy_id,1) = hvy_tmp(:,:,:,1:nc,hvy_id) - hvy_block(:,:,:,2*nc+1:3*nc,hvy_id)
                enddo
                ! compute norms
                call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:,1), tree_ID, "L2", norm(3*nc+1:4*nc))
                if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Diff spectral FD L2: ", norm(3*nc+1), " ---"
                call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:,1), tree_ID, "L1", norm(4*nc+1:5*nc))
                if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Diff spectral FD L1: ", norm(4*nc+1), " ---"
                call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:,1), tree_ID, "Linfty", norm(5*nc+1:6*nc))
                if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Diff spectral FD Linfty: ", norm(5*nc+1), " ---"

                call append_t_file('multigrid-compare.t', (/dble(i_cycle), norm(1), norm(2), norm(3), norm(4), norm(5), norm(6)/))
            else
                call append_t_file('multigrid-compare.t', (/dble(i_cycle), norm(1), norm(2), norm(3)/))
            endif
        endif

    enddo

    ! delete all non-leaf blocks with daughters as we for now do not have any use for them
    call prune_fulltree2leafs(params, tree_ID)
    
    ! save file under new name
    write(fname, '(A)') "u_00000000000.h5"
    call saveHDF5_tree(fname, 0.0_rk, it, 1, params, hvy_tmp, tree_ID )

    ! save file under new name
    write(fname, '(A)') "res_00000000000.h5"
    call saveHDF5_tree(fname, 0.0_rk, it, 1, params, hvy_tmp(:,:,:,nc+1:size(hvy_tmp,4),:), tree_ID )

    ! save file under new name
    write(fname, '(A)') "b-out_00000000000.h5"
    call saveHDF5_tree(fname, 0.0_rk, it, 1, params, hvy_block, tree_ID )

    t_block = MPI_Wtime()
    call fft_destroy(params)
    call toc( "fft destroy", 10100, MPI_Wtime()-t_block )

    call summarize_profiling( WABBIT_COMM )

end subroutine proto_GS_multigrid