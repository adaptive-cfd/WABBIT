!> \brief this is a prototype - it tries out the Gauss-Seidel Multigrid to solve the poisson equation
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
    integer(kind=ik)                   :: tree_ID=1, ic, nc, i_cycle, it, tc_length, mpierr, g(1:3)

    character(len=cshort)              :: fname, cycle_type
    logical                            :: exist_u, exist_uFD
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
    call get_command_argument(2, file_b)

    ! does the user need help?
    if (file_b=='--help' .or. file_b=='--h' .or. file_b=='-h') then
        if (params%rank==0) then
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Wabbit prototype: Multi-Grid Gauss-Seidel solver"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Read in a data field (2D or 3D) b and solves the system laplacian u = b"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') "./wabbit-post --proto-GS-multigrid b.h5 [u.h5] [u_FD.h5] WAVELET LAPLACE-ORDER IT_CYCLE IT_GS"
            write(*, '(A)') "    b.h5 = RHS field to be solved"
            write(*, '(A)') "    u.h5 = solution field for comparison (optional)"
            write(*, '(A)') "    u_FD.h5 = solution field with FD accuracy for comparison (optional)"
            write(*, '(A)') "    WAVELET = CDFXY"
            write(*, '(A)') "    LAPLACE-ORDER = 2, 4, 6, 8 for the order of the Laplace operator"
            write(*, '(A)') "    IT_CYCLE = Number of cycles to perform"
            write(*, '(A)') "    IT_GS = Number of GS-sweeps for upwards iterations"
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
    call get_command_argument(4 + merge(1,0, exist_u) + merge(1,0, exist_uFD), params%laplacian_order)
    ! set number of cycles
    call get_command_argument(5 + merge(1,0, exist_u) + merge(1,0, exist_uFD), cycle_type)
    read (cycle_type, *) params%laplacian_cycle_it
    call get_command_argument(6 + merge(1,0, exist_u) + merge(1,0, exist_uFD), cycle_type)
    read (cycle_type, *) params%laplacian_GS_it

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file_b, lgt_n(tree_ID), time, it, domain, Bs, tc_length, params%dim, &
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

    ! set only one variable for now
    params%n_eqn = 1 + merge(1,0, exist_u) + merge(1,0, exist_uFD)
    nc = 1

    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component(1) = "0"
    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component = 1
    params%order_discretization = "FD_4th_central"
    params%laplacian_coarsest = "FFT"
    params%FFT_accuracy = "spectral"  ! FD or spectral

    Bs = params%Bs

    call setup_wavelet(params, params%g)
    call setup_laplacian_stencils(params, params%g)

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
            i_cycle = 3
        else
            call readHDF5vct_tree( (/file_b, file_u/), params, hvy_block, tree_ID)
            i_cycle = 2
        endif
    else
        call readHDF5vct_tree( (/file_b/), params, hvy_block, tree_ID)
        i_cycle = 1
    endif
    call componentWiseNorm_tree(params, hvy_block(:,:,:,1:i_cycle,:), tree_ID, "Mean", norm(1:i_cycle))
    if (params%rank == 0) write(*, '(A, 3(1x,es10.3))') "--- Mean values of read data: ", norm(1:i_cycle)
    do ic = 1,i_cycle
        do k_block = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k_block, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            hvy_block(:,:,:,ic,hvy_id) = hvy_block(:,:,:,ic,hvy_id) - norm(ic)
        enddo
    enddo

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

    ! For the Mehrstellenverfahren, we actually do solve the system Au = Bb
    ! So before doing anything, we apply the matrix B on b, this is a large tensorial matrix, but we only have to apply this operation once
    if (params%laplacian_order == "FD_6th_mehrstellen") then
        t_block = MPI_Wtime()
        do k_block = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k_block, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            call get_block_spacing_origin( params, lgt_ID, x0, dx )
            do ic = 1, nc
                ! apply the stencil to the RHS
                call GS_compute_Ax(params, hvy_block(:,:,:,ic,hvy_id), hvy_tmp(:,:,:,ic,hvy_id), dx, apply_B_RHS=.true.)
            enddo
            hvy_block(:,:,:,1:nc,hvy_id) = hvy_tmp(:,:,:,1:nc,hvy_id)
        enddo
        call toc( "RHS Preparation Bb", 10000, MPI_Wtime()-t_block )

        t_block = MPI_Wtime()
        call sync_ghosts_tree(params, hvy_block(:,:,:,1:nc,:), tree_ID)
        call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )
    endif

    ! init solution as zero
    do k_block = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        hvy_tmp(:,:,:,1:nc,hvy_id) = 0.0_rk
    enddo

    ! compute several v- or f-cycles
    do i_cycle = 1,params%laplacian_cycle_it

        call multigrid_vcycle(params, hvy_tmp(:,:,:,1:nc,:), hvy_block(:,:,:,1:nc,:), hvy_tmp(:,:,:,nc+1:size(hvy_tmp,4),:), tree_ID, verbose=.true.)

        ! laplacian is invariant to shifts of constant values
        ! our values are defined with zero mean for comparison
        ! as multigrid might accidently introduce a constant offset, we remove it
        call componentWiseNorm_tree(params, hvy_tmp(:,:,:,1:nc,:), tree_ID, "Mean", norm(1:nc))
        do k_block = 1, hvy_n(tree_ID)
            do ic = 1,nc
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

    ! for not saving the intermediate values, also change all save times to 0
    enddo

    ! delete all non-leaf blocks with daughters as we for now do not have any use for them
    call prune_fulltree2leafs(params, tree_ID)
    
    ! save file under new name
    write(fname, '(A, i4.4, A)') "u_00", 0, "00000.h5"
    call saveHDF5_tree(fname, dble(0), 0, 1, params, hvy_tmp, tree_ID )

    ! save file under new name
    write(fname, '(A, i4.4, A)') "res_00", 0, "00000.h5"
    call saveHDF5_tree(fname, dble(0), 0, 1, params, hvy_tmp(:,:,:,nc+1:size(hvy_tmp,4),:), tree_ID )

    ! save file under new name
    write(fname, '(A, I4.4, A)') "b-out_00", 0, "00000.h5"
    call saveHDF5_tree(fname, dble(0), 0, 1, params, hvy_block, tree_ID )

    ! ! for saving intermediate value, also change all save times to i_cycle
    !     ! save spectral u at new time position
    !     if (exist_u) then
    !         write(fname, '(A, I4.4, A)') "u-spectral_00", i_cycle, "00000.h5"
    !         call saveHDF5_tree(fname, dble(i_cycle), i_cycle, 1, params, hvy_block(:,:,:,nc+1:2*nc,:), tree_ID )
    !         if (exist_uFD) then
    !             write(fname, '(A, I4.4, A)') "u-FD_00", i_cycle, "00000.h5"
    !             call saveHDF5_tree(fname, dble(i_cycle), i_cycle, 1, params, hvy_block(:,:,:,2*nc+1:3*nc,:), tree_ID )
    !         endif
    !     endif
    !     call init_full_tree(params, tree_ID, set_ref=0)
    ! enddo

    t_block = MPI_Wtime()
    call fft_destroy(params)
    call toc( "fft destroy", 10100, MPI_Wtime()-t_block )

    call summarize_profiling( WABBIT_COMM )

end subroutine proto_GS_multigrid