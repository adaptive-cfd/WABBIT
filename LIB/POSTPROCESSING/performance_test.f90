subroutine performance_test(params)
    use mpi
    use module_helpers
    use module_MPI
    use module_params           ! global parameters
    use module_timing
    use module_mesh             ! mesh manipulation subroutines
    use module_time_step
    use module_unit_test
    use module_bridge_interface ! bridge implementation of wabbit
    use module_forestMetaData

    implicit none
    type (type_params), intent(inout)   :: params

    ! perform 20 time steps per mesh.
    integer(kind=ik)                    :: N_timesteps, N_grids, equi_level, Bs
    real(kind=rk)                       :: target_grid_density, grid_exponent

    integer(kind=ik)                    :: number_procs, ierr, rank
    real(kind=rk), allocatable          :: t0_timesteps(:)

    real(kind=rk), allocatable          :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable          :: hvy_mask(:, :, :, :, :)
    real(kind=rk), allocatable          :: hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable          :: hvy_tmp(:, :, :, :, :)

    real(kind=rk)                       :: time = 0.0_rk
    integer(kind=ik)                    :: iteration = 0
    character(len=cshort)               :: filename
    integer(kind=ik)                    :: k, Nblocks_rhs, Nblocks, it, lgt_n_tmp, j, a
    real(kind=rk)                       :: t0, dt, t4
    logical :: error_OOM, debug_setting

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    call disable_all_t_files_output()

    !---------------------------------------------------------------------------
    ! Initialize parameters,bridge and grid
    !---------------------------------------------------------------------------
    ! read in the parameter file to setup the case
    ! get the second command line argument: this should be the ini-file name
    call get_command_argument( 2, filename )

    ! does the user need help?
    if (filename=='--help' .or. filename=='--h' .or. filename=='-h') then
        if (params%rank==0) then
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " Wabbit performance test"
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " This tool performs performance benchmarks by creating grids"
            write(*,'(A)') " and executing multiple time steps with refinement and"
            write(*,'(A)') " coarsening operations to measure computational efficiency."
            write(*,'(A)') ""
            write(*,'(A)') " ./wabbit-post --performance-test PARAMS.ini [options]"
            write(*,'(A)') ""
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " Optional arguments:"
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " --N_timesteps=15           Number of time steps per grid (default: 15)"
            write(*,'(A)') " --N_grids=50               Number of different grids to test (default: 50)"
            write(*,'(A)') " --target_grid_density=0.11 Target grid density (default: 0.11)"
            write(*,'(A)') " --equi=J                   Test only equidistant grids from level 0 to J, this disables the target_grid_density and N_grids options (default: -1, meaning disabled)"
            write(*,'(A)') " --Bs=32                    Block size (overrides ini file)"
            write(*,'(A)') " --wavelet=CDF40            Wavelet type (overrides ini file)"
            write(*,'(A)') " --g=4                      Number of ghost nodes (overrides ini file)"
            write(*,'(A)') " --debug                    Enable all debug output"
            write(*,'(A)') "-----------------------------------------------------------"
        end if
        return
    endif

    ! Get optional command line arguments
    call get_cmd_arg( "--N_timesteps", N_timesteps, default=15 )
    call get_cmd_arg( "--N_grids", N_grids, default=50 )
    ! this default should be 1 / 2**dim * (2**dim-1) / 2**dim * 0.95
    call get_cmd_arg( "--target_grid_density", target_grid_density, default=0.11_rk )
    call get_cmd_arg( "--equi", equi_level, default=-1 )
    call get_cmd_arg( "--debug", debug_setting, default=.false.)

    ! If equi mode is set, override N_grids
    if (equi_level > 0) then
        N_grids = equi_level + 1
    endif

    ! Allocate timestep array based on N_timesteps
    allocate(t0_timesteps(1:N_timesteps))
    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, filename )

    if (equi_level > params%Jmax) then
        params%Jmax = equi_level
        write(*,'("Warning: equi level (",i0,") is higher than Jmax (",i0,"). Setting Jmax = equi level.")') equi_level, params%Jmax
    endif
    
    ! Override parameters if specified on command line
    call get_cmd_arg( "--Bs", Bs, default=-1 )
    call get_cmd_arg( "--wavelet", params%wavelet, default=params%wavelet )
    call get_cmd_arg( "--g", params%g, default=params%g )
    if (Bs > 0) then
        params%Bs(1:params%dim) = Bs
    endif
    call setup_wavelet(params, params%g)

    ! write bs, wavelet and g to user just to be save
    if (params%rank==0) then
        write(*,'("Using wavelet = ", A, ", g = ", i2, ", Bs = ", 3(i3, 1x))') params%wavelet, params%g, params%Bs(1:params%dim)
    end if
    
    ! initializes the communicator for Wabbit and creates a bridge if needed
    call initialize_communicator(params)
    ! have the pysics module read their own parameters
    call init_physics_modules( params, filename, params%N_mask_components )
    ! allocate memory for heavy, light, work and neighbor data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_work=hvy_work, hvy_mask=hvy_mask)

    ! reset the grid: all blocks are inactive and empty
    call reset_tree( params, .true., tree_ID=tree_ID_flow )

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    ! Test is as follows:
    ! 1 Create a random grid with given density
    ! 2 perform N_timesteps times: (for statistical results)
    !   2.1 Refine everywhere
    !   2.2 Evolve in time
    !   2.3 Coarsen everhwyere
    !

    if (equi_level > 0) then
        if (params%rank == 0) write(*,'(A, i0, A)') "Starting performance test in equi mode with levels 0 to ", equi_level, "."
    else
        ! compute exponent, so that for the first grid we have 1 block per proc
        ! take into consideration, that this is then still refined, so we divide it by 2^d
        grid_exponent = log(target_grid_density * dble(params%number_blocks)*(2.0_rk**params%dim)) / log(dble(N_grids))
    endif

    do a = 1, N_grids

        ! in random mode - set grid density
        if (equi_level <= 0) then
            params%max_grid_density = target_grid_density * (dble(a) / dble(N_grids))**grid_exponent
        endif

        call reset_all_timings()

        if (equi_level > 0) then
            ! In equi mode, create equidistant grid at specific level
            call createEquidistantGrid_tree( params, hvy_block, a-1, verbosity=.true., tree_ID=tree_ID_flow )
        else
            ! Random grid mode
            call createRandomGrid_tree( params, hvy_block, hvy_tmp, level_init=3, verbosity=.true., iterations=10, tree_ID=tree_ID_flow )

            ! check for debug_setting if every block has minimum 1 block, because otherwise gathering bugs out
            ! if that happens, just skip this grid
            if (debug_setting) then
                call MPI_Allreduce(hvy_n(tree_ID_flow), j, 1, MPI_INTEGER4, MPI_MIN, WABBIT_COMM, ierr)
                if (j == 0) then
                    if (params%rank==0) write(*,'(A, i0, A)') "Skipping grid ", a, " because at least one process has 0 blocks, which causes gathering to fail in debug mode."
                    cycle
                endif
            endif
        endif

        ! on the grid, set some random data
        do k = 1, hvy_n(tree_ID_flow)
            call random_data( hvy_block(:, :, :, :, hvy_active(k, tree_ID_flow) ) )
        enddo

        ! let's debug everything we do in case the user requests it.
        if (debug_setting) then
            params%debug_balanceLoad = .true.
            params%debug_poisson = .true.
            params%debug_pruned2full = .true.
            params%debug_refinement = .true.
            params%debug_sync = .true.
            params%debug_wavelet_decompose = .true.
            params%debug_wavelet_reconstruct = .true.
        endif

        do j = 1, N_timesteps
            
            t0 = MPI_wtime()
            ! refine everywhere
            t4 = MPI_wtime()
            if ( params%adapt_tree ) then
                call sync_ghosts_tree( params, hvy_block, tree_ID_flow )

                call refine_tree( params, hvy_block, "everywhere", tree_ID=tree_ID_flow, error_OOM=error_OOM )
            endif
            call toc( "TOPLEVEL: refinement", 10, MPI_wtime()-t4)
            Nblocks_rhs = lgt_n(tree_ID_flow)


            ! internal loop over time steps: if desired, we perform more than one time step
            ! before adapting the grid again
            do it = 1, params%N_dt_per_grid
                t4 = MPI_wtime()
                call timeStep_tree( time, dt, iteration, params, hvy_block, hvy_work, hvy_mask, hvy_tmp, tree_ID=tree_ID_flow )
                call toc( "TOPLEVEL: time stepper", 11, MPI_wtime()-t4)
            enddo

            ! Adapt mesh (coarsening where possible)
            t4 = MPI_wtime()
            if ( params%adapt_tree ) then
                ! actual coarsening
                call adapt_tree( time, params, hvy_block, tree_ID_flow, "everywhere", hvy_tmp )
            endif
            call toc( "TOPLEVEL: adapt mesh", 14, MPI_wtime()-t4)
            Nblocks = lgt_n(tree_ID_flow)

            t0_timesteps(j) = MPI_wtime() - t0
            call toc( "TOPLEVEL: TOTAL", 9, t0_timesteps(j) )
            
            if (params%rank==0) write(*,'(A, i3, A, i3, A, f8.4, A)') "   Step ", j, "/", N_timesteps, ": ", t0_timesteps(j), " s"
        enddo

        ! let's disable debugging because I don't want it to be polluted by random grid generation
        if (debug_setting) then
            params%debug_balanceLoad = .false.
            params%debug_poisson = .false.
            params%debug_pruned2full = .false.
            params%debug_refinement = .false.
            params%debug_sync = .false.
            params%debug_wavelet_decompose = .false.
            params%debug_wavelet_reconstruct = .false.
        endif

        if (params%rank==0) then
            write(*,*) "max_grid_density", params%max_grid_density
            write(*,'("Result on this grid: Nb=(",i9,"/",i9,") tcpu=",es12.4," Ncpu=",i6," Nb_allocated=",i9,&
            &" set_grid_density=",f12.4," Jmin=",i2," Jmax=",i2)') &
            Nblocks_rhs, Nblocks, &
            dble(params%number_procs)*sum(t0_timesteps) / dble(N_timesteps) / (dble(Nblocks_rhs)*product(params%Bs(1:params%dim)-1)), &
            params%number_procs, size(lgt_block, 1), params%max_grid_density, &
            minActiveLevel_tree(tree_ID_flow), &
            maxActiveLevel_tree(tree_ID_flow)
        endif

        call summarize_profiling( WABBIT_COMM )

    enddo


    call deallocate_forest(params, hvy_block, hvy_work, hvy_tmp)
end subroutine
