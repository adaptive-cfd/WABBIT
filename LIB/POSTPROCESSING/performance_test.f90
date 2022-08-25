subroutine performance_test(params)
    use mpi
    use module_helpers
    use module_MPI
    use module_params           ! global parameters
    use module_timing
    use module_initialization   ! init data module
    use module_mesh             ! mesh manipulation subroutines
    use module_time_step
    use module_unit_test
    use module_bridge_interface ! bridge implementation of wabbit
    use module_forestMetaData

    implicit none
    type (type_params), intent(inout)   :: params

    ! perform 20 time steps per mesh.
    integer, parameter                  :: N_timesteps = 15
    integer, parameter                  :: N_grids = 50
    real(kind=rk), parameter            :: target_grid_density = 0.11

    integer(kind=ik)                    :: number_procs, ierr, rank
    real(kind=rk)                       :: t0_timesteps(1:N_timesteps)

    real(kind=rk), allocatable          :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable          :: hvy_mask(:, :, :, :, :)
    real(kind=rk), allocatable          :: hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable          :: hvy_tmp(:, :, :, :, :)

    real(kind=rk)                       :: time = 0.0_rk
    integer(kind=ik)                    :: iteration = 0
    character(len=cshort)                   :: filename
    integer(kind=ik)                    :: k, Nblocks_rhs, Nblocks, it, lgt_n_tmp, j, a
    real(kind=rk)                       :: t0, dt, t4

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
    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, filename )
    ! initializes the communicator for Wabbit and creates a bridge if needed
    call initialize_communicator(params)
    ! have the pysics module read their own parameters
    call init_physics_modules( params, filename, params%N_mask_components )
    ! allocate memory for heavy, light, work and neighbor data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_work=hvy_work)


    ! reset the grid: all blocks are inactive and empty
    call reset_tree( params, .true., tree_ID=tree_ID_flow )

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )


    params%max_grid_density = target_grid_density / real(N_grids)

    do a = 1, N_grids

        call createRandomGrid_tree( params, hvy_block, hvy_tmp, Jmin=1, verbosity=.true., iterations=10, tree_ID=tree_ID_flow )

        ! on the grid, set some random data
        do k = 1, hvy_n(tree_ID_flow)
            call random_data( hvy_block(:, :, :, :, hvy_active(k, tree_ID_flow) ) )
        enddo

        call reset_all_timings()

        do j = 1, N_timesteps
            t0 = MPI_wtime()
            ! refine everywhere
            t4 = MPI_wtime()
            if ( params%adapt_tree ) then
                call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

                call refine_tree( params, hvy_block, "everywhere", tree_ID=tree_ID_flow )
            endif
            call toc( "TOPLEVEL: refinement", MPI_wtime()-t4)
            Nblocks_rhs = lgt_n(tree_ID_flow)


            ! internal loop over time steps: if desired, we perform more than one time step
            ! before adapting the grid again
            do it = 1, params%N_dt_per_grid
                t4 = MPI_wtime()
                call timeStep_tree( time, dt, iteration, params, hvy_block, hvy_work, hvy_mask, hvy_tmp, tree_ID=tree_ID_flow )
                call toc( "TOPLEVEL: time stepper", MPI_wtime()-t4)
            enddo

            ! Adapt mesh (coarsening where possible)
            t4 = MPI_wtime()
            if ( params%adapt_tree ) then
                ! actual coarsening
                call adapt_tree( time, params, hvy_block, tree_ID_flow, "everywhere", hvy_tmp )
            endif
            call toc( "TOPLEVEL: adapt mesh", MPI_wtime()-t4)
            Nblocks = lgt_n(tree_ID_flow)

            t0_timesteps(j) = MPI_wtime() - t0
        enddo

        if (params%rank==0) then
            write(*,*) "max_grid_density", params%max_grid_density
            write(*,'("Result on this grid: Nb=(",i9,"/",i9,") tcpu=",es12.4," Ncpu=",i6," Nb_allocated=",i9,&
            &" set_grid_density=",f12.4," Jmin=",i2," Jmax=",i2)') &
            Nblocks_rhs, Nblocks, &
            dble(params%number_procs)*sum(t0_timesteps) / dble(N_timesteps) / (dble(Nblocks_rhs)*product(params%Bs-1)), &
            params%number_procs, size(lgt_block, 1), params%max_grid_density, &
            minActiveLevel_tree(tree_ID_flow), &
            maxActiveLevel_tree(tree_ID_flow)
        endif

        call summarize_profiling( WABBIT_COMM )

        ! next grid will be denser
        params%max_grid_density = params%max_grid_density + target_grid_density / real(N_grids)
    enddo


    call deallocate_forest(params, hvy_block, hvy_work, hvy_tmp)
end subroutine
