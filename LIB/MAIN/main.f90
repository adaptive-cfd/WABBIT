!> \brief main program, init all data, start time loop, output on screen during program run
! ********************************************************************************************
!> \image html rhs.svg width=600
!> \image html rhs.eps

program main

    use mpi
    use module_helpers
    use module_MPI
    use module_params           ! global parameters
    use module_timing
    use module_initialization   ! init data module
    use module_mesh             ! mesh manipulation subroutines
    use module_IO               ! IO module
    use module_time_step        ! time step module
    use module_unit_test        ! unit test module
    use module_bridge_interface ! bridge implementation of wabbit
    use module_mask
    ! this module is the saving wrapper (e.g. save state vector or vorticity)
    ! it exists to disentangle module_forest and module_IO
    use module_saving

    implicit none

    integer(kind=ik)                    :: ierr           ! MPI error variable
    integer(kind=ik)                    :: rank           ! process rank
    integer(kind=ik)                    :: number_procs   ! number of processes
    real(kind=rk)                       :: t0, t1, t2     ! cpu time variables for running time calculation
    type (type_params)                  :: params         ! user defined parameter structure

    ! light data array (the grid metadata, block treecodes etc)
    ! line number = ( 1 + proc_rank ) * heavy_data_line_number
    ! column(1:max_treelevel): block treecode, treecode -1 => block is inactive
    ! column(max_treelevel + IDX_MESH_LVL): treecode length = mesh level
    ! column(max_treelevel + IDX_REFINE_STS):   refinement status (-1..coarsen / 0...no change / +1...refine)
    integer(kind=ik), allocatable       :: lgt_block(:, :)

    ! heavy data array (the actual state vector). Is synchronized.
    ! dim 1-3: x,y,z coord ( 1:Bs+2*g )
    ! dim 4: components ( 1:number_equations)
    ! dim 5: block id  ( 1:number_blocks )
    real(kind=rk), allocatable          :: hvy_block(:, :, :, :, :)
    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    real(kind=rk), allocatable          :: hvy_mask(:, :, :, :, :)
    !!!!!! => renaming: hvy_block -> hvy_state

    ! hvy work array: the slots for RHS evaluation (e.g. 5 for a RK4)
    ! no synchronization performed
    ! dim 1-3: x,y,z coord ( 1:Bs+2*g )
    ! dim 4: components ( 1:number_equations)
    ! dim 5: RHS slot (k1,k2 etc for RK4)
    ! dim 6: block id  ( 1:number_blocks )
    real(kind=rk), allocatable          :: hvy_work(:, :, :, :, :, :)
    !!!!!! => renaming: hvy_work -> hvy_rhs

    ! work array (e.g. for saving fields, state vector conversions etc)
    ! no synchronization performed
    ! dim 1-3: x,y,z coord ( 1:Bs+2*g )
    ! dim 4: components ( 1:number_equations)
    ! dim 5: block id  ( 1:number_blocks )
    ! This array can be used for work data.
    real(kind=rk), allocatable          :: hvy_tmp(:, :, :, :, :)
    !!!!!! => renaming: hvy_tmp -> hvy_work

    ! neighbor array (heavy data) -> number_lines   = number_blocks (correspond to heavy data id)
    !                             -> number_columns = 16 (...different neighbor relations:
    ! '__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW' )
    !         saved data -> -1 ... no neighbor
    !                    -> light data id in corresponding column
    integer(kind=ik), allocatable       :: hvy_neighbor(:,:)

    ! The following list contains the numerical treecode and the lightID for the active blocks
    ! in a sorted fashion. this is very important for finding blocks. usually, in the rest of the code,
    ! a treecode is an array and this is handy. for finding a block however, this is not true,
    ! here, having a single, unique number is a lot faster. these numbers are called numerical treecodes.
    integer(kind=tsize), allocatable    :: lgt_sortednumlist(:, :, :)
    ! list of active blocks (light data) for each tree
    integer(kind=ik), allocatable       :: lgt_active(:, :)
    ! number of active blocks (light data) for each tree
    integer(kind=ik), allocatable       :: lgt_n(:)
    ! list of active blocks (heavy data) for each tree
    integer(kind=ik), allocatable       :: hvy_active(:, :)
    ! number of active blocks (heavy data) for each tree
    integer(kind=ik), allocatable       :: hvy_n(:)
    ! time loop variables
    real(kind=rk)                       :: time, output_time
    integer(kind=ik)                    :: iteration
    ! filename of *.ini file used to read parameters
    character(len=cshort)                   :: filename
    integer(kind=ik)                    :: k, Nblocks_rhs, Nblocks, it, tree_N, lgt_n_tmp, mpicode
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, t4, tstart, dt
    ! decide if data is saved or not
    logical                             :: it_is_time_to_save_data=.false., test_failed, keep_running=.true.
    logical                             :: overwrite

    ! init time loop
    time          = 0.0_rk
    output_time   = 0.0_rk
    iteration     = 0

    ! init mpi
    call MPI_Init(ierr)
    ! determine process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    params%rank = rank
    ! determine process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)
    params%number_procs=number_procs
    ! output MPI status
    WABBIT_COMM=MPI_COMM_WORLD

    if (rank==0) then
        write(*,'(80("_"))')
        write(*, '("MPI: using ", i5, " processes")') params%number_procs
        write(*,'("MPI: code build with NON-blocking send/recv in transfer (block_xfer_nonblocking.f90)")')
    end if

    call print_command_line_arguments()

    ! start time
    sub_t0 = MPI_Wtime()
    call cpu_time(t0)
    tstart = MPI_wtime()

    !---------------------------------------------------------------------------
    ! Initialize parameters,bridge and grid
    !---------------------------------------------------------------------------
    ! read in the parameter file to setup the case
    ! get the second command line argument: this should be the ini-file name
    call get_command_argument( 1, filename )
    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, filename )
    ! initializes the communicator for Wabbit and creates a bridge if needed
    call initialize_communicator(params)
    ! have the pysics module read their own parameters
    call init_physics_modules( params, filename, params%N_mask_components )
    ! allocate memory for heavy, light, work and neighbor data
    call allocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
        hvy_active, lgt_sortednumlist, hvy_work, hvy_tmp, hvy_mask, hvy_n, lgt_n)

    ! reset the grid: all blocks are inactive and empty
    call reset_tree( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, &
    lgt_sortednumlist, .true., tree_ID=tree_ID_flow )

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    !---------------------------------------------------------------------------
    ! Unit tests
    !---------------------------------------------------------------------------
    if (params%test_treecode) then
       call unit_test_treecode( params )
    end if

    ! perform a convergence test on ghost node sync'ing
    if (params%test_ghost_nodes_synch) then
        call unit_test_ghost_nodes_synchronization( params, lgt_block, hvy_block, hvy_work, &
        hvy_tmp, hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist, hvy_n, lgt_n )
    endif

    call reset_tree( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, &
    lgt_sortednumlist, .true., tree_ID=tree_ID_flow )


    !---------------------------------------------------------------------------
    ! Initial condition
    !---------------------------------------------------------------------------
    ! On all blocks, set the initial condition (incl. synchronize ghosts)
    call set_initial_grid( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_n, hvy_n, lgt_sortednumlist, tree_ID_flow, params%adapt_inicond, time, iteration, hvy_mask, hvy_tmp )

    if ((.not. params%read_from_files .or. params%adapt_inicond).and.(time>=params%write_time_first)) then
        ! save initial condition to disk (unless we're reading from file and do not adapt,
        ! in which case this makes no sense)
        ! we need to sync ghost nodes in order to compute the vorticity, if it is used and stored.
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

        ! NOte new versions (>16/12/2017) call physics module routines call prepare_save_data. These
        ! routines create the fields to be stored in the work array hvy_work in the first 1:params%N_fields_saved
        ! slots. the state vector (hvy_block) is copied if desired.
        call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active, &
        lgt_n, lgt_sortednumlist, hvy_n, hvy_tmp, hvy_active, hvy_mask, hvy_neighbor, tree_ID_flow )
    end if

    ! decide whether the ascii log files are overwritten and re-initialized (on startup
    ! this is what we want) or just openend for writing (when resuming a backup, of course we
    ! do not want old data to be erased)
    overwrite = .false.
    if (rank==0 .and. iteration==0) then
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

    if (rank==0) then
        call Initialize_runtime_control_file()
    endif

    ! next write time for reloaded data
    if (params%write_method == 'fixed_time') then
        params%next_write_time = real(floor(time/params%write_time), kind=rk)*params%write_time + params%write_time
        params%next_stats_time = real(floor(time/params%tsave_stats), kind=rk)*params%tsave_stats + params%tsave_stats

        ! sometimes, rarely, floor can be tricky: say we resume a run at t=2.26 and write_time is 0.01. if the h5 file
        ! is exactly at 2.26 but maybe the last digit flips (machine precision) it happens very rarely that
        ! floor(2.26/0.01) = 225 (and not 226). Then, WABBIT misses the next write time (and produces no output anymore...)
        ! correct for this mistake:
        if (abs(params%next_write_time-time)<=1.0e-10_rk) params%next_write_time = params%next_write_time + params%write_time
        if (abs(params%next_stats_time-time)<=1.0e-10_rk) params%next_stats_time = params%next_stats_time + params%tsave_stats
    end if


    ! timing
    call toc( "init_data", MPI_wtime()-sub_t0 )

    !---------------------------------------------------------------------------
    ! main time loop
    !---------------------------------------------------------------------------
    if (rank==0) then
        write(*,*) "params%next_write_time=", params%next_write_time
        write(*,*) "params%next_stats_time=", params%next_stats_time
        write(*,*) ""
        write(*,*) "        --------------------------------------------------"
        write(*,*) "        | On your marks, ready, starting main time loop! |"
        write(*,*) "        --------------------------------------------------"
        write(*,*) ""
    endif

    keep_running = .true.
    do while ( time<params%time_max .and. iteration<params%nt .and. keep_running)
        t2 = MPI_wtime()

        !***********************************************************************
        ! check redundant nodes
        !***********************************************************************
        t4 = MPI_wtime()
        if (params%check_redundant_nodes) then
            ! run the internal test for the ghost nodes.
            call check_unique_origin(params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), &
            hvy_n(tree_ID_flow), test_failed)

            if (test_failed) then
                call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active, &
                lgt_n, lgt_sortednumlist, hvy_n, hvy_tmp, hvy_active, hvy_mask, hvy_neighbor, tree_ID_flow )
                call abort(111111,"Same origin of ghost nodes check failed - stopping.")
            endif
        endif
        call toc( "TOPLEVEL: check ghost nodes", MPI_wtime()-t4)


        !***********************************************************************
        ! MPI bridge (used e.g. for particle-fluid coupling)
        !***********************************************************************
        if (params%bridge_exists) then
            call send_lgt_data (lgt_block,lgt_active(:,tree_ID_flow),lgt_n(tree_ID_flow),params)
            call serve_data_request(lgt_block, hvy_block, hvy_tmp, &
            hvy_neighbor, hvy_active(:,tree_ID_flow), lgt_active(:,tree_ID_flow), &
            lgt_n(tree_ID_flow), hvy_n(tree_ID_flow), params)
        endif

        ! check if we still have enough memory left: for very large simulations
        ! we cannot affort to have them fail without the possibility to resume them
        if (not_enough_memory(params, lgt_block, lgt_active, lgt_n)) then
            ! yippieh, a goto statement. thats soooo 90s.
            goto 17
        endif

        !***********************************************************************
        ! refine everywhere
        !***********************************************************************
        t4 = MPI_wtime()
        if ( params%adapt_tree ) then
            ! synchronization before refinement (because the interpolation takes place on the extended blocks
            ! including the ghost nodes)
            ! Note: at this point the grid is rather coarse (fewer blocks), and the sync step is rather cheap.
            ! Snych'ing becomes much mor expensive one the grid is refined.
            call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

            ! refine the mesh. Note: afterwards, it can happen that two blocks on the same level differ
            ! in their redundant nodes, but the ghost node sync'ing later on will correct these mistakes.
            call refine_tree( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
            lgt_sortednumlist, hvy_active, hvy_n, "everywhere", tree_ID=tree_ID_flow )
        endif
        call toc( "TOPLEVEL: refinement", MPI_wtime()-t4)
        Nblocks_rhs = lgt_n(tree_ID_flow)


        !***********************************************************************
        ! evolve solution in time
        !***********************************************************************
        ! internal loop over time steps: if desired, we perform more than one time step
        ! before adapting the grid again. this can further reduce the overhead of adaptivity.
        ! Note: the non-linear terms can create finer scales than resolved on the grid. they
        ! are usually filtered by the coarsening/refinement round trip. So if you do more than one time step
        ! on the grid, consider using a filter.
        do it = 1, params%N_dt_per_grid
            !*******************************************************************
            ! advance in time (make one time step)
            !*******************************************************************
            t4 = MPI_wtime()
            call timeStep_tree( time, dt, iteration, params, lgt_block, hvy_block, hvy_work, hvy_mask, hvy_tmp, hvy_neighbor, &
            hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist, tree_ID_flow )
            call toc( "TOPLEVEL: time stepper", MPI_wtime()-t4)

            ! determine if it is time to save data
            it_is_time_to_save_data = .false.
            if ((params%write_method=='fixed_freq' .and. modulo(iteration, params%write_freq)==0) .or. &
                (params%write_method=='fixed_time' .and. abs(time - params%next_write_time)<1.0e-12_rk)) then
                it_is_time_to_save_data = .true.
            endif
            ! do not save any output before this time (so maybe revoke the previous decision)
            if (time<=params%write_time_first) then
                it_is_time_to_save_data = .false.
            endif
            ! save after walltime unit is not affected by write_time_first
            if ((MPI_wtime()-tstart) - params%walltime_last_write > params%walltime_write*3600.0_rk) then
                params%walltime_last_write = MPI_wtime()-tstart
                it_is_time_to_save_data = .true.
            endif
            ! it can rarely happen that not all proc arrive at the same time at the above condition, then some decide to
            ! save data and others do not. this is a rare but severe problem, to solve it, synchronize:
            call MPI_BCAST( it_is_time_to_save_data, 1, MPI_LOGICAL, 0, WABBIT_COMM, mpicode )

            !*******************************************************************
            ! filter
            !*******************************************************************
            t4 = MPI_wtime()
            if (params%filter_type /= "no_filter") then
                if (modulo(iteration, params%filter_freq) == 0 .and. params%filter_freq > 0 .or. it_is_time_to_save_data) then
                    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

                    call filter_wrapper(time, params, hvy_block, hvy_tmp, hvy_mask, lgt_block, hvy_active(:,tree_ID_flow), &
                    hvy_n(tree_ID_flow), hvy_neighbor)
                end if
            end if
            call toc( "TOPLEVEL: filter", MPI_wtime()-t4)

            !*******************************************************************
            ! statistics
            !*******************************************************************
            if ( (modulo(iteration, params%nsave_stats)==0).or.(abs(time - params%next_stats_time)<1e-12_rk) ) then
                ! we need to sync ghost nodes for some derived qtys, for sure
                call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow))

                call statistics_wrapper(time, dt, params, hvy_block, hvy_tmp, hvy_mask, lgt_block, &
                lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, hvy_neighbor)
            endif

            ! if multiple time steps are performed on the same grid, we have to be careful
            ! not to skip past saving time intervals. Therefore, if it is time to save, we
            ! interupt the inner loop prematurely.
            if (it_is_time_to_save_data) exit
        enddo

        !***********************************************************************
        ! Adapt mesh (coarsening where possible)
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
                call create_mask_tree(params, time, lgt_block, hvy_mask, hvy_tmp, &
                hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist )

                ! actual coarsening (including the mask function)
                call adapt_tree( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
                lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_ID_flow, params%coarsening_indicator, hvy_tmp, hvy_mask )
            else
                ! actual coarsening (no mask function is required)
                call adapt_tree( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
                lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_ID_flow, params%coarsening_indicator, hvy_tmp )
            endif
        endif
        call toc( "TOPLEVEL: adapt mesh", MPI_wtime()-t4)
        Nblocks = lgt_n(tree_ID_flow)

        !***********************************************************************
        ! Write fields to HDF5 file
        !***********************************************************************
        if (it_is_time_to_save_data) then
            ! we need to sync ghost nodes in order to compute the vorticity, if it is used and stored.
            call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow))

            ! NOTE new versions (>16/12/2017) call physics module routines call prepare_save_data. These
            ! routines create the fields to be stored in the work array hvy_tmp in the first 1:params%N_fields_saved
            ! slots. the state vector (hvy_block) is copied if desired.
            call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active, &
            lgt_n, lgt_sortednumlist, hvy_n, hvy_tmp, hvy_active, hvy_mask, hvy_neighbor, tree_ID_flow )

            output_time = time
            params%next_write_time = params%next_write_time + params%write_time
        endif

        t2 = MPI_wtime() - t2
        ! output on screen
        if (rank==0) then
            write(*, '("RUN: it=",i7,1x," time=",f16.9,1x,"t_cpu=",es12.4," Nb=(",i6,"/",i6,") Jmin=",i2," Jmax=",i2, " dt=",es8.1)') &
             iteration, time, t2, Nblocks_rhs, Nblocks, &
             minActiveLevel_tree( lgt_block, tree_ID_flow, lgt_active, lgt_n ), &
             maxActiveLevel_tree( lgt_block, tree_ID_flow, lgt_active, lgt_n ), dt

             ! prior to 11/04/2019, this file was called timesteps_info.t but it was missing some important
             ! information, so I renamed it when adding those (since post-scripts would no longer be compatible
             ! it made sense to me to change the name)
             call append_t_file( 'performance.t', (/time, dble(iteration), t2, dble(Nblocks_rhs), dble(Nblocks), &
             dble(minActiveLevel_tree( lgt_block, tree_ID_flow, lgt_active, lgt_n )), &
             dble(maxActiveLevel_tree( lgt_block, tree_ID_flow, lgt_active, lgt_n )), &
             dble(params%number_procs) /) )
        end if

        !***********************************************************************
        ! walltime limiter
        !***********************************************************************
        ! maximum walltime allowed for simulations (in hours). The run will be stopped if this duration
        ! is exceeded. This is useful on real clusters, where the walltime of a job is limited, and the
        ! system kills the job regardless of whether we're done or not. If WABBIT itself ends execution,
        ! a backup is written and you can resume the simulation right where it stopped
        if ( (rank==0) .and. (MPI_wtime()-tstart)/3600.0_rk >= params%walltime_max ) then
            if (rank==0) write(*,'("WE ARE OUT OF WALLTIME: STOP. ",g12.3,"h / ",g12.3,"h")') (MPI_wtime()-tstart)/3600.0_rk, params%walltime_max
            keep_running = .false.
        endif
        ! it can rarely happen that not all proc arrive at the same time at the above condition, then some decide to
        ! stop and others not. this is a rare but severe problem, to solve it, synchronize:
        call MPI_BCAST( keep_running, 1, MPI_LOGICAL, 0, WABBIT_COMM, mpicode )

        !***********************************************************************
        ! runtime control
        !***********************************************************************
        ! it happens quite often that one wants to end a simulation prematurely, but
        ! one also wants to be able to resume it. the usual "kill" on clusters will terminate
        ! wabbit immediately and not write a backup. Hence, it is possible to terminate
        ! wabbit by writing "save_stop" to "runtime_control"
        if ( runtime_control_stop() ) then
            if (rank==0) write(*,*) "WE RECVED THE STOP COMMAND: WRITE BACKUP; THEN BYEBYE"
            keep_running = .false.
        endif
    end do

    !***************************************************************************
    ! end of main time loop
    !***************************************************************************
17  if (rank==0) write(*,*) "This is the end of the main time loop!"


    !*******************************************************************
    ! statistics ( last time )
    !*******************************************************************
    if ( (modulo(iteration, params%nsave_stats)==0).or.(abs(time - params%next_stats_time)<1e-12_rk) ) then
        ! we need to sync ghost nodes for some derived qtys, for sure
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow))

        call statistics_wrapper(time, dt, params, hvy_block, hvy_tmp, hvy_mask, lgt_block, &
        lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, hvy_neighbor)
    endif


    ! close and flush all existings *.t files
    call close_all_t_files()

    ! save end field to disk, only if this data is not saved already
    if ( abs(output_time-time) > 1e-10_rk ) then
        ! we need to sync ghost nodes in order to compute the vorticity, if it is used and stored.
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow))

        ! filter before write out
        if ( params%filter_freq > 0 .and. params%filter_type/="no_filter") then
            call filter_wrapper(time, params, hvy_block, hvy_tmp, hvy_mask, lgt_block, hvy_active(:,tree_ID_flow), &
            hvy_n(tree_ID_flow), hvy_neighbor)
        end if

        ! Note new versions (>16/12/2017) call physics module routines call prepare_save_data. These
        ! routines create the fields to be stored in the work array hvy_tmp in the first 1:params%N_fields_saved
        ! slots. the state vector (hvy_block) is copied if desired.
        call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active, &
        lgt_n, lgt_sortednumlist, hvy_n, hvy_tmp, hvy_active, hvy_mask, hvy_neighbor, tree_ID_flow )
    end if

    ! MPI Barrier before program ends
    call MPI_Barrier(WABBIT_COMM, ierr)

    ! make a summary of the program parts, which have been profiled using toc(...)
    ! and print it to stdout
    call summarize_profiling( WABBIT_COMM )

    call deallocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active,&
        hvy_active, lgt_sortednumlist, hvy_work, hvy_tmp, hvy_n, lgt_n )

    ! computing time output on screen
    call cpu_time(t1)
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("END: cpu-time = ",f16.4, " s")')  t1-t0
    end if

    ! on HPC clusters, one is often confronted with a limited runtime (walltime limit)
    ! and has to resume expensive simulations quite often. this process is usually automatized
    ! but for this it is nice to have a quick-to-evaluate criterion to know if a run ended normally
    ! or with an error. So write an empty success file if the run ended normally
    if (rank==0 .and. .not. params%out_of_memory) then
        open (77, file='success', status='replace')
        close(77)
    endif

    ! end mpi
    call MPI_Finalize(ierr)

end program main
