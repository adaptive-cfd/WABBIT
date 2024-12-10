program main
    use mpi
    use module_wavelets
    use module_helpers
    use module_MPI
    use module_params           ! global parameters
    use module_timing
    use module_initialization   ! init data module
    use module_mesh             ! mesh manipulation subroutines
    use module_forestMetaData
    use module_time_step        ! time step module
    use module_unit_test        ! unit test module
    use module_bridge_interface ! bridge implementation of wabbit
    use module_t_files
    ! this module is the saving wrapper (e.g. save state vector or vorticity)
    ! it exists to disentangle module_forest and module_IO
    use module_saving

    implicit none

    integer(kind=ik)                    :: ierr           ! MPI error variable
    integer(kind=ik)                    :: rank           ! process rank
    integer(kind=ik)                    :: number_procs   ! number of processes
    real(kind=rk)                       :: t0, t1, t2     ! cpu time variables for running time calculation
    type (type_params)                  :: params         ! user defined parameter structure


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

    ! time loop variables
    real(kind=rk)                       :: time, output_time
    integer(kind=ik)                    :: iteration
    ! filename of *.ini file used to read parameters
    character(len=clong)                :: filename
    ! some variables, loop, debug information and more
    integer(kind=ik)                    :: k, Nblocks_rhs, Nblocks, it, lgt_n_tmp, mpicode, Jmin1, Jmax1, Jmin2, Jmax2
    character(len=2*clong)              :: write_statement  ! dynamically build the output lengths
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, t4, tstart, dt
    ! decide if data is saved or not
    logical                             :: it_is_time_to_save_data=.false., test_failed, keep_running=.true.
    logical                             :: overwrite

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


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
        write(*,'(20("─"), A21, 39("─"))') "   STARTING wabbit   "
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
    ! initialize wavelet
    call setup_wavelet(params)
    ! initializes the communicator for Wabbit and creates a bridge if needed
    call initialize_communicator(params)
    ! have the pysics module read their own parameters
    call init_physics_modules( params, filename, params%N_mask_components )
    ! allocate memory for heavy, light, work and neighbor data
    call allocate_forest(params, hvy_block, hvy_work, hvy_tmp, hvy_mask)

    ! reset the grid: all blocks are inactive and empty
    call reset_tree( params, .true., tree_ID=tree_ID_flow )

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    !---------------------------------------------------------------------------
    ! Unit tests
    !    - treecode test for neighbor finding
    !    - sync test for correct patch allocation
    !    - godes nodes sync test for correct prediction order
    !    - wavelet decomposition test for filter banks
    !    - refine-coarsen test for refinement and adaption processes
    !---------------------------------------------------------------------------
    if (any((/ params%test_treecode, params%test_ghost_nodes_synch, params%test_wavelet_decomposition/)) .and. params%rank==0) then
        write(*,'("     /¯¯¯\                       /¯¯¯\       /¯¯¯\       /¯¯¯\       /¯¯¯\  ")')
        write(*,'("    /     \     Unit Tests      /     \     /     \     /     \     /     \ ")')
        write(*,'("___/       \___________________/       \___/       \___/       \___/       \")')
        write(*, '("")')
    end if


    if (params%test_treecode) then
       call unit_test_treecode( params, hvy_block, hvy_work, hvy_tmp, tree_ID_flow, abort_on_fail=.true. )
    end if

    ! perform a convergence test on ghost node sync'ing
    if (params%test_ghost_nodes_synch) then
        call unit_test_Sync( params, hvy_block, hvy_work, hvy_tmp, tree_ID_flow, abort_on_fail=.true.)
        call unit_test_ghostSync( params, hvy_block, hvy_work, hvy_tmp, tree_ID_flow, abort_on_fail=.false. )
    endif

    ! check if the wavelet filter banks are okay.
    if (params%test_wavelet_decomposition) then
        call unit_test_waveletDecomposition( params, hvy_block, hvy_work, hvy_tmp, tree_ID_flow )
        call unit_test_refineCoarsen( params, hvy_block, hvy_work, hvy_tmp, tree_ID_flow )
        call unit_test_waveletDecomposition_invertibility( params, hvy_block, hvy_work, hvy_tmp, tree_ID_flow )
    endif

    ! timing
    call toc( "TOPLEVEL: unit tests", 17, MPI_wtime()-sub_t0 )
    sub_t0 = MPI_Wtime()

    !---------------------------------------------------------------------------
    ! Initial condition
    !---------------------------------------------------------------------------
    ! reset the grid: all blocks are inactive and empty (in case the unit tests leave us with some data)
    call reset_tree( params, .true., tree_ID=tree_ID_flow )
    ! On all blocks, set the initial condition (incl. synchronize ghosts)
    call setInitialCondition_tree( params, hvy_block, tree_ID_flow, params%adapt_inicond, time, iteration, hvy_mask, hvy_tmp, hvy_work=hvy_work)

    if ((.not. params%read_from_files .or. params%adapt_inicond).and.(time>=params%write_time_first)) then
        ! NOTE: new versions (>16/12/2017) call physics module routines call prepare_save_data. These
        ! routines create the fields to be stored in the work array hvy_work in the first 1:params%N_fields_saved
        ! slots. the state vector (hvy_block) is copied if desired.
        call save_data( iteration, time, params, hvy_block, hvy_tmp, hvy_mask, tree_ID_flow )
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
    call init_t_file('thresholding.t', overwrite)

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
    call toc( "TOPLEVEL: init_data", 15, MPI_wtime()-sub_t0 )
    
    !---------------------------------------------------------------------------
    ! main time loop
    !---------------------------------------------------------------------------
    if (rank==0) then
        write(*,*) "params%next_write_time=", params%next_write_time
        write(*,*) "params%next_stats_time=", params%next_stats_time
        write(*,*) ""
        write(*,'(10(" "), "╔", 48("═"), "╗")') 
        write(*,'(10(" "), A)') "║ On your marks, ready, starting main time loop! ║"
        write(*,'(10(" "), "╚", 48("═"), "╝")')
        write(*,*) ""
    endif

    keep_running = .true.
    do while ( time<params%time_max .and. iteration<params%nt .and. keep_running)
        t2 = MPI_wtime()

        !***********************************************************************
        ! MPI bridge (used e.g. for particle-fluid coupling)
        !***********************************************************************
        if (params%bridge_exists) then
            call abort(2408229, "those routines must be updated to work on the right TREEs")
            ! call send_lgt_data (params)
            ! call serve_data_request(hvy_block, hvy_tmp, params)
        endif

        ! check if we still have enough memory left: for very large simulations
        ! we cannot affort to have them fail without the possibility to resume them
        if (notEnoughMemoryToRefineEverywhere_tree(params, tree_ID_flow, time)) then
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
            ! Snych'ing becomes much more expensive one the grid is refined.
            call sync_ghosts_tree( params, hvy_block, tree_ID_flow )

            ! refine the mesh after refinement_indicator, usually "everywhere" or "significant"
            ! for "significant", the refinement flags from the last adapt_tree call are reused. This might not be given for the first iteration so we just skip this
            if (params%refinement_indicator == "significant" .and. iteration == 0) then
                call refine_tree( params, hvy_block, hvy_tmp, "everywhere", tree_ID=tree_ID_flow)
            else
                call refine_tree( params, hvy_block, hvy_tmp, params%refinement_indicator, tree_ID=tree_ID_flow )
            endif
        endif
        call toc( "TOPLEVEL: refinement", 10, MPI_wtime()-t4)
        Nblocks_rhs = lgt_n(tree_ID_flow)

        Jmin1 = minActiveLevel_tree(tree_ID_flow)
        Jmax1 = maxActiveLevel_tree(tree_ID_flow)

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
            call timeStep_tree( time, dt, iteration, params, hvy_block, hvy_work, hvy_mask, hvy_tmp, tree_ID_flow )
            call toc( "TOPLEVEL: time stepper", 11, MPI_wtime()-t4)

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
                    call sync_ghosts_tree( params, hvy_block, tree_ID_flow )

                    call filter_wrapper(time, params, hvy_block, hvy_tmp, hvy_mask, tree_ID_flow)
                end if
            end if
            call toc( "TOPLEVEL: filter", 12, MPI_wtime()-t4)

            !*******************************************************************
            ! statistics
            !*******************************************************************
            t4 = MPI_wtime()
            if ( (modulo(iteration, params%nsave_stats)==0).or.(abs(time - params%next_stats_time)<1e-12_rk) ) then
                ! we need to sync ghost nodes for some derived qtys, for sure
                call sync_ghosts_RHS_tree( params, hvy_block, tree_ID_flow )

                call statistics_wrapper(time, dt, params, hvy_block, hvy_tmp, hvy_mask, tree_ID_flow)
            endif
            call toc( "TOPLEVEL: statistics", 13, MPI_wtime()-t4)

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
            ! NOTE new versions (>16/12/2017) call physics module routines call prepare_save_data. These
            ! routines create the fields to be stored in the work array hvy_tmp in the first 1:params%N_fields_saved
            ! slots. the state vector (hvy_block) is copied if desired.
            call save_data( iteration, time, params, hvy_block, hvy_tmp, hvy_mask, tree_ID_flow )

            output_time = time
            params%next_write_time = params%next_write_time + params%write_time
        endif

        t2 = MPI_wtime() - t2
        ! output on screen
        if (rank==0) then
            write(*, '("RUN: it=",i7)', advance='no') iteration
            write(*, '(" time=",f16.9, " t_wall=",es10.3)', advance='no') time, t2
            write(*, '(" Nb=(",i6,"/",i6,")")', advance='no') Nblocks_rhs, Nblocks
            write(*, '(" J=(",i2,":",i2,"/",i2,":",i2, ")")', advance='no') Jmin1, Jmax1, minActiveLevel_tree(tree_ID_flow), maxActiveLevel_tree(tree_ID_flow)
            write(*, '(" dt=",es8.1," mem=",i3,"%")', advance='no') dt, nint((dble(Nblocks_rhs+lgt_n(tree_ID_mask))/dble(size(lgt_block,1)))*100.0_rk)
            write(*, '("")')  ! line break

             ! prior to 11/04/2019, this file was called timesteps_info.t but it was missing some important
             ! information, so I renamed it when adding those (since post-scripts would no longer be compatible
             ! it made sense to me to change the name)
             call append_t_file( 'performance.t', (/time, dble(iteration), t2, dble(Nblocks_rhs), dble(Nblocks), &
             dble(minActiveLevel_tree(tree_ID_flow)), &
             dble(maxActiveLevel_tree(tree_ID_flow)), &
             dble(params%number_procs), dble(size(lgt_block,1)) /) )
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
        ! To reduce IO, do this after 15 iterations in 2D but in every time step in 3D (timeSteps are longer in 3D)
        if ((modulo(iteration,15)==0 .and. params%dim==2).or.(params%dim==3)) then
        if ( runtime_control_stop() ) then
            if (rank==0) write(*,*) "WE RECVED THE STOP COMMAND: WRITE BACKUP; THEN BYEBYE"
            keep_running = .false.
        endif
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
        call sync_ghosts_RHS_tree( params, hvy_block, tree_ID_flow)

        call statistics_wrapper(time, dt, params, hvy_block, hvy_tmp, hvy_mask, tree_ID_flow)
    endif


    ! close and flush all existings *.t files
    call close_all_t_files()

    ! save end field to disk, only if this data is not saved already
    if ( abs(output_time-time) > 1e-10_rk ) then
        ! Note new versions (>16/12/2017) call physics module routines call prepare_save_data. These
        ! routines create the fields to be stored in the work array hvy_tmp in the first 1:params%N_fields_saved
        ! slots. the state vector (hvy_block) is copied if desired.
        call save_data( iteration, time, params, hvy_block, hvy_tmp, hvy_mask, tree_ID_flow )
    end if

    ! MPI Barrier before program ends
    call MPI_Barrier(WABBIT_COMM, ierr)

    ! make a summary of the program parts, which have been profiled using toc(...)
    ! and print it to stdout, safe total time just before to insert it as well
    call toc( "TOPLEVEL: TOTAL", 9, MPI_wtime()-t0)
    call summarize_profiling( WABBIT_COMM )

    call deallocate_forest(params, hvy_block, hvy_work, hvy_tmp)

    ! computing time output on screen
    call cpu_time(t1)
    if (rank==0) then
        write(*,'(80("─"))')
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
