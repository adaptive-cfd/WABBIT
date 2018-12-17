!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name    main.f90
!> \version 0.5
!> \author  msr
!
!> \brief main program, init all data, start time loop, output on screen during program run
!
!>
!!
!! = log ======================================================================================
!! \n
!! 04/11/16 - switch to v0.4 \n
!! 23/11/16 - use computing time array for simple performance tests \n
!! 07/12/16 - now uses heavy work data array \n
!! 25/01/17 - switch to 3D, v0.5
! ********************************************************************************************
!> \image html rhs.svg width=600
!> \image html rhs.eps

program main

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    use module_helpers
    use module_MPI
    ! global parameters
    use module_params
    ! debug module
    use module_debug
    ! init data module
    use module_initialization
    ! mesh manipulation subroutines
    use module_mesh
    ! IO module
    use module_IO
    ! time step module
    use module_time_step
    ! unit test module
    use module_unit_test
    ! bridge implementation of wabbit
    use module_bridge_interface

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0, t1, t2

    ! user defined parameter structure
    type (type_params)                  :: params

    ! light data array  -> line number = ( 1 + proc_rank ) * heavy_data_line_number
    !                   -> column(1:max_treelevel): block treecode, treecode -1 => block is inactive
    !                   -> column(max_treelevel + idx_mesh_lvl): treecode length = mesh level
    !                   -> column(max_treelevel + idx_refine_sts):   refinement status (-1..coarsen / 0...no change / +1...refine)
    integer(kind=ik), allocatable       :: lgt_block(:, :)

    !                   -> dim 1: x coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 2: y coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 3: z coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 4: components ( 1:number_equations)
    ! heavy data array  -> dim 5: block id  ( 1:number_blocks )
    real(kind=rk), allocatable          :: hvy_block(:, :, :, :, :)

    !                   -> dim 1: x coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 2: y coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 3: z coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 4: components ( 1:number_equations)
    !                   -> dim 5: RHS slot (k1,k2 etc for RK4)
    ! heavy work array  -> dim 6: block id  ( 1:number_blocks )
    real(kind=rk), allocatable          :: hvy_work(:, :, :, :, :, :)

    !                   -> dim 1: x coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 2: y coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 3: z coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 4: components ( 1:number_equations)
    ! heavy data array  -> dim 5: block id  ( 1:number_blocks )
    ! This array can be used for work data.
    real(kind=rk), allocatable          :: hvy_tmp(:, :, :, :, :)

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
    integer(kind=tsize), allocatable    :: lgt_sortednumlist(:,:)

    ! list of active blocks (light data)
    integer(kind=ik), allocatable       :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik)                    :: lgt_n

    ! list of active blocks (heavy data)
    integer(kind=ik), allocatable       :: hvy_active(:)
    ! number of active blocks (heavy data)
    integer(kind=ik)                    :: hvy_n

    integer(kind=ik), allocatable       :: blocks_per_rank(:)

    ! time loop variables
    real(kind=rk)                       :: time, output_time
    integer(kind=ik)                    :: iteration

    ! filename of *.ini file used to read parameters
    character(len=80)                   :: filename

    ! loop variable
    integer(kind=ik)                    :: k, max_neighbors, Nblocks_rhs, Nblocks, it

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, t4, tstart, dt
    ! decide if data is saved or not
    logical                             :: it_is_time_to_save_data, test_failed, keep_running=.true.
!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! init time loop
    time          = 0.0_rk
    output_time   = 0.0_rk
    iteration     = 0

!---------------------------------------------------------------------------------------------
! main body

    ! init mpi
    call MPI_Init(ierr)
    ! determine process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    params%rank = rank
    ! determine process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)
    params%number_procs=number_procs
    allocate(blocks_per_rank(1:number_procs))
    ! output MPI status
    WABBIT_COMM=MPI_COMM_WORLD

    if (rank==0) then
        write(*,'(80("_"))')
        write(*, '("MPI: using ", i5, " processes")') params%number_procs
#ifdef BLOCKINGSENDRECV
        write(*,'("MPI: code build with blocking send/recv in transfer (block_xfer_blocking.f90)")')
#else
        write(*,'("MPI: code build with NON-blocking send/recv in transfer (block_xfer_nonblocking.f90)")')
#endif
    end if


    ! start time
    sub_t0 = MPI_Wtime()
    call cpu_time(t0)
    tstart = MPI_wtime()


    ! unit test off
    params%unit_test    = .false.

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
    call init_physics_modules( params, filename )
    ! allocate memory for heavy, light, work and neighbor data
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
        hvy_active, lgt_sortednumlist, hvy_work, hvy_tmp)
    ! reset the grid: all blocks are inactive and empty
    call reset_grid( params, lgt_block, hvy_block, hvy_work, hvy_tmp, hvy_neighbor, lgt_active, &
         lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
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
        hvy_tmp, hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist )
    endif

    call reset_grid( params, lgt_block, hvy_block, hvy_work, hvy_tmp, hvy_neighbor, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )


    !---------------------------------------------------------------------------
    ! Initial condition
    !---------------------------------------------------------------------------
    ! On all blocks, set the initial condition (incl. synchronize ghosts)
    call set_initial_grid( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_n, hvy_n, lgt_sortednumlist, params%adapt_inicond, time, iteration, hvy_tmp )

    if (.not. params%read_from_files .or. params%adapt_inicond) then
        ! save initial condition to disk (unless we're reading from file and do not adapt,
        ! in which case this makes no sense)
        ! we need to sync ghost nodes in order to compute the vorticity, if it is used and stored.
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

        ! NOte new versions (>16/12/2017) call physics module routines call prepare_save_data. These
        ! routines create the fields to be stored in the work array hvy_work in the first 1:params%N_fields_saved
        ! slots. the state vector (hvy_block) is copied if desired.
        call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n, hvy_tmp, hvy_active )

    end if

    if (rank==0 .and. iteration==0) then
        open (77, file='meanflow.t', status='replace')
        close(77)
        open (77, file='forces.t', status='replace')
        close(77)
        open (77, file='e_kin.t', status='replace')
        close(77)
        open (77, file='enstrophy.t', status='replace')
        close(77)
        open (44, file='dt.t', status='replace')
        close(44)
        open (44, file='timesteps_info.t', status='replace')
        close(44)
        open (44, file='blocks_per_mpirank.t', status='replace')
        close(44)
        open (44, file='blocks_per_mpirank_rhs.t', status='replace')
        close(44)
        open (44, file='eps_norm.t', status='replace')
        close(44)
        open (44, file='div.t', status='replace')
        close(44)
        open (44, file='block_dist.dat', status='replace')
        close(44)
        open (44, file='mask_volume.t', status='replace')
        close(44)
        open (44, file='u_residual.t', status='replace')
        close(44)
        open (44, file='krylov_err.t', status='replace')
        close(44)
        open (44, file='umag.t', status='replace')
        write(44,'(5(A15,1x))') "%          time","u_max","c0","MachNumber","u_eigen"
        close(44)
        open (44, file='CFL.t', status='replace')
        write(44,'(4(A15,1x))') "%          time","CFL","CFL_nu","CFL_eta"
        close(44)
    endif
    call Initialize_runtime_control_file()

    ! next write time for reloaded data
    if (params%write_method == 'fixed_time') then
        params%next_write_time = floor(time/params%next_write_time)*params%next_write_time + params%next_write_time
        params%next_stats_time = floor(time/params%next_stats_time)*params%next_stats_time + params%next_stats_time
    end if

    ! max neighbor num
    !> \todo move max neighbor num to params struct
    if ( params%threeD_case ) then
        max_neighbors = 74 ! 3D
    else
        max_neighbors = 12 ! 2D
    end if

    ! timing
    call toc( params, "init_data", MPI_wtime()-sub_t0 )

    !---------------------------------------------------------------------------
    ! main time loop
    !---------------------------------------------------------------------------
    if (rank==0) write(*,*) "starting main time loop"
    keep_running = .true.

    do while ( time<params%time_max .and. iteration<params%nt .and. keep_running)
        t2 = MPI_wtime()

        !***********************************************************************
        ! check redundant nodes
        !***********************************************************************
        t4 = MPI_wtime()
        if (params%check_redundant_nodes) then
            ! run the internal test for the ghost nodes.
            call check_unique_origin(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, test_failed)

            if (test_failed) then
                call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n, hvy_tmp, hvy_active )
                call abort(111111,"Same origin of ghost nodes check failed - stopping.")
            endif
        endif
        call toc( params, "TOPLEVEL: check ghost nodes", MPI_wtime()-t4)


        !***********************************************************************
        ! MPI bridge (used e.g. for particle-fluid coupling)
        !***********************************************************************
        if (params%bridge_exists) then
            call send_lgt_data (lgt_block,lgt_active,lgt_n,params)
            call serve_data_request(lgt_block, hvy_block, hvy_tmp, hvy_neighbor, hvy_active, &
                                    lgt_active, lgt_n, hvy_n,params)
        endif


        !***********************************************************************
        ! refine everywhere
        !***********************************************************************
        t4 = MPI_wtime()
        if ( params%adapt_mesh ) then
            ! synchronization before refinement (because the interpolation takes place on the extended blocks
            ! including the ghost nodes)
            call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

            ! refine the mesh. afterwards, it can happen that two blocks on the same level differ in their redunant nodes.
            call refine_mesh( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, lgt_active, lgt_n, &
            lgt_sortednumlist, hvy_active, hvy_n, "everywhere" )
        endif
        call toc( params, "TOPLEVEL: refinement", MPI_wtime()-t4)
        Nblocks_rhs = lgt_n

        !***********************************************************************
        ! update grid quantities
        !***********************************************************************
        ! While the state vector and many work variables (such as the mask function for penalization)
        ! are explicitly time dependent, some other quantities are not. They are rather grid-dependent
        ! but need not to be updated in every RK or krylov substep. Hence, those quantities are updated
        ! after the mesh is changed (i.e. after refine_mesh) and then kept constant during the evolution
        ! time step.
        ! An example for such a quantity would be geometry factors on non-cartesian grids, but also the
        ! body of an insect in tethered (=fixed) flight. In the latter example, only the wings need to be
        ! generated at every time t. This example generalizes to any combination of stationary and moving
        ! obstacle, i.e. insect behind fractal tree.
        ! Updating those grid-depend quantities is a task for the physics modules: they should provide interfaces,
        ! if they require such qantities. In many cases, the grid_qtys are probably not used.
        ! Please note that in the current implementation, hvy_tmp also plays the role of a work array
        t4 = MPI_wtime()
        call update_grid_qyts( time, params, lgt_block, hvy_tmp, hvy_active, hvy_n )
        call toc( params, "TOPLEVEL: update_grid_qyts", MPI_wtime()-t4)


        !***********************************************************************
        ! evolve solution in time
        !***********************************************************************
        ! internal loop over time steps: if desired, we perform more than one time step
        ! before adapting the grid again. this can further reduce the overhead of adaptivity
        ! Note: the non-linear terms can create finer scales than resolved on the grid. they
        ! are usually filtered by the coarsening/refinement round trip. So if you do more than one time step
        ! on the grid, consider using a filter.
        do it = 1, params%N_dt_per_grid
            !*******************************************************************
            ! advance in time (make one time step)
            !*******************************************************************
            t4 = MPI_wtime()
            call time_stepper( time, dt, params, lgt_block, hvy_block, hvy_work, hvy_tmp, hvy_neighbor, &
            hvy_active, lgt_active, lgt_n, hvy_n )
            call toc( params, "TOPLEVEL: time stepper", MPI_wtime()-t4)
            iteration = iteration + 1

            ! determine if it is time to save data our not.
            it_is_time_to_save_data = .false.
            if ((params%write_method=='fixed_freq' .and. modulo(iteration, params%write_freq)==0) .or. &
                (params%write_method=='fixed_time' .and. abs(time - params%next_write_time)<1e-12_rk)) then
                it_is_time_to_save_data=.true.
            endif

            !*******************************************************************
            ! filter
            !*******************************************************************
            t4 = MPI_wtime()
            if (params%filter_type /= "no_filter") then
                if (modulo(iteration, params%filter_freq) == 0 .and. params%filter_freq > 0 .or. it_is_time_to_save_data) then
                    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

                    call filter_wrapper(time, params, hvy_block, hvy_tmp, lgt_block, hvy_active, hvy_n)
                end if
            end if
            call toc( params, "TOPLEVEL: filter", MPI_wtime()-t4)

            ! it is useful to save the number of blocks per rank into a log file.
            call blocks_per_mpirank( params, blocks_per_rank, hvy_n)
            if (rank==0) then
                 open(14,file='blocks_per_mpirank_rhs.t',status='unknown',position='append')
                 write (14,'(g15.8,1x,i6,1x,i6,1x,i3,1x,i3,1x,4096(i4,1x))') time, iteration, lgt_n, &
                 min_active_level( lgt_block, lgt_active, lgt_n ), &
                 max_active_level( lgt_block, lgt_active, lgt_n ), blocks_per_rank
                 close(14)
            end if

            !*******************************************************************
            ! statistics
            !*******************************************************************
            if ( (modulo(iteration, params%nsave_stats)==0).or.(abs(time - params%next_stats_time)<1e-12_rk) ) then
                ! we need to sync ghost nodes for some derived qtys, for sure
                call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

                call statistics_wrapper(time, dt, params, hvy_block, hvy_tmp, lgt_block, hvy_active, hvy_n)
                params%next_stats_time = params%next_stats_time + params%tsave_stats
            endif
        enddo

        !***********************************************************************
        ! Adapt mesh (coarsening where possible)
        !***********************************************************************
        t4 = MPI_wtime()
        ! adapt the mesh
        if ( params%adapt_mesh ) then
            call adapt_mesh( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
            lgt_n, lgt_sortednumlist, hvy_active, hvy_n, params%coarsening_indicator, hvy_tmp )
        endif
        call toc( params, "TOPLEVEL: adapt mesh", MPI_wtime()-t4)
        Nblocks = lgt_n

        !***********************************************************************
        ! Write fields to HDF5 file
        !***********************************************************************
        if (it_is_time_to_save_data) then
            ! we need to sync ghost nodes in order to compute the vorticity, if it is used and stored.
            call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

            ! NOTE new versions (>16/12/2017) call physics module routines call prepare_save_data. These
            ! routines create the fields to be stored in the work array hvy_tmp in the first 1:params%N_fields_saved
            ! slots. the state vector (hvy_block) is copied if desired.
            call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active, &
            lgt_n, hvy_n, hvy_tmp, hvy_active )

            output_time = time
            params%next_write_time = params%next_write_time + params%write_time
        endif

        ! at the end of a time step, we increase the total counters/timers for all measurements
        ! by what has been done in the last time step, then we flush the current timing to disk.
        call timing_next_timestep( params, iteration )

        ! it is useful to save the number of blocks per rank into a log file.
        call blocks_per_mpirank( params, blocks_per_rank, hvy_n)

        t2 = MPI_wtime() - t2
        ! output on screen
        if (rank==0) then
            write(*, '("RUN: it=",i7,1x," time=",f16.9,1x,"t_cpu=",es12.4," Nb=(",i6,"/",i6,") Jmin=",i2," Jmax=",i2)') &
             iteration, time, t2, Nblocks_rhs, Nblocks, min_active_level( lgt_block, lgt_active, lgt_n ), &
             max_active_level( lgt_block, lgt_active, lgt_n )

             open(14,file='timesteps_info.t',status='unknown',position='append')
             write (14,'(2(g15.8,1x),i6,1x,i5,1x,i2,1x,i2)') time, t2, iteration, lgt_n, min_active_level( lgt_block, lgt_active, lgt_n ), &
             max_active_level( lgt_block, lgt_active, lgt_n )
             close(14)

             open(14,file='blocks_per_mpirank.t',status='unknown',position='append')
             write (14,'(g15.8,1x,i6,1x,i6,1x,i3,1x,i3,1x,4096(i4,1x))') time, iteration, lgt_n, &
             min_active_level( lgt_block, lgt_active, lgt_n ), &
             max_active_level( lgt_block, lgt_active, lgt_n ), blocks_per_rank
             close(14)
        end if

        !***********************************************************************
        ! walltime limiter
        !***********************************************************************
        ! maximum walltime allowed for simulations (in hours). The run will be stopped if this duration
        ! is exceeded. This is useful on real clusters, where the walltime of a job is limited, and the
        ! system kills the job regardless of whether we're done or not. If WABBIT itself ends execution,
        ! a backup is written and you can resume the simulation right where it stopped
        if ( (MPI_wtime()-tstart)/3600.0_rk >= params%walltime_max ) then
            if (rank==0) write(*,*) "WE ARE OUT OF WALLTIME AND STOPPING NOW!"
            keep_running = .false.
        endif

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
    if (rank==0) write(*,*) "This is the end of the main time loop!"

    ! save end field to disk, only if timestep is not saved already
    if ( abs(output_time-time) > 1e-10_rk ) then
        ! we need to sync ghost nodes in order to compute the vorticity, if it is used and stored.
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

        ! filter before write out
        if ( params%filter_freq > 0 .and. params%filter_type/="no_filter") then
            call filter_wrapper(time, params, hvy_block, hvy_tmp, lgt_block, hvy_active, hvy_n)
        end if

        ! NOte new versions (>16/12/2017) call physics module routines call prepare_save_data. These
        ! routines create the fields to be stored in the work array hvy_tmp in the first 1:params%N_fields_saved
        ! slots. the state vector (hvy_block) is copied if desired.
        call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n, hvy_tmp, hvy_active )
    end if

    ! at the end of a time step, we increase the total counters/timers for all measurements
    ! by what has been done in the last time step, then we flush the current timing to disk.
    call timing_next_timestep( params, iteration )

    ! MPI Barrier before program ends
    call MPI_Barrier(WABBIT_COMM, ierr)

    ! make a summary of the program parts, which have been profiled using toc(...)
    ! and print it to stdout
    call summarize_profiling( params, WABBIT_COMM )

    call deallocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active,&
        hvy_active, lgt_sortednumlist, hvy_work, hvy_tmp )

    ! computing time output on screen
    call cpu_time(t1)
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("END: cpu-time = ",f16.4, " s")')  t1-t0
    end if

    deallocate(blocks_per_rank)
    ! end mpi
    call MPI_Finalize(ierr)

end program main
