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

    use mpi
    use module_helpers
    use module_MPI
    ! global parameters
    use module_params
    ! timing module
    use module_timing
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
    character(len=80)                   :: filename
    ! loop variable
    integer(kind=ik)                    :: k, Nblocks_rhs, Nblocks, it
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, t4, tstart, dt
    ! decide if data is saved or not
    logical                             :: it_is_time_to_save_data=.false., test_failed, keep_running=.true.


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
    ! output MPI status
    WABBIT_COMM=MPI_COMM_WORLD

    if (rank==0) then
        write(*,'(80("_"))')
        write(*, '("MPI: using ", i5, " processes")') params%number_procs
        write(*,'("MPI: code build with NON-blocking send/recv in transfer (block_xfer_nonblocking.f90)")')
    end if


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
    ! have the pysics module read their own parameters. They also decide how many grid-qtys
    ! they want
    call init_physics_modules( params, filename, params%N_mask_components )
    ! allocate memory for heavy, light, work and neighbor data
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
        hvy_active, lgt_sortednumlist, hvy_work, hvy_tmp, hvy_mask, hvy_n, lgt_n)
    !  the grid: all blocks are inactive and empty
    call reset_grid( params, lgt_block, hvy_block, hvy_work, hvy_tmp, hvy_neighbor, lgt_active(:,tree_ID_flow), &
         lgt_n(tree_ID_flow), hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow), .true. )
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
        hvy_tmp, hvy_neighbor, lgt_active(:,tree_ID_flow), hvy_active(:,tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow) )
    endif

    call reset_grid( params, lgt_block, hvy_block, hvy_work, hvy_tmp, hvy_neighbor, lgt_active(:,tree_ID_flow), &
    lgt_n(tree_ID_flow), hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow), .true. )


    !---------------------------------------------------------------------------
    ! Initial condition
    !---------------------------------------------------------------------------
    ! On all blocks, set the initial condition (incl. synchronize ghosts)
    call set_initial_grid( params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,tree_ID_flow), hvy_active(:,tree_ID_flow), &
    lgt_n(tree_ID_flow), hvy_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow), params%adapt_inicond, time, iteration, hvy_mask )

    if (.not. params%read_from_files .or. params%adapt_inicond) then
        ! save initial condition to disk (unless we're reading from file and do not adapt,
        ! in which case this makes no sense)
        ! we need to sync ghost nodes in order to compute the vorticity, if it is used and stored.
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

        ! NOte new versions (>16/12/2017) call physics module routines call prepare_save_data. These
        ! routines create the fields to be stored in the work array hvy_work in the first 1:params%N_fields_saved
        ! slots. the state vector (hvy_block) is copied if desired.
        call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active(:,tree_ID_flow), &
        lgt_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow), hvy_n(tree_ID_flow), hvy_tmp, hvy_active(:,tree_ID_flow), hvy_mask,&
        hvy_neighbor )

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
        open (44, file='performance.t', status='replace')
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
        open (14, file='moments.t', status='replace')
        close(14)
        open (14, file='aero_power.t', status='replace')
        close(14)
        open (14, file='forces_body.t', status='replace')
        close(14)
        open (14, file='moments_body.t', status='replace')
        close(14)
        open (14, file='forces_leftwing.t', status='replace')
        close(14)
        open (14, file='moments_leftwing.t', status='replace')
        close(14)
        open (14, file='forces_rightwing.t', status='replace')
        close(14)
        open (14, file='moments_rightwing.t', status='replace')
        close(14)
    endif

    if (rank==0) then
        call Initialize_runtime_control_file()
    endif

    ! next write time for reloaded data
    if (params%write_method == 'fixed_time') then
        params%next_write_time = floor(time/params%next_write_time)*params%next_write_time + params%next_write_time
        params%next_stats_time = floor(time/params%next_stats_time)*params%next_stats_time + params%next_stats_time
    end if


    ! timing
    call toc( "init_data", MPI_wtime()-sub_t0 )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! mask
    ! call read_field2tree(params, (/"chi_00.h5"/), 1, 2, tree_n, &
    ! lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    ! hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)
    !
    ! call prune_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    !     hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id=2)

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
            call check_unique_origin(params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), &
            hvy_n(tree_ID_flow), test_failed)

            if (test_failed) then
                call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active(:,tree_ID_flow), &
                lgt_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow), hvy_n(tree_ID_flow), hvy_tmp, &
                hvy_active(:,tree_ID_flow), hvy_mask, hvy_neighbor )
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


        !***********************************************************************
        ! refine everywhere
        !***********************************************************************
        t4 = MPI_wtime()
        if ( params%adapt_mesh ) then
            ! synchronization before refinement (because the interpolation takes place on the extended blocks
            ! including the ghost nodes)
            ! Note: at this point the grid is rather coarse (fewer blocks), and the sync step is rather cheap.
            ! Snych'ing becomes much mor expensive one the grid is refined.
            call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

            ! refine the mesh. Note: afterwards, it can happen that two blocks on the same level differ
            ! in their redundant nodes, but the ghost node sync'ing later on will correct these mistakes.
            call refine_mesh( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, lgt_active(:,tree_ID_flow), lgt_n(tree_ID_flow), &
            lgt_sortednumlist(:,:,tree_ID_flow), hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), "everywhere" )
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
            call time_stepper( time, dt, params, lgt_block, hvy_block, hvy_work, hvy_mask, hvy_neighbor, &
            hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), lgt_active(:,tree_ID_flow), &
            lgt_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow) )
            call toc( "TOPLEVEL: time stepper", MPI_wtime()-t4)
            iteration = iteration + 1

            ! determine if it is time to save data
            it_is_time_to_save_data = .false.
            if ((params%write_method=='fixed_freq' .and. modulo(iteration, params%write_freq)==0) .or. &
                (params%write_method=='fixed_time' .and. abs(time - params%next_write_time)<1e-12_rk)) then
                it_is_time_to_save_data= .true.
            endif

            !*******************************************************************
            ! filter
            !*******************************************************************
            t4 = MPI_wtime()
            if (params%filter_type /= "no_filter") then
                if (modulo(iteration, params%filter_freq) == 0 .and. params%filter_freq > 0 .or. it_is_time_to_save_data) then
                    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow) )

                    call filter_wrapper(time, params, hvy_block, hvy_tmp, lgt_block, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow))
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
                lgt_active(:,tree_ID_flow), lgt_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow), &
                hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), hvy_neighbor)

                params%next_stats_time = params%next_stats_time + params%tsave_stats
            endif
        enddo

        !***********************************************************************
        ! Adapt mesh (coarsening where possible)
        !***********************************************************************
        t4 = MPI_wtime()
        ! adapt the mesh
        if ( params%adapt_mesh ) then
            call adapt_mesh( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,tree_ID_flow), &
            lgt_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow), hvy_active(:,tree_ID_flow), &
            hvy_n(tree_ID_flow), params%coarsening_indicator, hvy_tmp, hvy_mask )
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
            call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active(:,tree_ID_flow), &
            lgt_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow),  hvy_n(tree_ID_flow), hvy_tmp, &
            hvy_active(:,tree_ID_flow), hvy_mask, hvy_neighbor )

            output_time = time
            params%next_write_time = params%next_write_time + params%write_time
        endif

        t2 = MPI_wtime() - t2
        ! output on screen
        if (rank==0) then
            write(*, '("RUN: it=",i7,1x," time=",f16.9,1x,"t_cpu=",es12.4," Nb=(",i6,"/",i6,") Jmin=",i2," Jmax=",i2)') &
             iteration, time, t2, Nblocks_rhs, Nblocks, &
             min_active_level( lgt_block, lgt_active(:,tree_ID_flow), lgt_n(tree_ID_flow) ), &
             max_active_level( lgt_block, lgt_active(:,tree_ID_flow), lgt_n(tree_ID_flow) )

             ! prior to 11/04/2019, this file was called timesteps_info.t but it was missing some important
             ! information, so I renamed it when adding those (since post-scripts would no longer be compatible
             ! it made sense to me to change the name)
             open(14,file='performance.t',status='unknown',position='append')
             write (14,'(g15.8,1x,i9,1x,g15.8,1x,i6,1x,i6,1x,i2,1x,i2,1x,i6)') time, iteration, t2, Nblocks_rhs, Nblocks, &
             min_active_level( lgt_block, lgt_active(:,tree_ID_flow), lgt_n(tree_ID_flow) ), &
             max_active_level( lgt_block, lgt_active(:,tree_ID_flow), lgt_n(tree_ID_flow) ), &
             params%number_procs
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

    ! save end field to disk, only if this data is not saved already
    if ( abs(output_time-time) > 1e-10_rk ) then
        ! we need to sync ghost nodes in order to compute the vorticity, if it is used and stored.
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow))

        ! filter before write out
        if ( params%filter_freq > 0 .and. params%filter_type/="no_filter") then
            call filter_wrapper(time, params, hvy_block, hvy_tmp, lgt_block, hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow))
        end if

        ! Note new versions (>16/12/2017) call physics module routines call prepare_save_data. These
        ! routines create the fields to be stored in the work array hvy_tmp in the first 1:params%N_fields_saved
        ! slots. the state vector (hvy_block) is copied if desired.
        call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active(:,tree_ID_flow), &
        lgt_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow), hvy_n(tree_ID_flow), &
        hvy_tmp, hvy_active(:,tree_ID_flow), hvy_mask, hvy_neighbor )
    end if

    ! MPI Barrier before program ends
    call MPI_Barrier(WABBIT_COMM, ierr)

    ! make a summary of the program parts, which have been profiled using toc(...)
    ! and print it to stdout
    call summarize_profiling( WABBIT_COMM )

    call deallocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active,&
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
    if (rank==0) then
        open (77, file='success', status='replace')
        close(77)
    endif

    ! end mpi
    call MPI_Finalize(ierr)

end program main
