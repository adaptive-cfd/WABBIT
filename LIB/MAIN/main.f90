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
    ! module of bridge to other mpi worlds
    use bridgefluid
    !>TODO is this still needed here ???????
    use module_ACM_new

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
    real(kind=rk)                       :: t0, t1

    ! user defined parameter structure
    type (type_params)                  :: params

    ! light data array  -> line number = ( 1 + proc_rank ) * heavy_data_line_number
    !                   -> column(1:max_treelevel): block treecode, treecode -1 => block is inactive
    !                   -> column(max_treelevel+1): treecode length = mesh level
    !                   -> column(max_treelevel+2):   refinement status (-1..coarsen / 0...no change / +1...refine)
    integer(kind=ik), allocatable       :: lgt_block(:, :)

    !                   -> dim 1: x coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 2: y coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 3: z coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 4: data type ( 1:number_data_fields)
    ! heavy data array  -> dim 5: block id  ( 1:number_blocks )
    real(kind=rk), allocatable          :: hvy_block(:, :, :, :, :)

    !                   -> dim 1: x coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 2: y coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 3: z coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 4: data type ( old data, k1, k2, k3, k4 )
    ! heavy work array  -> dim 5: block id  ( 1:number_blocks )
    real(kind=rk), allocatable          :: hvy_work(:, :, :, :, :)

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

    ! time loop variables
    real(kind=rk)                       :: time, output_time
    integer(kind=ik)                    :: iteration

    ! filename of *.ini file used to read parameters
    character(len=80)                   :: filename

    ! loop variable
    integer(kind=ik)                    :: k, max_neighbors

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0

    ! allocate com lists and com matrix here
    ! communication lists:
    ! dim 1: list elements
    ! dim 2: columns
    !                       1   rank of sender process
    !                       2   rank of receiver process
    !                       3   sender block heavy data id
    !                       4   receiver block heavy data id
    !                       5   sender block neighborhood to receiver (dirs id)
    !                       6   difference between sender-receiver level
    ! dim 3: receiver proc rank
    ! dim 4: synch stage
    integer(kind=ik), allocatable       :: com_lists(:, :, :, :)

    ! communications matrix:
    ! count the number of communications between procs
    ! row/column number encodes process rank + 1
    ! com matrix pos: position in send buffer
    ! dim 3: synch stage
    integer(kind=ik), allocatable       :: com_matrix(:,:,:)

    ! send/receive buffer, integer and real
    ! allocate in init substep not in synchronize subroutine, to avoid slow down when using
    ! large numbers of processes and blocks per process
    integer(kind=ik), allocatable       :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), allocatable          :: real_send_buffer(:,:), real_receive_buffer(:,:)

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
    call MPI_Comm_rank(WABBIT_COMM, rank, ierr)
    params%rank         = rank
    ! determine process number
    call MPI_Comm_size(WABBIT_COMM, number_procs, ierr)
    params%number_procs = number_procs
    ! output MPI status
    if (rank==0) then
        write(*,'(80("_"))')
        write(*, '("MPI: using ", i5, " processes")') params%number_procs
    end if


    ! start time
    sub_t0 = MPI_Wtime()
    call cpu_time(t0)


    ! unit test off
    params%unit_test    = .false.

    ! are we running in 2D or 3D mode? Check that from the command line call.
    call decide_if_running_2D_or_3D(params)

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
    call init_physics_modules( params, filename )
    ! allocate memory for heavy, light, work and neighbor data
    call allocate_grid( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )
    ! reset the grid: all blocks are inactive and empty
    call reset_grid( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
    ! initalize debugging ( this is mainly time measurements )
    call allocate_init_debugging( params )
    ! allocate communication arrays
    call allocate_com_arrays(params, com_lists, com_matrix)

    !---------------------------------------------------------------------------
    ! Unit tests
    !---------------------------------------------------------------------------
    if (params%test_treecode) then
       call unit_test_treecode( params )
    end if
    ! perform a convergence test on ghost node sync'ing
    ! I don't see a good reason to skip this test ever - I removed the condition here.
    call unit_test_ghost_nodes_synchronization( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist, com_lists, com_matrix, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )
    call reset_grid( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )

!    if (params%test_wavelet_comp) then
!        call unit_test_wavelet_compression( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active )
        ! reset the grid: all blocks are inactive and empty
!        call reset_grid( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
!    end if


    !---------------------------------------------------------------------------
    ! Initial condition
    !---------------------------------------------------------------------------
    ! On all blocks, set the initial condition
    call set_initial_grid( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_n, hvy_n, lgt_sortednumlist, params%adapt_inicond, com_lists, com_matrix, int_send_buffer, &
    int_receive_buffer, real_send_buffer, real_receive_buffer, time, iteration )

    if (params%initial_cond /= "read_from_files") then
        ! save initial condition to disk
        ! we need to sync ghost nodes in order to compute the vorticity, if it is used and stored.
        call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, .true., int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )

        ! NOte new versions (>16/12/2017) call physics module routines call prepare_save_data. These
        ! routines create the fields to be stored in the work array hvy_work in the first 1:params%N_fields_saved
        ! slots. the state vector (hvy_block) is copied if desired.
        call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n, hvy_work, hvy_active )
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

    do while ( time<params%time_max .and. iteration<params%nt)

        ! new iteration
        iteration = iteration + 1

        ! refine everywhere
        if ( params%adapt_mesh ) then
            call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, .true., int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )
            call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, "everywhere" )
        endif
     !+++++++++++ serve any data request from the other side +++++++++++++
        if (params%bridge_exists) then
            call send_lgt_data (lgt_block,lgt_active,lgt_n,params)
            call serve_data_request(lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n,params)    
        endif
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! advance in time
        call time_stepper( time, params, lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n, com_lists(1:hvy_n*max_neighbors,:,:,:), com_matrix, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )

        ! filter
        if (modulo(iteration, params%filter_freq) == 0 .and. params%filter_freq > 0 .and. params%filter_type/="no_filter") then
            call filter_block( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists(1:hvy_n*max_neighbors,:,:,:), com_matrix, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )
        end if

        ! adapt the mesh
        if ( params%adapt_mesh ) then
            call adapt_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, "threshold", com_lists, com_matrix, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )
        endif

        ! statistics
        if ( (modulo(iteration, params%nsave_stats)==0).or.(abs(time - params%next_stats_time)<1e-12_rk) ) then
          ! we need to sync ghost nodes for some derived qtys, for sure
          call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, .true., int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )

          ! TODO make this nicer
          if (iteration==1 ) then
            open (15, file='meanflow.t', status='replace')
            close(15)
          endif

          call statistics_wrapper(time, params, hvy_block, hvy_work, lgt_block, hvy_active, hvy_n)
          params%next_stats_time = params%next_stats_time + params%tsave_stats
        endif

        ! write data to disk
        if ( (params%write_method=='fixed_freq' .and. modulo(iteration, params%write_freq)==0).or.(params%write_method=='fixed_time' .and. abs(time - params%next_write_time)<1e-12_rk) ) then
          ! we need to sync ghost nodes in order to compute the vorticity, if it is used and stored.
          call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, .true., int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )

          ! NOTE new versions (>16/12/2017) call physics module routines call prepare_save_data. These
          ! routines create the fields to be stored in the work array hvy_work in the first 1:params%N_fields_saved
          ! slots. the state vector (hvy_block) is copied if desired.
          call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n, hvy_work, hvy_active )

          output_time = time
          params%next_write_time = params%next_write_time + params%write_time
        endif

        ! at the end of a time step, we increase the total counters/timers for all measurements
        ! by what has been done in the last time step, then we flush the current timing to disk.
        call timing_next_timestep( params, iteration )

        ! output on screen
        if (rank==0) then
            write(*, '("RUN: iteration=",i7,3x," time=",f16.9,3x," active blocks=",i7," Jmin=",i2," Jmax=",i2)') &
             iteration, time, lgt_n, min_active_level( lgt_block, lgt_active, lgt_n ), &
             max_active_level( lgt_block, lgt_active, lgt_n )

             open(14,file='timesteps_info.t',status='unknown',position='append')
             write (14,'((g15.8,1x),i6,1x,i5,1x,i2,1x,i2)') time, iteration, lgt_n, min_active_level( lgt_block, lgt_active, lgt_n ), &
             max_active_level( lgt_block, lgt_active, lgt_n )
             close(14)
        end if
    end do
    !---------------------------------------------------------------------------
    ! end of main time loop
    !---------------------------------------------------------------------------
    if (rank==0) write(*,*) "This is the end of the main time loop!"

    ! save end field to disk, only if timestep is not saved allready
    if ( abs(output_time-time) > 1e-10_rk ) then
      ! we need to sync ghost nodes in order to compute the vorticity, if it is used and stored.
      call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, .true., int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )

      ! NOte new versions (>16/12/2017) call physics module routines call prepare_save_data. These
      ! routines create the fields to be stored in the work array hvy_work in the first 1:params%N_fields_saved
      ! slots. the state vector (hvy_block) is copied if desired.
      call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n, hvy_work, hvy_active )
    end if

    ! at the end of a time step, we increase the total counters/timers for all measurements
    ! by what has been done in the last time step, then we flush the current timing to disk.
    call timing_next_timestep( params, iteration )

    ! MPI Barrier before program ends
    call MPI_Barrier(WABBIT_COMM, ierr)

    ! debug info output
    if ( params%debug ) then
        ! sum times
        debug%comp_time(:,2) = 0.0_rk
        call MPI_Allreduce(debug%comp_time(:,4), debug%comp_time(:,2), size(debug%comp_time,1), MPI_REAL8, MPI_SUM, WABBIT_COMM, ierr)
        ! MPI Barrier before program ends
        call MPI_Barrier(WABBIT_COMM, ierr)

        ! average times
        debug%comp_time(:,2) = debug%comp_time(:,2) / params%number_procs
        ! standard deviation
        debug%comp_time(:,3) = 0.0_rk
        debug%comp_time(:,4) = (debug%comp_time(:,4) - debug%comp_time(:,2))**2.0_rk
        call MPI_Allreduce(debug%comp_time(:,4), debug%comp_time(:,3), size(debug%comp_time,1), MPI_REAL8, MPI_SUM, WABBIT_COMM, ierr)
        ! MPI Barrier before program ends
        call MPI_Barrier(WABBIT_COMM, ierr)

        if (params%number_procs == 1) then
            debug%comp_time(:,3) = 0.0_rk
        else
            debug%comp_time(:,3) = sqrt(debug%comp_time(:,3) / ( params%number_procs - 1 ))
        end if

        ! output
        if (rank==0) then
            write(*,'(80("_"))')
            write(*, '("debug times (average value +- standard deviation) :")')
            k = 1
            do while ( debug%name_comp_time(k) /= "---" )

                ! write name
                write(*, '(a)', advance='no') debug%name_comp_time(k)
                ! write average time
                write(*, '(2x,f12.3)', advance='no') debug%comp_time(k,2)
                ! write standard deviation
                write(*, '(2x,f12.3)', advance='no') debug%comp_time(k,3)
                ! next line
                write(*,*)
                ! loop variable
                k = k + 1

            end do

            write(*,'(80("_"))')
            write(*, '("sum: ", 2x,f12.3)', advance='yes') sum(debug%comp_time(:,2))

        end if

    end if

    ! computing time output on screen
    call cpu_time(t1)
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("END: cpu-time = ",f16.4, " s")')  t1-t0
    end if

    ! end mpi
    call MPI_Finalize(ierr)

end program main
