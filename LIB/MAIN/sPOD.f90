program sPOD

    use mpi
    use module_helpers
    use module_forest
    use module_MPI
    use module_params
    use module_debug
    use module_mesh
    use module_bridge_interface

    implicit none
    ! -----------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------------
    ! MPI data
    integer(kind=ik)   :: ierr, rank, number_procs
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0, t1, t2
    ! user defined parameter structure
    type (type_params)                  :: params
    ! light data array  -> line number = ( 1 + proc_rank ) * heavy_data_line_number
    !                   -> column(1:max_treelevel): block treecode, treecode -1 => block is inactive
    !                   -> column(max_treelevel + idx_mesh_lvl): treecode length = mesh level
    !                   -> column(max_treelevel + idx_refine_sts):   refinement status (-1..coarsen / 0...no change / +1...refine)
    integer(kind=ik), allocatable       :: lgt_block(:, :, :)
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
    integer(kind=ik), allocatable       :: lgt_active(:,:)
    ! number of active blocks (light data)
    integer(kind=ik)                    :: lgt_n
    ! list of active blocks (heavy data)
    integer(kind=ik), allocatable       :: hvy_active(:)
    ! number of active blocks (heavy data)
    integer(kind=ik)                    :: hvy_n
    integer(kind=ik), allocatable       :: blocks_per_rank(:)
    integer(kind=ik)                    :: iteration
    ! filename of *.ini file used to read parameters
    character(len=80)                   :: filename
    ! loop variable
    integer(kind=ik)                    :: k, Nblocks_rhs, Nblocks, it,tree_id
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, t4, tstart, dt, time
    ! decide if data is saved or not
    logical                             :: it_is_time_to_save_data=.false., test_failed, keep_running=.true.
    ! -----------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------------

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

    ! start time
    sub_t0 = MPI_Wtime()
    call cpu_time(t0)
    tstart = MPI_wtime()


    if ( rank == 0 ) then
        write(*,*) " HELLO here is sPODi your personal assistant..."
        write(*,*) " What can i do for you ?"
        write(*,*) ""
    end if
    !---------------------------------------------------------------------------
    ! Initialize parameters,bridge and grid
    !---------------------------------------------------------------------------
    ! read in the parameter file to setup the case
    ! get the second command line argument: this should be the ini-file name
    call get_command_argument( 1, filename )
    ! read ini-file and save parameters in struct
    call read_forest_params( params, filename )
    ! initializes the communicator for Wabbit and creates a bridge if needed
    call initialize_communicator(params)

    allocate(params%input_files(1))
    params%input_files = "rho_000000000000.h5"
    !call get_command_argument( 2, params%input_files(1) )

    tree_id=1
    call read_field2tree(params, params%input_files, size(params%input_files), &
                        lgt_block, hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)
    call read_field2tree(params, params%input_files, size(params%input_files), tree_id, &
                        lgt_block, lgt_active(:,tree_id), lgt_n, lgt_sortednumlist, &
                        hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)
    keep_running = .true.
    iteration     = 0

! #################################################################################
!                   MAIN LOOP
! #################################################################################
    call abort(23456,"stop here")
! #################################################################################
    if (rank==0) write(*,*) "Completed all tasks! Ending Program"
    ! MPI Barrier before program ends
    call MPI_Barrier(WABBIT_COMM, ierr)

    ! make a summary of the program parts, which have been profiled using toc(...)
    ! and print it to stdout
    call summarize_profiling( params, WABBIT_COMM )

    !call deallocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active,&
        !hvy_active, lgt_sortednumlist, hvy_work, hvy_tmp )

    ! computing time output on screen
    call cpu_time(t1)
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("END: cpu-time = ",f16.4, " s")')  t1-t0
    end if

    deallocate(blocks_per_rank)
    ! end mpi
    call MPI_Finalize(ierr)

end program sPOD
