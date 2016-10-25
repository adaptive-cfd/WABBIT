! ********************************
! WABBIT
! --------------------------------
!
! initialize all data (params, fields, blocks, ...)
!
! name: init_data.f90
! date: 25.10.2016
! author: msr, engels
! version: 0.3
!
! ********************************

subroutine init_data()

    use mpi
    use module_params
    use module_blocks
    use ini_files_parser

    implicit none

    type(inifile)       :: FILE
    character(len=80)   :: infile
    integer             :: read_logical, rank, ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! get the first command line argument
    call get_command_argument(1,infile)

    ! read the file, only process 0 create output on screen
    if (rank==0) then
        call read_ini_file(FILE, infile, .true.)
    else
        call read_ini_file(FILE, infile, .false.)
    end if

    !***************************************************************************
    ! BLOCK STRUCTURE PARAMETERS
    !***************************************************************************
    call read_param(FILE,'Blocks','size_domain',blocks_params%size_domain, 513 )
    call read_param(FILE,'Blocks','blocksize',blocks_params%size_block, 17 )
    call read_param(FILE,'Blocks','number_data_fields',blocks_params%number_data_fields, 1)
    call read_param(FILE,'Blocks','ghosts',blocks_params%number_ghost_nodes, 4 )

    !***************************************************************************
    ! GENERAL PARAMETERS (STRUCT "FILE")
    !***************************************************************************
    ! time loop parameter
    call read_param(FILE,'Time','time_max',params%time_max, 200.0_rk )
    call read_param(FILE,'Time','CFL',params%CFL, 0.5_rk )
    ! output write frequency
    call read_param(FILE,'Time','write_freq',params%write_freq, 25 )
    ! convective velocity
    call read_param(FILE,'Physics','u0',params%u0, (/1.0_rk, 0.5_rk/) )
    ! diffusion coeffcient
    call read_param(FILE,'Physics','nu',params%nu, 0.0_rk )
    ! domain size
    call read_param(FILE,'Physics','Lx',params%Lx, 256.0_rk )
    call read_param(FILE,'Physics','Ly',params%Ly, 256.0_rk )
    ! initial condition
    call read_param(FILE,'Physics','inicond',params%inicond, "gauss_blob" )

    call read_param(FILE,'Physics','boundary',params%boundary, "---" )

    ! eps for coarsen and refine the block
    call read_param(FILE,'Blocks','eps',params%eps, 1e-3_rk )
    ! set treelevel
    call read_param(FILE,'Blocks','max_treelevel',params%max_treelevel, 6 )
    call read_param(FILE,'Blocks','min_treelevel',params%min_treelevel, 1 )
    ! order of predictor for refinement
    call read_param(FILE,'Blocks','order_predictor',params%order_predictor, "multiresolution_4th" )
    ! discretization order
    call read_param(FILE,'Discretization','order_discretization',params%order_discretization, "FD_4th_central_optimized" )
    ! switch to turn on|off mesh refinement
    call read_param(FILE,'Blocks','adapt_mesh',read_logical, 1 )
    if (read_logical==1) then
        blocks_params%adapt_mesh = .true.
    else
        blocks_params%adapt_mesh = .false.
    end if

    ! read number of maximal blocks for memory allocation
    ! default value for number of max_blocks is: 4^(maxlevel-1) + 1, +1 needed for coarsening if all blocks at start on max_treelevel
    call read_param(FILE,'Blocks','number_max_blocks',blocks_params%number_max_blocks, 4**(params%max_treelevel)+1 )
    call read_param(FILE,'Blocks','number_max_blocks_data',blocks_params%number_max_blocks_data, blocks_params%number_max_blocks )

    ! output on screen
    if (rank==0) then
        write(*,'(80("-"))')
        write(*,*) "INITIALIZATION PHASE"
        write(*,*) "we use a maximum number of blocks of:", blocks_params%number_max_blocks
        write(*,*) "nghosts:", blocks_params%number_ghost_nodes
        write(*,'(80("-"))')
    end if

    ! allocate the individual block's memory
    call allocate_block_memory()
    ! initial data field
    call initial_condition_dense_field()
    ! clean up
    call clean_ini_file(FILE)

end subroutine init_data
