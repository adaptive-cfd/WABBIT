! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: init_data.f90
! version: 0.4
! author: msr
!
! initialize all data: read params from ini file, allocate memory, initialize starting condition
! and decompose start matrix into block data
!
! input:    - parameter array
!           - light data array
!           - heavy data array
!           - neighbor data array
!           - light and heavy active block list
! output:   - filled user defined data structure for global params
!           - initialized light and heavy data arrays
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4, now run complete initialization within these subroutine and return
!            initialized block data to main program
! ********************************************************************************************
subroutine init_data(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(out)                 :: params

    ! light data array
    integer(kind=ik), allocatable, intent(out)      :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk), allocatable, intent(out)         :: hvy_block(:, :, :, :)
    ! neighbor array (heavy data) -> number_lines = number_blocks * 16 (...different neighbor relations:
    ! '__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW' )
    !         saved data -> -1 ... no neighbor
    !                    -> light data line number (id)
    integer(kind=ik), allocatable, intent(out)      :: hvy_neighbor(:,:)
    ! list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)      :: lgt_active(:)
    ! list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)      :: hvy_active(:)

    ! MPI error variable
    integer(kind=ik)                                :: ierr
    ! process rank
    integer(kind=ik)                                :: rank
    ! number of processes
    integer(kind=ik)                                :: number_procs

    ! inifile name
    character(len=80)                               :: filename
    ! inifile structure
    type(inifile)                                   :: FILE

    ! auxiliary variable for reading logicals
    integer(kind=ik)                                :: read_logical
    ! allocation error variabel
    integer(kind=ik)                                :: allocate_error

    ! initial data field
    real(kind=rk), allocatable                      :: phi(:, :)

    ! loop variable
    integer(kind=ik)                                :: k

    ! cpu time variables for running time calculation
    real(kind=rk)                                   :: sub_t0, sub_t1

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    ! determinate process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)

    ! get the first command line argument
    call get_command_argument(1, filename)

    ! read the file, only process 0 should create output on screen
    if (rank==0) then
        call read_ini_file(FILE, filename, .true.)
    else
        call read_ini_file(FILE, filename, .false.)
    end if

    !***************************************************************************
    ! read BLOCK parameters
    !
    ! read number_domain_nodes
    call read_param(FILE, 'Blocks', 'number_domain_nodes', params%number_domain_nodes, 1 )
    ! read number_block_nodes
    call read_param(FILE, 'Blocks', 'number_block_nodes', params%number_block_nodes, 1 )
    ! read number_ghost_nodes
    call read_param(FILE, 'Blocks', 'number_ghost_nodes', params%number_ghost_nodes, 1 )
    ! read number_blocks
    call read_param(FILE, 'Blocks', 'number_blocks', params%number_blocks, 1 )
    if (params%number_blocks == -1) then
      if (rank==0) write(*,*) "blocks-per-rank is -1, so WABBIT decides automatically"
      params%number_blocks = ( (params%number_domain_nodes-1) / (params%number_block_nodes-1) )**2 / number_procs
    endif

    ! read number_data_fields
    call read_param(FILE, 'Blocks', 'number_data_fields', params%number_data_fields, 1 )
    ! set number of fields in heavy data
    params%number_fields = params%number_data_fields + 6
    ! read threshold value
    call read_param(FILE, 'Blocks', 'eps', params%eps, 1e-3_rk )
    ! read treelevel bounds
    call read_param(FILE, 'Blocks', 'max_treelevel', params%max_treelevel, 5 )
    call read_param(FILE, 'Blocks', 'min_treelevel', params%min_treelevel, 1 )
    ! read switch to turn on|off mesh refinement
    call read_param(FILE, 'Blocks', 'adapt_mesh', read_logical, 1 )
    if ( read_logical == 1 ) then
        params%adapt_mesh = .true.
    else
        params%adapt_mesh = .false.
    end if
    ! block distribution
    call read_param(FILE, 'Blocks', 'block_dist', params%block_distribution, "---" )

    !***************************************************************************
    ! read TIME parameters
    !
    ! read time_max
    call read_param(FILE, 'Time', 'time_max', params%time_max, 1.0_rk )
    ! read CFL number
    call read_param(FILE, 'Time', 'CFL', params%CFL, 0.5_rk )
    ! read output write frequency
    call read_param(FILE, 'Time', 'write_freq', params%write_freq, 25 )

    !***************************************************************************
    ! read PHYSICS parameters
    !
    ! first: allocate memory in params structure (need 2*data_fields for velocity
    ! and 1*data_fields for diffusion coefficient)
    allocate( params%u0( 2*params%number_data_fields ), stat=allocate_error )
    allocate( params%nu( params%number_data_fields ), stat=allocate_error )
    ! read velocity
    call read_param(FILE, 'Physics', 'u0', params%u0, params%u0 )
    ! read diffusion
    call read_param(FILE, 'Physics', 'nu', params%nu, params%nu )
    ! domain size
    call read_param(FILE, 'Physics', 'Lx', params%Lx, 256.0_rk )
    call read_param(FILE, 'Physics', 'Ly', params%Ly, 256.0_rk )
    ! initial condition
    call read_param(FILE, 'Physics', 'initial_cond', params%initial_cond, "---" )

    !***************************************************************************
    ! read DISCRETIZATION parameters
    !
    ! discretization order
    call read_param(FILE, 'Discretization', 'order_discretization', params%order_discretization, "---" )
    ! order of predictor for refinement
    call read_param(FILE, 'Discretization', 'order_predictor', params%order_predictor, "---" )
    ! boundary condition
    call read_param(FILE, 'Discretization', 'boundary_cond', params%boundary_cond, "---" )

    !***************************************************************************
    ! read DEBUG parameters
    !
    ! discretization order
    call read_param(FILE, 'Debug', 'debug', read_logical, 1 )
    if ( read_logical == 1 ) then
        params%debug = .true.
    else
        params%debug = .false.
    end if

    !***************************************************************************
    ! allocate light/heavy data, initialize start field and write block data
    !
    ! allocate block_list
    call allocate_block_list( lgt_block, params%number_blocks, params%max_treelevel )
    ! allocate heavy data
    call allocate_block_data( hvy_block, params%number_blocks, params%number_block_nodes, params%number_ghost_nodes, params%number_fields )

    ! initial data field
    select case( params%initial_cond )
        case ("gauss-blob","gauss_blob")
              call inicond_gauss_blob(phi, params%number_domain_nodes, params%Lx, params%Ly)

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: initial condition is unknown"
            write(*,*) params%initial_cond
            stop

    end select

    ! decompose init field phi to block data
    ! first: init light and heavy data for datafield 1, create starting block distribution
    call initial_block_distribution( params, lgt_block, hvy_block, phi )

    ! second: write heavy data for other datafields
    do k = 3, params%number_fields
        !block_data( :, :, k, : ) = 0.0_rk
        hvy_block( :, :, k, : ) = hvy_block( :, :, 2, : )
    end do

    ! allocate active list
    allocate( lgt_active( size(lgt_block, 1) ), stat=allocate_error )
    allocate( hvy_active( size(hvy_block, 4) ), stat=allocate_error )

    ! ------------------------------------------------------------------------------------------------------
    ! init neighbor data array
    allocate( hvy_neighbor( params%number_blocks, 16 ), stat=allocate_error )
    hvy_neighbor = -1

    ! ------------------------------------------------------------------------------------------------------
    ! init debug data
    ! note: fix size of time measurements array
    if ( params%debug ) then
        ! allocate array for time measurements - data
        allocate( debug%comp_time( 20, 3 ), stat=allocate_error )
        debug%comp_time = 0.0_rk
        ! allocate array for time measurements - names
        allocate( debug%name_comp_time( 20 ), stat=allocate_error )
        debug%name_comp_time = "---"
    end if

    ! clean up
    call clean_ini_file(FILE)
    deallocate( phi, stat=allocate_error )

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "init_data" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "init_data"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0

    end if

end subroutine init_data
