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
! output:   - filled user defined data structure for global params
!           - initialized light and heavy data arrays
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4, now run complete initialization within these subroutine and return
!            initialized block data to main program
! ********************************************************************************************

subroutine init_data(params, block_list, block_data, neighbor_list)

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! ini file parser module
    use module_ini_files_parser

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(out)                 :: params

    ! light data array
    integer(kind=ik), allocatable, intent(out)      :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), allocatable, intent(out)         :: block_data(:, :, :, :)
    ! neighbor array (heavy data) -> number_lines = number_blocks * 16 (...different neighbor relations:
    ! '__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW' )
    !         saved data -> -1 ... no neighbor
    !                    -> light data line number (id)
    integer(kind=ik), allocatable, intent(out)      :: neighbor_list(:)

    ! MPI error variable
    integer(kind=ik)                                :: ierr
    ! process rank
    integer(kind=ik)                                :: rank

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

!---------------------------------------------------------------------------------------------
! interfaces

    interface
        subroutine allocate_block_list(block_list, number_blocks, max_treelevel)
            use module_params
            integer(kind=ik), allocatable, intent(out)  :: block_list(:, :)
            integer(kind=ik), intent(in)                :: number_blocks
            integer(kind=ik), intent(in)                :: max_treelevel
        end subroutine allocate_block_list

        subroutine allocate_block_data(block_data, number_blocks, Bs, g, dF)
            use module_params
            real(kind=rk), allocatable, intent(out)     :: block_data(:, :, :, :)
            integer(kind=ik), intent(in)                :: number_blocks
            integer(kind=ik), intent(in)                :: Bs, g, dF
        end subroutine allocate_block_data

        subroutine inicond_gauss_blob(phi, Ds, Lx, Ly)
            use module_params
            real(kind=rk), allocatable, intent(out)     :: phi(:, :)
            integer(kind=ik), intent(in)                :: Ds
            real(kind=rk), intent(in)                   :: Lx, Ly
        end subroutine inicond_gauss_blob

        subroutine initial_block_distribution(params, block_list, block_data, phi)
            use module_params
            type (type_params), intent(in)              :: params
            integer(kind=ik), intent(inout)             :: block_list(:, :)
            real(kind=rk), intent(inout)                :: block_data(:, :, :, :)
            real(kind=rk), intent(in)                   :: phi(:, :)
        end subroutine initial_block_distribution
    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

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
    ! allocate light/heavy data, initialize start field and write block data
    !
    ! allocate block_list
    call allocate_block_list( block_list, params%number_blocks, params%max_treelevel )
    ! allocate heavy data
    call allocate_block_data( block_data, params%number_blocks, params%number_block_nodes, params%number_ghost_nodes, params%number_fields )

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
    call initial_block_distribution( params, block_list, block_data, phi )

    ! second: write heavy data for other datafields
    do k = 3, params%number_fields
        block_data( :, :, k, : ) = 0.0_rk
    end do

    ! ------------------------------------------------------------------------------------------------------
    ! init neighbor data array
    allocate( neighbor_list( params%number_blocks*16 ), stat=allocate_error )
    neighbor_list = -1

    ! clean up
    call clean_ini_file(FILE)
    deallocate( phi, stat=allocate_error )

end subroutine init_data
