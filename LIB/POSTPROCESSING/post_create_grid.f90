subroutine post_create_grid(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_helpers
    use module_forestMetaData

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=cshort) :: grid_type, fname_out
    integer(kind=ik) :: level, iterations, tree_ID, Bs_input(3), Jmax_input, Jmin_input, dim_input
    real(kind=rk) :: max_grid_density_input
    real(kind=rk), allocatable :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable :: hvy_tmp(:, :, :, :, :)
    real(kind=rk) :: time

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(2, grid_type)

    ! does the user need help?
    if (grid_type=='--help' .or. grid_type=='--h' .or. grid_type=='-h') then
        if (params%rank==0) then
            write(*,*) "--------------------------------------------------------------------------------"
            write(*,*) "Create and save a grid to file"
            write(*,*) "--------------------------------------------------------------------------------"
            write(*,*) ""
            write(*,*) "Usage:"
            write(*,*) ""
            write(*,*) "  Create an equidistant grid:"
            write(*,*) "    ./wabbit-post --create-grid equidistant OUTPUT.h5 --Bs=16 --Jmax=5 [--Jmin=1] [--level=4] [--dim=2]"
            write(*,*) ""
            write(*,*) "  Create a random grid:"
            write(*,*) "    ./wabbit-post --create-grid random OUTPUT.h5 --Bs=16 --Jmax=5 [--Jmin=1] [--level=3] [--iterations=10] [--dim=2] [--max-grid-density=0.5]"
            write(*,*) ""
            write(*,*) "--------------------------------------------------------------------------------"
            write(*,*) "Parameters:"
            write(*,*) ""
            write(*,*) "  grid_type: 'equidistant' or 'random'"
            write(*,*) "      Type of grid to create"
            write(*,*) ""
            write(*,*) "  OUTPUT.h5:"
            write(*,*) "      Output filename for the generated grid"
            write(*,*) ""
            write(*,*) "  --Bs=N (required)"
            write(*,*) "      Block size (number of grid points per block per dimension)"
            write(*,*) ""
            write(*,*) "  --Jmax=N (required)"
            write(*,*) "      Maximum refinement level"
            write(*,*) ""
            write(*,*) "  --Jmin=N (default: 1)"
            write(*,*) "      Minimum refinement level"
            write(*,*) ""
            write(*,*) "  --dim=N (default: 2)"
            write(*,*) "      Dimension of the grid (2 or 3)"
            write(*,*) ""
            write(*,*) "  --level=N (default: 1 for equidistant, 3 for random)"
            write(*,*) "      For equidistant: level of the grid (all blocks on this level)"
            write(*,*) "      For random: initial level before random refinement/coarsening"
            write(*,*) ""
            write(*,*) "  --iterations=N (default: 10, only for random grids)"
            write(*,*) "      Number of refinement/coarsening iterations for random grid generation"
            write(*,*) ""
            write(*,*) "  --max-grid-density=X (default: 0.8, only for random grids)"
            write(*,*) "      Maximum grid density (fraction of total available blocks)"
            write(*,*) ""
            write(*,*) "--------------------------------------------------------------------------------"
        end if
        return
    endif

    ! Check if grid_type is valid
    if (grid_type /= 'equidistant' .and. grid_type /= 'random') then
        if (params%rank == 0) then
            write(*,*) "ERROR: Invalid grid type: ", trim(adjustl(grid_type))
            write(*,*) "Must be 'equidistant' or 'random'"
            write(*,*) "Use --help for usage information"
        endif
        call abort(123456, "Invalid grid type for --create-grid")
    endif

    ! Read output filename
    call get_command_argument(3, fname_out)

    !-----------------------------------------------------------------------------------------------------
    ! Read parameters from command line
    !-----------------------------------------------------------------------------------------------------
    ! Required parameters
    call get_cmd_arg("--Bs", Bs_input(1), default=-1_ik)
    call get_cmd_arg("--Jmax", Jmax_input, default=-1_ik)

    ! Check that required parameters were provided
    if (Bs_input(1) == -1 .or. Jmax_input == -1) then
        if (params%rank == 0) then
            write(*,*) "ERROR: Missing required parameters!"
            write(*,*) "Required: --Bs=N --Jmax=N"
            write(*,*) "Use --help for usage information"
        endif
        call abort(123457, "Missing required parameters for --create-grid")
    endif

    ! Optional parameters
    call get_cmd_arg("--Jmin", Jmin_input, default=1_ik)
    call get_cmd_arg("--dim", dim_input, default=2_ik)
    call get_cmd_arg("--max-grid-density", max_grid_density_input, default=0.8_rk)

    ! Set up parameters
    Bs_input(2) = Bs_input(1)
    if (dim_input == 3) then
        Bs_input(3) = Bs_input(1)
    else
        Bs_input(3) = 1
    endif
    
    params%wavelet = "CDF20"
    params%Bs = Bs_input
    params%Jmax = Jmax_input
    params%Jmin = Jmin_input
    params%dim = dim_input
    params%n_eqn = 1
    params%domain_size = (/ 1.0_rk, 1.0_rk, 1.0_rk /)
    params%number_blocks = 2**(params%dim * params%Jmax + 2)  ! allocate enough memory
    params%forest_size = 1  ! we only need one tree
    params%block_distribution = "sfc_hilbert"
    params%time_step_method = 'none'
    params%order_predictor = "multiresolution_4th"
    params%physics_type = "ACM-new"
    params%max_grid_density = max_grid_density_input

    if (params%rank == 0) then
        write(*,*) "--------------------------------------------------------------------------------"
        write(*,*) "Creating ", trim(adjustl(grid_type)), " grid"
        write(*,*) "Output file:    ", trim(adjustl(fname_out))
        write(*,'("Block size:     Bs=",i3)') params%Bs(1)
        write(*,'("Ghost nodes:    g=",i2)') params%g
        write(*,'("Min level:      Jmin=",i2)') params%Jmin
        write(*,'("Max level:      Jmax=",i2)') params%Jmax
        write(*,'("Dimension:      dim=",i1)') params%dim
        write(*,*) "--------------------------------------------------------------------------------"
    endif

    !-----------------------------------------------------------------------------------------------------
    ! Allocate memory
    !-----------------------------------------------------------------------------------------------------
    call setup_wavelet(params, params%g)
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)
    call init_ghost_nodes(params)
    call reset_forest(params)

    tree_ID = 1
    time = 0.0_rk

    !-----------------------------------------------------------------------------------------------------
    ! Create grid based on type
    !-----------------------------------------------------------------------------------------------------
    select case(grid_type)
    case('equidistant')
        ! Get level from command line or use 1 as default
        call get_cmd_arg("--level", level, default=1_ik)

        if (params%rank == 0) then
            write(*,'("Creating equidistant grid at level J=",i2)') level
        endif

        call createEquidistantGrid_tree(params, hvy_block, level, .true., tree_ID)

    case('random')
        ! Get parameters for random grid
        call get_cmd_arg("--level", level, default=3_ik)
        call get_cmd_arg("--iterations", iterations, default=10_ik)

        if (params%rank == 0) then
            write(*,'("Creating random grid with initial level J=",i2," and ",i3," iterations")') level, iterations
        endif

        call createRandomGrid_tree(params, hvy_block, hvy_tmp, level_init=level, &
                                   verbosity=.true., iterations=iterations, tree_ID=tree_ID)
    end select

    !-----------------------------------------------------------------------------------------------------
    ! Save grid to file
    !-----------------------------------------------------------------------------------------------------
    if (params%rank == 0) then
        write(*,*) "--------------------------------------------------------------------------------"
        write(*,'("Saving grid with ",i7," active blocks to file: ", A)') lgt_n(tree_ID), trim(adjustl(fname_out))
    endif

    ! Initialize block data to zero (or some simple pattern)
    ! Users can modify this to set specific initial data
    hvy_block(:,:,:,:,:) = 0.0_rk

    call saveHDF5_tree(fname_out, time, 1_ik, 1, params, hvy_block, tree_ID)

    if (params%rank == 0) then
        write(*,*) "Grid saved successfully!"
        write(*,*) "--------------------------------------------------------------------------------"
    endif

end subroutine post_create_grid
