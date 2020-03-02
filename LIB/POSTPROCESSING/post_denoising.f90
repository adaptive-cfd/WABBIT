subroutine post_denoising(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_mpi
    use module_initialization
    use module_helpers

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=80)      :: indicator="threshold-state-vector", file_in, file_out
    real(kind=rk)          :: time, eps=-1.0_rk
    integer(kind=ik), allocatable           :: lgt_block(:, :)
    real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :)
    integer(kind=ik), allocatable           :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable           :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable        :: lgt_sortednumlist(:,:)
    integer(kind=ik)                        :: max_neighbors, level, k, tc_length
    integer(hid_t)                          :: file_id
    character(len=80)                       :: order
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik)                        :: treecode_size, number_dense_blocks, i, l, dim, hvy_n, lgt_n, iteration
    !-----------------------------------------------------------------------------------------------------

    call get_command_argument(2, file_in)
    if (file_in == '--help' .or. file_in == '--h') then
        if ( params%rank==0 ) then
            write(*,*) "--------------------------------------------------------------"
            write(*,*) "                DENOISING "
            write(*,*) "--------------------------------------------------------------"
            write(*,*) " --eps-norm=L2"
            write(*,*) " --wavelet=CDF44"
            write(*,*) " --input=data_00000.h5 "
            write(*,*) " --output=data_00000.h5 "
            write(*,*) "-------------------------------------------------------------"
            write(*,*)
        end if
        return
    end if

    !----------------------------------
    ! read parameters
    !----------------------------------
    call get_cmd_arg_str( "--eps-norm", params%eps_norm, default="L2" )
    call get_cmd_arg_str( "--input", file_in, default="" )
    call get_cmd_arg_str( "--output", file_out, default="" )
    call get_cmd_arg_str( "--indicator", indicator, default="threshold-state-vector" )
    call get_cmd_arg_str( "--wavelet", order, default="CDF44" )


    ! Check parameters for correct inputs:
    if (order == "CDF20") then
        params%harten_multiresolution = .true.
        params%order_predictor = "multiresolution_2nd"
        params%n_ghosts = 2_ik

    elseif (order == "CDF40") then
        params%harten_multiresolution = .true.
        params%order_predictor = "multiresolution_4th"
        params%n_ghosts = 4_ik

    elseif (order == "CDF44") then
        params%harten_multiresolution = .false.
        params%wavelet_transform_type = 'biorthogonal'
        params%order_predictor = "multiresolution_4th"
        params%wavelet='CDF4,4'
        params%n_ghosts = 6_ik
    else
        call abort(20030202, "The --order parameter is not correctly set [CDF40, CDF20, CDF44]")

    end if


    params%coarsening_indicator = indicator
    params%forest_size = 1
    params%n_eqn = 1
    allocate(params%field_names(params%n_eqn))
    allocate(params%threshold_state_vector_component(params%n_eqn))
    params%eps_normalized = .true.
    params%threshold_state_vector_component = .true.
    params%block_distribution = "sfc_hilbert"
    allocate(params%butcher_tableau(1,1))

    !-------------------------------------------
    ! check and find common params in all h5-files
    !-------------------------------------------
    call read_attributes(file_in, lgt_n, time, iteration, params%domain_size, params%Bs, params%max_treelevel, params%dim)

    ! no refinement is made in this postprocessing tool; we therefore allocate about
    ! the number of blocks in the file (and not much more than that)
    params%number_blocks = ceiling(  real(lgt_n)/real(params%number_procs) )

    !----------------------------------
    ! allocate data and reset grid
    !----------------------------------
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist)

    ! reset the grid: all blocks are inactive and empty
    call reset_tree( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true., tree_ID=1 )

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    call read_mesh(file_in, params, lgt_n, hvy_n, lgt_block)
    call read_field(file_in, 1, params, hvy_block, hvy_n)

    call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_ID=1)
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )


    call write_field(file_out, time, iteration, 1, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n, hvy_active )

end subroutine post_denoising
