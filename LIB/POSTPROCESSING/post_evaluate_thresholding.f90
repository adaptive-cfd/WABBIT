
subroutine post_evaluate_thresholding(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: fname, indicator, operator, order
    real(kind=rk)                      :: time
    integer(kind=ik)                   :: iteration, k, lgtID, tc_length, g
    integer(kind=ik), dimension(3)     :: Bs

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :), hvy_mask(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvyID
    logical                            :: bool_dump

    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    integer(kind=ik)                   :: nwork, level, Jmin_active, Jmax_active

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(1, operator)
    call get_command_argument(2, fname)

    ! does the user need help?
    if (operator=='--help' .or. operator=='--h') then
        if (params%rank==0) then
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Wabbit postprocessing: evaluate wavelet thresholding"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " evaluates wavelet thresholding for an existing field."
            write(*, '(A)') " This tool does, however, not know how to work with masks or multiple fields yet."
            write(*, '(A)') ""
            write(*, '(A)') " ./wabbit-post --evaluate-wavelet-thresholding ux_000.h5 --eps=1.0e-3 --wavelet=CDF40"
            write(*, '(A)') ""
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " --eps=1.0e-3"
            write(*, '(A)') " --eps-norm=Linfty"
            write(*, '(A)') " --eps-normalized=1"
            write(*, '(A)') " --wavelet=CDF44"
            write(*, '(A)') " --indicator=threshold-state-vector"
            write(*, '(A)') " --dump=0                              - print everything to command line"
            write(*, '(A)') "-----------------------------------------------------------"
        end if
        return
    endif

    call get_cmd_arg( "--eps-normalized", params%eps_normalized, default=.true. )
    call get_cmd_arg( "--eps-norm", params%eps_norm, default="L2" )
    call get_cmd_arg( "--eps", params%eps, default=-1.0_rk )
    call get_cmd_arg( "--indicator", indicator, default="threshold-state-vector" )
    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF44" )
    call get_cmd_arg( "--dump", bool_dump, default=.false. )

    ! initialize wavelet transform
    ! also, set number of ghost nodes params%G to minimal value for this wavelet
    call setup_wavelet(params, params%g, params%g_rhs)

    if (params%eps < 0.0_rk) then
        call abort(2303191,"You must specify the threshold value --eps")
    endif

    params%coarsening_indicator = "threshold-state-vector"
    params%forest_size = 1
    params%n_eqn = 1

    call check_file_exists(trim(fname))


    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(fname, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)


    params%Jmax = tc_length
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs
    params%threshold_mask = .false.
    params%force_maxlevel_dealiasing = .false.
    allocate(params%butcher_tableau(1,1))

    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component(1) = "0"

    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component = .true.

    Bs = params%Bs
    g  = params%g


    ! no refinement is made in this postprocessing tool; we therefore allocate about
    ! the number of blocks in the file (and not much more than that)
    params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) )

    nwork = 1

    ! ! have the pysics module read their own parameters - used for mask
    ! call init_physics_modules( params, filename, params%N_mask_components )

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_mask=hvy_mask, neqn_hvy_tmp=nwork)

    ! read input data
    call readHDF5vct_tree( (/fname/), params, hvy_block, tree_ID)

    call sync_ghosts_tree( params, hvy_block, tree_ID)


    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)

    call coarseningIndicator_tree( time, params, hvy_block, hvy_tmp, tree_ID, params%coarsening_indicator, &
        ignore_maxlevel=.true., input_is_WD=.false.)

    ! print to command line but only if requested
    if (params%rank == 0 .and. bool_dump) then
        do k = 1, lgt_n(tree_ID)
            lgtID = lgt_active(k, tree_ID)
            call lgt2hvy(hvyID, lgtID, params%rank, params%number_blocks)
            write(*,'(i6,1x,i2,1x,i2,4x)') lgtID,  lgt_block( lgtID, IDX_MESH_LVL), lgt_block( lgtID, IDX_REFINE_STS )
        enddo
    endif


    call saveHDF5_tree(fname, time, iteration, 1, params, hvy_block, tree_ID )
end subroutine
