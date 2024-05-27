
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

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvyID

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
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Wabbit postprocessing: evaluate wavelet thresholding"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " evaluates wavelet thresholding for an existing field"
            write(*,*) ""
            write(*,*) " ./wabbit-post --evaluate-wavelet-thresholding ux_000.h5 --eps=1.0e-3 --wavelet=CDF40"
            write(*,*) ""
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " --eps=1.0e-3"
            write(*,*) " --eps-norm=Linfty"
            write(*,*) " --wavelet=CDF44"
            write(*,*) " --eps-normalized=1"
            write(*,*) "-----------------------------------------------------------"
        end if
        return
    endif

    ! because we print the details, theyre heavy, so some parallel issuses, easier on monoproc
    if (params%number_procs>1) call abort(9991911, "serial routine, use NCPU=1")

    call get_cmd_arg( "--eps-normalized", params%eps_normalized, default=.true. )
    call get_cmd_arg( "--eps-norm", params%eps_norm, default="L2" )
    call get_cmd_arg( "--eps", params%eps, default=-1.0_rk )
    call get_cmd_arg( "--indicator", indicator, default="threshold-state-vector" )
    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF40" )

    ! initialize wavelet transform
    ! also, set number of ghost nodes params%G to minimal value for this wavelet
    call setup_wavelet(params, params%g)

    call setup_wavelet(params)

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

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, neqn_hvy_tmp=nwork)

    ! read input data
    call readHDF5vct_tree( (/fname/), params, hvy_block, tree_ID)

    call sync_ghosts_all( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID))


    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)

    do level = Jmax_active, Jmin_active, -1
        write(*,*) level
        call coarseningIndicator_level( time, params, level, hvy_block, hvy_tmp, tree_ID, params%coarsening_indicator, iteration, ignore_maxlevel=.true., input_is_WD=.false.)
    enddo



    if (params%rank == 0) then
        do k = 1, lgt_n(tree_ID)
            lgtID = lgt_active(k, tree_ID)
            call lgt2hvy(hvyID, lgtID, params%rank, params%number_blocks)
            write(*,'(i6,1x,i2,1x,i2,4x)') lgtID,  lgt_block( lgtID, IDX_MESH_LVL), lgt_block( lgtID, IDX_REFINE_STS )
        enddo
    endif

    call saveHDF5_tree(fname, time, iteration, 1, params, hvy_block, tree_ID )
end subroutine
