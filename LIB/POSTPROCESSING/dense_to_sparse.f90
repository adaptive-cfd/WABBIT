!> \brief postprocessing for sparsing data from a dense wabbit field
!-----------------------------------------------------------------------------------------------------

subroutine dense_to_sparse(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_initialization
    use module_helpers
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)       :: params
    character(len=cshort)                   :: indicator="threshold-state-vector", file_in, args
    character(len=cshort)                   :: tail_string
    real(kind=rk)                           :: time, eps=-1.0_rk
    integer(kind=ik)                        :: iteration
    character(len=cshort), allocatable      :: file_out(:)
    real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable              :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                        :: max_neighbors, level, k, tc_length, lgt_n_tmp
    integer(kind=ik), dimension(3)          :: Bs
    integer(hid_t)                          :: file_id
    character(len=2)                        :: level_in
    character(len=cshort)                   :: order
    real(kind=rk), dimension(3)             :: domain
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik)                        :: treecode_size, number_dense_blocks, i, l, dim


    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    call get_command_argument(2, file_in)
    if (file_in == '--help' .or. file_in == '--h') then
        if ( params%rank==0 ) then
            write(*,*) "--------------------------------------------------------------"
            write(*,*) "                DENSE TO SPARSE "
            write(*,*) "--------------------------------------------------------------"
            write(*,*) "postprocessing subroutine sparse a mesh with a given detail treshold"
            write(*,*) " "
            write(*,*) "Command:"
            write(*,*) "./wabbit-post --dense-to-sparse "
            write(*,*) "-------------------------------------------------------------"
            write(*,*) " Parameters: "
            write(*,*) "  --eps-normalized="
            write(*,*) "  --eps-norm="
            write(*,*) "  --eps="
            write(*,*) "  --indicator="
            write(*,*) "  --order="
            write(*,*) "  --files="
            write(*,*) "-------------------------------------------------------------"
            write(*,*)
        end if
        return
    end if

    !----------------------------------
    ! read parameters
    !----------------------------------
    call get_cmd_arg( "--eps-normalized", params%eps_normalized, default=.true. )
    call get_cmd_arg( "--eps-norm", params%eps_norm, default="L2" )
    call get_cmd_arg( "--eps", params%eps, default=-1.0_rk )
    call get_cmd_arg( "--indicator", indicator, default="threshold-state-vector" )
    ! --order and --wavelet are synonyms, but --wavelet overules (--order is deprecated)
    call get_cmd_arg( "--order", order, default="CDF40" )
    call get_cmd_arg( "--wavelet", params%wavelet, default=order )
    call get_cmd_arg( "--files", params%input_files )

    if (params%eps < 0.0_rk) then
        call abort(2303191,"You must specify the threshold value --eps")
    endif

    ! initialize wavelet transform
    ! also, set number of ghost nodes params%G to minimal value for this wavelet
    call setup_wavelet(params, params%g)

    params%coarsening_indicator = indicator
    params%forest_size = 1
    params%physics_type = "postprocessing"

    params%n_eqn = size(params%input_files)
    allocate(params%field_names(params%n_eqn))
    allocate(file_out(params%n_eqn))
    allocate(params%threshold_state_vector_component(params%n_eqn))
    params%threshold_state_vector_component = .true.

    !-------------------------------------------
    ! check and find common params in all h5-files
    !-------------------------------------------
    call read_attributes(params%input_files(1), lgt_n_tmp, time, iteration, params%domain_size, &
    params%Bs,params%Jmax, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    do i = 1, params%n_eqn
        file_in = params%input_files(i)
        call check_file_exists(trim(file_in))
        call read_attributes(file_in, lgt_n_tmp, time, iteration, domain, Bs, level, dim)

        params%Jmin = 1
        params%Jmax = max(params%Jmax, level) ! find the maximal level of all snapshot

        if (any(params%Bs .ne. Bs)) call abort( 203192, " Block size is not consistent ")
        if (params%dim .ne. dim) call abort(243191,"Dimensions do not agree!")
        if ( abs(sum(params%domain_size(1:dim) - domain(1:dim))) > 1e-14 ) call abort( 203195, "Domain size is not consistent ")

        ! Concatenate "sparse" with filename
        params%input_files(i) = trim(file_in)
        file_out(i) = trim(file_in)
    end do

    ! in postprocessing, it is important to be sure that the parameter struct is correctly filled:
    ! most variables are unfortunately not automatically set to reasonable values. In simulations,
    ! the ini files parser takes care of that (by the passed default arguments). But in postprocessing
    ! we do not read an ini file, so defaults may not be set.
    allocate(params%butcher_tableau(1,1))
    params%block_distribution="sfc_hilbert"

    ! read attributes from file. This is especially important for the number of
    ! blocks the file contains: this will be the number of active blocks right
    ! after reading.
    if (params%dim==3) then
        ! how many blocks do we need for the desired level?
        number_dense_blocks = 8_ik**level
        max_neighbors = 74
    else
        number_dense_blocks = 4_ik**level
        max_neighbors = 16
    end if

    if (params%rank==0) then
        write(*,'(80("-"))')
        write(*,*) "Wabbit dense-to-sparse."
        do i = 1, params%n_eqn
            write(*,'(A20,1x,A80)') "Reading file:", params%input_files(i)
            write(*,'(A20,1x,A80)') "Writing to file:", file_out(i)
        end do
        write(*,'(A20,1x,A80)') "Predictor used:", params%order_predictor
        write(*,'(A20,1x,A8)') "Wavelets used:", params%wavelet
        write(*,'(A20,1x,es9.3)') "eps:", params%eps
        write(*,'(A20,1x,A80)')"wavelet normalization:", params%eps_norm
        write(*,'(A20,1x,A80)')"indicator:", params%coarsening_indicator
        write(*,'(80("-"))')
    endif

    ! is lgt_n > number_dense_blocks (downsampling)? if true, allocate lgt_n blocks
    !> \todo change that for 3d case
    params%number_blocks = ceiling( 4.0*dble(max(lgt_n_tmp, number_dense_blocks)) / dble(params%number_procs) )

    if (params%rank==0) then
        write(*,'("Data dimension: ",i1,"D")') params%dim
        write(*,'("File contains Nb=",i6," blocks of size Bs=",i4," x ",i4," x ",i4)') lgt_n_tmp, Bs(1),Bs(2),Bs(3)
        write(*,'("Domain size is ",3(g12.4,1x))') domain
        write(*,'("Time=",g12.4," it=",i9)') time, iteration
        write(*,'("Length of treecodes in file=",i3," in memory=",i3)') level, params%Jmax
        write(*,'("NCPU=",i6)') params%number_procs
        write(*,'("File   Nb=",i6," blocks")') lgt_n_tmp
        write(*,'("Memory Nb=",i6)') params%number_blocks
        write(*,'("Dense  Nb=",i6)') number_dense_blocks
    endif

    !----------------------------------
    ! allocate data and reset grid
    !----------------------------------
    call allocate_forest(params, hvy_block, hvy_work=hvy_work, hvy_tmp=hvy_tmp)

    ! reset the grid: all blocks are inactive and empty
    call reset_tree( params, .true., tree_ID=tree_ID_flow )

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    !----------------------------------
    ! READ Grid and coarse if possible
    !----------------------------------
    params%adapt_tree = .true.
    params%adapt_inicond = .true.
    params%read_from_files = .true.
    params%force_maxlevel_dealiasing = .false.
    params%threshold_mask = .false.
    call setInitialCondition_tree( params, hvy_block, tree_ID_flow, params%adapt_inicond, time, iteration, hvy_tmp=hvy_tmp )

    !----------------------------------
    ! Write sparse files
    !----------------------------------
    do i = 1, params%n_eqn
        call saveHDF5_tree(file_out(i), time, iteration, i, params, hvy_block, tree_ID_flow)
    enddo

    call deallocate_forest(params, hvy_block, hvy_work, hvy_tmp=hvy_tmp)
end subroutine dense_to_sparse
