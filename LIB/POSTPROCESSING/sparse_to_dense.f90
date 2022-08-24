subroutine sparse_to_dense(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_mpi
    use module_globals
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)  :: file_in
    character(len=cshort)  :: file_out
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvy_id

    integer(kind=ik)                        :: max_neighbors, level, k, tc_length
    integer(kind=ik), dimension(3)          :: Bs
    integer(hid_t)                          :: file_id
    character(len=2)                        :: level_in, order
    character(len=cshort)                       :: operator
    real(kind=rk), dimension(3)             :: domain
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik)                        :: number_dense_blocks

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    call get_command_argument(2, file_in)
    call get_command_argument(3, file_out)

    if (file_in == '--help' .or. file_in == '--h') then
        if ( params%rank==0 ) then
            write(*,*) "--------------------------------------------------------------"
            write(*,*) "                SPARSE to DENSE "
            write(*,*) "--------------------------------------------------------------"
            write(*,*) "postprocessing subroutine to refine/coarse mesh to a uniform"
            write(*,*) "grid (up and downsampling ensured)."
            write(*,*) "Command:"
            write(*,*) "./wabbit-post --sparse-to-dense source.h5 target.h5 [target_treelevel order-predictor(2 or 4)]"
            write(*,*) "-------------------------------------------------------------"
            write(*,*) "Optional Inputs: "
            write(*,*) "  1. target_treelevel = number specifying the desired treelevel"
            write(*,*) "  (default is the max treelevel of the source file) "
            write(*,*) "  2. order-predictor = consistency order or the predictor stencil"
            write(*,*) "  (default is preditor order 4) "
            write(*,*)
            write(*,*)
        end if
        return
    end if


    ! get values from command line (filename and level for interpolation)
    call check_file_exists(trim(file_in))
    call read_attributes(file_in, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, &
    params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    if (len_trim(file_out)==0) then
        call abort(0909191,"You must specify a name for the target! See --sparse-to-dense --help")
    endif

    ! check if optional arguments are specified
    if (command_argument_count()<4) then
        ! set defaults
        order = "4"
        level = tc_length
    else
        call get_command_argument(4, level_in)
        read(level_in,*) level
        call get_command_argument(5, order)
    endif

    write(*,*) "order=", order

    call get_cmd_arg( "--operator", operator, default="sparse-to-dense")

    if (order == "4") then
        params%wavelet_transform_type = 'harten-multiresolution'
        params%order_predictor = "multiresolution_4th"
        params%n_ghosts = 4_ik
        params%wavelet = "CDF4,0"
    elseif (order == "44") then
        params%wavelet_transform_type = "biorthogonal"
        params%order_predictor = "multiresolution_4th"
        params%n_ghosts = 6_ik
        params%wavelet = "CDF4,4"
    elseif (order == "2") then
        params%wavelet_transform_type = "harten-multiresolution"
        params%order_predictor = "multiresolution_2nd"
        params%n_ghosts = 2_ik
        params%wavelet = "CDF2,0"
    else
        call abort(392,"ERROR: chosen predictor order invalid or not (yet) implemented. choose between 4 (multiresolution_4th) and 2 (multiresolution_2nd)")
    end if

    ! in postprocessing, it is important to be sure that the parameter struct is correctly filled:
    ! most variables are unfortunately not automatically set to reasonable values. In simulations,
    ! the ini files parser takes care of that (by the passed default arguments). But in postprocessing
    ! we do not read an ini file, so defaults may not be set.
    allocate(params%butcher_tableau(1,1))
    ! we read only one datafield in this routine
    params%n_eqn  = 1
    params%block_distribution="sfc_hilbert"

    if (params%dim==3) then
        ! how many blocks do we need for the desired level?
        number_dense_blocks = 8_ik**level
        max_neighbors = 74
    else
        number_dense_blocks = 4_ik**level
        max_neighbors = 12
    end if

    if (params%rank==0) then
        write(*,'(80("-"))')
        write(*,*) "Wabbit sparse-to-dense. Will read a wabbit field and return a"
        write(*,*) "full grid with all blocks at the chosen level."
        write(*,'(A20,1x,A80)') "Reading file:", file_in
        write(*,'(A20,1x,A80)') "Writing to file:", file_out
        write(*,'(A20,1x,A80)') "Predictor used:", params%order_predictor
        write(*,'(A20,1x,i3," => ",i9," Blocks")') "Target level:", level, number_dense_blocks

        write(*,'(A40,1x,A40)') "params%wavelet_transform_type", params%wavelet_transform_type
        write(*,'(A40,1x,A40)') "params%order_predictor", params%order_predictor
        write(*,'(A40,1x,A40)') "params%wavelet", params%wavelet
        write(*,'(A40,1x,i2)') "params%n_ghosts", params%n_ghosts
        write(*,'(80("-"))')
    endif

    ! set max_treelevel for allocation of hvy_block
    params%max_treelevel = max(level, tc_length)
    params%min_treelevel = level
    params%Bs = Bs
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)

    ! is lgt_n > number_dense_blocks (downsampling)? if true, allocate lgt_n blocks
    !> \todo change that for 3d case
    params%number_blocks = ceiling( 4.0*dble(max(lgt_n(tree_ID), number_dense_blocks)) / dble(params%number_procs) )

    if (params%rank==0) then
        write(*,'("Data dimension: ",i1,"D")') params%dim
        write(*,'("File contains Nb=",i6," blocks of size Bs=",i4," x ",i4," x ",i4)') lgt_n(tree_ID), Bs(1),Bs(2),Bs(3)
        write(*,'("Domain size is ",3(g12.4,1x))') domain
        write(*,'("Time=",g12.4," it=",i9)') time, iteration
        write(*,'("Length of treecodes in file=",i3," in memory=",i3)') tc_length, params%max_treelevel
        write(*,'("   NCPU=",i6)') params%number_procs
        write(*,'("File   Nb=",i6," blocks")') lgt_n(tree_ID)
        write(*,'("Memory Nb=",i6)') params%number_blocks
        write(*,'("Dense  Nb=",i6)') number_dense_blocks
    endif

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)

    ! read input data
    call readHDF5vct_tree( (/file_in/), params, hvy_block, tree_ID)

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call updateMetadata_tree(params, tree_ID)

    ! balance the load
    call balanceLoad_tree(params, hvy_block, tree_ID)

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

    if (operator=="sparse-to-dense") then
        call toEquidistant_tree(params, hvy_block, hvy_tmp, tree_ID, target_level=level)

    elseif (operator=="refine-coarsen") then
        write(*,*) "starting at", lgt_n(tree_ID)

        call refine_tree( params, hvy_block, "everywhere", tree_ID )

        write(*,*) "refined to", lgt_n(tree_ID)

        params%threshold_mask = .false.
        params%coarsening_indicator = "threshold-state-vector"
        params%physics_type = "ConvDiff-new"
        params%ghost_nodes_redundant_point_coarseWins = .false.
        params%iter_ghosts = .false.
        params%eps_normalized = .false.
        params%force_maxlevel_dealiasing = .false.
        params%min_treelevel = 1

        call adapt_tree( time, params, hvy_block, tree_ID_flow, "everywhere", hvy_tmp )

        write(*,*) "coarsened to", lgt_n(tree_ID)
    endif

    call saveHDF5_tree(file_out, time, iteration, 1, params, hvy_block, tree_ID)

    if (params%rank==0 ) then
        write(*,'("Wrote data of input-file: ",A," now on uniform grid (level",i3, ") to file: ",A)') &
        trim(adjustl(file_in)), level, trim(adjustl(file_out))
        write(*,'("Minlevel:", i3," Maxlevel:", i3, " (should be identical now)")') &
        minActiveLevel_tree( tree_ID ),&
        maxActiveLevel_tree( tree_ID )
    end if

    call deallocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)
end subroutine sparse_to_dense
