!##############################################################
!> This function reads a set of data into a specified tree.
!> The set of data can involve multiple files. Each file should
!> contain a quantity at the same snapshot (i.e. same iteration/time).
!> If no data has been read in before. This function will allocate
!> the heavy data for you.
subroutine read_field2tree(params, fnames, N_files, tree_id, tree_n, &
    lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)
    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params           !< params structure
    integer(kind=ik), intent(in)      :: N_files     !< number of fields/quantities to read
    character(len=*), intent(in)      :: fnames(N_files)  !< file names
    integer(kind=ik), intent(inout)   :: tree_n       !< number of active trees
    integer(kind=ik), intent(in)      :: tree_id     !< number of the tree
    integer(kind=ik), allocatable, intent(inout)   :: lgt_n(:),hvy_n(:)   !< number of active light and heavy blocks
    integer(kind=ik), ALLOCATABLE, intent(inout)   :: lgt_block(:, : )  !< light data array
    real(kind=rk), ALLOCATABLE, intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    integer(kind=ik), ALLOCATABLE, intent(inout)   :: hvy_neighbor(:,:)!< neighbor array
    integer(kind=ik), ALLOCATABLE, intent(inout)   :: lgt_active(:,:), hvy_active(:,:) !< active lists
    integer(kind=tsize), ALLOCATABLE, intent(inout):: lgt_sortednumlist(:,:,:)
    real(kind=rk), ALLOCATABLE, intent(inout)      :: hvy_tmp(:, :, :, :, :) ! used for saving, filtering, and helper qtys
    !-------------------------------- ---------------------------------
    integer(kind=ik)  :: iteration, dF, tc_length, dim, i, lgt_n_tmp, &
    rank, level, Bs(3), fsize, number_dense_blocks
    real(kind=rk)     :: time, domain(3)

    ! set MPI parameter
    rank  = params%rank
    fsize = params%forest_size

    if (tree_id <= tree_n) then
        ! the tree altready exists: to overwrite it, we first delete the existing one
        call delete_tree(params, lgt_block, lgt_active, lgt_n, tree_id)

        call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
        lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)
    endif

    ! From any file (here: 1st in list), we read the essential parameters of the tree:
    ! most importantly Bs, dim, length of treecodes. NOTE: all trees in the forest must
    ! have the same Bs.
    call read_attributes(fnames(1), lgt_n_tmp, time, iteration, domain, Bs, tc_length, dim)

    ! Check if one tree already exists in the forest:
    ! If it doesnt initialize some important parameters like the Block size Bs
    ! and the spatial dimension of the data.
    ! If it does chFeck if the parameters are suitable for the forest
    if (allocated(hvy_block)) then
        if ( any(params%Bs /= Bs) ) then
            ! NOTE: all trees in the forest must have the same Bs.
            call abort(100119, "ERROR: Trying to read a tree (grid) with Bs different from existing trees in forest!")
        endif
        if ( params%max_treelevel < tc_length) call abort(10119,"maximal treelevel incompatible")
    else
        params%max_treelevel = tc_length
        params%dim = dim
        params%n_eqn = N_files
        params%N_fields_saved = N_files
        params%Bs = Bs
        params%domain_size = domain

        ! we have to allocate grid if this routine is called for the first time
        call allocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
        hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp, hvy_n=hvy_n, lgt_n=lgt_n)

        hvy_neighbor = -1
        hvy_n = 0
        lgt_n = 0 ! reset number of active light blocks
        tree_n = 0 ! reset number of trees in forest
    endif



    ! read treecode from first input file
    call read_tree(fnames, N_files, params, lgt_n_tmp, lgt_block, hvy_block, hvy_tmp, tree_id)

    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)

    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:,tree_id),&
    lgt_n(tree_id), lgt_sortednumlist(:,:,tree_id), hvy_active(:,tree_id) , hvy_n(tree_id) )

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_id), hvy_n(tree_id) )

end subroutine read_field2tree
!##############################################################



!-------------------------------------------------------------------------------
! Read a tree (=grid) in the forest. A tree can consist of several datafields (think
! of the state vector) that all live on the same grid.
!
!
!-------------------------------------------------------------------------------
subroutine read_tree(fnames, N_files, params, lgt_n, lgt_block, hvy_block, hvy_tmp, tree_id_optional)
    implicit none

    !-------------------------------- ---------------------------------
    integer(kind=ik), intent(in)   :: N_files         !< number of files to read in this tree (=grid)
    character(len=*), intent(in)   :: fnames(N_files) !< list of files to be read
    type (type_params), intent(in) :: params          !< user defined parameter structure
    ! NOTE: we pass the routine lgt_n, which must be previously read from one of the files.
    ! This is done using the read_attributes subroutine. NOTE: This can possibly be merged here
    ! and mustnt neccessarily be done in the caller
    integer(kind=ik), intent(inout) :: lgt_n !< number of active blocks (heavy and light data)
    integer(kind=ik), intent(inout) :: lgt_block(:,:) !< light data array

    integer(kind=ik), optional, intent(in) :: tree_id_optional !< index of the tree you want to save the data in
    real(kind=rk), intent(inout) :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    real(kind=rk), intent(inout) :: hvy_tmp(:, :, :, :, :) !< heavy data array - block data
    !-------------------------------- ---------------------------------
    integer(kind=ik) :: k, N, rank, number_procs, ierr, treecode_size, tree_id, Bs(3), g, dF
    integer(kind=ik) :: ubounds(2), lbounds(2), blocks_per_rank_list(0:params%number_procs-1)
    integer(kind=ik) :: free_hvy_id, free_lgt_id, my_hvy_n
    integer(hsize_t) :: dims_treecode(2)
    integer(kind=ik), dimension(:,:), allocatable :: block_treecode
    integer(hid_t)        :: file_id
    character(len = 80) :: fname

    if (present(tree_id_optional)) then
        tree_id = tree_id_optional
    else
        tree_id = 1
    endif

    ! set MPI parameters
    rank         = params%rank
    number_procs = params%number_procs
    N            = params%number_blocks
    fname        = fnames(1) ! take the first of the N_files to read the lgt_data from
    Bs           = params%Bs
    g            = params%n_ghosts

    if (params%rank==0) write(*,*) "read_tree tries to read lgt_n=", lgt_n, "from file ", trim(adjustl(fname))

    ! open the file
    call check_file_exists(fname)
    call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)


    !-----------------------------------------------------------------------------
    ! Step 0: define how many blocks per CPU
    !-----------------------------------------------------------------------------
    call set_desired_num_blocks_per_rank(params, blocks_per_rank_list, lgt_n)

    ! number of active blocks on my process
    my_hvy_n = blocks_per_rank_list(rank)


    !-----------------------------------------------------------------------------
    ! Step 1: read light data (treecodes)
    !-----------------------------------------------------------------------------
    ! check what Jmax was saved in file (check length of treecode in file)
    call get_size_datafield(2, file_id, "block_treecode", dims_treecode)

    ! compare treecode lengths
    if (dims_treecode(1)>params%max_treelevel) then
        ! treecode in input file is greater than the new one, abort and output on screen
        ! NOTE this can be made working if not all levels in the file are actually used (e.g. level_max=17
        ! but active level=4). On the other hand, that appears to be rare.
        call abort(73947887, "ERROR: Treecode in file is longer than what is set in INI file.")
    end if

    allocate( block_treecode(1:dims_treecode(1), 1:my_hvy_n) )
    block_treecode = -1

    ! tell the hdf5 wrapper what part of the global [ n_active x max_treelevel + IDX_REFINE_STS]
    ! array we want to hold, so that all CPU can read from the same file simultaneously
    ! (note zero-based offset):
    lbounds = (/0, sum(blocks_per_rank_list(0:rank-1))/)
    ubounds = (/int(dims_treecode(1),4)-1, lbounds(2) + my_hvy_n - 1/)
    call read_dset_mpi_hdf5_2D(file_id, "block_treecode", lbounds, ubounds, block_treecode)

    ! close file and HDF5 library
    call close_file_hdf5(file_id)

    !-----------------------------------------------------------------------------
    ! Step 2: read heavy data from files into hvy_tmp
    !-----------------------------------------------------------------------------
    do dF = 1, N_files
        ! read data from file
        call check_file_exists(trim(fnames(dF)))
        call read_field(fnames(dF), dF, params, hvy_tmp, my_hvy_n )
    end do

    !-----------------------------------------------------------------------------
    ! Step 3: assemble both hvy and lgt data in tree
    !-----------------------------------------------------------------------------
    do k = 1, my_hvy_n
        call get_free_local_light_id( params, rank, lgt_block, free_lgt_id)

        call lgt_id_to_hvy_id( free_hvy_id, free_lgt_id, rank, N )
        ! copy treecode
        lgt_block(free_lgt_id, 1:dims_treecode(1)) = block_treecode(1:dims_treecode(1), k)
        ! set mesh level
        lgt_block(free_lgt_id, params%max_treelevel+IDX_MESH_LVL) = treecode_size(block_treecode(:,k), size(block_treecode,1))
        ! set refinement status
        lgt_block(free_lgt_id, params%max_treelevel+IDX_REFINE_STS) = 0
        ! set number of the tree
        lgt_block(free_lgt_id, params%max_treelevel+IDX_TREE_ID) = tree_id
        ! copy actual data
        do dF = 1, params%n_eqn
            hvy_block( :, :, :, df, free_hvy_id ) = hvy_tmp(:, :, :, df, k)
        end do
    end do

    if ( rank == 0 ) then
        write(*,'("Stored in Tree_id: ",i3)') tree_id
        write(*,'("Nblocks=",i6," (on all cpus)")') lgt_n
    end if

    ! synchronize light data. This is necessary as all CPUs above created their blocks locally.
    ! As they all pass the same do loops, the counter array blocks_per_rank_list does not have to
    ! be synced. However, the light data has to.
    call synchronize_lgt_data( params, lgt_block, refinement_status_only=.false. )

    deallocate(block_treecode)

    ! it is useful to print out the information on active levels in the file
    ! to get an idea how it looks like and if the desired dense level is larger
    ! or smaller
    if (params%rank==0) then
        write(*,'("In the file we just read, Jmin=",i3," Jmax=",i3)') min_active_level( lgt_block ), &
        max_active_level( lgt_block )
    endif

end subroutine read_tree
!##############################################


!##############################################################
!> save a specific field with specified tree_id
!> \TODO This function will replace write_field in the future
subroutine write_tree_field(fname, params, lgt_block, lgt_active, hvy_block, &
    lgt_n, hvy_n, hvy_active, dF, tree_id, time, iteration )

    implicit none

    !-----------------------------------------------------------------
    character(len=*), intent(in) :: fname  !< filename
    type (type_params), intent(in) :: params
    integer(kind=ik), intent(in)   :: lgt_block(:, :) !< ligh block data
    integer(kind=ik), intent(in)   :: lgt_active(:,:) !< list of active blocks for each tree
    integer(kind=ik), intent(in)   :: lgt_n(:), hvy_n(:) !< length of active lists
    integer(kind=ik), intent(in)   :: hvy_active(:,:) !< list of active hvy blocks
    real(kind=rk), intent(in)      :: hvy_block(:, :, :, :, :)
    !--------------------
    ! optional parameter
    !--------------------
    integer(kind=ik), optional, intent(in) :: tree_id !< id of the tree (default: 1)
    real(kind=rk), optional, intent(in)    :: time !< time loop parameters (default: 0.0)
    integer(kind=ik), optional, intent(in) :: iteration !< iteration of the solver (default: 1)
    integer(kind=ik), optional,intent(in)  :: dF !< datafield number (default: 1)
    !-----------------------------------------------------------------
    integer (kind=ik) ::tree_hvy_n, treeid, it, dataField
    real (kind=rk) :: t

    ! check the optional params
    if (present(dF)) then
        dataField = dF
    else
        dataField = 1
    endif
    if (present(tree_id)) then
        treeid = tree_id
    else
        treeid = 1
    endif
    if (present(time)) then
        t = time
    else
        t = 0.0_rk
    endif
    if (present(iteration)) then
        it = iteration
    else
        it = 0
    endif

    ! write the data
    call write_field(fname, t, it, dataField, params, lgt_block, hvy_block, &
    lgt_active(:,treeid), lgt_n(treeid), hvy_n(treeid), hvy_active(:, treeid))
end subroutine
!##############################################################
