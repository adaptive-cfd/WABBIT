!-------------------------------------------------------------------------------
! Read a tree (=grid) in the forest. A tree can consist of several datafields (think
! of the state vector) that all live on the same grid. They are stored in individual
! HDF5 files.
!
! Note: all files are supposed to share the same grid (components of a vector). No
! test is made if that is true: if the grids are different, garbage will be read.
!
!-------------------------------------------------------------------------------
subroutine readHDF5vct_tree(fnames, params, hvy_block, tree_ID, time, iteration, verbosity)
    implicit none

    character(len=*), intent(in)   :: fnames(1:)      !< list of files to be read (=number of vector components)
    type (type_params), intent(in) :: params          !< user defined parameter structure
    integer(kind=ik), intent(in)   :: tree_ID         !< index of the tree you want to save the data in
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :) !< heavy data array - data are stored herein

    logical, intent(in), optional  :: verbosity       !< if verbosity==True generates log output
    real(kind=rk), intent(out), optional   :: time
    integer(kind=ik), intent(out), optional   :: iteration

    real(kind=rk), allocatable :: hvy_buffer(:, :, :, :, :) !< buffer for reading from HDF5 file

    integer(kind=ik) :: N_files         !< number of files to read in this tree (=number of vector components)
    integer(kind=ik) :: NblocksFile
    integer(kind=ik) :: k, N, rank, number_procs, ierr, treecode_size, Bs(3), g, dF, i
    integer(kind=ik) :: ubounds(2), lbounds(2), blocks_per_rank_list(0:params%number_procs-1)
    integer(kind=ik), dimension(4) :: ubounds3D, lbounds3D       ! offset variables
    integer(kind=ik), dimension(3) :: ubounds2D, lbounds2D
    integer(kind=ik) :: free_hvy_id, free_lgt_id, my_hvy_n, version(1), datarank, Bs_file(1:3)=0
    integer(hsize_t) :: dims_treecode(2)
    integer(kind=ik), dimension(:,:), allocatable :: block_treecode
    integer(hid_t)        :: file_id
    integer(kind=hsize_t) :: size_field(1:4)
    logical :: verbose = .true.

    if (present(verbosity)) verbose=verbosity


    rank         = params%rank
    number_procs = params%number_procs
    N            = params%number_blocks
    Bs           = params%Bs
    g            = params%n_ghosts
    NblocksFile  = getNumberBlocksH5File(fnames(1))
    N_files      = size(fnames)

    if (tree_ID <= tree_N) then
        ! the tree already exists: to overwrite it, we first delete the existing one
        call delete_tree(params, tree_ID)
        call createActiveSortedLists_forest(params)
    endif

    !-----------------------------------------------------------------------------
    ! Header information (useful for debugging if something goes wrong)
    !-----------------------------------------------------------------------------
    call open_file_hdf5( trim(adjustl(fnames(1))), file_id, .false.)
    ! datarank is rank of block array (3D rank=4, 2D rank=3)
    call get_rank_datafield(file_id, "blocks", datarank)
    ! the size of the array is Bs x Bs x Bs x Nb
    call get_size_datafield(datarank, file_id, "blocks", size_field(1:datarank))
    ! copy first 3 entries to Bs
    Bs_file(1:datarank-1) = size_field(1:datarank-1)
    ! Files created using newGhostNodes branch (after 08 04 2020) contain a version number
    call read_attribute( file_id, "blocks", "version", version)
    ! read time stamp (if desired)
    if (present(time)) then
        call read_attribute(file_id, "blocks", "time", time)
    endif
    ! read iteration (if desired)
    if (present(iteration)) then
        call read_attribute(file_id, "blocks", "iteration", iteration)
    endif
    ! done. (We'll open it again!)
    call close_file_hdf5(file_id)


    if (params%rank==0 .and. verbose) then
        do i = 1, N_files
            write(*,'("READING: Filename #",i2," = ",A)') i, trim(adjustl(fnames(i)))
        enddo
        write(*,'("READING: VERSION of file(s)= ",i9)') version(1)
        write(*,'("READING: Expected Nblocks  = ",i9," (on all cpus)")') NblocksFile
        write(*,'("READING: datarank          = ",i1," ---> ",i1,"D data")') datarank, datarank-1
        write(*,'("READING: Bs_file           = ",3(i3,1x),"Bs_memory = ",3(i3,1x))') Bs_file, params%Bs
        write(*,'("READING: Stored in tree_ID = ",i3)') tree_ID
    endif

    ! check if at least all files have same number of blocks (which still does not mean
    ! that the grids are identical)
    do i = 1, N_files
        if ( NblocksFile /= getNumberBlocksH5File(fnames(i)) ) then
            call abort(2408221, "You are trying to read files with different number of blocks &
            &into one vector. You probably have the wrong files.")
        endif
    enddo

    ! will the file fit in the allocated memory?
    if (NblocksFile > size(hvy_block,5)) then
        call abort(2408223, "The file we read is larger than the allocated array!")
    endif


    !----------------------------------------------------------------------------------
    ! version check
    !----------------------------------------------------------------------------------
    ! Files created using the newGhostNodes branch (after 08 April 2020) contain a version number.
    ! if the number is not found, version=0.
    if (version(1) == 20200408 .and. rank==0) then
        write(*,*) "--------------------------------------------------------------------"
        write(*,*) "-----WARNING----------WARNING----------WARNING----------WARNING-----"
        write(*,*) "--------------------------------------------------------------------"
        write(*,*) "The file we are trying to read is generated with an intermediate version"
        write(*,*) "of wabbit (after 08 April 2020). In this the file, the grid"
        write(*,*) "definition does not include a redundant point, i.e., a block is defined"
        write(*,*) "with spacing dx = L*2^-J / Bs. This definition was a dead-end, as it lead"
        write(*,*) "to instabilities and other problems. Current versions of WABBIT include a redundant point again."
        write(*,*) ""
        write(*,*) "The newGhostNodes branch still STORED the redundant point for visualization."
        write(*,*) "This was simply the first ghost node. If Bs was odd, this lead to an even number"
        write(*,*) "of points, and this cannot be read with present code versions."
        write(*,*) ""
        write(*,*) "A workaround must be done in preprocessing: upsampling to equidistant resolution and"
        write(*,*) "re-gridding is required. I am truely sorry for this."
        write(*,*) ""
        write(*,*) "If bs was even, the resulting data size is odd, and the file can be read. note however"
        write(*,*) "that the resolution changes slightly, and results cannot be perfectly identical to what would"
        write(*,*) "have been obtained with the newGhostNodes branch"
        write(*,*) "--------------------------------------------------------------------"
        write(*,*) "-----WARNING----------WARNING----------WARNING----------WARNING-----"
        write(*,*) "--------------------------------------------------------------------"
    endif



    !-----------------------------------------------------------------------------
    ! Step 0: define how many blocks per CPU
    !-----------------------------------------------------------------------------
    ! this list contains (on each mpirank) the number of blocks for each mpirank. note
    ! zero indexing as required by MPI
    ! set list to the average value
    blocks_per_rank_list(:) = NblocksFile / number_procs

    ! as this does not necessarily work out, distribute remaining blocks on the first CPUs
    if (mod(NblocksFile, number_procs) > 0) then
        blocks_per_rank_list(0:mod(NblocksFile, number_procs)-1) = &
        blocks_per_rank_list(0:mod(NblocksFile, number_procs)-1) + 1
    end if

    ! number of active blocks on my process (= the amount of blocks to be read)
    my_hvy_n = blocks_per_rank_list(rank)


    !-----------------------------------------------------------------------------
    ! Step 1: read light data (treecodes)
    !-----------------------------------------------------------------------------
    ! open the file
    call check_file_exists(fnames(1))
    call open_file_hdf5(fnames(1), file_id, .false.)

    ! check what Jmax was saved in file (check length of treecode in file)
    call get_size_datafield(2, file_id, "block_treecode", dims_treecode)

    ! compare treecode lengths
    if (dims_treecode(1) > params%max_treelevel) then
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

    ! actual reading of treecodes (= the description of the tree)
    call read_dset_mpi_hdf5_2D(file_id, "block_treecode", lbounds, ubounds, block_treecode)

    ! close file and HDF5 library
    call close_file_hdf5(file_id)


    !-----------------------------------------------------------------------------
    ! Step 2: read heavy data from files into the buffer
    !-----------------------------------------------------------------------------

    if ( params%dim == 3 ) then
        allocate( hvy_buffer(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g, 1:N_files, 1:my_hvy_n) )
        ! tell the hdf5 wrapper what part of the global [Bsx x Bsy x Bsz x hvy_n]
        ! array we want to hold, so that all CPU can read from the same file simultaneously
        ! (note zero-based offset):
        lbounds3D = (/0,0,0,sum(blocks_per_rank_list(0:rank-1))/)
        ubounds3D = (/Bs(1)-1,Bs(2)-1,Bs(3)-1,lbounds3D(4)+my_hvy_n-1/)

        ! 3D data case
        ! actual reading of file (into hvy_buffer! not hvy_work)
        do dF = 1, N_files
            call open_file_hdf5( trim(adjustl(fnames(dF))), file_id, .false.)

            if (rank==0 .and. verbose) then
                write(*,'("READING: Reading datafield ",i2," from file ",A)') dF, trim(adjustl(fnames(dF)))
            end if

            call read_dset_mpi_hdf5_4D(file_id, "blocks", lbounds3D, ubounds3D, &
            hvy_buffer(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, 1:my_hvy_n))

            ! close file and HDF5 library
            call close_file_hdf5(file_id)
        enddo
    else
        allocate( hvy_buffer(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1, 1:N_files, 1:my_hvy_n) )
        ! tell the hdf5 wrapper what part of the global [Bsx x Bsy x 1 x hvy_n]
        ! array we want to hold, so that all CPU can read from the same file simultaneously
        ! (note zero-based offset):
        lbounds2D = (/0,0,sum(blocks_per_rank_list(0:rank-1))/)
        ubounds2D = (/Bs(1)-1,Bs(2)-1,lbounds2D(3)+my_hvy_n-1/)

        ! 2D data case
        ! actual reading of file (into hvy_buffer! not hvy_work)
        do dF = 1, N_files
            call open_file_hdf5( trim(adjustl(fnames(dF))), file_id, .false.)

            if (rank==0 .and. verbose) then
                write(*,'("READING: Reading datafield ",i2," from file ",A)') dF, trim(adjustl(fnames(dF)))
            end if

            call read_dset_mpi_hdf5_3D(file_id, "blocks", lbounds2D, ubounds2D, &
            hvy_buffer(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, dF, 1:my_hvy_n))

            ! close file and HDF5 library
            call close_file_hdf5(file_id)
        enddo
    endif


    !-----------------------------------------------------------------------------
    ! Step 3: assemble both hvy and lgt data in tree
    !-----------------------------------------------------------------------------
    do k = 1, my_hvy_n
        call get_free_local_light_id( params, rank, free_lgt_id)

        call lgt2hvy( free_hvy_id, free_lgt_id, rank, N )
        ! copy treecode
        lgt_block(free_lgt_id, 1:dims_treecode(1)) = block_treecode(1:dims_treecode(1), k)
        ! set mesh level
        lgt_block(free_lgt_id, params%max_treelevel+IDX_MESH_LVL) = treecode_size(block_treecode(:,k), size(block_treecode,1))
        ! set refinement status
        lgt_block(free_lgt_id, params%max_treelevel+IDX_REFINE_STS) = 0
        ! set number of the tree
        lgt_block(free_lgt_id, params%max_treelevel+IDX_TREE_ID) = tree_id
        ! copy actual data
        do dF = 1, N_files
            hvy_block( :, :, :, dF, free_hvy_id ) = hvy_buffer(:, :, :, dF, k)
        end do
    end do


    ! synchronize light data. This is necessary as all CPUs above created their blocks locally.
    ! As they all pass the same do loops, the counter array blocks_per_rank_list does not have to
    ! be synced. However, the light data has to.
    call synchronize_lgt_data( params, refinement_status_only=.false. )

    ! it is good practice that this routine returns a working forest, i.e., all meta
    ! data is updated.
    call createActiveSortedLists_forest(params)
    call updateNeighbors_tree(params, tree_ID)

    call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

    deallocate(hvy_buffer)
    deallocate(block_treecode)

    ! it is useful to print out the information on active levels in the file
    ! to get an idea how it looks like and if the desired dense level is larger
    ! or smaller
    if (params%rank==0 .and. verbose ) then
        write(*,'("READING: In the file(s) we just read, Jmin=",i3," Jmax=",i3)') minActiveLevel_tree(tree_ID), &
        maxActiveLevel_tree(tree_ID)
    endif

end subroutine readHDF5vct_tree
