!> \brief read mesh properties of a field saved in a hdf5-file
!! input:    - parameter array
!!           - name of the file we want to read from
!! output:   - light block array
!!           - number of active blocks (light and heavy)
! ********************************************************************************************

subroutine read_mesh(fname, params, lgt_n, hvy_n, lgt_block, tree_id_optional )
    implicit none

    character(len=*), intent(in)      :: fname                                  !> file name
    type (type_params), intent(in)    :: params                                 !> user defined parameter structure
    integer(kind=ik), intent(inout)   :: hvy_n, lgt_n                           !> number of active blocks (heavy and light data)
    integer(kind=ik), intent(inout)   :: lgt_block(:,:)                         !> light data array
    integer(kind=ik), optional, intent(in)   :: tree_id_optional                !> index of the tree you want to save the field in
    integer(kind=ik), dimension(:,:), allocatable :: block_treecode             ! treecode array
    integer(hid_t)        :: file_id
    integer(kind=ik)      :: rank, number_procs
    integer(kind=ik)      :: g, datarank                                        ! grid parameter
    integer(kind=ik), dimension(3)    :: Bs
    ! offset variables
    integer(kind=ik)      :: ubounds(2), lbounds(2), Bs_file(1:3)
    integer(kind=ik)      :: blocks_per_rank_list(0:params%number_procs-1)
    integer(kind=ik)      :: lgt_id, k, tree_id
    integer(kind=ik)      :: ierr
    integer(kind=ik)      :: treecode_size, tmp(1), version(1)
    integer(hsize_t)      :: dims_treecode(2)
    integer(kind=hsize_t) :: size_field(1:4)

    if (present(tree_id_optional)) then
        tree_id=tree_id_optional
    else
        tree_id=1
    endif

    rank         = params%rank
    number_procs = params%number_procs
    Bs   = params%Bs
    g    = params%n_ghosts
    Bs_file      = 1


    call check_file_exists(fname)
    ! open the file
    call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)

    ! how much blocks are we going to read?
    call read_attribute(file_id, "blocks", "total_number_blocks", tmp)
    lgt_n = tmp(1)
    call get_rank_datafield(file_id, "blocks", datarank)
    call get_size_datafield(datarank, file_id, "blocks", size_field(1:datarank))
    ! Files created using newGhostNodes branch (after 08 04 2020) contain a version number
    call read_attribute( file_id, "blocks", "version", version)
    Bs_file(1:datarank-1) = size_field(1:datarank-1)

    !---------------------------------------------------------------------------
    ! Header, useful for finding errors
    !---------------------------------------------------------------------------
    if ( rank == 0 ) then
        write(*,'(80("~"))')
        write(*,'(A)') "READING: initializing grid from file... (NOT reading the actual data!)"
        write(*,'(80("~"))')
        write(*,'("Filename         = ",A)') trim(adjustl(fname))
        write(*,'("VERSION of file  = ",i9)') version(1)
        write(*,'("Expected Nblocks = ",i9," (on all cpus)")') lgt_n
        write(*,'("datarank         = ",i1," ---> ",i1,"D data")') datarank, datarank-1
        write(*,'("Bs_file          = ",3(i3,1x))') Bs_file
        write(*,'("Bs_memory        = ",3(i3,1x))') params%Bs
    end if

    !----------------------------------------------------------------------------------
    ! version check
    !----------------------------------------------------------------------------------
    ! Files created using the newGhostNodes branch (after 08 April 2020) contain a version number.
    ! if the number is not found, version=0.
    call read_attribute( file_id, "blocks", "version", version)

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


    ! If the file has the wrong dimensions, we cannot read it
    if (maxval(Bs(1:datarank-1)-Bs_file(1:datarank-1)) > 0) then
        call abort(20030218, "ERROR:read_mesh:We try to load a file which has the wrong block size!")
    endif

    ! do we have enough memory?
    if (lgt_n > size(lgt_block,1)) then
        call abort(20030219, "ERROR:read_mesh:We try to load a file which will not fit in the memory!")
    endif

    if (lgt_n <= 0) then
        call abort(20030219, "We try to read an empty data file! Something is wrong, Nb<=0")
    end if

    ! Nblocks per CPU
    ! this list contains (on each mpirank) the number of blocks for each mpirank. note
    ! zero indexing as required by MPI
    ! set list to the average value
    blocks_per_rank_list(:) = lgt_n / number_procs

    ! as this does not necessarily work out, distribute remaining blocks on the first CPUs
    if (mod(lgt_n, number_procs) > 0) then
        blocks_per_rank_list(0:mod(lgt_n, number_procs)-1) = &
            blocks_per_rank_list(0:mod(lgt_n, number_procs)-1) + 1
    end if

    ! some error control -> did we loose blocks? should never happen.
    if ( sum(blocks_per_rank_list) /= lgt_n) then
        call abort(1028,"ERROR: while reading from file, we seem to have gained/lost some blocks during distribution...")
    end if

    ! number of active blocks on my process
    hvy_n = blocks_per_rank_list(rank)

    ! check what max_level was saved in file (check dimensions of treecode in file)
    call get_size_datafield(2, file_id, "block_treecode", dims_treecode)

    ! compare dimensions
    if (dims_treecode(1)>params%max_treelevel) then
        ! treecode in input file is greater than the new one, abort and output on screen
        ! NOTE this can be made working if not all levels in the file are actually used (e.g. level_max=17
        ! but active level=4). On the other hand, that appears to be rare.
        call abort(73947887, "ERROR: Treecode in file is longer than what is set in INI file.")
    end if

    allocate(block_treecode(1:dims_treecode(1), 1:hvy_n))
    block_treecode = -1

    ! tell the hdf5 wrapper what part of the global [ n_active x max_treelevel + IDX_REFINE_STS]
    ! array we want to hold, so that all CPU can read from the same file simultaneously
    ! (note zero-based offset):
    lbounds = (/0, sum(blocks_per_rank_list(0:rank-1))/)
    ubounds = (/int(dims_treecode(1),4)-1, lbounds(2) + hvy_n - 1/)
    call read_dset_mpi_hdf5_2D(file_id, "block_treecode", lbounds, ubounds, block_treecode)

    ! close file and HDF5 library
    call close_file_hdf5(file_id)

    ! this resetting is expensive but we do it only once:
    lgt_block = -1

     do k = 1, hvy_n
        call hvy2lgt( lgt_id, k, rank, params%number_blocks )
        ! copy treecode
        lgt_block(lgt_id, 1:dims_treecode(1)) = block_treecode(1:dims_treecode(1), k)
        ! set mesh level
        lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL) = treecode_size(block_treecode(:,k), size(block_treecode,1))
        ! set refinement status
        lgt_block(lgt_id, params%max_treelevel+IDX_REFINE_STS) = 0
        ! set number of the tree
        lgt_block(lgt_id, params%max_treelevel+IDX_TREE_ID) = tree_id
    end do

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
        write(*,*) "Done reading MESH! (NOT the actual data!)"
        write(*,'(80("~"))')
    endif

end subroutine read_mesh
