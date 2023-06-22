subroutine saveHDF5_tree(fname, time, iteration, dF, params, hvy_block, tree_ID, no_sync)
    implicit none

    character(len=*), intent(in)        :: fname                    !> file name
    real(kind=rk), intent(in)           :: time                     !> time loop parameters
    integer(kind=ik), intent(in)        :: iteration
    integer(kind=ik), intent(in)        :: dF                       !> datafield number
    type (type_params), intent(in)      :: params                   !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :) !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID
    logical, optional, intent(in)       :: no_sync

    integer(kind=ik)                    :: rank, lgt_rank                       ! process rank
    integer(kind=ik)                    :: k, hvy_id, l, lgt_id, status         ! loop variable
    integer(kind=ik)                    :: g, dim                               ! grid parameter
    integer(kind=ik), dimension(3)      :: Bs
    ! block data buffer, need for compact data storage
    real(kind=rk), allocatable          :: myblockbuffer(:,:,:,:)
    ! coordinates and spacing arrays
    real(kind=rk), allocatable          :: coords_origin(:,:), coords_spacing(:,:)
    ! treecode array
    integer(kind=ik), allocatable       :: block_treecode(:,:)
    integer(hid_t)                      :: file_id
    ! offset variables
    integer(kind=ik), dimension(1:4)    :: ubounds3D, lbounds3D
    integer(kind=ik), dimension(1:3)    :: ubounds2D, lbounds2D, periodic_BC, symmetry_BC

    character(len=cshort)               :: arg
    type(INIFILE)                       :: FILE

    logical, parameter                  :: save_ghosts = .false.

    ! procs per rank array
    integer, dimension(:), allocatable  :: actual_blocks_per_proc

    ! spacing and origin (new)
    real(kind=rk) :: xx0(1:3) , ddx(1:3), sparsity_Jcurrent, sparsity_Jmax
    integer(kind=ik), allocatable :: procs(:), lgt_ids(:), refinement_status(:)
    logical :: no_sync2
    integer(kind=ik) :: Jmin_active, Jmax_active

    no_sync2 = .false.
    if (present(no_sync)) no_sync2 = no_sync

    ! uniqueGrid modification
    if (.not. no_sync2) then
        ! because when saving pruned trees, sync is not possible...
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, &
        hvy_active(:,tree_ID), hvy_n(tree_ID) )
    endif
    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)

    rank = params%rank
    Bs   = params%Bs
    g    = params%g
    dim  = params%dim
    periodic_BC = 0_ik
    symmetry_BC = 0_ik

    where (params%periodic_BC)
        periodic_BC = 1_ik
    end where

    where (params%symmetry_BC)
        symmetry_BC = 1_ik
    end where

    ! to know our position in the last index of the 4D output array, we need to
    ! know how many blocks all procs have
    allocate(actual_blocks_per_proc( 0:params%number_procs-1 ))
    if (save_ghosts) then
        allocate(myblockbuffer( 1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g, 1:hvy_n(tree_ID) ), stat=status)
    else
        ! uniqueGrid modification: we save the first ghost node in addition to the internal
        ! points for visualization
        allocate(myblockbuffer( 1:Bs(1)+1, 1:Bs(2)+1, 1:Bs(3)+1, 1:hvy_n(tree_ID) ), stat=status)
    endif

    if (status /= 0) then
        call abort(2510191, "IO: sorry, but buffer allocation failed! At least the weather is clearing up. Go outside.")
    endif

    allocate(coords_spacing(1:3, 1:hvy_n(tree_ID)) )
    allocate(coords_origin(1:3, 1:hvy_n(tree_ID)))
    allocate(procs(1:hvy_n(tree_ID)))
    allocate(lgt_ids(1:hvy_n(tree_ID)))
    allocate(refinement_status(1:hvy_n(tree_ID)))
    procs = rank
    allocate(block_treecode(1:params%Jmax, 1:hvy_n(tree_ID)))

    coords_origin = 7.0e6_rk

    if (lgt_n(tree_ID) < 1 ) then
        write(*,*) "Trying to save tree_ID=", tree_ID, " but tree is empty!"
        call abort(291019, "you try to save an empty mesh.")
    endif

#ifdef DEV
    ! first: check if field contains NaNs
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        if (block_contains_NaN(hvy_block(:,:,:,dF,hvy_id))) call abort(0201, "ERROR: Field "//get_dsetname(fname)//" contains NaNs!! We should not save this...")
    end do
#endif

    ! output on screen
    if (rank == 0) then
        sparsity_Jcurrent = dble(lgt_n(tree_ID)) / dble(2**maxActiveLevel_tree(tree_ID))**dim
        sparsity_Jmax = dble(lgt_n(tree_ID)) / dble(2**params%Jmax)**dim

        write(*,'("IO: saving HDF5 file t=", f15.8," file = ",A," J=(",i2," | ",i2,")", " Nblocks=",i7," sparsity=(",f5.1,"% / ",f5.1,"%)")') &
        time, trim(adjustl(fname)), Jmin_active, Jmax_active, lgt_n(tree_ID), 100.0*sparsity_Jcurrent, 100.0*sparsity_Jmax
    endif

    ! we need to know how many blocks each rank actually holds, and all procs need to
    ! know that distribution for all other procs in order to know what portion of the array
    ! they must write to.
    call blocks_per_mpirank( params, actual_blocks_per_proc, tree_ID )

    ! fill blocks buffer (we cannot use the bvy_block array as it is not contiguous, i.e.
    ! it may contain holes)
    if ( params%dim == 3 ) then

        ! tell the hdf5 wrapper what part of the global [bsx x bsy x bsz x n_active]
        ! array we hold, so that all CPU can write to the same file simultaneously
        ! (note zero-based offset):
        lbounds3D = (/1, 1, 1, sum(actual_blocks_per_proc(0:rank-1))+1/) - 1
        ! ubounds3D = (/Bs(1), Bs(2), Bs(3), lbounds3D(4)+hvy_n(tree_ID)/) - 1
        ! uniqueGrid modification:
        ubounds3D = (/Bs(1)+1, Bs(2)+1, Bs(3)+1, lbounds3D(4)+hvy_n(tree_ID)/) - 1
        if (save_ghosts) ubounds3D = (/Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, lbounds3D(4)+hvy_n(tree_ID)/) - 1

    else

        ! tell the hdf5 wrapper what part of the global [bsx x bsy x n_active]
        ! array we hold, so that all CPU can write to the same file simultaneously
        ! (note zero-based offset):
        lbounds2D = (/1, 1, sum(actual_blocks_per_proc(0:rank-1))+1/) - 1
        ! ubounds2D = (/Bs(1), Bs(2), lbounds2D(3)+hvy_n(tree_ID)/) - 1
        ! uniqueGrid modification:
        ubounds2D = (/Bs(1)+1, Bs(2)+1, lbounds2D(3)+hvy_n(tree_ID)/) - 1
        if (save_ghosts) ubounds2D = (/Bs(1)+2*g, Bs(2)+2*g, lbounds2D(3)+hvy_n(tree_ID)/) - 1

    endif

    l = 1
    ! loop over all active light block IDs, check if it is mine, if so, copy the block to the buffer
    do k = 1, lgt_n(tree_ID)
        lgt_id = lgt_active(k, tree_ID)
        ! calculate proc rank from light data line number
        call lgt2proc( lgt_rank, lgt_id, params%number_blocks )
        ! calculate heavy block id corresponding to light id
        call lgt2hvy( hvy_id, lgt_id, rank, params%number_blocks )

        ! if I own this block, I copy it to the buffer.
        ! also extract block coordinate origin and spacing
        if (lgt_rank == rank) then
            ! compute block spacing and origin from the treecode
            call get_block_spacing_origin( params, lgt_id, xx0, ddx )


            if ( params%dim == 3 ) then
                ! 3D
                if (save_ghosts) then
                    myblockbuffer(:,:,:,l) = hvy_block( :, :, :, dF, hvy_id)
                else
                    ! myblockbuffer(:,:,:,l) = hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, hvy_id)
                    ! uniqueGrid modification:
                    myblockbuffer(:,:,:,l) = hvy_block( g+1:Bs(1)+g+1, g+1:Bs(2)+g+1, g+1:Bs(3)+g+1, dF, hvy_id)
                endif

                ! note reverse ordering (paraview uses C style, we fortran...life can be hard)
                coords_origin(1,l) = xx0(3)
                coords_origin(2,l) = xx0(2)
                coords_origin(3,l) = xx0(1)
                coords_spacing(1,l) = ddx(3)
                coords_spacing(2,l) = ddx(2)
                coords_spacing(3,l) = ddx(1)

                if (save_ghosts) then
                    coords_origin(1,l) = xx0(3) -dble(g)*ddx(3)
                    coords_origin(2,l) = xx0(2) -dble(g)*ddx(2)
                    coords_origin(3,l) = xx0(1) -dble(g)*ddx(1)
                endif

                ! copy treecode (we'll save it to file as well)
                block_treecode(:,l) = lgt_block( lgt_id, 1:params%Jmax )
            else
                ! 2D
                if (save_ghosts) then
                    myblockbuffer(:,:,1,l) = hvy_block(:, :, 1, dF, hvy_id)
                else
                    ! myblockbuffer(:,:,1,l) = hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, 1, dF, hvy_id)
                    ! uniqueGrid modification:
                    myblockbuffer(:,:,1,l) = hvy_block( g+1:Bs(1)+g+1, g+1:Bs(2)+g+1, 1, dF, hvy_id)
                endif

                ! note reverse ordering (paraview uses C style, we fortran...life can be hard)
                coords_origin(1,l) = xx0(2)
                coords_origin(2,l) = xx0(1)
                coords_spacing(1,l) = ddx(2)
                coords_spacing(2,l) = ddx(1)

                if (save_ghosts) then
                    coords_origin(1,l) = xx0(2) -dble(g)*ddx(1)
                    coords_origin(2,l) = xx0(1) -dble(g)*ddx(1)
                endif
                ! copy treecode (we'll save it to file as well)
                block_treecode(:,l) = lgt_block( lgt_id, 1:params%Jmax )
            endif


            refinement_status(l) = lgt_block( lgt_id, params%Jmax+IDX_REFINE_STS )
            lgt_ids(l) = lgt_id

            ! next block
            l = l + 1
        endif

    end do

    ! open the file
    call open_file_hdf5( trim(adjustl(fname)), file_id, .true.)

    ! write heavy block data to disk
    if ( params%dim == 3 ) then
        ! 3D data case
        call write_dset_mpi_hdf5_4D(file_id, "blocks", lbounds3D, ubounds3D, myblockbuffer)
        call write_attribute(file_id, "blocks", "domain-size", (/params%domain_size(1), params%domain_size(2), params%domain_size(3)/))
        call write_attribute(file_id, "blocks", "periodic_BC", periodic_BC )
        call write_attribute(file_id, "blocks", "symmetry_BC", symmetry_BC )
        call write_dset_mpi_hdf5_2D(file_id, "coords_origin", (/0,lbounds3D(4)/), (/2,ubounds3D(4)/), coords_origin)
        call write_dset_mpi_hdf5_2D(file_id, "coords_spacing", (/0,lbounds3D(4)/), (/2,ubounds3D(4)/), coords_spacing)
        call write_dset_mpi_hdf5_2D(file_id, "block_treecode", (/0,lbounds3D(4)/), (/params%Jmax-1,ubounds3D(4)/), block_treecode)
        call write_int_dset_mpi_hdf5_1D(file_id, "procs", (/lbounds3D(4)/), (/ubounds3D(4)/), procs)
        call write_int_dset_mpi_hdf5_1D(file_id, "refinement_status", (/lbounds3D(4)/), (/ubounds3D(4)/), refinement_status)
        call write_int_dset_mpi_hdf5_1D(file_id, "lgt_ids", (/lbounds3D(4)/), (/ubounds3D(4)/), lgt_ids)
    else
        ! 2D data case
        call write_dset_mpi_hdf5_3D(file_id, "blocks", lbounds2D, ubounds2D, myblockbuffer(:,:,1,:))
        call write_attribute(file_id, "blocks", "domain-size", (/params%domain_size(1), params%domain_size(2)/))
        call write_attribute(file_id, "blocks", "periodic_BC", periodic_BC )
        call write_attribute(file_id, "blocks", "symmetry_BC", symmetry_BC )
        call write_dset_mpi_hdf5_2D(file_id, "coords_origin", (/0,lbounds2D(3)/), (/1,ubounds2D(3)/), coords_origin(1:2,:))
        call write_dset_mpi_hdf5_2D(file_id, "coords_spacing", (/0,lbounds2D(3)/), (/1,ubounds2D(3)/), coords_spacing(1:2,:))
        call write_dset_mpi_hdf5_2D(file_id, "block_treecode", (/0,lbounds2D(3)/), (/params%Jmax-1,ubounds2D(3)/), block_treecode)
        call write_int_dset_mpi_hdf5_1D(file_id, "procs", (/lbounds2D(3)/), (/ubounds2D(3)/), procs)
        call write_int_dset_mpi_hdf5_1D(file_id, "refinement_status", (/lbounds2D(3)/), (/ubounds2D(3)/), refinement_status)
        call write_int_dset_mpi_hdf5_1D(file_id, "lgt_ids", (/lbounds2D(3)/), (/ubounds2D(3)/), lgt_ids)
    endif


    ! add additional annotations
    call write_attribute(file_id, "blocks", "version", (/20231602/)) ! this is used to distinguish wabbit file formats
    call write_attribute(file_id, "blocks", "block-size", Bs)
    call write_attribute(file_id, "blocks", "time", (/time/))
    call write_attribute(file_id, "blocks", "iteration", (/iteration/))
    call write_attribute(file_id, "blocks", "total_number_blocks", (/lgt_n(tree_ID)/))

    ! close file and HDF5 library
    call close_file_hdf5(file_id)

    ! clean up
    deallocate(actual_blocks_per_proc)
    deallocate(myblockbuffer)
    deallocate(coords_origin)
    deallocate(coords_spacing)
    deallocate(block_treecode, procs, refinement_status, lgt_ids)

    ! check if we find a *.ini file name in the command line call
    ! if we do, read it, and append it to the HDF5 file. this way, data
    ! and parameters are always together. thomas, 16/02/2019

    ! if (params%rank==0) then
    !     call open_file_hdf5_serial( trim(adjustl(fname)), file_id, .true.)
    !     do k = 1, COMMAND_ARGUMENT_COUNT()
    !         call get_command_argument( k, arg )
    !         if (index(arg, ".ini") /= 0) then
    !             ! found the ini file.
    !             call read_ini_file(FILE, arg, .false., remove_comments=.false.)
    !             call write_string_dset_hdf5(file_id, "params", FILE%PARAMS, maxcolumns)
    !             call write_attribute(file_id, "params", "filename", arg)
    !         end if
    !     enddo
    !     call close_file_hdf5_serial(file_id)
    ! endif
end subroutine saveHDF5_tree



!-------------------------------------------------------------------------------
! Read a tree (=grid) in the forest. A tree can consist of several datafields (think
! of the state vector) that all live on the same grid. They are stored in individual
! HDF5 files.
!
! Note: all files are supposed to share the same grid (components of a vector). No
! test is made if that is true: if the grids are different, garbage will be read.
!
!-------------------------------------------------------------------------------
subroutine readHDF5vct_tree(fnames, params, hvy_block, tree_ID, time, iteration, verbosity, synchronize_ghosts)
    implicit none

    character(len=*), intent(in)   :: fnames(1:)      !< list of files to be read (=number of vector components)
    type (type_params), intent(in) :: params          !< user defined parameter structure
    integer(kind=ik), intent(in)   :: tree_ID         !< index of the tree you want to save the data in
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :) !< heavy data array - data are stored herein

    logical, intent(in), optional  :: verbosity       !< if verbosity==True generates log output
    logical, intent(in), optional  :: synchronize_ghosts
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
    g            = params%g
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
    if (NblocksFile > size(hvy_block,5)*number_procs) then
        call abort(2408223, "The file we read is larger than the allocated array!")
    endif


    !----------------------------------------------------------------------------------
    ! version check
    !----------------------------------------------------------------------------------
    ! Files created using the newGhostNodes branch (after 08 April 2020) contain a version number.
    ! if the number is not found, version=0.
    !
    ! VERSION   DESCRIPTION
    ! =======   ===========
    !        0  redundantGrid created by all WABBIT versions before 1st attempt to use uniqueGrid
    ! 20200408  uniqueGrid created by the intermediate
    ! 20200902  redundantGrid after we first abandonned the uniqueGrid
    ! 20231602  Current: uniqueGrid and with biorthogonal wavelets: equivalent to 20200408
    !
    if ((version(1) < 20231602) .and. (version(1)>0)) then
        if (rank == 0) then
            write(*,*) "--------------------------------------------------------------------"
            write(*,*) "-----ERROR------------ERROR------------ERROR------------ERROR-------"
            write(*,*) "--------------------------------------------------------------------"
            write(*,*) "The file we are trying to read is generated with an old version"
            write(*,*) "of wabbit (before 2023). In this the file, the grid"
            write(*,*) "definition includes a redundant point, i.e., a block is defined"
            write(*,*) "with spacing dx = L*2^-J / Bs-1. This definition was replaced by uniqueGrid."
        endif
        call abort(2304061, "Wrong WABBIT file format.")
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
    if (dims_treecode(1) > params%Jmax) then
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
        ! NOTE: uniqueGrid stores the FIRST GHOST NODE as well as the interior points,
        ! important for visualization. Hence, Bs+1 points are read.
        allocate( hvy_buffer(1:Bs(1)+1, 1:Bs(2)+1, 1:Bs(3)+1, 1:N_files, 1:my_hvy_n) )
        ! tell the hdf5 wrapper what part of the global [Bsx x Bsy x Bsz x hvy_n]
        ! array we want to hold, so that all CPU can read from the same file simultaneously
        ! (note zero-based offset):
        lbounds3D = (/0, 0, 0, sum(blocks_per_rank_list(0:rank-1))/)
        ubounds3D = (/Bs(1), Bs(2), Bs(3), lbounds3D(4)+my_hvy_n-1/)

        ! 3D data case
        ! actual reading of file (into hvy_buffer! not hvy_work)
        do dF = 1, N_files
            call open_file_hdf5( trim(adjustl(fnames(dF))), file_id, .false.)

            if (rank==0 .and. verbose) then
                write(*,'("READING: Reading datafield ",i2," from file ",A)') dF, trim(adjustl(fnames(dF)))
            end if

            call read_dset_mpi_hdf5_4D(file_id, "blocks", lbounds3D, ubounds3D, &
            hvy_buffer(:, :, :, dF, 1:my_hvy_n))

            ! close file and HDF5 library
            call close_file_hdf5(file_id)
        enddo
    else
        ! NOTE: uniqueGrid stores the FIRST GHOST NODE as well as the interior points,
        ! important for visualization. Hence, Bs+1 points are read.
        allocate( hvy_buffer(1:Bs(1)+1, 1:Bs(2)+1, 1, 1:N_files, 1:my_hvy_n) )
        ! tell the hdf5 wrapper what part of the global [Bsx x Bsy x 1 x hvy_n]
        ! array we want to hold, so that all CPU can read from the same file simultaneously
        ! (note zero-based offset):
        lbounds2D = (/0,0,sum(blocks_per_rank_list(0:rank-1))/)
        ubounds2D = (/Bs(1),Bs(2),lbounds2D(3)+my_hvy_n-1/)

        ! 2D data case
        ! actual reading of file (into hvy_buffer! not hvy_work)
        do dF = 1, N_files
            call open_file_hdf5( trim(adjustl(fnames(dF))), file_id, .false.)

            if (rank==0 .and. verbose) then
                write(*,'("READING: Reading datafield ",i2," from file ",A)') dF, trim(adjustl(fnames(dF)))
            end if

            call read_dset_mpi_hdf5_3D(file_id, "blocks", lbounds2D, ubounds2D, &
            hvy_buffer(:, :, 1, dF, 1:my_hvy_n))

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
        lgt_block(free_lgt_id, params%Jmax+IDX_MESH_LVL) = treecode_size(block_treecode(:,k), size(block_treecode,1))
        ! set refinement status
        lgt_block(free_lgt_id, params%Jmax+IDX_REFINE_STS) = 0
        ! set number of the tree
        lgt_block(free_lgt_id, params%Jmax+IDX_TREE_ID) = tree_id
        ! copy actual data (form buffer to actual data array)
        do dF = 1, N_files
            if (params%dim == 3) then
                ! NOTE: uniqueGrid stores the FIRST GHOST NODE as well as the interior points,
                ! important for visualization. Hence, Bs+1 points are read.
                hvy_block( g+1:Bs(1)+g+1, g+1:Bs(2)+g+1, g+1:Bs(3)+g+1, dF, free_hvy_id ) = hvy_buffer(:, :, :, dF, k)
            else
                hvy_block( g+1:Bs(1)+g+1, g+1:Bs(2)+g+1, :, dF, free_hvy_id ) = hvy_buffer(:, :, :, dF, k)
            endif
        end do
    end do

    deallocate(hvy_buffer)
    deallocate(block_treecode)


    ! synchronize light data. This is necessary as all CPUs above created their blocks locally.
    ! As they all pass the same do loops, the counter array blocks_per_rank_list does not have to
    ! be synced. However, the light data has to.
    call synchronize_lgt_data( params, refinement_status_only=.false. )

    ! it is good practice that this routine returns a working forest, i.e., all meta
    ! data is updated.
    call createActiveSortedLists_forest(params)
    call updateNeighbors_tree(params, tree_ID)

    if (present(synchronize_ghosts)) then
        if (synchronize_ghosts) then
            call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )
        endif
    else
        call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )
    endif

    ! it is useful to print out the information on active levels in the file
    ! to get an idea how it looks like and if the desired dense level is larger
    ! or smaller
    if (params%rank==0 .and. verbose ) then
        write(*,'("READING: In the file(s) we just read, Jmin=",i3," Jmax=",i3)') minActiveLevel_tree(tree_ID), &
        maxActiveLevel_tree(tree_ID)
    endif

end subroutine readHDF5vct_tree


!> \brief read attributes saved in a hdf5-file

subroutine read_attributes(fname, nBlocksFile, time, iteration, domain, Bs, tc_length, dim, &
    periodic_BC, symmetry_BC, verbosity)

    implicit none
    !> file name
    character(len=*), intent(in)                  :: fname
    !> number of active blocks (required to allocate light data, prior to reading)
    integer(kind=ik), intent(out)                 :: nBlocksFile
    !> time (to be read from file)
    real(kind=rk), intent(out)                    :: time
    !> iteration (to be read from file)
    integer(kind=ik), intent(out)                 :: iteration
    !> blocksize in the file (required to allocate light data, prior to reading)
    integer(kind=ik), dimension(3), intent(out)   :: Bs
    !> length of treecodes in the file (required to allocate light data, prior to reading)
    integer(kind=ik), intent(out)                 :: tc_length
    !> data dimensionality (2 or 3)
    integer(kind=ik), intent(out)                 :: dim
    !> domain size
    real(kind=rk), dimension(3), intent(out)      :: domain
    logical, intent(in), optional                 :: verbosity !< if verbosity==True generates log output
    logical, intent(inout), optional              :: periodic_BC(1:3), symmetry_BC(1:3)

    integer(kind=ik), dimension(1)                :: iiteration, number_blocks, version
    real(kind=rk), dimension(1)                   :: ttime
    integer(hid_t)                                :: file_id
    integer                                       :: datarank
    integer(kind=ik)                              :: rank, ierr, tmp1(1:3), tmp2(1:3), Nb
    integer(kind=hsize_t)                         :: size_field(1:4)
    integer(hsize_t), dimension(2)                :: dims_treecode
    logical :: verbose = .true.

    if (present(verbosity)) verbose=verbosity

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call check_file_exists(fname)

    ! open the file
    call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)

    !----------------------------------------------------------------------------------
    ! version check
    !----------------------------------------------------------------------------------
    ! Files created using the newGhostNodes branch (after 08 April 2020) contain a version number.
    ! if the number is not found, version=0.
    call read_attribute( file_id, "blocks", "version", version)

    ! if (version(1) == 20200408 .and. rank==0) then
    !     write(*,*) "--------------------------------------------------------------------"
    !     write(*,*) "-----WARNING----------WARNING----------WARNING----------WARNING-----"
    !     write(*,*) "--------------------------------------------------------------------"
    !     write(*,*) "The file we are trying to read is generated with an intermediate version"
    !     write(*,*) "of wabbit (after 08 April 2020). In this the file, the grid"
    !     write(*,*) "definition does not include a redundant point, i.e., a block is defined"
    !     write(*,*) "with spacing dx = L*2^-J / Bs. This definition was a dead-end, as it lead"
    !     write(*,*) "to instabilities and other problems. Current versions of WABBIT include a redundant point again."
    !     write(*,*) ""
    !     write(*,*) "The newGhostNodes branch still STORED the redundant point for visualization."
    !     write(*,*) "This was simply the first ghost node. If Bs was odd, this lead to an even number"
    !     write(*,*) "of points, and this cannot be read with present code versions."
    !     write(*,*) ""
    !     write(*,*) "A workaround must be done in preprocessing: upsampling to equidistant resolution and"
    !     write(*,*) "re-gridding is required. I am truely sorry for this."
    !     write(*,*) ""
    !     write(*,*) "If bs was even, the resulting data size is odd, and the file can be read. note however"
    !     write(*,*) "that the resolution changes slightly, and results cannot be perfectly identical to what would"
    !     write(*,*) "have been obtained with the newGhostNodes branch"
    !     write(*,*) "--------------------------------------------------------------------"
    !     write(*,*) "-----WARNING----------WARNING----------WARNING----------WARNING-----"
    !     write(*,*) "--------------------------------------------------------------------"
    ! endif

    if (present(periodic_BC)) then
        call read_attribute(file_id, "blocks", "periodic_BC", tmp1, (/1_ik, 1_ik, 1_ik/))

        periodic_BC = .false.
        where (tmp1 > 0_ik)
            periodic_BC = .true.
        end where
    endif

    if (present(symmetry_BC)) then
        call read_attribute(file_id, "blocks", "symmetry_BC", tmp2, (/0_ik, 0_ik, 0_ik/))

        symmetry_BC = .false.
        where (tmp2 > 0_ik)
            symmetry_BC = .true.
        end where
    endif

    call read_attribute(file_id, "blocks", "domain-size", domain)
    call read_attribute(file_id, "blocks", "time", ttime)
    time = ttime(1)
    call read_attribute(file_id, "blocks", "iteration", iiteration)
    iteration = iiteration(1)
    call read_attribute(file_id, "blocks", "total_number_blocks", number_blocks)
    nBlocksFile = number_blocks(1)

    !---------------------------------------------------------------------------
    ! Number of blocks and blocksize
    !---------------------------------------------------------------------------
    ! check if we deal with 2D or 3D data

    call read_attribute(file_id, "blocks", "block-size", Bs)
    call get_rank_datafield(file_id, "blocks", datarank)

    if (datarank == 3) then
        ! 2D data
        call get_size_datafield( datarank, file_id, "blocks", size_field(1:datarank))
        Bs(1) = int( size_field(1), kind=ik)-1 !uniqueGrid modification
        Bs(2) = int( size_field(2), kind=ik)-1 !uniqueGrid modification
        Bs(3) = 1
        Nb = int( size_field(3), kind=ik)
        domain(3) = 0.0_rk
        dim = 2
    elseif (datarank == 4) then
        ! 3D data
        call get_size_datafield( datarank, file_id, "blocks", size_field(1:datarank))
        Bs(1) = int( size_field(1), kind=ik)-1 !uniqueGrid modification
        Bs(2) = int( size_field(2), kind=ik)-1 !uniqueGrid modification
        Bs(3) = int( size_field(3), kind=ik)-1 !uniqueGrid modification
        Nb = int( size_field(4), kind=ik)
        dim = 3

    else
        ! crazy data
        call abort(33321, "Datarank neither 2d nor 3d..that is unusual.")

    endif

    !---------------------------------------------------------------------------
    ! length of treecodes in file
    !---------------------------------------------------------------------------
    ! NOTE: we do store only the treecode, not the level or refinement status
    ! so the length of this array is indeed the treecode length, and not treecode_length+2
    call get_size_datafield(2, file_id, "block_treecode", dims_treecode)
    tc_length = int(dims_treecode(1), kind=ik)

    ! close file and HDF5 library
    call close_file_hdf5(file_id)

    if (rank == 0 .and. verbose) then
        write(*,'(80("~"))')
        write(*,*) "read_attributes.f90: Read important numbers from a file."
        write(*,'(80("~"))')
        write(*,*) "These numbers are used either for initializiation/allocation (in post-"
        write(*,*) "processing) or to check if a file we try to load matches the"
        write(*,*) "specification of the current simulation."
        write(*,*) "We read:"
        write(*,'(" file         = ",A)') trim(adjustl(fname))
        write(*,'(" Bs           = ",3(i3,1x))') Bs
        write(*,'(" domain       = ",3(g15.8,1x))') domain
        write(*,'(" tc_length    = ",i3)') tc_length
        write(*,'(" nBlocksFile  = ",i8)') nBlocksFile
        write(*,'(" dim          = ",i8)') dim
        if (present(periodic_BC)) write(*,'(" periodic_BC  = ",3(L1))') periodic_BC
        if (present(symmetry_BC)) write(*,'(" symmetry_BC  = ",3(L1))') symmetry_BC
        write(*,'(80("~"))')
    endif

    if (modulo(Bs(1),2) /= 0) then
        write(*,*) Bs
        call abort(202009021, "Blocksize Bs(1) is an odd number, which this code version cannot handle.")
    endif
    if (modulo(Bs(2),2) /= 0) then
        call abort(202009021, "Blocksize Bs(2) is an odd number, which this code version cannot handle.")
    endif
    if (modulo(Bs(3),2) /= 0 .and. dim==3) then
        call abort(202009021, "Blocksize Bs(3) is an odd number, which this code version cannot handle.")
    endif
end subroutine read_attributes



! returns the number of blocks stored in an HDF5 file
function getNumberBlocksH5File(fname)
    character(len=*), intent(in) :: fname
    integer(kind=ik) :: getNumberBlocksH5File

    integer(hid_t) :: file_id

    ! open the file
    call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)

    ! read number of blocks...
    call read_attribute(file_id, "blocks", "total_number_blocks", getNumberBlocksH5File)

    ! close file and HDF5 library
    call close_file_hdf5(file_id)
end function




!##############################################################
!> This function reads a set of data into a specified tree.
!> The set of data can involve multiple files. Each file should
!> contain a quantity at the same snapshot (i.e. same iteration/time :: the same grid!).
!> If no data has been read in before, this function will allocate
!> the heavy data for you.
subroutine read_field2tree(params, fnames, N_files, tree_ID, hvy_block, verbosity)
    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params           !< params structure
    integer(kind=ik), intent(in)      :: N_files     !< number of fields/quantities to read
    character(len=*), intent(in)      :: fnames(N_files)  !< file names
    integer(kind=ik), intent(in)      :: tree_ID     !< number of the tree
    real(kind=rk), ALLOCATABLE, intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    logical, intent(in), optional :: verbosity !< if verbosity==True generates log output
    !-------------------------------- ---------------------------------
    integer(kind=ik)  :: iteration, dF, tc_length, dim, i, lgt_n_tmp, &
    rank, level, Bs(3), fsize, number_dense_blocks
    real(kind=rk)     :: time, domain(3)
    logical :: verbose = .true.

    if (present(verbosity)) verbose=verbosity

    ! set MPI parameter
    rank  = params%rank
    fsize = params%forest_size

    if (tree_ID <= tree_n) then
        ! the tree already exists: to overwrite it, we first delete the existing one
        call delete_tree(params, tree_ID)
        call createActiveSortedLists_forest(params)
    endif

    ! From any file (here: 1st in list), we read the essential parameters of the tree:
    ! most importantly Bs, dim, length of treecodes. NOTE: all trees in the forest must
    ! have the same Bs.
    call read_attributes(fnames(1), lgt_n_tmp, time, iteration, domain, Bs, &
    tc_length, dim, verbosity=verbose, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    ! Check if one tree already exists in the forest:
    ! If it doesnt initialize some important parameters like the Block size Bs
    ! and the spatial dimension of the data.
    ! If it does chFeck if the parameters are suitable for the forest
    if (allocated(hvy_block)) then
        if ( any(params%Bs /= Bs) ) then
            ! NOTE: all trees in the forest must have the same Bs.
            call abort(100119, "ERROR: Trying to read a tree (grid) with Bs different from existing trees in forest!")
        endif

        if (params%Jmax < tc_length) call abort(10119, "maximal treelevel incompatible")
    else
        params%Jmax = tc_length
        params%dim = dim
        params%n_eqn = N_files
        params%N_fields_saved = N_files
        params%Bs = Bs
        params%domain_size = domain

        ! we have to allocate grid if this routine is called for the first time
        call allocate_forest(params, hvy_block)

        hvy_neighbor = -1
        hvy_n = 0
        lgt_n = 0 ! reset number of active light blocks
        tree_n = 0 ! reset number of trees in forest
    endif

    call readHDF5vct_tree(fnames, params, hvy_block, tree_id, verbosity=.true.)


    call createActiveSortedLists_forest(params)
    call updateNeighbors_tree(params, tree_ID)

    call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

end subroutine read_field2tree
