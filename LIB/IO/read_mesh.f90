!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name read_mesh.f90
!> \version 0.5
!> \author sm
!
!> \brief read mesh properties of a field saved in a hdf5-file
!
!>
!! input:
!!           - parameter array
!!           - name of the file we want to read from
!!
!! output:
!!           - light block array
!!           - number of active blocks (light and heavy)
!!
!!
!! = log ======================================================================================
!! \n
!! 29/09/17 - create
!
! ********************************************************************************************

subroutine read_mesh(fname, params, lgt_n, hvy_n, lgt_block, tree_id_optional )
    implicit none

    !> file name
    character(len=*), intent(in)      :: fname
    !> user defined parameter structure
    type (type_params), intent(in)    :: params
    !> number of active blocks (heavy and light data)
    integer(kind=ik), intent(inout)   :: hvy_n, lgt_n
    !> light data array
    integer(kind=ik), intent(inout)   :: lgt_block(:,:)
     !> index of the tree you want to save the field in
    integer(kind=ik), optional, intent(in)   :: tree_id_optional

    ! treecode array
    integer(kind=ik), dimension(:,:), allocatable :: block_treecode
    ! file id integer
    integer(hid_t)        :: file_id
    ! process rank, number of procs
    integer(kind=ik)      :: rank, number_procs
    ! grid parameter
    integer(kind=ik)      :: g, datarank
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
    Bs           = params%Bs
    g            = params%n_ghosts
    Bs_file      = 1


    call check_file_exists(fname)
    ! open the file
    call open_file_hdf5(trim(adjustl(fname)), file_id, .false.)

    ! how much blocks are we going to read?
    call read_attribute(file_id, "blocks", "total_number_blocks", tmp)
    lgt_n = tmp(1)

    call get_rank_datafield(file_id, "blocks", datarank)
    call get_size_datafield(datarank, file_id, "blocks", size_field(1:datarank))


    ! Files created after 08 04 2020 contain a version number. If we try to read older
    ! files, version=0 is returned.
    call read_attribute( file_id, "blocks", "version", version)


    if (version(1) == 20200408) then
        ! Newer files also contain the block size in addition:
        call read_attribute( file_id, "blocks", "block-size", Bs_file)
    else
        ! The block size is also coded as the size of the dataset. Note newer files
        ! save in fact ONE MORE point (Bs+1). The additional point is the first ghost node.
        ! This is done purely because in visualization, paraview does not understand
        ! that there is no missing data.

        ! this is an old file: the dimension of the array in the HDF5 file must be
        ! (Bs+1)
        Bs_file(1:datarank-1) = size_field(1:datarank-1) - 1

        write(*,*) "--------------------------------------------------------------------"
        write(*,*) "-----WARNING----------WARNING----------WARNING----------WARNING-----"
        write(*,*) "--------------------------------------------------------------------"
        write(*,*) "The file we are trying to read is generated with an older version"
        write(*,*) "of wabbit (prior to 08 April 2020). In this the file, the grid"
        write(*,*) "definition includes a redundant point, i.e., a block is defined"
        write(*,*) "with spacing dx = L*2^-J / (Bs-1). In new versions, the redundant"
        write(*,*) "point is removed, leaving us with the spacing dx = L*2^-J / Bs"
        write(*,*) "The file can be read; provided that the new Bs (in this version of the code)"
        write(*,*) "is the old one -1, this is an EVEN number, which the code supports as"
        write(*,*) "of 29 Apr 2020."
        write(*,*) "--------------------------------------------------------------------"
        write(*,*) "-----WARNING----------WARNING----------WARNING----------WARNING-----"
        write(*,*) "--------------------------------------------------------------------"
    endif



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
    endif

    ! Nblocks per CPU
    ! this list contains (on each mpirank) the number of blocks for each mpirank. note
    ! zero indexing as required by MPI
    ! set list to the average value
    blocks_per_rank_list(:) = lgt_n / number_procs

    ! as this does not necessarily work out, distribute remaining blocks on the first CPUs
    if (mod(lgt_n, number_procs) > 0) then
        blocks_per_rank_list(0:mod(lgt_n, number_procs)-1) = blocks_per_rank_list(0:mod(lgt_n, number_procs)-1) + 1
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
        call hvy_id_to_lgt_id( lgt_id, k, rank, params%number_blocks )
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
    if (rank==0) then
        write(*,'("In the file we just read, Jmin=",i3," Jmax=",i3)') min_active_level( lgt_block ), &
        max_active_level( lgt_block )
        write(*,*) "Done reading MESH! (NOT the actual data!)"
        write(*,'(80("~"))')
    endif

end subroutine read_mesh
