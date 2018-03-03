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

subroutine read_mesh(fname, params, lgt_n, hvy_n, lgt_block)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> file name
    character(len=*), intent(in)                  :: fname
    !> user defined parameter structure
    type (type_params), intent(in)                :: params
    !> number of active blocks (heavy and light data)
    integer(kind=ik), intent(inout)               :: hvy_n, lgt_n
    !> light data array
    integer(kind=ik), intent(inout)               :: lgt_block(:,:)

    ! file id integer
    integer(hid_t)                                :: file_id
    ! process rank, number of procs
    integer(kind=ik)                              :: rank, number_procs
    ! grid parameter
    integer(kind=ik)                              :: Bs, g
    ! offset variables
    integer(kind=ik), dimension(2)                :: ubounds, lbounds
    ! treecode array
    integer(kind=ik), dimension(:,:), allocatable :: block_treecode
    integer(kind=ik), dimension(:,:), allocatable :: my_lgt_block
    integer(kind=ik)                              :: blocks_per_rank_list(0:params%number_procs-1) 
    ! loop variables
    integer(kind=rk)                              :: lgt_id, k
    ! error variable
    integer(kind=ik)                              :: ierr
    integer(kind=ik)                              :: treecode_size
    integer(hsize_t), dimension(2)                :: dims_treecode
!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank         = params%rank
    number_procs = params%number_procs

    ! grid parameter
    Bs   = params%number_block_nodes
    g    = params%number_ghost_nodes
    
    lgt_id = 0
!---------------------------------------------------------------------------------------------
! main body

    call check_file_exists(fname)
    ! open the file
    call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)

    if ( (rank == 0) ) then
        write(*,'(80("_"))')
        write(*,'(A)') "READING: initializing grid from file..."
        write(*,'( "Nblocks=",i6," (on all cpus)")') lgt_n
        ! check if there is already some data on the grid
        if ( maxval(lgt_block(:,1))>=0 ) then
            write(*,'(A)') "ERROR: READ_MESH is called with NON_EMPTY DATA!!!!!"
        end if
    end if

    ! Nblocks per CPU
    ! this list contains (on each mpirank) the number of blocks for each mpirank. note
    ! zero indexing as required by MPI
    ! set list to the average value
    blocks_per_rank_list = lgt_n / number_procs

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
    if (dims_treecode(1)<=params%max_treelevel) then
        ! max_treelevels are the same or
        ! want to continue with greater max_treevel, read old treecode and fill up the rest with -1
    else
        ! treecode in input file is greater than the new one, abort and output on screen
        if (params%threeD_case) write(*,'("ERROR: max_treelevel is smaller than saved in file, this is not possible.",/ ,"max_treelevel in ini-file:",i4," in input-file:",i4)') params%max_treelevel, dims_treecode(1)
        call MPI_ABORT( MPI_COMM_WORLD, 10004, ierr)
    end if

    allocate(block_treecode(1:dims_treecode(1), 1:hvy_n))
    allocate (my_lgt_block(size(lgt_block,1), size(lgt_block,2)))
    my_lgt_block = -1
    block_treecode = -1
    ! tell the hdf5 wrapper what part of the global [ n_active x max_treelevel + 2]
    ! array we want to hold, so that all CPU can read from the same file simultaneously
    ! (note zero-based offset):
    lbounds = (/0, sum(blocks_per_rank_list(0:rank-1))/)
    ubounds = (/int(dims_treecode(1),4)-1, lbounds(2) + hvy_n - 1/)
    call read_dset_mpi_hdf5_2D(file_id, "block_treecode", lbounds, ubounds, block_treecode)
    
        
    ! close file and HDF5 library
    call close_file_hdf5(file_id)
     do k=1, hvy_n
        call hvy_id_to_lgt_id( lgt_id, k, rank, params%number_blocks )
        ! copy treecode
        my_lgt_block(lgt_id,1:dims_treecode(1)) = block_treecode(1:dims_treecode(1),k)
        ! set mesh level
        my_lgt_block(lgt_id, params%max_treelevel+1) = treecode_size(block_treecode(:,k), dims_treecode(1))
        ! set refinement status 
        my_lgt_block(lgt_id, params%max_treelevel+2) = 0
    end do
    ! synchronize light data. This is necessary as all CPUs above created their blocks locally.
    ! As they all pass the same do loops, the counter array blocks_per_rank_list does not have to
    ! be synced. However, the light data has to.
    lgt_block = -1
    call MPI_Allreduce(my_lgt_block, lgt_block, size(lgt_block,1)*size(lgt_block,2), MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ierr)

    deallocate(my_lgt_block)
    deallocate(block_treecode)

end subroutine read_mesh
