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
    character(len=*), intent(in)      :: fname
    !> user defined parameter structure
    type (type_params), intent(in)    :: params
    !> number of active blocks (heavy and light data)
    integer(kind=ik), intent(inout)   :: hvy_n, lgt_n
    !> light data array
    integer(kind=ik), intent(inout)   :: lgt_block(:,:)

    ! treecode array
    integer(kind=ik), dimension(:,:), allocatable :: block_treecode
    ! file id integer
    integer(hid_t)        :: file_id
    ! process rank, number of procs
    integer(kind=ik)      :: rank, number_procs
    ! grid parameter
    integer(kind=ik)      :: Bs, g
    ! offset variables
    integer(kind=ik)      :: ubounds(2), lbounds(2)
    integer(kind=ik)      :: blocks_per_rank_list(0:params%number_procs-1)
    integer(kind=ik)      :: lgt_id, k
    integer(kind=ik)      :: ierr
    integer(kind=ik)      :: treecode_size
    integer(hsize_t)      :: dims_treecode(2)

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank         = params%rank
    number_procs = params%number_procs
    ! grid parameter
    Bs   = params%Bs
    g    = params%n_ghosts

!---------------------------------------------------------------------------------------------
! main body

    call check_file_exists(fname)
    ! open the file
    call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)

    if ( rank == 0 ) then
        write(*,'(80("_"))')
        write(*,'(A)') "READING: initializing grid from file..."
        write(*,'("Filename: ",A)') trim(adjustl(fname))
        ! NOTE: we pass the routine lgt_n, which must be previsouly read from the
        ! file. This is done using read_attributes. This can possibly be merged here
        ! and mustnt be done in the caller
        write(*,'("Expected Nblocks=",i6," (on all cpus)")') lgt_n
        ! check if there is already some data on the grid
        if ( maxval(lgt_block(:,1))>=0 ) then
            write(*,'(A)') "ERROR: READ_MESH is called with NON_EMPTY DATA!!!!!"
        end if
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

    ! tell the hdf5 wrapper what part of the global [ n_active x max_treelevel + idx_refine_sts]
    ! array we want to hold, so that all CPU can read from the same file simultaneously
    ! (note zero-based offset):
    lbounds = (/0, sum(blocks_per_rank_list(0:rank-1))/)
    ubounds = (/int(dims_treecode(1),4)-1, lbounds(2) + hvy_n - 1/)
    call read_dset_mpi_hdf5_2D(file_id, "block_treecode", lbounds, ubounds, block_treecode)

    ! close file and HDF5 library
    call close_file_hdf5(file_id)

    ! this is expensive but we do it only once:
    lgt_block = -1

     do k = 1, hvy_n
        call hvy_id_to_lgt_id( lgt_id, k, rank, params%number_blocks )
        ! copy treecode
        lgt_block(lgt_id, 1:dims_treecode(1)) = block_treecode(1:dims_treecode(1), k)
        ! set mesh level
        lgt_block(lgt_id, params%max_treelevel+idx_mesh_lvl) = treecode_size(block_treecode(:,k), size(block_treecode,1))
        ! set refinement status
        lgt_block(lgt_id, params%max_treelevel+idx_refine_sts) = 0
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
    endif

end subroutine read_mesh
