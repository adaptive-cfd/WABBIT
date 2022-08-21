!> \brief read data of a single datafield dF at iteration and time t
! ********************************************************************************************

subroutine read_field(fname, dF, params, hvy_block, hvy_n, tree_ID)
    implicit none

    character(len=*), intent(in)        :: fname                      !> file name
    integer(kind=ik), intent(in)        :: dF                         !> datafield number
    type (type_params), intent(in)      :: params                     !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    integer(kind=ik), intent(in)        :: hvy_n(:)                   !> number of heavy and light active blocks
    integer(kind=ik), intent(in)        :: tree_ID

    integer(hid_t)                      :: file_id
    integer(kind=ik)                    :: rank                       ! process rank
    integer(kind=ik)                    :: g, k                       ! grid parameter
    integer(kind=ik), dimension(3)      :: Bs
    integer(kind=ik), dimension(4)      :: ubounds3D, lbounds3D       ! offset variables
    integer(kind=ik), dimension(3)      :: ubounds2D, lbounds2D

    ! procs per rank array
    integer, dimension(:), allocatable  :: actual_blocks_per_proc

    rank = params%rank
    g    = params%n_ghosts
    Bs   = params%Bs
    allocate(actual_blocks_per_proc( 0:params%number_procs-1 ))

    call check_file_exists(fname)
    ! open the file
    call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)
    call blocks_per_mpirank( params, actual_blocks_per_proc, hvy_n(tree_ID) )
    if ( params%dim == 3 ) then

        ! tell the hdf5 wrapper what part of the global [Bsx x Bsy x Bsz x hvy_n]
        ! array we want to hold, so that all CPU can read from the same file simultaneously
        ! (note zero-based offset):
        lbounds3D = (/0,0,0,sum(actual_blocks_per_proc(0:rank-1))/)
        ubounds3D = (/Bs(1)-1,Bs(2)-1,Bs(3)-1,lbounds3D(4)+hvy_n(tree_ID)-1/)

    else

        ! tell the hdf5 wrapper what part of the global [Bsx x Bsy x 1 x hvy_n]
        ! array we want to hold, so that all CPU can read from the same file simultaneously
        ! (note zero-based offset):
        lbounds2D = (/0,0,sum(actual_blocks_per_proc(0:rank-1))/)
        ubounds2D = (/Bs(1)-1,Bs(2)-1,lbounds2D(3)+hvy_n(tree_ID)-1/)

    endif

    ! print a message
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("READING: Reading datafield ",i2," from file ",A)') dF,&
        trim(adjustl(fname))
    end if

    ! actual reading of file
    if ( params%dim == 3 ) then
        ! 3D data case
        call read_dset_mpi_hdf5_4D(file_id, "blocks", lbounds3D, ubounds3D, &
            hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, 1:hvy_n(tree_ID)))
    else
        ! 2D data case
        call read_dset_mpi_hdf5_3D(file_id, "blocks", lbounds2D, ubounds2D, &
            hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, dF, 1:hvy_n(tree_ID)))
    end if

    ! close file and HDF5 library
    call close_file_hdf5(file_id)
    ! check if field contains NaNs
    do k = 1, hvy_n(tree_ID)
        if ( params%dim == 3 ) then
            if (block_contains_NaN(hvy_block(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,dF,k))) &
                call abort(0200, "ERROR: Saved field "//get_dsetname(fname)//" contains NaNs!! I don't want to read from this file!")
        else
            if (block_contains_NaN(hvy_block(g+1:Bs(1)+g,g+1:Bs(2)+g,:,dF,k))) &
                call abort(0200, "ERROR: Saved field "//get_dsetname(fname)//" contains NaNs!! I don't want to read from this file!")
        end if
    end do
end subroutine read_field
