!> \brief write data of a single datafield dF at timestep iteration and time t
!! input:    - time loop parameter
!!           - datafield number
!!           - parameter array
!!           - light data array
!!           - heavy data array
! ********************************************************************************************

subroutine saveHDF5_tree(fname, time, iteration, dF, params, hvy_block, tree_ID)
    implicit none

    character(len=*), intent(in)        :: fname                                !> file name
    real(kind=rk), intent(in)           :: time                                 !> time loop parameters
    integer(kind=ik), intent(in)        :: iteration
    integer(kind=ik), intent(in)        :: dF                                   !> datafield number
    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    real(kind=rk), intent(in)           :: hvy_block(:, :, :, :, :)             !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID

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

    logical, parameter                  :: save_ghosts = .False.

    ! procs per rank array
    integer, dimension(:), allocatable  :: actual_blocks_per_proc

    ! spacing and origin (new)
    real(kind=rk) :: xx0(1:3) , ddx(1:3), sparsity_Jcurrent, sparsity_Jmax
    integer(kind=ik), allocatable :: procs(:), lgt_ids(:), refinement_status(:)


    rank = params%rank
    Bs   = params%Bs
    g    = params%n_ghosts
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
        allocate(myblockbuffer( 1:Bs(1), 1:Bs(2), 1:Bs(3), 1:hvy_n(tree_ID) ), stat=status)
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
    allocate(block_treecode(1:params%max_treelevel, 1:hvy_n(tree_ID)))

    coords_origin = 7.0e6_rk

    if (lgt_n(tree_ID) < 1 ) call abort(291019, "you try to save an empty mesh.")


    ! first: check if field contains NaNs
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        if (block_contains_NaN(hvy_block(:,:,:,dF,hvy_id))) call abort(0201, "ERROR: Field "//get_dsetname(fname)//" contains NaNs!! We should not save this...")
    end do

    ! output on screen
    if (rank == 0) then
        sparsity_Jcurrent = dble(lgt_n(tree_ID)) / dble(2**maxActiveLevel_tree(tree_ID))**dim
        sparsity_Jmax = dble(lgt_n(tree_ID)) / dble(2**params%max_treelevel)**dim

        write(*,'("IO: writing data for time = ", f15.8," file = ",A," Nblocks=",i7," sparsity=(",f5.1,"% / ",f5.1,"%)")') &
        time, trim(adjustl(fname)), lgt_n(tree_ID), 100.0*sparsity_Jcurrent, 100.0*sparsity_Jmax
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
        ubounds3D = (/Bs(1), Bs(2), Bs(3), lbounds3D(4)+hvy_n(tree_ID)/) - 1
        if (save_ghosts) ubounds3D = (/Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, lbounds3D(4)+hvy_n(tree_ID)/) - 1

    else

        ! tell the hdf5 wrapper what part of the global [bsx x bsy x n_active]
        ! array we hold, so that all CPU can write to the same file simultaneously
        ! (note zero-based offset):
        lbounds2D = (/1, 1, sum(actual_blocks_per_proc(0:rank-1))+1/) - 1
        ubounds2D = (/Bs(1), Bs(2), lbounds2D(3)+hvy_n(tree_ID)/) - 1
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
                    myblockbuffer(:,:,:,l) = hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, dF, hvy_id)
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
                block_treecode(:,l) = lgt_block( lgt_id, 1:params%max_treelevel )
            else
                ! 2D
                if (save_ghosts) then
                    myblockbuffer(:,:,1,l) = hvy_block(:, :, 1, dF, hvy_id)
                else
                    myblockbuffer(:,:,1,l) = hvy_block( g+1:Bs(1)+g, g+1:Bs(2)+g, 1, dF, hvy_id)
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
                block_treecode(:,l) = lgt_block( lgt_id, 1:params%max_treelevel )
            endif


            refinement_status(l) = lgt_block( lgt_id, params%max_treelevel+IDX_REFINE_STS )
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
        call write_dset_mpi_hdf5_2D(file_id, "block_treecode", (/0,lbounds3D(4)/), (/params%max_treelevel-1,ubounds3D(4)/), block_treecode)
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
        call write_dset_mpi_hdf5_2D(file_id, "block_treecode", (/0,lbounds2D(3)/), (/params%max_treelevel-1,ubounds2D(3)/), block_treecode)
        call write_int_dset_mpi_hdf5_1D(file_id, "procs", (/lbounds2D(3)/), (/ubounds2D(3)/), procs)
        call write_int_dset_mpi_hdf5_1D(file_id, "refinement_status", (/lbounds2D(3)/), (/ubounds2D(3)/), refinement_status)
        call write_int_dset_mpi_hdf5_1D(file_id, "lgt_ids", (/lbounds2D(3)/), (/ubounds2D(3)/), lgt_ids)
    endif


    ! add additional annotations
    call write_attribute(file_id, "blocks", "version", (/20200902/)) ! this is used to distinguish wabbit file formats
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
