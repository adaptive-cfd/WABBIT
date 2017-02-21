

subroutine write_field( fname, time, iteration, dF, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_n)

  ! modules
  use hdf5
  use module_mpi
  ! variables

  implicit none
  character(len=*), intent(in)        :: fname
  ! time loop parameters
  real(kind=rk), intent(in)           :: time
  integer(kind=ik), intent(in)        :: iteration
  ! datafield number
  integer(kind=ik), intent(in)        :: dF
  ! user defined parameter structure
  type (type_params), intent(in)      :: params
  ! light data array
  integer(kind=ik), intent(in)        :: lgt_block(:, :)
  ! heavy data array - block data
  real(kind=rk), intent(in)           :: hvy_block(:, :, :, :, :)
  ! heavy data array - neifghbor data
  integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
  ! list of active blocks (light data)
  integer(kind=ik), intent(in)        :: lgt_active(:)
  ! number of active blocks (light data)
  integer(kind=ik), intent(in)        :: lgt_n, hvy_n

  ! process rank
  integer(kind=ik)                    :: myrank, lgt_rank
  ! loop variable
  integer(kind=ik)                    :: k, i, hvy_id, l
  ! grid parameter
  integer(kind=ik)                    :: Bs, g
  real(kind=rk), allocatable :: myblockbuffer3D(:,:,:,:),myblockbuffer2D(:,:,:)
  real(kind=rk), allocatable :: coords_origin(:,:), coords_spacing(:,:)
  integer(hid_t):: file_id
  integer :: error
  integer,dimension(1:4) :: ubounds3D,lbounds3D
  integer,dimension(1:3) :: ubounds2D,lbounds2D
  integer, dimension(:), allocatable :: actual_blocks_per_proc

  ! set MPI parameters
  myrank = params%rank
  Bs = params%number_block_nodes
  g  = params%number_ghost_nodes

  if (myrank == 0) then
    write(*,'("IO: writing data for time = ", f15.8," file = ",A," active blocks=",i5)') time, trim(adjustl(fname)), lgt_n
  endif

  ! to know our position in the last index of the 4D output array, we need to
  ! know how many blocks all procs have
  allocate(actual_blocks_per_proc(0:params%number_procs-1))
  call blocks_per_mpirank( params, actual_blocks_per_proc, hvy_n )

  ! fill blocks buffer (we cannot use the bvy_block array as it is not contiguous, i.e.
  ! it may contain holes)
  if ( params%threeD_case ) then
    ! 3D
    allocate (myblockbuffer3D(1:bs,1:bs,1:bs,1:hvy_n), coords_spacing(1:3,1:hvy_n), coords_origin(1:3,1:hvy_n))
    ! tell the hdf5 wrapper what part of the global [bs x bs x bs x n_active]
    ! array we hold, so that all CPU can write to the same file simultaneously
    ! (note zero-based offset):
    lbounds3D = (/1,1,1,sum(actual_blocks_per_proc(0:myrank-1))+1/) - 1
    ubounds3D = (/bs-1,bs-1,bs-1,lbounds3D(4)+hvy_n-1/)
  else
    ! 2D
    allocate (myblockbuffer2D(1:bs,1:bs,1:hvy_n), coords_spacing(1:2,1:hvy_n), coords_origin(1:2,1:hvy_n))
    ! tell the hdf5 wrapper what part of the global [bs x bs x bs x n_active]
    ! array we hold, so that all CPU can write to the same file simultaneously
    ! (note zero-based offset):
    lbounds2D = (/1,1,sum(actual_blocks_per_proc(0:myrank-1))+1/) - 1
    ubounds2D = (/bs-1,bs-1,lbounds2D(3)+hvy_n-1/)
  endif

  l = 1
  ! loop over all active light block IDs, check if it is mine, if so, copy the block to the buffer
  do k = 1, lgt_n
      ! calculate proc rank from light data line number
      call lgt_id_to_proc_rank( lgt_rank, lgt_active(k), params%number_blocks )
      ! calculate heavy block id corresponding to light id
      call lgt_id_to_hvy_id( hvy_id, lgt_active(k), myrank, params%number_blocks )
      ! if I own this block, I copy it to the buffer.
      ! also extract block coordinate origin and spacing
      if (lgt_rank == myrank) then
        if ( params%threeD_case ) then
          ! 3D
          myblockbuffer3D(:,:,:,l) = hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, hvy_id)
          coords_origin(1:3,l)   = hvy_block( 1:3, 1, 1, 1, hvy_id)
          coords_spacing(1:3,l)  = abs(hvy_block( 1:3, 2, 1, 1, hvy_id) - hvy_block( 1:3, 1, 1, 1, hvy_id) )
        else
          ! 2D
          myblockbuffer2D(:,:,l) = hvy_block( g+1:Bs+g, g+1:Bs+g, 1, dF, hvy_id)
          coords_origin(1:2,l)   = hvy_block( 1:2, 1, 1, 1, hvy_id)
          coords_spacing(1:2,l)  = abs(hvy_block( 1:2, 2, 1, 1, hvy_id) - hvy_block( 1:2, 1, 1, 1, hvy_id) )
        endif
        l = l + 1
      endif
  end do

  ! open the file
  call open_file_hdf5( trim(adjustl(fname)), file_id, .true.)

  ! write heavy block data to disk
  if ( params%threeD_case ) then
    ! 3D data case
    call write_dset_mpi_hdf5_4D(file_id, "blocks", lbounds3D, ubounds3D, myblockbuffer3D)
    call write_attribute(file_id, "blocks", "domain-size", (/params%Lx, params%Ly, params%Lz/))
    call write_dset_mpi_hdf5_2D(file_id, "coords_origin", (/0,lbounds3D(4)/), (/2,ubounds3D(4)/), coords_origin)
    call write_dset_mpi_hdf5_2D(file_id, "coords_spacing", (/0,lbounds3D(4)/), (/2,ubounds3D(4)/), coords_spacing)
  else
    ! 2D data case
    call write_dset_mpi_hdf5_3D(file_id, "blocks", lbounds2D, ubounds2D, myblockbuffer2D)
    call write_attribute(file_id, "blocks", "domain-size", (/params%Lx, params%Ly/))
    call write_dset_mpi_hdf5_2D(file_id, "coords_origin", (/0,lbounds2D(3)/), (/1,ubounds2D(3)/), coords_origin)
    call write_dset_mpi_hdf5_2D(file_id, "coords_spacing", (/0,lbounds2D(3)/), (/1,ubounds2D(3)/), coords_spacing)
  endif

  ! add aditional annotations
  call write_attribute(file_id, "blocks", "time", (/time/))
  call write_attribute(file_id, "blocks", "iteration", (/iteration/))

  ! close file and HDF5 library
  call close_file_hdf5(file_id)

  if (allocated(myblockbuffer3D)) deallocate(myblockbuffer3D)
  if (allocated(myblockbuffer2D)) deallocate(myblockbuffer2D)

  deallocate(actual_blocks_per_proc, coords_origin, coords_spacing)
end subroutine write_field
