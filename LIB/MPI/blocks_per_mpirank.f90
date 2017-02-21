subroutine blocks_per_mpirank( params, actual_blocks_per_proc, hvy_n )
  implicit none
  ! user defined parameter structure
  type (type_params), intent(in)      :: params
  integer, intent(in) :: hvy_n
  integer, dimension(0:params%number_procs-1), intent(out) :: actual_blocks_per_proc
  integer, dimension(0:params%number_procs-1) :: tmp_actual_blocks_per_proc

  integer :: error, myrank

  myrank = params%rank

  tmp_actual_blocks_per_proc = -1
  tmp_actual_blocks_per_proc(myrank) = hvy_n
  ! call MPI_Allreduce(my_block_list, lgt_block, size(lgt_block,1)*size(lgt_block,2), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(tmp_actual_blocks_per_proc, actual_blocks_per_proc, params%number_procs, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, error)
end subroutine
