subroutine blocks_per_mpirank( params, actual_blocks_per_proc, tree_ID )
    use module_forestMetaData
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    integer, intent(in) :: tree_ID
    integer, dimension(0:params%number_procs-1), intent(out) :: actual_blocks_per_proc

    integer :: error

    actual_blocks_per_proc = -1
    actual_blocks_per_proc(params%rank) = hvy_n(tree_ID)
    call MPI_ALLREDUCE(MPI_IN_PLACE, actual_blocks_per_proc, params%number_procs, MPI_INTEGER, MPI_MAX, WABBIT_COMM, error)
end subroutine
