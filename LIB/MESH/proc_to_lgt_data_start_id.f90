!> \brief return start index on light data corresponding to proc rank

subroutine proc_to_lgt_data_start_id( lgt_start, rank, N )
    use module_params

    implicit none

    integer(kind=ik), intent(out)       :: lgt_start      !> light data start index
    integer(kind=ik), intent(in)        :: rank
    integer(kind=ik), intent(in)        :: N              !> number of blocks per proc

    !> Note zero-based MPIRANK so we add the 1 here:
    lgt_start = rank*N + 1
    
end subroutine proc_to_lgt_data_start_id
