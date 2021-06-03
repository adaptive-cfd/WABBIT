!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name proc_to_lgt_data_start_id.f90
!> \version 0.4
!> \author msr
!
!> \brief return start index on light data corresponding to proc rank
!
subroutine proc_to_lgt_data_start_id( lgt_start, rank, N )
    use module_params

    implicit none

    !> light data start index
    integer(kind=ik), intent(out)       :: lgt_start

    integer(kind=ik), intent(in)        :: rank

    !> number of blocks per proc
    integer(kind=ik), intent(in)        :: N

    !> Note zero-based MPIRANK so we add the 1 here:
    lgt_start = rank*N + 1
end subroutine proc_to_lgt_data_start_id
