! ********************************
! WABBIT
! --------------------------------
!
! synchronize ghosts nodes
!
! name: synchronize_internal_ghosts.f90
! date: 26.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine synchronize_ghosts()

    use mpi
    use module_params
    use module_blocks
    use module_interpolation

    implicit none

    ! communication list and plan
    integer(kind=ik), dimension(:,:), allocatable   :: com_list
    integer(kind=ik), dimension(2000, 2)            :: com_plan

    integer(kind=ik)                                :: i, k, rank, ierr, n_proc, n_com, allocate_error

    ! allocate com list, maximal size = max number blocks (all blocks on maxlevel) * 8 neighbors
    allocate( com_list( blocks_params%number_max_blocks*16 , 7), stat=allocate_error )

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, n_proc, ierr)

    ! reset com-list, com_plan
    com_list = -99
    com_plan = -99

    ! create com_list
    call create_com_list(com_list, n_com, blocks_params%number_max_blocks*16)

    ! sort com_list, create com_plan
    call sort_com_list(com_list, com_plan, n_proc, n_com, blocks_params%number_max_blocks*16)

!    ! test output
!    i = 1
!        if (rank == 0) then
!            do while ( com_list(i,1) /= -99 )
!                write(*,'(a,i3.1,a,i3.1,a,i3.1,a,i3.1,a,i3.1,a,i3.1,a)', advance='no') "com-id: ", com_list(i,1), ", rank: ", com_list(i,2), ", neighbor-rank: ", com_list(i,3), ", block: ", com_list(i,4), ", neighbor-block: ", com_list(i,5), ", dir: "
!                write(*, '(a)', advance='no') blocks(com_list(i,4))%neighbor_dir(com_list(i,6))
!                write(*, '(a)', advance='no') " "
!                write(*, '(a)', advance='no') blocks(com_list(i,5))%neighbor_dir(com_list(i,7))
!                write(*,*)
!                i = i + 1
!            end do
!            write(*,'(80("#"))')
!    end if
!    i = 1
!        if (rank == 0) then
!            do while ( com_plan(i,1) /= -99 )
!                write(*,'(a,i3.1,x,i3.1)') "com-plan: ", com_plan(i, 1), com_plan(i, 2)
!                i = i + 1
!            end do
!            write(*,'(80("#"))')
!    end if
!    call MPI_Barrier(MPI_COMM_WORLD, ierr)
!    stop

! reset all ghost nodes for testing
do i = 1, blocks_params%number_max_blocks_data
    if ( blocks_data(i)%block_id /= -1 ) then
        blocks_data(i)%data_fields(1)%data_(1:blocks_params%number_ghost_nodes, :) = 9.0e9_rk
        blocks_data(i)%data_fields(1)%data_(blocks_params%number_ghost_nodes+blocks_params%size_block+1:2*blocks_params%number_ghost_nodes+blocks_params%size_block, :) = 9.0e9_rk
        blocks_data(i)%data_fields(1)%data_(:, 1:blocks_params%number_ghost_nodes) = 9.0e9_rk
        blocks_data(i)%data_fields(1)%data_(:, blocks_params%number_ghost_nodes+blocks_params%size_block+1:2*blocks_params%number_ghost_nodes+blocks_params%size_block) = 9.0e9_rk
    end if
end do

    ! synchronize ghost nodes
    i       = 1
    k       = i
    ! loop over com_plan
    do while ( com_plan(i, 1) /= -99 )

        if ( (com_list(k, 2) == rank) .or. (com_list(k, 3) == rank) ) then

            ! proc has to send/receive data
            call send_receive_data( k, com_plan(i, 1), com_plan(i, 2), com_list, blocks_params%number_max_blocks*8*2 )

            ! next step in com_plan
            if (com_plan(i, 1) == 1) then
                ! internal com
                k = k + com_plan(i, 2)
                i = i + 1
            else
                ! external com
                k = k + 2*com_plan(i, 2)
                i = i + 2
            end if

        else

            ! nothing to do, go to next communication
            if (com_plan(i, 1) == 1) then
                ! internal com
                k = k + com_plan(i, 2)
                i = i + 1
            else
                ! external com
                k = k + 2*com_plan(i, 2)
                i = i + 2
            end if

        end if

    end do

    deallocate( com_list, stat=allocate_error )

end subroutine synchronize_ghosts
