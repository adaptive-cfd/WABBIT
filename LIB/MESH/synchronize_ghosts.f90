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
    integer(kind=ik), dimension(200, 7) :: com_list
    integer(kind=ik), dimension(200, 2) :: com_plan

    integer(kind=ik)                    :: i, k, rank, ierr, n_proc, n_com

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, n_proc, ierr)

    ! reset com-list, com_plan
    com_list = -99
    com_plan = -99

    ! create com_list
    call create_com_list(com_list, n_com)

    ! sort com_list, create com_plan
    call sort_com_list(com_list, com_plan, n_proc, n_com)

    ! synchronize ghost nodes
    i       = 1
    k       = i
    ! loop over com_plan
    do while ( com_plan(i, 1) /= -99 )

        if ( (com_list(k, 2) == rank) .or. (com_list(k, 3) == rank) ) then

            ! proc has to send/receive data
            call send_receive_data( k, com_plan(i, 1), com_plan(i, 2), com_list )
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

end subroutine synchronize_ghosts
