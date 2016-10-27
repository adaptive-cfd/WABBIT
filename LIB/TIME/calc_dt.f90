! ********************************
! WABBIT
! --------------------------------
!
! calculate timestep
!
! name: calc_dt.f90
! date: 26.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine calc_dt(dt)

    use mpi
    use module_params
    use module_blocks

    implicit none

    real(kind=rk), intent(inout) 	:: dt

    real(kind=rk)                   :: dx, dt_loc
    integer                         :: k, N, rank, ierr, n_procs, tag
    integer                         :: status(MPI_status_size)

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, n_procs, ierr)

    tag = 0

    N = blocks_params%number_max_blocks_data

    ! first: calculate local time step
    ! minimal spacing
    dx = 9.0e9_rk

    ! loop over all blocks
    do k = 1, N
        if ( blocks_data(k)%block_id /= -1 ) then
            dx = min(dx, blocks_data(k)%dx )
        end if
    end do

    ! time step
    dt_loc = params%CFL * dx / norm2(params%u0)

    ! second: calculate minimum over all local time steps
    ! proc 0 collect all steps and broadcast them after min operation
    if (rank==0) then
        dt = dt_loc
    end if

    do k = 1, n_procs-1

        ! collect local time steps
        if (rank==0) then
            call MPI_Recv(dt_loc, 1, MPI_REAL8, k, tag, MPI_COMM_WORLD, status, ierr)
            dt = min(dt, dt_loc)
        elseif (rank==k) then
            call MPI_Send( dt_loc, 1, MPI_REAL8, 0, tag, MPI_COMM_WORLD, ierr)
        end if

    end do

    ! broadcast min time step
    call MPI_Bcast(dt, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

end subroutine calc_dt
