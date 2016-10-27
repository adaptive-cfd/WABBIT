! ********************************
! WABBIT
! --------------------------------
!
! main program, time loop
!
! name: main.f90
! date: 25.10.2016
! author: msr
! version: 0.3
!
! ********************************

program main

    use mpi
    use module_params
    use module_blocks

    implicit none

    ! time loop variables
    real(kind=rk) 	                    :: time
    integer(kind=ik)	                :: iteration, active_blocks
    ! cpu time variables
    real(kind=rk)                       :: t0, t1
    ! MPI variables
    integer(kind=ik)                    :: ierr, rank, n_proc

    ! initialize local variables
    time          = 0.0_rk
    iteration     = 0
    active_blocks = 0

    ! init mpi
    call MPI_Init(ierr)

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, n_proc, ierr)

    ! cpu time start
    call cpu_time(t0)

    ! initializing data
    call init_data()

    ! create block tree
    call matrix_to_block_tree()

    ! update neighbor relations
    call update_neighbors()

    ! output MPI status
    if (rank==0) then
        write(*, '("MPI: using ", i7, " processes")') n_proc
        write(*,'(80("-"))')
    end if

    ! save start field to disk
    call save_data(iteration, time)

    ! main time loop
    do while ( time < params%time_max )

        iteration = iteration + 1

        ! refine every block to create the safety zone
!        if (blocks_params%adapt_mesh) call refine_everywhere()

        ! update the neighbor relations
        call update_neighbors()

        ! advance in time
        call time_step_RK4(time)

        ! adapt the mesh
!        if (blocks_params%adapt_mesh) call adapt_mesh()

        ! write data to disk
        if (modulo(iteration, params%write_freq) == 0) then
          call save_data(iteration, time)
        endif

        ! output on screen
        if (rank==0) then
            call block_count(active_blocks)
            write(*, '("iteration=",i5,3x," time=",f10.6,3x," active blocks=",i7)') iteration, time, active_blocks
            write(*,'(80("-"))')
        end if

   end do

    ! save end field to disk
    call save_data(iteration, time)

    ! cpu time calculation and writing
    call cpu_time(t1)
    if (rank==0) then
        write(*,'(a,f10.6)') "cpu-time = ", t1-t0
    end if

    ! end mpi
    call MPI_Finalize(ierr)

end program main
