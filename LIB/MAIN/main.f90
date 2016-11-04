! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: main.f90
! version: 0.4
! author: msr
!
! main program, init all data, start time loop, output on screen during program run
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4
! ********************************************************************************************

program main

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0, t1

    ! user defined parameter structure
    type (type_params)                  :: params

    ! light data array  -> line number = ( 1 + proc_rank ) * heavy_data_line_number
    !                   -> column(1:max_treelevel): block treecode, treecode -1 => block is inactive
    !                   -> column(max_treelevel+1):   refinement status (-1..coarsen / 0...no change / +1...refine)
    integer(kind=ik), allocatable       :: block_list(:, :)

    ! heavy data array  -> dim 1: block id  ( 1:number_blocks )
    !                   -> dim 2: x coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 3: y coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 4: data type ( 1:number_data_fields, data_old, k1, k2, k3,
    !                                       k4 [for runge kutta] )
    real(kind=rk), allocatable          :: block_data(:, :, :, :)

!---------------------------------------------------------------------------------------------
! interfaces

    interface
        subroutine init_data(params, block_list, block_data)
            use module_params
            type (type_params), intent(out)             :: params
            integer(kind=ik), allocatable, intent(out)  :: block_list(:, :)
            real(kind=rk), allocatable, intent(out)     :: block_data(:, :, :, :)
        end subroutine init_data
    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! init mpi
    call MPI_Init(ierr)
    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    ! determinate process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)

    ! cpu start time
    call cpu_time(t0)

    ! output MPI status
    if (rank==0) then
        write(*,'(80("_"))')
        write(*, '("MPI: using ", i5, " processes")') number_procs
    end if

    ! initializing data
    call init_data( params, block_list, block_data )

    ! cpu end time and output on screen
    call cpu_time(t1)
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("END: cpu-time = ",f10.6, " s")')  t1-t0
    end if

    ! end mpi
    call MPI_Finalize(ierr)


!    use module_params
!    use module_blocks
!
!    implicit none
!
!    ! time loop variables
!    real(kind=rk) 	                    :: time
!    integer(kind=ik)	                :: iteration, active_blocks
!    ! cpu time variables
!    real(kind=rk)                       :: t0, t1
!    ! MPI variables
!    integer(kind=ik)                    :: ierr, rank, n_proc
!
!    ! initialize local variables
!    time          = 0.0_rk
!    iteration     = 0
!    active_blocks = 0
!
!    ! init mpi
!    call MPI_Init(ierr)
!    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
!    call MPI_Comm_size(MPI_COMM_WORLD, n_proc, ierr)
!
!    ! cpu time start
!    call cpu_time(t0)
!
!    ! initializing data
!    call init_data()
!
!    ! create block tree
!    call matrix_to_block_tree()
!
!    ! update neighbor relations
!    call update_neighbors()
!
!    ! output MPI status
!    if (rank==0) then
!        write(*, '("MPI: using ", i7, " processes")') n_proc
!        write(*,'(80("-"))')
!    end if
!
!    ! save start field to disk
!    call save_data(iteration, time)
!
!    ! main time loop
!    do while ( time < params%time_max )
!
!        iteration = iteration + 1
!
!        ! refine every block to create the safety zone
!        if (blocks_params%adapt_mesh) call refine_everywhere()
!
!        ! update the neighbor relations
!        call update_neighbors()
!
!        ! advance in time
!        call time_step_RK4(time)
!
!        ! adapt the mesh
!        if (blocks_params%adapt_mesh) call adapt_mesh()
!
!        ! write data to disk
!        if (modulo(iteration, params%write_freq) == 0) then
!          call save_data(iteration, time)
!        endif
!
!        ! output on screen
!        if (rank==0) then
!            call block_count(active_blocks)
!            write(*, '("iteration=",i5,3x," time=",f10.6,3x," active blocks=",i7)') iteration, time, active_blocks
!            write(*,'(80("-"))')
!        end if
!
!    end do
!
!    ! save end field to disk
!    call save_data(iteration, time)
!
!    ! cpu time calculation and writing
!    call cpu_time(t1)
!    if (rank==0) then
!        write(*,'(a,f10.6)') "cpu-time = ", t1-t0
!    end if
!
!    ! end mpi
!    call MPI_Finalize(ierr)

end program main
