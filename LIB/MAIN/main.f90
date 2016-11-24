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
! 23/11/16 - use computing time array for simple performance tests
! ********************************************************************************************

program main

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! debug module
    use module_debug
    ! init data module
    use module_init
    ! mesh manipulation subroutines
    use module_mesh

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
    !                   -> column(max_treelevel+1): treecode length = mesh level
    !                   -> column(max_treelevel+2):   refinement status (-1..coarsen / 0...no change / +1...refine)
    integer(kind=ik), allocatable       :: lgt_block(:, :)

    !                   -> dim 1: x coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 2: y coord   ( 1:number_block_nodes+2*number_ghost_nodes )
    !                   -> dim 3: data type ( field_1, 2:number_data_fields+1, data_old, k1, k2, k3,
    !                                       k4 [for runge kutta] )
    ! heavy data array  -> dim 4: block id  ( 1:number_blocks )
    !           field_1 (to save mixed data):   line 1: x coordinates
    !                                           line 2: y coordinates
    real(kind=rk), allocatable          :: hvy_block(:, :, :, :)

    ! neighbor array (heavy data) -> number_lines = number_blocks * 16 (...different neighbor relations:
    ! '__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW' )
    !         saved data -> -1 ... no neighbor
    !                    -> light data line number (id)
    integer(kind=ik), allocatable       :: hvy_neighbor(:)

    ! list of active blocks (light data)
    integer(kind=ik), allocatable       :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik)                    :: lgt_n

    ! list of active blocks (heavy data)
    integer(kind=ik), allocatable       :: hvy_active(:)
    ! number of active blocks (heavy data)
    integer(kind=ik)                    :: hvy_n

    ! time loop variables
    real(kind=rk)                       :: time
    integer(kind=ik)                    :: iteration

    ! loop variable
    integer(kind=ik)                    :: k, l

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! init time loop
    time          = 0.0_rk
    iteration     = 0

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
        open (15, file='load_balancing.t', status='replace')
        close(15)
        open (15, file='blocks_per_rank.t', status='replace')
        close(15)
    end if

    ! initializing data
    call init_data( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active )

    ! create lists of active blocks (light and heavy data)
    call create_lgt_active_list( lgt_block, lgt_active, lgt_n )
    call create_hvy_active_list( lgt_block, hvy_active, hvy_n )

    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

!    ! save start data
!    call save_data( iteration, time, params, block_list, block_data, neighbor_list )
!
!    ! main time loop
!    do while ( time < params%time_max )
!
!        iteration = iteration + 1
!
!        ! start time
!        sub_t0 = MPI_Wtime()
!        ! refine every block to create the safety zone
!        if ( params%adapt_mesh ) call refine_everywhere( params, block_list, block_data )
!        ! end time
!        sub_t1 = MPI_Wtime()
!        if ( params%debug ) then
!            debug%name_comp_time(2) = "refine_everywhere"
!            debug%comp_time(rank+1, 2) = sub_t1 - sub_t0
!        end if
!
!        ! start time
!        sub_t0 = MPI_Wtime()
!        ! update neighbor relations
!        call update_neighbors( block_list, neighbor_list, params%number_blocks, params%max_treelevel )
!        ! end time
!        sub_t1 = MPI_Wtime()
!        if ( params%debug ) then
!            debug%name_comp_time(3) = "update_neighbors"
!            debug%comp_time(rank+1, 3) = sub_t1 - sub_t0
!        end if
!
!        ! start time
!        sub_t0 = MPI_Wtime()
!        ! advance in time
!        call time_step_RK4( time, params, block_list, block_data, neighbor_list )
!        ! end time
!        sub_t1 = MPI_Wtime()
!        if ( params%debug ) then
!            debug%name_comp_time(4) = "time_step_RK4"
!            debug%comp_time(rank+1, 4) = sub_t1 - sub_t0
!        end if
!
!        ! adapt the mesh
!        if ( params%adapt_mesh ) call adapt_mesh( params, block_list, block_data, neighbor_list, debug )
!
!        ! output on screen
!        if (rank==0) then
!            write(*,'(80("-"))')
!            call block_count(block_list, block_number)
!            write(*, '("RUN: iteration=",i5,3x," time=",f10.6,3x," active blocks=",i7)') iteration, time, block_number
!
!        end if
!
!        ! write data to disk
!        if (modulo(iteration, params%write_freq) == 0) then
!          call save_data( iteration, time, params, block_list, block_data, neighbor_list )
!        endif
!
        ! output computing time for every proc
        if ( params%debug ) then
            do k = 1, number_procs
                if ( k-1 == rank ) then
                    write(*,'(80("."))')
                    write(*,'("RUN: computing time details for rank = ",i3)')  rank
                    l = 1
                    ! time array is used
                    do while ( debug%name_comp_time(l) /= "---" )
                        write(*,'(a, " : ", f10.6, " s")')  debug%name_comp_time(l), debug%comp_time(k, l)
                        l = l + 1
                    end do

                end if
                call MPI_Barrier(MPI_COMM_WORLD, ierr)
            end do
        end if

!    end do
!
!    ! save end field to disk
!    call save_data( iteration, time, params, block_list, block_data, neighbor_list )

    ! MPI Barrier before program ends
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! computing time output on screen
    call cpu_time(t1)
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("END: cpu-time = ",f10.6, " s")')  t1-t0
    end if

    ! end mpi
    call MPI_Finalize(ierr)

end program main
