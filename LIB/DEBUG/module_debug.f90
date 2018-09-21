
!===========================================================================
!> Module for all debug subroutines
!-------------------------------------------------------------------------!!!!
!> \details
!>
!> \version 0.4
!> \author msr
!!
!!
!> \date 29.11.16 - create
!! \date 30.04.19 - check if debug arrays are allocated
!!
!!
!! \todo  union with debug data structure
! *****************************************************************************
module module_debug

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! I ave added the moduöe treelib here because I want the redundant node check to
    ! spit out the numerical treecode of blocks.
    use module_treelib
    ! global parameters
    use module_params
    use module_interpolation
!---------------------------------------------------------------------------------------------
! variables
    implicit none

    PRIVATE

    !**********************************************************************************************
    ! These are the important routines that are visible to other modules:
    !**********************************************************************************************
    PUBLIC :: write_debug_times, toc, timing_next_timestep, summarize_profiling, write_block_distribution, &
              write_neighbors, write_lgt_data
    !*********************************************************************************************
    !> global user defined debug structure
    type type_debug

        ! computing time measurement array
        ! row number: id corresponding to names list
        ! column 1: number of subroutine calls for one time loop
        ! column 2: sum (time) of all subroutine calls for one time loop
        ! column 3: number of subroutine calls over complete program
        ! column 4: sum (time) of all subroutine calls over complete program
        real(kind=rk), dimension(:,:), allocatable          :: comp_time

        ! names of time measurements
        ! row number: id
        ! column: name
        character(len=40), dimension(:), allocatable        :: name_comp_time

    end type type_debug

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (type_debug), save                                 :: debug
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


contains

    ! lgt_block synchronization
    include "check_lgt_block_synchronization.f90"


!========================================================================================
    !> \brief write time measurements
    subroutine write_debug_times( iteration, params )
        implicit none
        !-----------------------------------------------------------
        integer(kind=ik), intent(in)     :: iteration !< iteration
        type (type_params), intent(in)   :: params !< user defined parameter structure
        !------------------------------------------------------------
        ! process rank
        integer(kind=ik)      :: k,file_size
        integer(kind=ik),save ::counter =0
        ! file name
        character(len=80)     :: fname

        write( fname,'(i5.5, "times.dat")') params%rank

        ! filesize in byte
        INQUIRE(FILE=fname, SIZE=file_size)
        ! if file is larger then 100MB or this function is called for the first time
        ! we replace the old file by a new file, otherwise we append the iterations
        if ( file_size/1e6> 100 .or. counter == 0 ) then
          open(unit=99,file=fname, status='replace')
        else
          open(unit=99,file=fname,  status='old', position='append', action='write')
        end if
        ! write file header
        write(99,'(80("_"))')
        write(99, '(42x, "calls", 2x, "sum", 5x, "time", 6x, "sum")', advance='no')
        write(99,*)

        ! write times
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! write name
            write(99, '(a)', advance='no') debug%name_comp_time(k)
            ! write number of calls
            write(99, '(2x,i5)', advance='no') int(debug%comp_time(k,1))
            ! write global number of calls
            write(99, '(2x,i7)', advance='no') int(debug%comp_time(k,3))
            ! write time
            write(99, '(2x,f12.6)', advance='no') debug%comp_time(k,2)
            ! write global time
            write(99, '(2x,f12.6)', advance='no') debug%comp_time(k,4)
            ! next line
            write(99,*)
            ! loop variable
            k = k + 1
        end do

        write(99,'(80("-"))')
        write(99, '("iteration: ", i7)', advance='no') iteration
        write(99,*)
        ! close file
        close(unit=99)

    end subroutine write_debug_times
!========================================================================================


!========================================================================================
    !> \brief write current block distribution to file
    subroutine write_block_distribution( dist_list, params )
        implicit none

       !---------------------------------------------------------------
        integer(kind=ik), intent(in)        :: dist_list(:)   !< iteration
        type (type_params), intent(in)      :: params         !< user defined parameter structure
        !---------------------------------------------------------------

        ! file existence variable
        logical                             :: file_exists
        ! file IO error variable, counter, and process rank
        integer(kind=ik)                    :: io_error, k,rank

        rank = params%rank
        ! check file existence, if not create file
        inquire(file="block_dist.dat", exist=file_exists)

        if ( rank == 0 ) then
            if (file_exists) then
                ! open for append
                open(unit=99,file="block_dist.dat",status='old', position="append", action='write', iostat=io_error)
            else
                ! first opening
                open(unit=99,file="block_dist.dat",status='new',action='write', iostat=io_error)
            end if
        end if

        ! write data
        if (rank == 0) then
            do k = 1, size(dist_list)
                write(99, '(i3,1x)', advance='no') dist_list(k)
            end do
            ! next line
            write(99,*)
            ! close file
            close(unit=99)
        end if

    end subroutine write_block_distribution
!========================================================================================



!========================================================================================
    !> allocates the arrays for profiling subroutines and parts of the code
    subroutine allocate_init_debugging(params)
      implicit none
        type (type_params), intent(in)          :: params

      ! note: fix size of time measurements array
      if ( params%debug ) then

          ! allocate array for time measurements - data
          allocate( debug%comp_time( 150, 4 )  )
          ! reset times
          debug%comp_time = 0.0_rk
          ! allocate array for time measurements - names
          allocate( debug%name_comp_time( 150 )  )
          ! reset names
          debug%name_comp_time = "---"
      end if

    end subroutine allocate_init_debugging
!========================================================================================



!========================================================================================
    !> For a given NAME, increase the function call counter by one and store the
    !> elapsed time in the global arrays.
    subroutine toc( params, name, t_elapsed_this, call_counter )
        implicit none
        type (type_params), intent(in) :: params
        character(len=*), intent(in) :: name
        real(kind=rk), intent(in) :: t_elapsed_this
        integer, optional, intent(in) :: call_counter

        integer :: k

        ! write time
        if ( params%debug ) then

            ! check if allocate_init_debbuging was called before
            if (.not. allocated(debug%name_comp_time)) then
                call allocate_init_debugging(params)
            endif

            ! find free or corresponding line
            k = 1
            do while ( debug%name_comp_time(k) /= "---" )
                ! entry for current subroutine exists
                if ( debug%name_comp_time(k) == name ) exit
                k = k + 1
            end do
            ! write time
            debug%name_comp_time(k) = name
            if (present(call_counter)) then
                debug%comp_time(k, 1)   = debug%comp_time(k, 1) + real( call_counter, kind=rk)
            else
                debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1.0_rk
            endif
            debug%comp_time(k, 2)   = debug%comp_time(k, 2) + t_elapsed_this
        end if

    end subroutine toc
!========================================================================================



!========================================================================================
    !> at the end of a time step, we increase the total counters/timers for all measurements
    !> by what has been done in the last time step, then we flush the current timing to disk.
    subroutine timing_next_timestep( params, iteration )
        implicit none
        type (type_params), intent(in) :: params
        integer(kind=ik), intent(in) :: iteration

        ! debug info
        if ( params%debug ) then
            ! sum and reset times and calls
            debug%comp_time(:,3) = debug%comp_time(:,3) + debug%comp_time(:,1)
            debug%comp_time(:,4) = debug%comp_time(:,4) + debug%comp_time(:,2)
            ! write debug infos to file
            call write_debug_times( iteration, params )
            ! reset loop values
            debug%comp_time(:,1) = 0.0_rk
            debug%comp_time(:,2) = 0.0_rk
        end if

    end subroutine timing_next_timestep
!========================================================================================

!========================================================================================
    !> This function summarizes the profile of the Simulation.
    !> It should be called on the end of the program, when the statistics of
    !> the profiled functions is large.
    !> \details
    !> The function displays the total sum of the cpu time spend
    !> in the profiled part of your program and its standard deviation in a tabel.
    subroutine summarize_profiling( params, comm )
        implicit none
        !---------------------------------------
        type (type_params), intent(in) :: params
        integer         , intent(in)   :: comm        !< MPI communicator
        !---------------------------------------
        integer :: rank,k,number_procs,ierr

        call MPI_Comm_rank(comm, rank, ierr)
        call MPI_Comm_size(comm, number_procs, ierr)
        ! debug info output
        if ( params%debug ) then
            ! sum times
            debug%comp_time(:,2) = 0.0_rk
            call MPI_Allreduce(debug%comp_time(:,4), debug%comp_time(:,2), size(debug%comp_time,1), &
                                MPI_REAL8, MPI_SUM,  comm, ierr)
            ! MPI Barrier before program ends
            call MPI_Barrier(comm, ierr)

            ! average times
            debug%comp_time(:,2) = debug%comp_time(:,2) / number_procs
            ! standard deviation
            debug%comp_time(:,3) = 0.0_rk
            debug%comp_time(:,4) = (debug%comp_time(:,4) - debug%comp_time(:,2))**2.0_rk
            call MPI_Allreduce( debug%comp_time(:,4), debug%comp_time(:,3), &
                                size(debug%comp_time,1), MPI_REAL8, MPI_SUM, comm, ierr)
            ! MPI Barrier before program ends
            call MPI_Barrier(comm, ierr)

            if (number_procs == 1) then
                debug%comp_time(:,3) = 0.0_rk
            else
                debug%comp_time(:,3) = sqrt(debug%comp_time(:,3) / ( number_procs - 1 ))
            end if

            ! output
            if (rank==0) then
                write(*,'(80("_"))')
                write(*, '("time (average value +- standard deviation) :")')
                k = 1
                do while ( debug%name_comp_time(k) /= "---" )
                    ! write name
                    write(*, '(a)', advance='no') debug%name_comp_time(k)
                    ! write average time
                    write(*, '(2x,f12.3)', advance='no') debug%comp_time(k,2)
                    ! write standard deviation
                    write(*, '(2x,f12.3)', advance='no') debug%comp_time(k,3)
                    ! next line
                    write(*,*)
                    ! loop variable
                    k = k + 1
                end do
                write(*,'(80("_"))')
                write(*, '("sum: ", 2x,f12.3)', advance='yes') sum(debug%comp_time(:,2))
            end if
        end if

    end subroutine summarize_profiling
!========================================================================================



end module module_debug
