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
module module_timing

!---------------------------------------------------------------------------------------------
! modules

    use mpi

!---------------------------------------------------------------------------------------------
! variables
    implicit none

    PRIVATE

    !**********************************************************************************************
    ! These are the important routines that are visible to other modules:
    !**********************************************************************************************
    PUBLIC :: write_times, toc, timing_next_timestep, summarize_profiling, setup_indiv_timings
    !*********************************************************************************************

    ! precision in this module
    integer, parameter, public      :: rk_timing=selected_real_kind(8)
    
    !> global user defined arrays

    ! computing time measurement array
    ! row number: id corresponding to names list
    ! column 1: number of subroutine calls for one time loop
    ! column 2: sum (time) of all subroutine calls for one time loop
    ! column 3: number of subroutine calls over complete program
    ! column 4: sum (time) of all subroutine calls over complete program
    real(kind=rk_timing), dimension(:,:), allocatable          :: comp_time

    ! names of time measurements
    ! row number: id
    ! column: name
    character(len=40), dimension(:), allocatable        :: name_comp_time

    logical :: write_individual_timings = .false.
    logical :: setup_completed = .false.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !type (type_timing), save :: timing
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!========================================================================================
    !> \brief reads flag of write_individual_timings from params
    subroutine  setup_indiv_timings( write_indiv_timings )

      implicit none
      ! value from params
      logical, intent(in) :: write_indiv_timings

      if ( .not. setup_completed ) then
          write_individual_timings = write_indiv_timings
      endif

      setup_completed = .true.

    end subroutine setup_indiv_timings

!========================================================================================
    !> \brief write time measurements
    subroutine write_times( iteration )
        implicit none
        !-----------------------------------------------------------
        integer, intent(in)     :: iteration !< iteration
        !------------------------------------------------------------
        ! process rank
        integer      :: k,file_size, proc_rank, ierr
        integer,save ::counter =0
        ! file name
        character(len=80)     :: fname

        call MPI_comm_rank(MPI_COMM_WORLD, proc_rank, ierr)
        write( fname,'(i5.5, "times.dat")') proc_rank

        ! filesize in byte
        INQUIRE(FILE=fname, SIZE=file_size)
        ! if file is larger than 100MB or this function is called for the first time
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
        do while (  name_comp_time(k) /= "---" )
            ! write name
            write(99, '(a)', advance='no')  name_comp_time(k)
            ! write number of calls
            write(99, '(2x,i5)', advance='no') int( comp_time(k,1))
            ! write global number of calls
            write(99, '(2x,i7)', advance='no') int( comp_time(k,3))
            ! write time
            write(99, '(2x,f12.6)', advance='no')  comp_time(k,2)
            ! write global time
            write(99, '(2x,f12.6)', advance='no')  comp_time(k,4)
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

    end subroutine write_times
!========================================================================================



!========================================================================================
    !> allocates the arrays for profiling subroutines and parts of the code
    subroutine allocate_init_timing

      implicit none

      ! note: fix size of time measurements array
      ! allocate array for time measurements - data
      allocate(  comp_time( 150, 4 )  )
      ! reset times
      comp_time = 0.0_rk_timing
      ! allocate array for time measurements - names
      allocate(  name_comp_time( 150 )  )
      ! reset names
      name_comp_time = "---"

    end subroutine allocate_init_timing
!========================================================================================



!========================================================================================
    !> For a given NAME, increase the function call counter by one and store the
    !> elapsed time in the global arrays.
    subroutine toc( name, t_elapsed_this, call_counter )
        implicit none
        character(len=*), intent(in) :: name
        real(kind=rk_timing), intent(in) :: t_elapsed_this
        integer, optional, intent(in) :: call_counter

        integer :: k

        ! write time
        ! check if allocate_init_debbuging was called before
        if (.not. allocated( name_comp_time)) then
            call allocate_init_timing
        endif

        ! find free or corresponding line
        k = 1
        do while (  name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if (  name_comp_time(k) == name ) exit
            k = k + 1
        end do
        ! write time
         name_comp_time(k) = name
        if (present(call_counter)) then
             comp_time(k, 1)   =  comp_time(k, 1) + real( call_counter, kind=rk_timing)
        else
             comp_time(k, 1)   =  comp_time(k, 1) + 1.0_rk_timing
        endif
         comp_time(k, 2)   =  comp_time(k, 2) + t_elapsed_this

    end subroutine toc
!========================================================================================



!========================================================================================
    !> at the end of a time step, we increase the total counters/timers for all measurements
    !> by what has been done in the last time step, then we flush the current timing to disk.
    subroutine timing_next_timestep( iteration )! params , ...
        implicit none
        !type (type_params), intent(in) :: params
        integer, intent(in) :: iteration

        ! debug info
        ! sum and reset times and calls
         comp_time(:,3) =  comp_time(:,3) +  comp_time(:,1)
         comp_time(:,4) =  comp_time(:,4) +  comp_time(:,2)
        ! write debug infos to file
        if (write_individual_timings) call write_times( iteration )
        ! reset loop values
         comp_time(:,1) = 0.0_rk_timing
         comp_time(:,2) = 0.0_rk_timing

    end subroutine timing_next_timestep
!========================================================================================

!========================================================================================
    !> This function summarizes the profile of the Simulation.
    !> It should be called on the end of the program, when the statistics of
    !> the profiled functions is large.
    !> \details
    !> The function displays the total sum of the cpu time spend
    !> in the profiled part of your program and its standard deviation in a tabel.
    subroutine summarize_profiling( comm )
        implicit none
        !---------------------------------------
        integer         , intent(in)   :: comm        !< MPI communicator
        !---------------------------------------
        integer :: rank,k,number_procs,ierr

        call MPI_Comm_rank(comm, rank, ierr)
        call MPI_Comm_size(comm, number_procs, ierr)
        ! debug info output
          ! sum times
           comp_time(:,2) = 0.0_rk_timing
          call MPI_Allreduce( comp_time(:,4),  comp_time(:,2), size( comp_time,1), &
                              MPI_REAL8, MPI_SUM,  comm, ierr)
          ! MPI Barrier before program ends
          call MPI_Barrier(comm, ierr)

          ! average times
           comp_time(:,2) =  comp_time(:,2) / number_procs
          ! standard deviation
           comp_time(:,3) = 0.0_rk_timing
           comp_time(:,4) = ( comp_time(:,4) -  comp_time(:,2))**2.0_rk_timing
          call MPI_Allreduce(  comp_time(:,4),  comp_time(:,3), &
                              size( comp_time,1), MPI_REAL8, MPI_SUM, comm, ierr)
          ! MPI Barrier before program ends
          call MPI_Barrier(comm, ierr)

          if (number_procs == 1) then
               comp_time(:,3) = 0.0_rk_timing
          else
               comp_time(:,3) = sqrt( comp_time(:,3) / ( number_procs - 1 ))
          end if

          ! output
          if (rank==0) then
              write(*,'(80("_"))')
              write(*, '("time (average value +- standard deviation) :")')
              k = 1
              do while (  name_comp_time(k) /= "---" )
                  ! write name
                  write(*, '(a)', advance='no')  name_comp_time(k)
                  ! write average time
                  write(*, '(2x,f12.3)', advance='no')  comp_time(k,2)
                  ! write standard deviation
                  write(*, '(2x,f12.3)', advance='no')  comp_time(k,3)
                  ! next line
                  write(*,*)
                  ! loop variable
                  k = k + 1
              end do
              write(*,'(80("_"))')
              write(*, '("sum: ", 2x,f12.3)', advance='yes') sum( comp_time(:,2))
          end if

    end subroutine summarize_profiling
!========================================================================================



end module module_timing
