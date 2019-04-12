!===========================================================================
!> Module to perform runtime timings using MPI_WTIME
! *****************************************************************************
module module_timing
    use mpi

    ! variables
    implicit none

    ! everything is private unless explicitly marked public
    PRIVATE

    !**********************************************************************************************
    ! These are the important routines that are visible to other modules:
    !**********************************************************************************************
    PUBLIC :: toc, summarize_profiling
    !*********************************************************************************************

    ! precision of reals in this module. we use this separated precision in order for
    ! the timing module to be completely independent of the rest of the code (and hence to be
    ! reusable in other projects)
    integer, parameter :: dp = selected_real_kind(8)

    ! we provide at most this many slots for timing:
    integer, PARAMETER :: MAX_TIMING_SLOTS = 250

    ! array of time measurements and call counters
    real(kind=dp), dimension(:,:), allocatable :: comp_time

    ! each timing slot gets a name so it can be easily identified
    character(len=100), dimension(:), allocatable :: name_comp_time

contains


    !> For a given NAME, increase the function call counter by one and store the
    !> elapsed time in the global arrays.
    subroutine toc( name, t_elapsed_this, call_counter )
        implicit none
        character(len=*), intent(in)  :: name
        real(kind=dp), intent(in)     :: t_elapsed_this
        integer, optional, intent(in) :: call_counter

        integer :: k

        ! check if allocate_init_debbuging was called before
        if (.not. allocated( name_comp_time)) then
            ! note: fix size of time measurements array
            ! allocate array for time measurements - data
            allocate(  comp_time( MAX_TIMING_SLOTS, 2 )  )
            ! reset times
            comp_time = 0.0_dp
            ! allocate array for time measurements - names
            allocate(  name_comp_time( MAX_TIMING_SLOTS )  )
            ! reset names
            name_comp_time = "---"
        endif

        ! find a free or the corresponding slot in the array:
        k = 1
        do while (  name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if (  name_comp_time(k) == name ) exit
            k = k + 1
        end do

        ! write time
        name_comp_time(k) = name
        if (present(call_counter)) then
            ! increase by the number given in argument
            comp_time(k, 1)   =  comp_time(k, 1) + real( call_counter, kind=dp)
        else
            ! increase by one
            comp_time(k, 1)   =  comp_time(k, 1) + 1.0_dp
        endif
        comp_time(k, 2)   =  comp_time(k, 2) + t_elapsed_this

    end subroutine toc


    !> This function summarizes the profile of the Simulation.
    !> It should be called on the end of the program, when the statistics of
    !> the profiled functions is large.
    !> \details
    !> The function displays the total sum of the cpu time spend
    !> in the profiled part of your program and its standard deviation in a tabel.
    subroutine summarize_profiling( comm )
        implicit none
        !---------------------------------------
        !< MPI communicator
        integer, intent(in)   :: comm
        !---------------------------------------
        integer :: rank,k,number_procs,ierr
        real(kind=dp), dimension(:), allocatable :: avg, std

        call MPI_Comm_rank(comm, rank, ierr)
        call MPI_Comm_size(comm, number_procs, ierr)

        if (.not. allocated(comp_time)) then
            if (rank==0) write(*,*) "Timing information not available (no toc used?)"
            return
        endif

        allocate(avg(1:size(comp_time,1)))
        allocate(std(1:size(comp_time,1)))

        ! sum times (over all mpi processes) for all slots
        call MPI_Allreduce( comp_time(:,2), avg, size(comp_time,1), MPI_REAL8, MPI_SUM,  comm, ierr)

        ! average times (over all mpi processes) for all slots
        avg = avg / dble(number_procs)

        ! standard deviation
        std = (comp_time(:,2) -  avg)**2.0_dp
        call MPI_Allreduce( MPI_IN_PLACE, std, size(comp_time,1), MPI_REAL8, MPI_SUM, comm, ierr)

        if (number_procs == 1) then
            std = 0.0_dp
        else
            std = sqrt( std / dble(number_procs - 1 ))
        end if

        ! output
        if (rank==0) then
            write(*,'(80("_"))')
            write(*, '("time (average value +- standard deviation) :")')
            k = 1
            do while (  name_comp_time(k) /= "---" )
                write(*,'(A100, 2(2x, f12.3))') name_comp_time(k), avg(k), std(k)
                k = k + 1
            end do
            write(*,'(80("_"))')
        end if

        ! MPI Barrier to be sure to see the above write statements
        call MPI_Barrier(comm, ierr)

        deallocate(avg)
        deallocate(std)

    end subroutine summarize_profiling

end module module_timing
