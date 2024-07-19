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
    PUBLIC :: toc, summarize_profiling, reset_all_timings
    !*********************************************************************************************

    ! precision of reals in this module. we use this separated precision in order for
    ! the timing module to be completely independent of the rest of the code (and hence to be
    ! reusable in other projects)
    integer, parameter :: dp = selected_real_kind(8)

    ! we provide at most this many slots for timing:
    integer, PARAMETER :: MAX_TIMING_SLOTS = 250

    ! array of time measurements and call counters
    real(kind=dp), dimension(:,:), allocatable :: comp_time

    ! each timing slot gets a name and number so it can be easily identified
    ! name is to identify it later by name, they will only be checked for uniqueness for DEV
    ! num is to uniquely address one timing, they should not overlap for different names
    character(len=80), dimension(:), allocatable :: name_comp_time
    integer, dimension(:, :), allocatable :: num_comp_time

contains

    !> Resets all timings to 0
    subroutine reset_all_timings()
        implicit none
        if (allocated(comp_time)) comp_time = 0.0_dp
    end subroutine

    ! List of codes that I somewhat arbitrarily set just to have some kind of sorting for toc calls
    !    9-  15 TOPLEVEL
    !   20-  24 timestep
    !   30-  33 RHS_WRAPPER
    !   50- 100 Fundamental functions
    !   50-  55    createActiveSortedLists
    !        59    updateMetadata_tree
    !   60-  63    synchronize_lgt_data
    !   70-  74    xfer_block_data
    !   80-  85    sync ghosts
    !   90-  92    balanceLoad_tree
    !  100- 200 adapt functions
    !  100- 115    adapt_tree
    !  120- 124    coarseningIndicator
    !  130- 131    ensureGradedness
    !  140- 146    refine_tree
    !  150- 158    coarseExtensionUpdate_tree
    !  200-1000 Other
    !  250- 258    forest
    !  350- 354    module_MOR
    ! 1000-XXXX Miscellaneous
    ! 1001-1002    Commented block based toc in coarsening indicator
    ! 1010-1015    Old coarse extension

    !> For a given NAME, increase the function call counter by one and store the
    !> elapsed time in the global arrays.
    subroutine toc( name, num, t_elapsed_this, call_counter )
        implicit none
        character(len=*), intent(in)  :: name
        integer, intent(in)  :: num
        real(kind=dp), intent(in)     :: t_elapsed_this
        integer, optional, intent(in) :: call_counter

        integer :: k, ierr

        ! check if allocate_init_debbuging was called before
        if (.not. allocated( name_comp_time)) then
            ! note: fix size of time measurements array
            ! allocate array for time measurements - data
            allocate(  comp_time( MAX_TIMING_SLOTS, 3 )  ) ! first is call counter, second sum, third max
            ! reset times
            comp_time = 0.0_dp
            ! allocate array for time measurements - names and nums
            allocate(  name_comp_time( MAX_TIMING_SLOTS )  )
            allocate(  num_comp_time( MAX_TIMING_SLOTS, 2 )  )  ! second slot is for sorting later on
            ! reset names and num
            name_comp_time = "---"
            num_comp_time(:, 1) = -1
        endif

        ! find a free or the corresponding slot in the array:
        k = 1
        do while (  num_comp_time(k, 1) /= -1 )
            ! entry for current subroutine exists
            if (  num_comp_time(k, 1) == num ) exit
            k = k + 1
        end do

        ! DEV: check if name is also identical. String comparison so not done for deployment
#ifdef DEV
        if (name_comp_time(k) /= name .and. name_comp_time(k) /= "---") then
            write(*, '("DEV WARNING: Measurement ", a, " and ", a, " have conflicting unique number ", i0)') name_comp_time(k), name, num
        endif
#endif

        ! write time
        name_comp_time(k) = name
        num_comp_time(k, 1) = num
        if (present(call_counter)) then
            ! increase by the number given in argument
            comp_time(k, 1)   =  comp_time(k, 1) + real( call_counter, kind=dp)
        else
            ! increase by one
            comp_time(k, 1)   =  comp_time(k, 1) + 1.0_dp
        endif
        comp_time(k, 2)   =  comp_time(k, 2) + t_elapsed_this
        comp_time(k, 3)   =  max(comp_time(k, 3), t_elapsed_this)
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
        integer :: rank, k_timings, k_sorted, number_procs,ierr
        real(kind=dp), dimension(:), allocatable :: avg, max, std

        call MPI_Comm_rank(comm, rank, ierr)
        call MPI_Comm_size(comm, number_procs, ierr)

        if (.not. allocated(comp_time)) then
            if (rank==0) write(*,*) "Timing information not available (no toc used?)"
            return
        endif

        allocate(avg(1:MAX_TIMING_SLOTS))
        allocate(max(1:MAX_TIMING_SLOTS))
        allocate(std(1:MAX_TIMING_SLOTS))

        ! sum times (over all mpi processes) for all slots
        call MPI_Allreduce( comp_time(:,2), avg, MAX_TIMING_SLOTS, MPI_REAL8, MPI_SUM,  comm, ierr)
        ! average times (over all mpi processes) for all slots
        avg = avg / dble(number_procs)

        ! max
        call MPI_Allreduce( comp_time(:,3), max, MAX_TIMING_SLOTS, MPI_REAL8, MPI_MAX,  comm, ierr)

        ! standard deviation
        std = (comp_time(:,2) -  avg)**2.0_dp
        call MPI_Allreduce( MPI_IN_PLACE, std, MAX_TIMING_SLOTS, MPI_REAL8, MPI_SUM, comm, ierr)

        if (number_procs == 1) then
            std = 0.0_dp
        else
            std = sqrt( std / dble(number_procs - 1 ))
        end if
        ! make it relative to avg and to percent
        do k_timings = 1, MAX_TIMING_SLOTS
            ! check if we divide by 0 and skip that as it is a zero-entry
            if (avg(k_timings) > 1e-6) then
                std(k_timings) = std(k_timings) / avg(k_timings) * 100
            ! else  ! if the call is really small then just disable averaging
            !     std(k_timings) = -1
            endif
        enddo

        ! write indices as unique ids into second entry so that we can retrieve it for the other arrays
        ! CONTINUE HERE
        do k_timings = 1, MAX_TIMING_SLOTS
            num_comp_time(k_timings, 2) = k_timings
        enddo

        ! sort array after the number we associated to it so that it looks nicely
        call interchange_sort_timing(num_comp_time, 1, MAX_TIMING_SLOTS)

        ! output, with DEV we output the unique ID as well to have an overview
        if (rank==0) then
            write(*,'(108("_"))')
#ifdef DEV
            write(*, '(A60, 2(A14), 3(A10))') "Timing name" // repeat(" ", 30), "AVG (s)", "MAX (s)", "STD (%)", "COUNT", "Unique ID"
#else
            write(*, '(A60, 2(A14), 2(A10))') "Timing name" // repeat(" ", 30), "AVG (s)", "MAX (s)", "STD (%)", "COUNT"
#endif
            write(*,'(108("_"))')
            do k_timings = 1, MAX_TIMING_SLOTS
                ! get the index in the arrays from the sorted list
                k_sorted = num_comp_time(k_timings, 2)
                if (num_comp_time(k_timings, 1) /= -1) then
#ifdef DEV
                write(*,'(A60, 2(2x, f12.3), 2x, f8.3, 2(2x, i8))') name_comp_time(k_sorted), avg(k_sorted), max(k_sorted), std(k_sorted), int(comp_time(k_sorted, 1)), num_comp_time(k_timings, 1)
#else
                write(*,'(A60, 2(2x, f12.3), 2x, f8.3, 2x, i8)') name_comp_time(k_sorted), avg(k_sorted), max(k_sorted), std(k_sorted), int(comp_time(k_sorted, 1))
#endif
                endif
            end do
            write(*,'(108("_"))')
        end if

        ! MPI Barrier to be sure to see the above write statements
        call MPI_Barrier(comm, ierr)

        deallocate(avg)
        deallocate(std)

    end subroutine summarize_profiling


    !> \brief Interchange algorithm sorting after one number, copy from quicksort but here to sort out the importing things
    !> \details This algorithm sorts the array a from position first to position last.
    subroutine interchange_sort_timing(a, left_end, right_end)
        implicit none
        integer, intent(inout) ::  a(:,:)
        integer :: left_end, right_end

        integer :: i, j
        integer, dimension(size(a, 2)) :: temp

        do i = left_end, right_end - 1
            do j = i+1, right_end
            if (a(j, 1) < a(i, 1)) then
                temp = a(i,:)
                a(i,:) = a(j,:)
                a(j,:) = temp
                end if
            end do
        end do

    end subroutine interchange_sort_timing

end module module_timing
