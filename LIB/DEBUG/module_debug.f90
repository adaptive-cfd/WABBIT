
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

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    ! lgt_block synchronization
    include "check_lgt_block_synchronization.f90"

    ! check ghost nodes
    include "check_ghost_nodes.f90"

    ! write time measurements
    include "write_debug_times.f90"

    ! write block distribution
    include "write_block_distribution.f90"

    ! use this to initalize the memory for debugging:
    include "allocate_init_debugging.f90"

    ! For a given NAME, increase the function call counter by one and store the
    ! elapsed time in the global arrays.
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
                call abort(5946,'ERROR [module_debug]: debug arrays are not allocated yet')
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


    ! at the end of a time step, we increase the total counters/timers for all measurements
    ! by what has been done in the last time step, then we flush the current timing to disk.
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

end module module_debug
