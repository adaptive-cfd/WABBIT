! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: ensure_completeness.f90
! version: 0.4
! author: msr, engels
!
! sets refinement status to -2 for all sister blocks, if coarsening is possible
!
! input:    - light data array
! output:   - light data array
!
! = log ======================================================================================
!
! 10/11/16 - switch to v0.4
! 05/04/17 - works for 2D and 3D data and uses readable find_sisters routine.
! ********************************************************************************************

subroutine ensure_completeness( params, lgt_block, lgt_active, lgt_n )

!---------------------------------------------------------------------------------------------
! modules


!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    ! list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n

    ! max treelevel
    integer(kind=ik)                    :: max_treelevel
    ! loop variables
    integer(kind=ik)                    :: k, l, N_sisters, status
    ! sister ids
    integer(kind=ik), allocatable       :: id(:)
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    max_treelevel = params%max_treelevel
    if (params%threeD_case) then
      N_sisters = 7
    else
      N_sisters = 3
    end if
    ! allocate the array which will hold the light ids of the sisters
    allocate( id(1:N_sisters) )
    id = -1

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! loop over all active blocks
    do k = 1, lgt_n

        ! if the block wants to coarsen, check refinement status of sister blocks
        if ( lgt_block( lgt_active(k), max_treelevel+2 ) == -1) then
            ! find sister IDs of the block we're looking at. If a sister is not found, -1
            ! is returned in the array id
            call find_sisters( params, lgt_active(k), id, lgt_block, lgt_active, lgt_n )

            ! if all sisters exists, then the array should not contain values smaller
            ! zero (-1 would mean not found)
            if ( minval(id) > 0 ) then

                ! now loop over all sisters, check if they also want to coarsen and have status -1
                ! only if all sisters agree to coarsen, they can all be merged into their mother block.
                status = -1
                do l = 1, N_sisters
                  status = max( status, lgt_block(id(l), max_treelevel+2) )
                end do

                ! if all agree, the status is -1, and we can indeed coarsen, set +2
                if ( status == -1 ) then
                  do l = 1, N_sisters
                    lgt_block( id(l), max_treelevel+2 )  = -2
                  end do
                end if

            end if
        end if
    end do

    deallocate( id )

    ! NOTE: this routine runs redundantly on all procs - they all do the same. That means
    ! by consequence that the light data here does not have to be synchronized.


    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "ensure_completeness" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "ensure_completeness"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine ensure_completeness
