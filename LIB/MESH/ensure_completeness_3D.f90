! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: ensure_completeness_3D.f90
! version: 0.5
! author: msr
!
! sets refinement status to -2 for all sister blocks, if coarsening is possible
!
! input:    - light data array
! output:   - light data array
!
! = log ======================================================================================
!
! 03/02/17 - create
!
! ********************************************************************************************

subroutine ensure_completeness_3D( params, lgt_block, lgt_active, lgt_n )

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
    integer(kind=ik)                    :: k, l, i, lgt_id

    ! sister ids
    integer(kind=ik)                    :: id(7)

    ! treecode variable
    integer(kind=ik)                    :: treecode(params%max_treelevel)
    ! block level
    integer(kind=ik)                    :: level

    ! exists variable
    logical                             :: exists

    ! process rank
    integer(kind=ik)                    :: rank

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    max_treelevel = params%max_treelevel

    i  = 0
    id = -1

!---------------------------------------------------------------------------------------------
! main body

    ! set MPI parameter
    rank         = params%rank

    ! start time
    sub_t0 = MPI_Wtime()

    ! loop over all active blocks
    do k = 1, lgt_n

        ! if the block want to coarsen, check refinement status of sister blocks
        if ( lgt_block( lgt_active(k), max_treelevel+2 ) == -1) then

            ! data from current block
            treecode = lgt_block( lgt_active(k), 1:max_treelevel )
            level    = lgt_block( lgt_active(k), max_treelevel+1 )

            ! reset id list
            id = -1

            ! get sister id's
            i = 0
            do l = 1, 8
                ! sister treecode differs only on last element
                if ( lgt_block( lgt_active(k), level ) /= l-1) then

                    i               = i + 1
                    treecode(level) = l-1
                    ! find block id, use exists subroutine
                    call does_block_exist(treecode, lgt_block, max_treelevel, exists, lgt_id, lgt_active, lgt_n)
                    ! block exists
                    if (exists) then
                        id(i) = lgt_id
                    else
                        ! sister does not exists, nothing to do
                    end if

                end if
            end do

            ! if all sisters exists
            if ( (id(1)>0) .and. (id(2)>0) .and. (id(3)>0) .and. (id(4)>0) .and. (id(5)>0) .and. (id(6)>0) .and. (id(7)>0) ) then

                ! if all sister blocks want to coarsen, then coarsening is allowed
                if ( ( lgt_block( id(1), max_treelevel+2 ) == -1) &
                        .and. ( lgt_block( id(2), max_treelevel+2 ) == -1) &
                        .and. ( lgt_block( id(3), max_treelevel+2 ) == -1) &
                        .and. ( lgt_block( id(4), max_treelevel+2 ) == -1) &
                        .and. ( lgt_block( id(5), max_treelevel+2 ) == -1) &
                        .and. ( lgt_block( id(6), max_treelevel+2 ) == -1) &
                        .and. ( lgt_block( id(7), max_treelevel+2 ) == -1) ) then

                    lgt_block( lgt_active(k), max_treelevel+2 )     = -2
                    lgt_block( id(1), max_treelevel+2 )             = -2
                    lgt_block( id(2), max_treelevel+2 )             = -2
                    lgt_block( id(3), max_treelevel+2 )             = -2
                    lgt_block( id(4), max_treelevel+2 )             = -2
                    lgt_block( id(5), max_treelevel+2 )             = -2
                    lgt_block( id(6), max_treelevel+2 )             = -2
                    lgt_block( id(7), max_treelevel+2 )             = -2
                end if

            end if

        end if

    end do

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

end subroutine ensure_completeness_3D
