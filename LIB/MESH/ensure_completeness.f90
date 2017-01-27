! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: ensure_completeness.f90
! version: 0.4
! author: msr
!
! sets refinement status to -2 for all sister blocks, if coarsening is possible
!
! input:    - light data array
! output:   - light data array
!
! = log ======================================================================================
!
! 10/11/16 - switch to v0.4
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

    ! light data array for blocks with coarsening status
    ! max size: number of active blocks
    !integer(kind=ik)                    :: coarse_blocks(lgt_n, size(lgt_block,1))
    ! number of blocks with coarsening status, loop variable
    !integer(kind=ik)                    :: n_coarse_blocks


    ! max treelevel
    integer(kind=ik)                    :: max_treelevel

    ! loop variables
    integer(kind=ik)                    :: k, l, i, lgt_id

    ! sister ids
    integer(kind=ik)                    :: id(3)

    ! treecode variable
    integer(kind=ik)                    :: treecode(params%max_treelevel)
    ! block level
    integer(kind=ik)                    :: level

    ! exists variable
    logical                             :: exists

    ! MPI error variable
    integer(kind=ik)                    :: ierr
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

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! start time
    sub_t0 = MPI_Wtime()

!    ! first: create list of blocks with refinement status -1
!    !-----------------------------------------------------------------------------------------
!    ! reset counter for blocks with coarsening status
!    n_coarse_blocks = 0
!    ! loop over all active blocks
!    do k = 1, lgt_n
!        ! block want to coarsen
!        if ( lgt_block( lgt_active(k), max_treelevel+2 ) == -1) then
!
!            ! increase counter
!            n_coarse_blocks = n_coarse_blocks + 1
!
!            ! save light data
!            ! light data id
!            coarse_blocks( n_coarse_blocks, 1) = lgt_active(k)
!            ! treecode
!            coarse_blocks( n_coarse_blocks, 2:max_treelevel+1) = lgt_block( lgt_active(k), 1:max_treelevel )
!
!        end if
!    end do
!
!    ! second: search sister blocks in coarse block list, if all 4 sister blocks are present,
!    ! set status in light data array
!    !-----------------------------------------------------------------------------------------
!
!    ! loop over blocks with coarsening status
!    do k = 1, n_coarse_blocks
!
!        ! check light data id
!        ! note: id -1 means, blocks was deleted in previous loop
!        if ( coarse_blocks( k, 1) /= -1 ) then
!
!            ! search sister blocks
!            !---------------------
!            ! reset id list
!            id = -1
!            ! reset sister count
!            i = 0
!            ! loop over remaining coarse block list
!            do l = k, n_coarse_blocks
!
!
!
!            end do
!
!        end if
!
!    end do

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
            do l = 1, 4
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
            if ( (id(1)>0) .and. (id(2)>0) .and. (id(3)>0) ) then

                ! if all sister blocks want to coarsen, then coarsening is allowed
                if ( ( lgt_block( id(1), max_treelevel+2 ) == -1) .and. ( lgt_block( id(2), max_treelevel+2 ) == -1) .and. ( lgt_block( id(3), max_treelevel+2 ) == -1) ) then
                    lgt_block( lgt_active(k), max_treelevel+2 )      = -2
                    lgt_block( id(1), max_treelevel+2 )  = -2
                    lgt_block( id(2), max_treelevel+2 )  = -2
                    lgt_block( id(3), max_treelevel+2 )  = -2
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

end subroutine ensure_completeness
