!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name create_external_com_list.f90
!> \version 0.5
!> \author msr
!
!> \brief create com matrix and com list for external communications
!!  check mesh level and synch stages
!
!>
!! input:    
!!           - params, light
!!
!! output:   
!!           - com matrix
!!           - com lists for external synchronization
!!
!! = log ======================================================================================
!! \n
!! 19/05/17 - create
!
! ********************************************************************************************

subroutine create_external_com_list(  params, lgt_block, hvy_neighbor, hvy_active, hvy_n, com_matrix, com_lists, synch_stage )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    !> communication lists: 
    !!       - dim 1: list elements 
    !!       - dim 2: columns
    !!                      - 1   rank of sender process
    !!                      - 2   rank of receiver process
    !!                      - 3   sender block heavy data id
    !!                      - 4   receiver block heavy data id
    !!                      - 5   sender block neighborhood to receiver (dirs id)
    !!                      - 6   difference between sender-receiver level
    !!       - dim 3: receiver proc rank
    integer(kind=ik), intent(inout)     :: com_lists(:, :, :)

    !> communications matrix: \n
    !! count the number of communications between procs
    !! row/column number encodes process rank + 1
    integer(kind=ik), intent(inout)     :: com_matrix(:,:)

    !> stage of synchronization
    integer(kind=ik), intent(in)        :: synch_stage

    ! loop variables
    integer(kind=ik)                    :: k, N, i, lgt_id, hvy_id, neighbor_num, neighborhood

    ! grid parameter
    integer(kind=ik)                    :: g, Bs

    ! process rank
    integer(kind=ik)                    :: rank, neighbor_rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! neighbor light data id
    integer(kind=ik)                    :: neighbor_light_id
    ! difference between current block and neighbor block level
    integer(kind=ik)                    :: level_diff

    ! receiver lists: receiver_pos [position in rank and count list],
    ! receiver_rank [proc rank list], receiver_count [number of communications to receiver]
    ! receiver_N [number of neighbor procs]
    integer(kind=ik), allocatable       :: receiver_pos(:), receiver_rank(:), receiver_count(:)
    integer(kind=ik)                    :: receiver_N

    ! synch switch
    logical                             :: synch

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! set MPI parameter
    rank            = params%rank
    number_procs    = params%number_procs

    ! receiver lists
    allocate( receiver_pos( number_procs) )
    allocate( receiver_rank( number_procs) )
    allocate( receiver_count( number_procs) )

    ! reset
    receiver_pos    =  0
    receiver_rank   = -1
    receiver_count  =  0
    receiver_N      =  0

    synch = .false.

!---------------------------------------------------------------------------------------------
! main body

    ! set loop number for 2D/3D case
    neighbor_num = size(hvy_neighbor, 2)

    ! loop over active heavy data
    do k = 1, hvy_n

        ! loop over all neighbors
        do i = 1, neighbor_num
            ! neighbor exists
            if ( hvy_neighbor( hvy_active(k), i ) /= -1 ) then

                ! neighbor light data id
                neighbor_light_id = hvy_neighbor( hvy_active(k), i )
                ! calculate light id
                call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )
                ! calculate the difference between block levels
                ! here: sender (me) - receiver (neighbor)
                ! if +1 -> sender on finer level
                level_diff = lgt_block( lgt_id, params%max_treelevel+1 ) - lgt_block( neighbor_light_id, params%max_treelevel+1 )
                !level_diff = lgt_block( neighbor_light_id, params%max_treelevel+1 ) - lgt_block( lgt_id, params%max_treelevel+1 )

                ! proof if neighbor internal or external
                call lgt_id_to_proc_rank( neighbor_rank, neighbor_light_id, N )

                ! check synch stage
                ! stage 1: level +1
                ! stage 2: level 0
                ! stage 3: level -1
                ! stage 4: special
                synch = .false.
                if ( (synch_stage == 1) .and. (level_diff == 1) ) then
                    synch = .true.
                end if
                if ( (synch_stage == 2) .and. (level_diff == 0) ) then
                    synch = .true.
                end if
                if ( (synch_stage == 3) .and. (level_diff == -1) ) then
                    synch = .true.
                end if

                if ( (synch_stage == 4) .and. (level_diff == 0) ) then
                    ! neighborhood NE
                    if ( i == 5 ) then
                        if ( (hvy_neighbor( hvy_active(k), 9) /= -1) .or. (hvy_neighbor( hvy_active(k), 13) /= -1) ) then
                            synch = .true.
                        end if
                    end if
                    ! neighborhood NW
                    if ( i == 6 ) then
                        if ( (hvy_neighbor( hvy_active(k), 10) /= -1) .or. (hvy_neighbor( hvy_active(k), 15) /= -1) ) then
                            synch = .true.
                        end if
                    end if
                    ! neighborhood SE
                    if ( i == 7 ) then
                        if ( (hvy_neighbor( hvy_active(k), 11) /= -1) .or. (hvy_neighbor( hvy_active(k), 14) /= -1) ) then
                            synch = .true.
                        end if
                    end if
                    ! neighborhood SW
                    if ( i == 8 ) then
                        if ( (hvy_neighbor( hvy_active(k), 12) /= -1) .or. (hvy_neighbor( hvy_active(k), 16) /= -1) ) then
                            synch = .true.
                        end if
                    end if
                end if

                if (synch) then

                    neighborhood = i !0

                    if ( rank == neighbor_rank ) then
                        ! internal neighbor
                        ! write communications matrix
                        !com_matrix(rank+1, rank+1) = com_matrix(rank+1, rank+1) + 1

                    else
                        ! external neighbor

                        ! neighbor heavy id
                        call lgt_id_to_hvy_id( hvy_id, neighbor_light_id, neighbor_rank, N )

                        ! check neighbor proc rank
                        if ( receiver_pos(neighbor_rank+1) == 0 ) then

                            ! first communication with neighbor proc
                            ! -------------------------------------------
                            ! set list position, increase number of neighbor procs by 1
                            receiver_N                      = receiver_N + 1
                            ! save list pos
                            receiver_pos(neighbor_rank+1)   = receiver_N
                            ! save neighbor rank
                            receiver_rank(receiver_N)       = neighbor_rank
                            ! count communications - here: first one
                            receiver_count(receiver_N)      = 1

                            ! external neighbor -> new com_lists entry (first entry)
                            com_lists( 1 , 1, neighbor_rank+1)  = rank
                            com_lists( 1 , 2, neighbor_rank+1)  = neighbor_rank
                            com_lists( 1 , 3, neighbor_rank+1)  = hvy_active(k)
                            com_lists( 1 , 4, neighbor_rank+1)  = hvy_id
                            com_lists( 1 , 5, neighbor_rank+1)  = neighborhood !i
                            com_lists( 1 , 6, neighbor_rank+1)  = level_diff

                        else

                            ! additional communication with neighbor proc
                            ! -------------------------------------------
                            ! count communications - +1
                            receiver_count( receiver_pos(neighbor_rank+1) )  = receiver_count( receiver_pos(neighbor_rank+1) ) + 1

                            ! external neighbor -> new com_lists entry
                            com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 1, neighbor_rank+1)  = rank
                            com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 2, neighbor_rank+1)  = neighbor_rank
                            com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 3, neighbor_rank+1)  = hvy_active(k)
                            com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 4, neighbor_rank+1)  = hvy_id
                            com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 5, neighbor_rank+1)  = neighborhood !i
                            com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 6, neighbor_rank+1)  = level_diff

                        end if

                    end if

                end if

            end if
        end do

    end do

    ! write my com matrix, loop over number of receiver procs, write counted communications
    do k = 1, receiver_N
        ! write matrix
        !com_matrix( rank+1, receiver_rank(k)+1 ) = receiver_count( receiver_pos( receiver_rank(k)+1 ) )
        com_matrix( receiver_rank(k)+1, rank+1 ) = receiver_count( receiver_pos( receiver_rank(k)+1 ) )
    end do

    ! clean up
    deallocate( receiver_pos )
    deallocate( receiver_rank )
    deallocate( receiver_count )

end subroutine create_external_com_list
