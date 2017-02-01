! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: synchronize_internal_nodes.f90
! version: 0.5
! author: msr
!
! synchronize internal ghosts nodes, create com matrix and com list for external communication
!
! input:    - params, light and heavy data
! output:   - heavy data
!           - com matrix
!           - com lists for external synchronization
!
! = log ======================================================================================
!
! 12/01/17 - create from old synchronize ghost routine
! 31/01/17 - switch to 3D, v0.5
!
! ********************************************************************************************

subroutine synchronize_internal_nodes(  params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_matrix, com_lists )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    ! heavy data array - neifghbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    ! list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    ! number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! communication lists:
    ! dim 1: list elements
    ! dim 2: columns
    !                       1   rank of sender process
    !                       2   rank of receiver process
    !                       3   sender block heavy data id
    !                       4   receiver block heavy data id
    !                       5   sender block neighborhood to receiver (dirs id)
    !                       6   difference between sender-receiver level
    ! dim 3: receiver proc rank
    integer(kind=ik), intent(inout)     :: com_lists(:, :, :)

    ! communications matrix:
    ! count the number of communications between procs
    ! row/column number encodes process rank + 1
    integer(kind=ik), intent(inout)     :: com_matrix(:,:)

    ! loop variables
    integer(kind=ik)                    :: k, N, i, lgt_id, hvy_id, neighbor_num

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

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

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
    allocate( receiver_pos( number_procs), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    allocate( receiver_rank( number_procs), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    allocate( receiver_count( number_procs), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    ! reset
    receiver_pos    =  0
    receiver_rank   = -1
    receiver_count  =  0
    receiver_N      =  0

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
                level_diff = lgt_block( lgt_id, params%max_treelevel+1 ) - lgt_block( neighbor_light_id, params%max_treelevel+1 )

                ! proof if neighbor internal or external
                call lgt_id_to_proc_rank( neighbor_rank, neighbor_light_id, N )

                if ( rank == neighbor_rank ) then
                    ! calculate internal heavy id
                    call lgt_id_to_hvy_id( hvy_id, neighbor_light_id, rank, N )

                    ! internal neighbor -> copy ghost nodes
                    if ( params%threeD_case ) then
                        ! 3D:
                        call copy_ghost_nodes_3D( params, hvy_block, hvy_active(k), hvy_id, i, level_diff )
                    else
                        ! 2D:
                        call copy_ghost_nodes_2D( params, hvy_block(:, :, 1, :, :), hvy_active(k), hvy_id, i, level_diff )
                    end if

                    ! write communications matrix
                    com_matrix(rank+1, rank+1) = com_matrix(rank+1, rank+1) + 1

                else
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
                        com_lists( 1 , 5, neighbor_rank+1)  = i
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
                        com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 5, neighbor_rank+1)  = i
                        com_lists( receiver_count( receiver_pos(neighbor_rank+1) ) , 6, neighbor_rank+1)  = level_diff

                    end if

                end if

            end if
        end do

    end do

    ! write my com matrix, loop over number of receiver procs, write counted communications
    do k = 1, receiver_N
        ! write matrix
        com_matrix( rank+1, receiver_rank(k)+1 ) = receiver_count( receiver_pos( receiver_rank(k)+1 ) )
    end do

    ! clean up
    deallocate( receiver_pos, stat=allocate_error )
    deallocate( receiver_rank, stat=allocate_error )
    deallocate( receiver_count, stat=allocate_error )

end subroutine synchronize_internal_nodes
