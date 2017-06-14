!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name check_redundant_nodes.f90
!> \version 0.5
!> \author msr
!
!> \brief check all redundant nodes, if there is a difference between nodes: stop
!> wabbit and write current state to file
!
!>
!! input:    - params, light and heavy data \n
!! output:   - heavy data array
!!
!> \details
!! = log ======================================================================================
!! \n
!! 09/05/17 - create
!
! ********************************************************************************************

subroutine check_redundant_nodes(  params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, stop_status )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! status of nodes check: if true: stops program
    logical, intent(inout)              :: stop_status

    ! loop variables
    integer(kind=ik)                    :: N, i, j

    ! grid parameter
    integer(kind=ik)                    :: g, Bs

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

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
    integer(kind=ik), allocatable       :: com_lists(:, :, :)

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! communications matrix:
    ! count the number of communications between procs
    ! row/column number encodes process rank + 1
    ! com matrix pos: position in send buffer
    integer(kind=ik), allocatable       :: com_matrix(:,:), my_com_matrix(:,:)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! set MPI parameter
    rank         = params%rank
    number_procs = params%number_procs

    ! allocate local com_lists
    allocate( com_lists( N*16, 6, number_procs), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    ! allocate com matrix
    allocate( com_matrix(number_procs, number_procs), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    allocate( my_com_matrix(number_procs, number_procs), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    ! reset com-list, com_plan, com matrix, receiver lists
    com_lists       = -1

    com_matrix      =  0
    my_com_matrix   =  0

    ! reset status
    stop_status = .false.

!---------------------------------------------------------------------------------------------
! main body

    ! ----------------------------------------------------------------------------------------
    ! first: check internal nodes, create com_list for external communications

    ! copy internal nodes and create com_matrix/com_lists for external communications
    call check_internal_nodes( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, my_com_matrix, com_lists, stop_status )

    ! synchronize com matrix
    call MPI_Allreduce(my_com_matrix, com_matrix, number_procs*number_procs, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! symmetrize com matrix - choose max com number
    do i = 1, number_procs
        do j = i+1, number_procs
            com_matrix(i,j) = max( com_matrix(i,j), com_matrix(j,i) )
            com_matrix(j,i) = com_matrix(i,j)
        end do
    end do

    ! call external nodes synchronization
    !call check_external_nodes(  params, hvy_block, com_lists, com_matrix, stop_status )

    ! clean up
    deallocate( com_lists, stat=allocate_error )
    deallocate( com_matrix, stat=allocate_error )
    deallocate( my_com_matrix, stat=allocate_error )

end subroutine check_redundant_nodes

! check internal nodes

subroutine check_internal_nodes(  params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_matrix, com_lists, stop_status )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neifghbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! status of nodes check: if true: stops program
    logical, intent(inout)              :: stop_status

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
                            !call copy_ghost_nodes_3D( params, hvy_block, hvy_active(k), hvy_id, i, level_diff )
                        else
                            ! 2D:
                            call check_ghost_nodes_2D( params, hvy_block(:, :, 1, :, :), hvy_active(k), hvy_id, i, level_diff, stop_status )
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

end subroutine check_internal_nodes

! check_ghost_nodes_2D

subroutine check_ghost_nodes_2D( params, hvy_block, sender_id, receiver_id, neighborhood, level_diff, stop_status)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :)
    !> heavy data id's
    integer(kind=ik), intent(in)                    :: sender_id, receiver_id
    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: level_diff

    ! status of nodes check: if true: stops program
    logical, intent(inout)                          :: stop_status

    ! grid parameter
    integer(kind=ik)                                :: Bs, g
    ! loop variables
    integer(kind=ik)                                :: dF

    ! difference between sender/receiver nodes
    real(kind=rk)                                   :: diff_norm
    ! error threshold
    real(kind=rk)                                   :: eps

!---------------------------------------------------------------------------------------------
! interfaces

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    ! set error threshold
    eps = 1e-12_rk
    !eps = 1e-16_rk

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    select case(neighborhood)
        ! '__N'
        case(1)
            if ( level_diff == 0 ) then
                ! sender/receiver on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( Bs+g, g+1:Bs+g, dF, receiver_id ) - &
                                           hvy_block( g+1, g+1:Bs+g, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (__N)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do
            end if

        ! '__E'
        case(2)
            if ( level_diff == 0 ) then
                ! sender/receiver on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1:Bs+g, g+1, dF, receiver_id ) - &
                                           hvy_block( g+1:Bs+g, Bs+g, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (__E)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do
            end if

        ! '__S'
        case(3)
            if ( level_diff == 0 ) then
                ! sender/receiver on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1, g+1:Bs+g, dF, receiver_id ) - &
                                           hvy_block( Bs+g, g+1:Bs+g, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (__S)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do
            end if

        ! '__W'
        case(4)
            if ( level_diff == 0 ) then
                ! sender/receiver on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1:Bs+g, Bs+g, dF, receiver_id ) - &
                                           hvy_block( g+1:Bs+g, g+1, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (__W)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do
            end if

        ! '_NE'
        case(5)
            ! loop over all datafields
            do dF = 2, params%number_data_fields+1
                diff_norm = sqrt((( hvy_block( Bs+g, g+1, dF, receiver_id ) - &
                                       hvy_block( g+1, Bs+g, dF, sender_id ) ))**2 )
                ! check error
                if ( diff_norm > eps ) then
                    write(*,*) "ERROR: difference in redundant nodes (_NE)"
                    ! save data
                    ! write proc rank into block data only for first error case
                    if (stop_status) then
                        ! do nothing
                    else
                        hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                        hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                        stop_status = .true.
                    end if
                end if
            end do

        ! '_NW'
        case(6)
            ! loop over all datafields
            do dF = 2, params%number_data_fields+1
                diff_norm = sqrt((( hvy_block( Bs+g, Bs+g, dF, receiver_id ) - &
                                       hvy_block( g+1, g+1, dF, sender_id ) ))**2 )
                ! check error
                if ( diff_norm > eps ) then
                    write(*,*) "ERROR: difference in redundant nodes (_NW)"
                    ! save data
                    ! write proc rank into block data only for first error case
                    if (stop_status) then
                        ! do nothing
                    else
                        hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                        hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                        stop_status = .true.
                    end if
                end if
            end do

        ! '_SE'
        case(7)
            ! loop over all datafields
            do dF = 2, params%number_data_fields+1
                diff_norm = sqrt((( hvy_block( g+1, g+1, dF, receiver_id ) - &
                                       hvy_block( Bs+g, Bs+g, dF, sender_id ) ))**2 )
                ! check error
                if ( diff_norm > eps ) then
                    write(*,*) "ERROR: difference in redundant nodes (_SE)"
                    ! save data
                    ! write proc rank into block data only for first error case
                    if (stop_status) then
                        ! do nothing
                    else
                        hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                        hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                        stop_status = .true.
                    end if
                end if
            end do

        ! '_SW'
        case(8)
            ! loop over all datafields
            do dF = 2, params%number_data_fields+1
                diff_norm = sqrt((( hvy_block( g+1, Bs+g, dF, receiver_id ) - &
                                       hvy_block( Bs+g, g+1, dF, sender_id ) ))**2 )
                ! check error
                if ( diff_norm > eps ) then
                    write(*,*) "ERROR: difference in redundant nodes (_SW)"
                    ! save data
                    ! write proc rank into block data only for first error case
                    if (stop_status) then
                        ! do nothing
                    else
                        hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                        hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                        stop_status = .true.
                    end if
                end if
            end do

        ! 'NNE'
        case(9)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( Bs+g, g+1:Bs+g:2, dF, receiver_id ) - &
                                       hvy_block( g+1, (Bs+1)/2+g:Bs+g, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (NNE)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( Bs+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) - &
                                       hvy_block( g+1, g+1:Bs+g:2, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (NNE)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do
            end if

        ! 'NNW'
        case(10)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( Bs+g, g+1:Bs+g:2, dF, receiver_id ) - &
                                       hvy_block( g+1, g+1:(Bs+1)/2+g, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (NNW)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( Bs+g, g+1:g+(Bs+1)/2, dF, receiver_id ) - &
                                       hvy_block( g+1, g+1:Bs+g:2, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (NNW)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do
            end if

        ! 'SSE'
        case(11)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1, g+1:Bs+g:2, dF, receiver_id ) - &
                                       hvy_block( Bs+g, (Bs+1)/2+g:Bs+g, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (SSE)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1, g+(Bs+1)/2:Bs+g, dF, receiver_id ) - &
                                       hvy_block( Bs+g, g+1:Bs+g:2, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (SSE)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do
            end if

        ! 'SSW'
        case(12)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1, g+1:Bs+g:2, dF, receiver_id ) - &
                                       hvy_block( Bs+g, g+1:(Bs+1)/2+g, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (SSW)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1, g+1:g+(Bs+1)/2, dF, receiver_id ) - &
                                       hvy_block( Bs+g, g+1:Bs+g:2, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (SSW)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            end if

        ! 'ENE'
        case(13)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1:Bs+g:2, g+1, dF, receiver_id ) - &
                                       hvy_block( g+1:(Bs+1)/2+g, Bs+g, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (ENE)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1:g+(Bs+1)/2, g+1, dF, receiver_id ) - &
                                       hvy_block( g+1:Bs+g:2, Bs+g, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (ENE)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do
            end if

        ! 'ESE'
        case(14)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1:Bs+g:2, g+1, dF, receiver_id ) - &
                                       hvy_block( g+(Bs+1)/2:Bs+g, Bs+g, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (ESE)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+(Bs+1)/2:Bs+g, g+1, dF, receiver_id ) - &
                                       hvy_block( g+1:Bs+g:2, Bs+g, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (ESE)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do
            end if

        ! 'WNW'
        case(15)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1:Bs+g:2, Bs+g, dF, receiver_id ) - &
                                       hvy_block( g+1:(Bs+1)/2+g, g+1, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (WNW)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1:g+(Bs+1)/2, Bs+g, dF, receiver_id ) - &
                                       hvy_block( g+1:Bs+g:2, g+1, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (WNW)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do
            end if

        ! 'WSW'
        case(16)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1:Bs+g:2, Bs+g, dF, receiver_id ) - &
                                       hvy_block( g+(Bs+1)/2:Bs+g, g+1, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (WNW)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+(Bs+1)/2:Bs+g, Bs+g, dF, receiver_id ) - &
                                       hvy_block( g+1:Bs+g:2, g+1, dF, sender_id ) ))**2 )
                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes (WNW)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, receiver_id ) = real( neighborhood, kind=rk)
                            hvy_block( :, :, dF, sender_id ) = real( neighborhood, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do
            end if

    end select

end subroutine check_ghost_nodes_2D

! --------------------------------------------------------------------------------------

subroutine check_external_nodes(  params, hvy_block, com_lists, com_matrix, stop_status )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)

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
    integer(kind=ik), intent(inout)        :: com_lists(:, :, :)

    ! communications matrix:
    ! count the number of communications between procs
    ! row/column number encodes process rank + 1
    ! com matrix pos: position in send buffer
    integer(kind=ik), intent(inout)        :: com_matrix(:,:)

    integer(kind=ik), allocatable       :: com_matrix_pos(:,:)

    ! loop variables
    integer(kind=ik)                    :: k, N, i, j

    ! grid parameter
    integer(kind=ik)                    :: g, Bs

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! send/receive buffer, integer and real
    integer(kind=ik), allocatable       :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), allocatable          :: real_send_buffer(:,:), real_receive_buffer(:,:)

    ! length of buffer array and column number in buffer, use for readability
    integer(kind=ik)                    :: int_N, real_N, buffer_pos

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! number of communications, number of neighboring procs
    integer(kind=ik)                    :: my_n_com, n_com, n_procs

    ! status of nodes check: if true: stops program
    logical, intent(inout)              :: stop_status

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! set MPI parameter
    rank         = params%rank
    number_procs = params%number_procs

    allocate( com_matrix_pos(number_procs, number_procs), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    com_matrix_pos  =  0

!---------------------------------------------------------------------------------------------
! main body

    ! max number of communications and neighboring procs - use for buffer allocation
    ! every proc loop over com matrix line
    call max_com_num_2( my_n_com, n_procs, com_matrix(rank+1,:), rank )

    ! synchronize max com number, because if not:
    ! RMA buffer displacement is not fixed, so we need synchronization there
    call MPI_Allreduce(my_n_com, n_com, 1, MPI_INTEGER4, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! for proc without neighbors: set n_procs to 1
    ! so we allocate arrays with second dimension=1
    if (n_procs==0) n_procs = 1

    ! next steps only for more than two procs
    if ( number_procs > 1 ) then

        ! ----------------------------------------------------------------------------------------
        ! second: allocate memory for send/receive buffer
        ! buffer size:
        !               number of columns: number of neighboring procs
        !               number of lines:
        !                   int buffer  - max number of communications  * 3 + 1 (length of real buffer)
        !                   real buffer - max number of communications * (Bs+g) * g * number of datafields
        !                   for 3D: real_buffer * Bs
        allocate( int_send_buffer( n_com * 3 + 1, n_procs ), stat=allocate_error )
        !call check_allocation(allocate_error)
        if ( allocate_error /= 0 ) then
            write(*,'(80("_"))')
            write(*,*) "ERROR: memory allocation fails"
            stop
        end if
        allocate( int_receive_buffer( n_com * 3 + 1, n_procs ), stat=allocate_error )
        !call check_allocation(allocate_error)
        if ( allocate_error /= 0 ) then
            write(*,'(80("_"))')
            write(*,*) "ERROR: memory allocation fails"
            stop
        end if

        if ( params%threeD_case ) then
            ! 3D:
            allocate( real_receive_buffer( n_com * (Bs+g) * Bs * params%number_data_fields, n_procs ), stat=allocate_error )
            !call check_allocation(allocate_error)
            if ( allocate_error /= 0 ) then
                write(*,'(80("_"))')
                write(*,*) "ERROR: memory allocation fails"
                stop
            end if
            allocate( real_send_buffer( n_com * (Bs+g) * Bs * params%number_data_fields, n_procs ), stat=allocate_error )
            !call check_allocation(allocate_error)
            if ( allocate_error /= 0 ) then
                write(*,'(80("_"))')
                write(*,*) "ERROR: memory allocation fails"
                stop
            end if
        else
            ! 2D:
            allocate( real_receive_buffer( n_com * (Bs+g) * params%number_data_fields, n_procs ), stat=allocate_error )
            !call check_allocation(allocate_error)
            if ( allocate_error /= 0 ) then
                write(*,'(80("_"))')
                write(*,*) "ERROR: memory allocation fails"
                stop
            end if
            allocate( real_send_buffer( n_com * (Bs+g) * params%number_data_fields, n_procs ), stat=allocate_error )
            !call check_allocation(allocate_error)
            if ( allocate_error /= 0 ) then
                write(*,'(80("_"))')
                write(*,*) "ERROR: memory allocation fails"
                stop
            end if
        end if

        ! ----------------------------------------------------------------------------------------
        ! third: fill send buffer
        ! int buffer:  store receiver block id, neighborhood and level difference (in order of neighbor proc rank, use com matrix)
        ! real buffer: store block data (in order of neighbor proc rank, use com matrix)
        ! first element of int buffer = length of real buffer (buffer_i)

        ! fill send buffer and position communication matrix
        call fill_send_buffer_for_redundant_check( params, hvy_block, com_lists, com_matrix(rank+1,:), rank, int_send_buffer, real_send_buffer )

        ! calculate position matrix: position is column in send buffer, so simply count the number of communications
        ! loop over all com_matrix elements
        do i = 1, size(com_matrix_pos,1)
            ! new line, means new proc: reset counter
            k = 1
            ! loop over communications
            do j = 1, size(com_matrix_pos,1)
                ! found external communication
                if ( (com_matrix(i,j) /= 0) .and. (i /= j) ) then
                    ! save com position
                    com_matrix_pos(i,j) = k
                    ! increase counter
                    k = k + 1

                end if
            end do
        end do

        ! ----------------------------------------------------------------------------------------
        ! fourth: get data for receive buffer

        ! communicate, fill receive buffer
        call isend_irecv_data_check_redundant( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix, com_matrix_pos )

        ! ----------------------------------------------------------------------------------------
        ! fifth: write receive buffer to heavy data

        ! loop over corresponding com matrix line
        do k = 1, number_procs

            ! received data from proc k-1
            if ( ( com_matrix(rank+1, k) > 0 ) .and. ( (rank+1) /= k ) ) then

                ! set buffer position and calculate length if integer/real buffer
                buffer_pos = com_matrix_pos(rank+1, k)
                int_N  = com_matrix(rank+1, k) * 3 + 1
                real_N = int_receive_buffer( 1, buffer_pos )

                ! read received data
                if ( params%threeD_case ) then
                    ! 3D:
                    !call write_receive_buffer_3D(params, int_receive_buffer(2:int_N, buffer_pos), real_receive_buffer(1:real_N, buffer_pos), hvy_block )
                else
                    ! 2D:
                    call write_receive_buffer_2D_check_redundant(params, int_receive_buffer(2:int_N, buffer_pos), real_receive_buffer(1:real_N, buffer_pos), hvy_block(:, :, 1, :, :), stop_status )
                end if

            end if

        end do

    end if

    ! clean up
    deallocate( com_matrix_pos, stat=allocate_error )

    deallocate( int_send_buffer, stat=allocate_error )
    deallocate( int_receive_buffer, stat=allocate_error )
    deallocate( real_send_buffer, stat=allocate_error )
    deallocate( real_receive_buffer, stat=allocate_error )

end subroutine check_external_nodes

!----------------------------------------------------------------------------------
subroutine max_com_num_2( N, P, com_matrix_line, rank )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> number of coms, number of procs
    integer(kind=ik), intent(out)                   :: N, P

    !> com matrix line
    integer(kind=ik), intent(in)                    :: com_matrix_line(:)

    !> proc rank
    integer(kind=ik), intent(in)                    :: rank

    ! loop variable
    integer(kind=ik)                                :: k

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset N
    N = 0
    ! reste P
    P = 0

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all line elements
    do k = 1, size(com_matrix_line,1)

        ! communication to other proc, do not count internal communications
        if ( (com_matrix_line(k) /= 0) .and. (k /= rank+1) ) then

            ! max number
            N = max(N, com_matrix_line(k))
            ! add one proc
            P = P + 1

        end if

    end do

end subroutine max_com_num_2

! ============================================================================================

subroutine fill_send_buffer_for_redundant_check( params, hvy_block, com_lists, com_matrix_line, rank, int_send_buffer, real_send_buffer )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params

    !> heavy data array - block data
    real(kind=rk), intent(in)                       :: hvy_block(:, :, :, :, :)

    !> communication lists:
    integer(kind=ik), intent(in)                    :: com_lists(:, :, :)

    !> com matrix line
    integer(kind=ik), intent(in)                    :: com_matrix_line(:)

    !> proc rank
    integer(kind=ik), intent(in)                    :: rank

    !> integer send buffer
    integer(kind=ik), intent(inout)                 :: int_send_buffer(:,:)
    !> real send buffer
    real(kind=rk), intent(inout)                    :: real_send_buffer(:,:)

    ! loop variable
    integer(kind=ik)                                :: k, i

    ! column number of send buffer, position in integer buffer
    integer(kind=ik)                                :: column_pos, int_pos

    ! send buffer for one proc
    real(kind=rk), allocatable                      :: proc_send_buffer(:)

    ! allocation error variable
    integer(kind=ik)                                :: allocate_error

    ! index of send buffer, return from create_send_buffer subroutine
    integer(kind=ik)                                :: buffer_i

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset column number
    column_pos = 1

    ! allocate proc send buffer, size = line size of real send buffer
    allocate( proc_send_buffer( size(real_send_buffer,1) ), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all line elements
    do k = 1, size(com_matrix_line,1)

        ! communication to other proc, do not work with internal communications
        if ( (com_matrix_line(k) /= 0) .and. (k /= rank+1) ) then

            ! first: real data
            ! ----------------

            ! write real send buffer for proc k
            if ( params%threeD_case ) then
                ! 3D:
                !call create_send_buffer_3D(params, hvy_block, com_lists( 1:com_matrix_line(k), :, k), com_matrix_line(k), proc_send_buffer, buffer_i)
            else
                ! 2D:
                call create_send_buffer_2D_check_redundant(params, hvy_block(:, :, 1, :, :), com_lists( 1:com_matrix_line(k), :, k), com_matrix_line(k), proc_send_buffer, buffer_i)
            end if

            ! real buffer entry
            real_send_buffer( 1 : buffer_i, column_pos ) = proc_send_buffer( 1 : buffer_i )

            ! second: integer data
            ! --------------------

            ! save real buffer length
            int_send_buffer(1, column_pos) = buffer_i

            ! reset position
            int_pos = 2

            ! loop over all communications to this proc
            do i = 1, com_matrix_line(k)

                ! int buffer entry: neighbor block id, neighborhood, level difference
                int_send_buffer( int_pos  , column_pos ) = com_lists( i, 4, k)
                int_send_buffer( int_pos+1, column_pos ) = com_lists( i, 5, k)
                int_send_buffer( int_pos+2, column_pos ) = com_lists( i, 6, k)
                ! increase int buffer position
                int_pos = int_pos + 3

            end do

            ! third: increase column number
            ! ------------------------------------
            column_pos = column_pos + 1

        end if

    end do

    ! clean up
    deallocate( proc_send_buffer, stat=allocate_error )

end subroutine fill_send_buffer_for_redundant_check

! ============================================================================================

subroutine create_send_buffer_2D_check_redundant(params, hvy_block, com_list, com_number, send_buff, buffer_i)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> heavy data array - block data
    real(kind=rk), intent(in)                       :: hvy_block(:, :, :, :)

    !> com list
    integer(kind=ik), intent(in)                    :: com_list(:, :)
    !> com_list id, number of communications
    integer(kind=ik), intent(in)                    :: com_number

    !> send buffer
    real(kind=rk), intent(out)                      :: send_buff(:)

    !> buffer index
    integer(kind=ik), intent(out)                   :: buffer_i

    ! grid parameter
    integer(kind=ik)                                :: Bs, g

    ! com list elements
    integer(kind=ik)                                :: my_block, neighbor_block, my_dir, level_diff

    ! loop variable
    integer(kind=ik)                                :: k, dF
!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    buffer_i         = 1

!---------------------------------------------------------------------------------------------
! main body

    ! fill send buffer
    do k = 1 , com_number

        my_block        = com_list( k, 3 )
        neighbor_block  = com_list( k, 4 )
        my_dir          = com_list( k, 5 )
        level_diff      = com_list( k, 6 )

        select case(my_dir)
            ! '__N'
            case(1)
                do dF = 2, params%number_data_fields+1
                    send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1, g+1:Bs+g, dF, my_block )
                    buffer_i                            = buffer_i + Bs
                end do

            ! '__E'
            case(2)
                do dF = 2, params%number_data_fields+1
                    send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1:Bs+g, Bs+g, dF, my_block )
                    buffer_i                            = buffer_i + Bs
                end do

            ! '__S'
            case(3)
                do dF = 2, params%number_data_fields+1
                    send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( Bs+g, g+1:Bs+g, dF, my_block )
                    buffer_i                            = buffer_i + Bs
                end do

            ! '__W'
            case(4)
                do dF = 2, params%number_data_fields+1
                    send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1:Bs+g, g+1, dF, my_block )
                    buffer_i                            = buffer_i + Bs
                end do

            ! '_NE'
            case(5)
                do dF = 2, params%number_data_fields+1
                    send_buff(buffer_i)                 = hvy_block( g+1, Bs+g, dF, my_block )
                    buffer_i                            = buffer_i + 1
                end do

            ! '_NW'
            case(6)
                do dF = 2, params%number_data_fields+1
                    send_buff(buffer_i)                 = hvy_block( g+1, g+1, dF, my_block )
                    buffer_i                            = buffer_i + 1
                end do

            ! '_SE'
            case(7)
                do dF = 2, params%number_data_fields+1
                    send_buff(buffer_i)                 = hvy_block( Bs+g, Bs+g, dF, my_block )
                    buffer_i                            = buffer_i + 1
                end do

            ! '_SW'
            case(8)
                do dF = 2, params%number_data_fields+1
                    send_buff(buffer_i)                 = hvy_block( Bs+g, g+1, dF, my_block )
                    buffer_i                            = buffer_i + 1
                end do

            ! 'NNE'
            case(9)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1, (Bs+1)/2+g:Bs+g, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1, g+1:Bs+g:2, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    end if
                end do

            ! 'NNW'
            case(10)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1, g+1:(Bs+1)/2+g, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1, g+1:Bs+g:2, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    end if
                end do

            ! 'SSE'
            case(11)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( Bs+g, (Bs+1)/2+g:Bs+g, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( Bs+g, g+1:Bs+g:2, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    end if
                end do

            ! 'SSW'
            case(12)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( Bs+g, g+1:(Bs+1)/2+g, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( Bs+g, g+1:Bs+g:2, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    end if
                end do

            ! 'ENE'
            case(13)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:(Bs+1)/2+g, Bs+g, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, Bs+g, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    end if
                end do

            ! 'ESE'
            case(14)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( (Bs+1)/2+g:Bs+g, Bs+g, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    elseif ( level_diff == 1 ) then
                         ! sender on higher level
                         ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, Bs+g, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    end if
                end do

            ! 'WNW'
            case(15)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:(Bs+1)/2+g, g+1, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+1, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    end if
                end do

            ! 'WSW'
            case(16)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( (Bs+1)/2+g:Bs+g, g+1, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! send data
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+1, dF, my_block )
                        buffer_i                                    = buffer_i + (Bs+1)/2

                    end if
                end do

        end select
    end do

    ! decrease buffer index (for using in other subroutines)
    buffer_i = buffer_i - 1

end subroutine create_send_buffer_2D_check_redundant

! ============================================================================================

subroutine isend_irecv_data_check_redundant( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix, com_matrix_pos )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> send/receive buffer, integer and real
    integer(kind=ik), intent(in)        :: int_send_buffer(:,:)
    integer(kind=ik), intent(out)       :: int_receive_buffer(:,:)

    real(kind=rk), intent(in)           :: real_send_buffer(:,:)
    real(kind=rk), intent(out)          :: real_receive_buffer(:,:)

    !> communications matrix: neighboring proc rank
    !> com matrix pos: position in send buffer
    integer(kind=ik), intent(in)        :: com_matrix(:,:), com_matrix_pos(:,:)

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs
    ! MPI status
    !integer                             :: status(MPI_status_size)

    ! MPI message tag
    integer(kind=ik)                    :: tag
    ! MPI request
    integer(kind=ik)                    :: send_request(size(com_matrix,1)), recv_request(size(com_matrix,1))

    ! column number of send buffer, column number of receive buffer, real data buffer length
    integer(kind=ik)                    :: send_pos, receive_pos, real_pos

    ! loop variable
    integer(kind=ik)                    :: k, i

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank            = params%rank
    number_procs    = params%number_procs

    ! set message tag
    tag = 0

!---------------------------------------------------------------------------------------------
! main body

    ! ----------------------------------------------------------------------------------------
    ! first: integer data

    ! reset communication counter
    i = 0

    ! reset request arrays
    recv_request = MPI_REQUEST_NULL
    send_request = MPI_REQUEST_NULL

    ! loop over corresponding com matrix line
    do k = 1, number_procs

        ! communication between proc rank and proc k-1
        if ( ( com_matrix(rank+1, k) > 0 ) .and. ( (rank+1) /= k ) ) then

            ! test output
            !write(*, '( "rank ", i3, " send/receive to/from rank " , i3)') rank, k-1

            ! increase communication counter
            i = i + 1

            tag = rank+1+k

            ! receive buffer column number, read from position matrix
            receive_pos = com_matrix_pos(rank+1, k)

            ! receive data
            call MPI_Irecv( int_receive_buffer(1, receive_pos), size(int_receive_buffer,1), MPI_INTEGER4, k-1, tag, MPI_COMM_WORLD, recv_request(i), ierr)

            ! send buffer column number, read position matrix
            send_pos = com_matrix_pos(rank+1, k)

            ! send data
            call MPI_Isend( int_send_buffer(1, send_pos), size(int_send_buffer,1), MPI_INTEGER4, k-1, tag, MPI_COMM_WORLD, send_request(i), ierr)

        end if

    end do

    ! synchronize non-blocking communications
    ! note: single status variable do not work with all compilers, so use MPI_STATUSES_IGNORE instead
    if (i>0) then
        call MPI_Waitall( i, send_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
        call MPI_Waitall( i, recv_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
    end if

    ! ----------------------------------------------------------------------------------------
    ! second: real data

    ! reset communication couter
    i = 0

    ! reset request arrays
    recv_request = MPI_REQUEST_NULL
    send_request = MPI_REQUEST_NULL

    ! loop over corresponding com matrix line
    do k = 1, number_procs

        ! communication between proc rank and proc k-1
        if ( ( com_matrix(rank+1, k) > 0 ) .and. ( (rank+1) /= k ) ) then

            ! increase communication counter
            i = i + 1

            tag = number_procs*10*(rank+1+k)

            ! receive buffer column number, read from position matrix
            receive_pos = com_matrix_pos(rank+1, k)

            ! real buffer length
            real_pos = int_receive_buffer(1, receive_pos)

            ! receive data
            call MPI_Irecv( real_receive_buffer(1, receive_pos), real_pos, MPI_REAL8, k-1, tag, MPI_COMM_WORLD, recv_request(i), ierr)
            !call MPI_Irecv( real_receive_buffer(1, receive_pos), real_pos, MPI_DOUBLE_PRECISION, k-1, tag, MPI_COMM_WORLD, recv_request(i), ierr)

            ! send buffer column number, read position matrix
            send_pos = com_matrix_pos(rank+1, k)

            ! real buffer length
            real_pos = int_send_buffer(1, send_pos)

            ! send data
            call MPI_Isend( real_send_buffer(1, send_pos), real_pos, MPI_REAL8, k-1, tag, MPI_COMM_WORLD, send_request(i), ierr)
            !call MPI_Isend( real_send_buffer(1, send_pos), real_pos, MPI_DOUBLE_PRECISION, k-1, tag, MPI_COMM_WORLD, send_request(i), ierr)

        end if

    end do

    ! synchronize non-blocking communications
    if (i>0) then
        call MPI_Waitall( i, send_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
        call MPI_Waitall( i, recv_request(1:i), MPI_STATUSES_IGNORE, ierr) !status, ierr)
    end if

end subroutine isend_irecv_data_check_redundant

! ============================================================================================

subroutine write_receive_buffer_2D_check_redundant(params, int_buffer, recv_buff, hvy_block, stop_status)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> send buffer
    integer(kind=ik), intent(in)                    :: int_buffer(:)
    !> send buffer
    real(kind=rk), intent(in)                       :: recv_buff(:)

    !> heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :)

    ! buffer index
    integer(kind=ik)                                :: buffer_i

    ! grid parameter
    integer(kind=ik)                                :: Bs, g

    ! com list elements
    integer(kind=ik)                                :: my_block, my_dir, level_diff

    ! loop variable
    integer(kind=ik)                                :: k, dF

    ! status of nodes check: if true: stops program
    logical, intent(inout)                          :: stop_status

    ! difference between sender/receiver nodes
    real(kind=rk)                                   :: diff_norm
    ! error threshold
    real(kind=rk)                                   :: eps

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    buffer_i         = 1

    ! set error threshold
    eps = 1e-12_rk

!---------------------------------------------------------------------------------------------
! main body

    ! write received data in block data
    do k = 1, size(int_buffer,1), 3

        my_block        = int_buffer( k )
        my_dir          = int_buffer( k+1 )
        level_diff      = int_buffer( k+2 )

        select case(my_dir)
            ! '__N'
            case(1)
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( Bs+g, g+1:Bs+g, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+Bs-1) ))**2 )

                    buffer_i = buffer_i + Bs

                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes, external (__N)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            ! '__E'
            case(2)
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1:Bs+g, g+1, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+Bs-1) ))**2 )

                    buffer_i = buffer_i + Bs

                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes, external (__E)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            ! '__S'
            case(3)
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1, g+1:Bs+g, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+Bs-1) ))**2 )

                    buffer_i = buffer_i + Bs

                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes, external (__S)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            ! '__W'
            case(4)
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt(sum(( hvy_block( g+1:Bs+g, Bs+g, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+Bs-1) ))**2 )

                    buffer_i = buffer_i + Bs

                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes, external (__W)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            ! '_NE'
            case(5)
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt((( hvy_block( Bs+g, g+1, dF, my_block ) - &
                                           recv_buff(buffer_i) ))**2 )

                    buffer_i = buffer_i + 1

                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes, external (_NE)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            ! '_NW'
            case(6)
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt((( hvy_block( Bs+g, Bs+g, dF, my_block ) - &
                                           recv_buff(buffer_i) ))**2 )

                    buffer_i = buffer_i + 1

                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes, external (_NW)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            ! '_SE'
            case(7)
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt((( hvy_block( g+1, g+1, dF, my_block ) - &
                                           recv_buff(buffer_i) ))**2 )

                    buffer_i = buffer_i + 1

                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes, external (_SE)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            ! '_SW'
            case(8)
                do dF = 2, params%number_data_fields+1
                    diff_norm = sqrt((( hvy_block( g+1, Bs+g, dF, my_block ) - &
                                           recv_buff(buffer_i) ))**2 )

                    buffer_i = buffer_i + 1

                    ! check error
                    if ( diff_norm > eps ) then
                        write(*,*) "ERROR: difference in redundant nodes, external (_SW)"
                        ! save data
                        ! write proc rank into block data only for first error case
                        if (stop_status) then
                            ! do nothing
                        else
                            hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                            stop_status = .true.
                        end if
                    end if
                end do

            ! 'NNE'
            case(9)
                if ( level_diff == -1 ) then
                    ! sender on lower level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( Bs+g, g+1:Bs+g:2, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (NNE)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do

                elseif ( level_diff == 1 ) then
                    ! sender on higher level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( Bs+g, g+(Bs+1)/2:Bs+g, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (NNE)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do
                end if

            ! 'NNW'
            case(10)
                if ( level_diff == -1 ) then
                    ! sender on lower level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( Bs+g, g+1:Bs+g:2, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (NNW)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do

                elseif ( level_diff == 1 ) then
                    ! sender on higher level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( Bs+g, g+1:g+(Bs+1)/2, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (NNW)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do
                end if

            ! 'SSE'
            case(11)
                if ( level_diff == -1 ) then
                    ! sender on lower level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( g+1, g+1:Bs+g:2, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (SSE)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do

                elseif ( level_diff == 1 ) then
                    ! sender on higher level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( g+1, g+(Bs+1)/2:Bs+g, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (SSE)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do
                end if

            ! 'SSW'
            case(12)
                if ( level_diff == -1 ) then
                    ! sender on lower level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( g+1, g+1:Bs+g:2, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (SSW)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do

                elseif ( level_diff == 1 ) then
                    ! sender on higher level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( g+1, g+1:g+(Bs+1)/2, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (SSW)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do
                end if

            ! 'ENE'
            case(13)
                if ( level_diff == -1 ) then
                    ! sender on lower level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( g+1:Bs+g:2, g+1, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (ENE)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do

                elseif ( level_diff == 1 ) then
                    ! sender on higher level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( g+1:g+(Bs+1)/2, g+1, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (ENE)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do
                end if

            ! 'ESE'
            case(14)
                if ( level_diff == -1 ) then
                    ! sender on lower level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( g+1:Bs+g:2, g+1, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (ESE)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do

                elseif ( level_diff == 1 ) then
                    ! sender on higher level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( g+(Bs+1)/2:Bs+g, g+1, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (ESE)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do
                end if

            ! 'WNW'
            case(15)
                if ( level_diff == -1 ) then
                    ! sender on lower level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( g+1:Bs+g:2, Bs+g, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (WNW)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do

                elseif ( level_diff == 1 ) then
                    ! sender on higher level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( g+1:g+(Bs+1)/2, Bs+g, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (WNW)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do
                end if

            ! 'WSW'
            case(16)
                if ( level_diff == -1 ) then
                    ! sender on lower level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( g+1:Bs+g:2, Bs+g, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (WSW)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do

                elseif ( level_diff == 1 ) then
                    ! sender on higher level
                    ! loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        diff_norm = sqrt(sum(( hvy_block( g+(Bs+1)/2:Bs+g, Bs+g, dF, my_block ) - &
                                           recv_buff(buffer_i:buffer_i+(Bs+1)/2-1) ))**2 )

                        buffer_i = buffer_i + (Bs+1)/2

                        ! check error
                        if ( diff_norm > eps ) then
                            write(*,*) "ERROR: difference in redundant nodes, external (WSW)"
                            ! save data
                            ! write proc rank into block data only for first error case
                            if (stop_status) then
                                ! do nothing
                            else
                                hvy_block( :, :, dF, my_block ) = real( my_dir, kind=rk)
                                stop_status = .true.
                            end if
                        end if
                    end do
                end if

        end select

    end do

end subroutine write_receive_buffer_2D_check_redundant
