!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name synchronize_ghosts.f90
!> \version 0.5
!> \author msr
!
!> \brief synchronize ghosts nodes
!! synchronization is done in 4 stages:
!!          stage 1: finer blocks give data to coarser blocks (restriction)
!!          stage 2: blocks on same level synch
!!          stage 3: coarser blocks synch with finer blocks (prediction)
!!          stage 4: additional step to correct interpolation errors (needed for stability)
!!
!! subroutine structure:
!! ---------------------
!!      1. create external communication list:
!!              loop over active heavy blocks, if neighbor needs ghost nodes and neighbor is located
!!              on other proc: new com list entry
!!      2. synchronize com matrix:
!!              needed for buffer allocation and communication order
!!      3. fill send buffer:
!!              sender proc loop over his com list and gather all data for all neighbor procs/blocks
!!      4. send/receive data
!!      5. write own and received ghost nodes into own blocks
!!              - loop over own active heavy data (own block has to receive data)
!!              - if neighbor exists: check level difference and synch stage
!!                note: here level diff is neighbor (sender) - me (receiver)
!!              - if something to synchronize: switch neighborhood, because, neighborhood is defined from
!!                senders point of view, but we loop actually over receiver blocks
!!              - if neighbor is internal: copy data, note: correct sender/receiver block ids!
!!              - external neighbor: search for neighbor data in received buffer
!!                note: reiceived buffer includes all data, and write routine can handle this, but
!!                at this point we work onm one specific block, maybe later we can switch back to old
!!                internal/external handling
!!  \todo TODO: rework ghost nodes writing routine to avoid receive buffer searching
!!  stage 4 handling:
!!      - external neighbor: check condition in com list creation - then send/receive data as before
!!      - internal neighbor: use new copy_redundant_nodes subroutine
!!  \todo TODO: if possible use same copy ghost nodes routine for all 4 stages
!!
!>
!! input:    - params, light and heavy data \n
!! output:   - heavy data array
!
!> \details
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4 \n
!! 06/01/17 - use RMA to synchronize data \n
!! 31/01/17 - switch to 3D, v0.5 \n
!! 12/04/17 - redundant ghost nodes workaround
!! 19/05/17 - switch to new synchronization routine (correct redundant nodes handling)
!! 16/06/17 - allocate all send/receive buffer in ini step -> huge performance boost
!
! ********************************************************************************************

subroutine synchronize_ghosts(  params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, grid_changed, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )

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

    ! grid stay fixed between two synch calls, so use old com_lists and com_matrix
    logical, intent(in)                 :: grid_changed

    ! loop variables
    integer(kind=ik)                    :: k, N, i, j, neighbor_num, synch_stage, hvy_id, lgt_id, neighbor_light_id, neighborhood, com_number

    ! synch switch
    logical                             :: synch

    ! difference between current block and neighbor block level
    integer(kind=ik)                    :: level_diff

    ! grid parameter
    integer(kind=ik)                    :: g, Bs

    ! process rank
    integer(kind=ik)                    :: rank, neighbor_rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! communication lists:
    integer(kind=ik), intent(inout)     :: com_lists(:, :, :, :)

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1, time_sum

    ! communications matrix:
    integer(kind=ik), intent(inout)     :: com_matrix(:,:,:)

    ! communications matrix:
    ! count the number of communications between procs
    ! row/column number encodes process rank + 1
    ! com matrix pos: position in send buffer
    integer(kind=ik), allocatable       :: com_matrix_pos(:,:)

    ! send/receive buffer, integer and real
    integer(kind=ik), intent(inout)      :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), intent(inout)         :: real_send_buffer(:,:), real_receive_buffer(:,:)

    ! number of communications, number of neighboring procs
    integer(kind=ik)                     :: my_n_com, n_procs

    ! indexes of buffer array and column number in buffer, use for readability
    integer(kind=ik)                     :: int_start, real_start, buffer_pos, real_N

    ! variable for non-uniform mesh correction: remove redundant node between fine->coarse blocks
    integer(kind=ik)                     :: rmv_redundant

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! start time
    sub_t0 = MPI_Wtime()

    time_sum = 0.0_rk

    N = params%number_blocks

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    rmv_redundant = 0

    ! set non-uniform mesh correction
    if ( params%non_uniform_mesh_correction ) then
        rmv_redundant = 1
    else
        rmv_redundant = 0
    end if

    ! set MPI parameter
    rank         = params%rank
    number_procs = params%number_procs

    ! allocate com position matrix
    allocate( com_matrix_pos(number_procs, number_procs) )

!    ! reset ghost nodes for all active blocks
!    if ( params%debug ) then
!        call reset_ghost_nodes(  params, hvy_block, hvy_active, hvy_n )
!    end if

!---------------------------------------------------------------------------------------------
! main body

    ! reset com-list, com_plan, com matrix, receiver lists
    if (grid_changed) then
        com_lists       = -99
        com_matrix      =  0
    end if

    ! synchronize in three stages
    ! stage 1: neighbors with level diff 0, -1
    ! stage 2: neighbors with level diff +1
    ! stage 3: redundant node correction

    do synch_stage = 1, 4

        ! reset com-list, com_plan, com matrix, receiver lists
        com_matrix_pos  =  0

        ! end time
        sub_t1 = MPI_Wtime()
        time_sum = time_sum + (sub_t1 - sub_t0)
        ! start time
        sub_t0 = MPI_Wtime()

        ! next steps only if the grid has changed
        if (grid_changed) then

            ! ----------------------------------------------------------------------------------------
            ! first: create com matrix and com list for external communications
            ! output depends on synch stage
            ! ----------------------------------------------------------------------------------------
            call create_external_com_list(  params, lgt_block, hvy_neighbor, hvy_active, hvy_n, com_matrix(:,:,synch_stage), com_lists(:,:,:,synch_stage), synch_stage )

            ! symmetrize com matrix - choose max com number
            do i = 1, number_procs
                do j = i+1, number_procs
                    com_matrix(i,j,synch_stage) = max( com_matrix(i,j,synch_stage), com_matrix(j,i,synch_stage) )
                    com_matrix(j,i,synch_stage) = com_matrix(i,j,synch_stage)
                end do
            end do

        end if

!        ! save com matrix
!        if ( params%debug ) then
!            call write_com_matrix( com_matrix )
!        end if

        ! ----------------------------------------------------------------------------------------
        ! third: fill send buffer
        ! max number of communications and neighboring procs - use for buffer allocation
        ! every proc loop over com matrix line
        ! ----------------------------------------------------------------------------------------
        call max_com_num( my_n_com, n_procs, com_matrix(rank+1,:,synch_stage), rank )

        ! for proc without neighbors: set n_procs to 1
        ! so we allocate arrays with second dimension=1
        if (n_procs==0) n_procs = 1

        ! end time
        sub_t1 = MPI_Wtime()
        ! write time
        if ( params%debug ) then
            ! find free or corresponding line
            k = 1
            do while ( debug%name_comp_time(k) /= "---" )
                ! entry for current subroutine exists
                if ( debug%name_comp_time(k) == "synch. ghosts - com list/matrix" ) exit
                k = k + 1
            end do
            ! write time
            debug%name_comp_time(k) = "synch. ghosts - com list/matrix"
            debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
            debug%comp_time(k, 2)   = debug%comp_time(k, 2) + (sub_t1 - sub_t0)
        end if

        ! start time
        sub_t0 = MPI_Wtime()

        ! next steps only for more than two procs
        if ( number_procs > 1 ) then

            ! ----------------------------------------------------------------------------------------
            ! fill send buffer
            ! int buffer:  store receiver block id, neighborhood and level difference (in order of neighbor proc rank, use com matrix)
            ! real buffer: store block data (in order of neighbor proc rank, use com matrix)
            ! first element of int buffer = length of real buffer (buffer_i)

            ! fill send buffer and position communication matrix
            call fill_send_buffer( params, hvy_block, com_lists(:,:,:,synch_stage), com_matrix(rank+1,:,synch_stage), rank, int_send_buffer, real_send_buffer )

            ! end time
            sub_t1 = MPI_Wtime()
            ! write time
            if ( params%debug ) then
                ! find free or corresponding line
                k = 1
                do while ( debug%name_comp_time(k) /= "---" )
                    ! entry for current subroutine exists
                    if ( debug%name_comp_time(k) == "synch. ghosts - fill send buffer" ) exit
                    k = k + 1
                end do
                ! write time
                debug%name_comp_time(k) = "synch. ghosts - fill send buffer"
                debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
                debug%comp_time(k, 2)   = debug%comp_time(k, 2) + (sub_t1 - sub_t0)
            end if

            ! start time
            sub_t0 = MPI_Wtime()

            ! ----------------------------------------------------------------------------------------
            ! fourth: send/receive data
            ! calculate position matrix: position is column in send buffer, so simply count the number of communications
            ! loop over all com_matrix elements
            ! ----------------------------------------------------------------------------------------
            do i = 1, size(com_matrix_pos,1)
                ! new line, means new proc: reset counter
                k = 1
                ! loop over communications
                do j = 1, size(com_matrix_pos,1)
                    ! found external communication
                    if ( (com_matrix(i,j,synch_stage) /= 0) .and. (i /= j) ) then
                        ! save com position
                        com_matrix_pos(i,j) = k
                        ! increase counter
                        k = k + 1
                    end if
                end do
            end do

            ! communicate, fill receive buffer
            call fill_receive_buffer( params, int_send_buffer, real_send_buffer, int_receive_buffer, real_receive_buffer, com_matrix(:,:,synch_stage), com_matrix_pos  )

        end if

        ! end time
        sub_t1 = MPI_Wtime()
        ! write time
        if ( params%debug ) then
            ! find free or corresponding line
            k = 1
            do while ( debug%name_comp_time(k) /= "---" )
                ! entry for current subroutine exists
                if ( debug%name_comp_time(k) == "synch. ghosts - communication" ) exit
                k = k + 1
            end do
            ! write time
            debug%name_comp_time(k) = "synch. ghosts - communication"
            debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
            debug%comp_time(k, 2)   = debug%comp_time(k, 2) + (sub_t1 - sub_t0)
        end if

        ! start time
        sub_t0 = MPI_Wtime()

        ! ----------------------------------------------------------------------------------------
        ! fifth: write ghost nodes
        ! set loop number for 2D/3D case
        ! ----------------------------------------------------------------------------------------
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
                    ! note: level diff is calculated sender-receiver, here the receiver looks at his neighbor
                    level_diff = lgt_block( neighbor_light_id, params%max_treelevel+1 ) - lgt_block( lgt_id, params%max_treelevel+1 )

                    ! proof if neighbor internal or external
                    call lgt_id_to_proc_rank( neighbor_rank, neighbor_light_id, N )

                    if ( (synch_stage == 1) .or. (synch_stage == 2) .or. (synch_stage == 3) ) then

                        synch = .false.
                        ! check synch stage and level difference
                        if ( (level_diff == 1) .and. (synch_stage == 1) ) then
                            synch = .true.
                        elseif ( (level_diff == 0) .and. (synch_stage == 2) ) then
                            synch = .true.
                        elseif ( (level_diff == -1) .and. (synch_stage == 3) ) then
                            synch = .true.
                        end if

                        if (synch) then

                            ! neighborhood
                            neighborhood = 0
                            select case (i)
                                case(1)
                                    neighborhood = 3
                                case(2)
                                    neighborhood = 4
                                case(3)
                                    neighborhood = 1
                                case(4)
                                    neighborhood = 2
                                case(5)
                                    neighborhood = 8
                                case(6)
                                    neighborhood = 7
                                case(7)
                                    neighborhood = 6
                                case(8)
                                    neighborhood = 5
                                case(9)
                                    neighborhood = 11
                                case(10)
                                    neighborhood = 12
                                case(11)
                                    neighborhood = 9
                                case(12)
                                    neighborhood = 10
                                case(13)
                                    neighborhood = 15
                                case(14)
                                    neighborhood = 16
                                case(15)
                                    neighborhood = 13
                                case(16)
                                    neighborhood = 14
                            end select

                            ! sender block on same or higher level
                            if ( rank == neighbor_rank ) then
                                ! calculate internal heavy id
                                call lgt_id_to_hvy_id( hvy_id, neighbor_light_id, rank, N )

                                ! internal neighbor
                                if ( params%threeD_case ) then
                                    ! 3D:
                                    !call copy_ghost_nodes_3D( params, hvy_block, hvy_active(k), hvy_id, i, level_diff )
                                else
                                    ! 2D:
                                    call copy_ghost_nodes_2D( params, hvy_block(:, :, 1, :, :), hvy_id, hvy_active(k), neighborhood, level_diff )
                                end if

                            else

                                ! external neighbor
                                ! neighbor heavy id
                                call lgt_id_to_hvy_id( hvy_id, neighbor_light_id, neighbor_rank, N )

                                ! position in receiver buffer array
                                buffer_pos = com_matrix_pos(rank+1, neighbor_rank+1)
                                ! length of real buffer
                                real_N = int_receive_buffer( 1, buffer_pos )

                                ! find data in receive buffer
                                ! loop over int receive buffer, note: first number = length of real receive buffer
                                ! int buffer has neighbor heavy id on first (start at 2) position of 3 integer numbers per neighborhood
                                j = 3
                                real_start = 1
                                int_start = 2

                                ! calculate com number
                                do com_number = 2, com_matrix(rank+1, neighbor_rank+1, synch_stage) * 3 + 1, 3
                                    if ( int_receive_buffer( com_number, buffer_pos ) == -99 ) exit
                                end do
                                com_number = com_number - 3

                                do while ( j < com_number + 2 )

                                    ! check heavy id, note: receiver heavy id is hvy_active(k)
                                    ! sender (neighbor) heavy id is hvy_id
                                    ! second check: neighborhood
                                    if ( (int_receive_buffer( j-1, buffer_pos ) == hvy_active(k) ) .and. (int_receive_buffer( j, buffer_pos ) == neighborhood) ) then
                                        ! own block found
                                        int_start = j-1
                                        exit

                                    else

                                        ! increase real buffer start index
                                        ! neighborhood on j+1 position
                                        if ( params%threeD_case ) then
                                            ! to do
                                        else
                                            ! 2D
                                            !level_diff = int_receive_buffer( j+1, buffer_pos )

                                            select case ( int_receive_buffer( j, buffer_pos ) )
                                                case(1:4)
                                                    real_start = real_start + Bs*g*params%number_data_fields
                                                case(5:8)
                                                    if ( int_receive_buffer( j+1, buffer_pos ) == 1 ) then
                                                        ! case +1
                                                        real_start = real_start + (g+rmv_redundant)*(g+rmv_redundant)*params%number_data_fields
                                                    elseif ( int_receive_buffer( j+1, buffer_pos ) == -1 ) then
                                                        ! case -1
                                                        real_start = real_start + g*g*params%number_data_fields
                                                    else
                                                        ! case 0
                                                        real_start = real_start + (g+rmv_redundant)*(g+rmv_redundant)*params%number_data_fields
                                                    end if
                                                case(9:16)
                                                    if ( int_receive_buffer( j+1, buffer_pos ) /= -1 ) then
                                                        ! case +1
                                                        real_start = real_start + (g+rmv_redundant)*((Bs+1)/2)*params%number_data_fields
                                                    else
                                                        ! case -1
                                                        real_start = real_start + (Bs+g)*g*params%number_data_fields
                                                    end if
                                            end select
                                        end if

                                        ! next block
                                        j = j + 3

                                    end if

                                end do

                                ! error case
                                if ( real_start > real_N ) then
                                   write(*,'(80("-"))')
                                   write(*, '("ERROR: block ", i3 , " not found in receive buffer")') hvy_active(k)
                                   stop
                                end if

                                ! write external ghost nodes
                                ! note: write_receive_buffer subroutine also works on complete receive buffer
                                ! use start and end index to work only with exactly one neighbor
                                if ( params%threeD_case ) then
                                    ! 3D:
                                    !call write_receive_buffer_3D(params, int_receive_buffer(2:int_N, buffer_pos), real_receive_buffer(1:real_N, buffer_pos), hvy_block )
                                else
                                    ! 2D:
                                    call write_receive_buffer_2D(params, int_receive_buffer(int_start:int_start+2, buffer_pos), real_receive_buffer(real_start:real_N, buffer_pos), hvy_block(:, :, 1, :, :) )
                                end if

                            end if
                        end if

                    elseif (synch_stage == 4) then
                        ! third stage: correction stage
                        if ( rank == neighbor_rank ) then

                            ! neighborhood
                            neighborhood = 0
                            select case (i)
                                case(5)
                                    neighborhood = 8
                                case(6)
                                    neighborhood = 7
                                case(7)
                                    neighborhood = 6
                                case(8)
                                    neighborhood = 5
                            end select

                            ! calculate internal heavy id
                            call lgt_id_to_hvy_id( hvy_id, neighbor_light_id, rank, N )

                            ! stage 3
                            if ( params%threeD_case ) then
                                ! 3D:
                            else
                                ! 2D:
                                call copy_redundant_nodes_2D( params, hvy_block(:, :, 1, :, :), hvy_id, hvy_active(k), neighborhood, level_diff, hvy_neighbor )
                            end if

                        else
                            ! external neighbor

                            ! neighbor heavy id
                            call lgt_id_to_hvy_id( hvy_id, neighbor_light_id, neighbor_rank, N )

                            ! position in receiver buffer array
                            buffer_pos = com_matrix_pos(rank+1, neighbor_rank+1)

                            ! next steps only if there is at least one communication
                            ! means: buffer_pos > 0

                            if ( buffer_pos > 0 ) then

                                ! length of real buffer
                                real_N = int_receive_buffer( 1, buffer_pos )

                                ! find data in receive buffer
                                ! loop over int receive buffer, note: first number = length of real receive buffer
                                ! int buffer has neighbor heavy id on first (start at 2) position of 3 integer numbers per neighborhood
                                j = 3
                                real_start = 1
                                int_start = 2

                                ! calculate com number
                                if ( int_receive_buffer( 1, buffer_pos ) /= -99 ) then
                                    ! receive buffer is not empty
                                    do com_number = 2, com_matrix(rank+1, neighbor_rank+1, synch_stage) * 3 + 1, 3
                                        if ( int_receive_buffer( com_number, buffer_pos ) == -99 ) exit
                                    end do
                                    com_number = com_number - 3
                                else
                                    ! empty receive buffer
                                    com_number = -99
                                end if

                                neighborhood = 0

                                do while (j < com_number + 2)

                                    ! calculate neighborhood
                                    ! note: sender stored receiver neighborhood, but receiver has sender on
                                    ! sender neighborhood in hvy_neighbor list
                                    neighborhood = 0
                                    select case (i)
                                        case(5)
                                            neighborhood = 8
                                        case(6)
                                            neighborhood = 7
                                        case(7)
                                            neighborhood = 6
                                        case(8)
                                            neighborhood = 5
                                    end select

                                    ! check heavy id, note: receiver heavy id is hvy_active(k)
                                    ! sender (neighbor) heavy id is hvy_id
                                    ! second check: neighborhood
                                    if ( (int_receive_buffer( j-1, buffer_pos ) == hvy_active(k) ) .and. (int_receive_buffer( j, buffer_pos ) == neighborhood) ) then
                                        ! own block found
                                        int_start = j-1
                                        exit

                                    else

                                        ! increase real buffer start index
                                        ! neighborhood on j+1 position
                                        if ( params%threeD_case ) then
                                            ! to do
                                        else
                                            ! 2D
                                            select case ( int_receive_buffer( j, buffer_pos ) )
                                                case(5:8)
                                                    ! case 0
                                                    real_start = real_start + (g+rmv_redundant)*(g+rmv_redundant)*params%number_data_fields
                                            end select
                                        end if

                                        ! next block
                                        j = j + 3

                                    end if

                                end do

                                if ( neighborhood /= 0 ) then

                                    if ( real_start < real_N ) then

                                        ! write external ghost nodes
                                        ! note: write_receive_buffer subroutine also works on complete receive buffer
                                        ! use start and end index to work only with exactly one neighbor
                                        if ( params%threeD_case ) then
                                            ! 3D:
                                            !call write_receive_buffer_3D(params, int_receive_buffer(2:int_N, buffer_pos), real_receive_buffer(1:real_N, buffer_pos), hvy_block )
                                        else
                                            ! 2D:
                                            call write_receive_buffer_2D(params, int_receive_buffer(int_start:int_start+2, buffer_pos), real_receive_buffer(real_start:real_N, buffer_pos), hvy_block(:, :, 1, :, :) )
                                        end if

                                    end if

                                end if

                            end if

                        end if

                    end if
                end if
            end do
        end do

    end do

    ! clean up
    deallocate( com_matrix_pos )

    ! end time
    sub_t1 = MPI_Wtime()
    time_sum = time_sum + (sub_t1 - sub_t0)
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "synch. ghosts" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "synch. ghosts"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + time_sum
    end if

end subroutine synchronize_ghosts
