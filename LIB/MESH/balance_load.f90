!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name balance_load_2D.f90
!> \version 0.4
!> \author msr
!
!> \brief balance the load
!
!> \image html balancing.svg "Load balancing" width=400
!
!> \details
!! input:    - params, light and heavy data, neighbor data, lists of active blocks \n
!! output:   - light and heavy data arrays
!! \n
!> = log ======================================================================================
!!\n
!! 08/11/16    - switch to v0.4 \n
!! 16/11/2016  - Avoid some communication by more carefully distributing the excess blocks \n
!! 05/12/2016  - add space filling curve distribution \n
!
!> \image html load_balancing.svg width=500
! ********************************************************************************************

subroutine balance_load( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! send/receive buffer, note: size is equal to block data array, because if a block want to send all his data
    real(kind=rk), allocatable          :: buffer_data( :, :, :, :, : )
    integer(kind=ik), allocatable       :: buffer_light( : )
    ! light id start
    integer(kind=ik)                    :: my_light_start
    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank, proc_dist_id, proc_data_id
    ! number of processes
    integer(kind=ik)                    :: number_procs
    ! MPI message tag
    integer(kind=ik)                    :: tag
    ! MPI status
    integer                             :: status(MPI_status_size)
    ! distribution type
    character(len=80)                   :: distribution
    ! block distribution lists
    integer(kind=ik), allocatable       :: opt_dist_list(:), dist_list(:), friends(:,:), affinity(:)
    ! loop variables
    integer(kind=ik)                    :: k, N, num_blocks, l, com_i, com_N, heavy_id, sfc_id
    ! com plan
    integer(kind=ik), allocatable       :: com_plan(:,:)
    ! size of data array
    integer(kind=ik)                    :: data_size
    ! free light/heavy data id
    integer(kind=ik)                    :: lgt_free_id, hvy_free_id, lgt_id
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0
    ! space filling curve list
    integer(kind=ik), allocatable       :: sfc_com_list(:,:), sfc_sorted_list(:,:)
    ! hilbert code
    integer(kind=ik)                    :: hilbertcode(params%max_treelevel)
    ! communicator
    integer(kind=ik)                    :: WABBIT_COMM

!---------------------------------------------------------------------------------------------
! variables initialization

    ! start time
    t0 = MPI_Wtime()

    tag = 0

    distribution    = params%block_distribution

    ! MPI_parameters
    rank = params%rank
    number_procs = params%number_procs
    WABBIT_COMM  = params%WABBIT_COMM
    ! allocate block to proc lists
    allocate( opt_dist_list(1:number_procs), dist_list(1:number_procs))

    ! allocate com plan, maximal number of communications: every proc send every other proc something
    allocate( com_plan( number_procs*(number_procs-1), 3 ) )
    com_plan = -1

    ! allocate sfc com list, maximal number of communications is when every proc wants to send all of his blocks
    allocate( sfc_com_list( number_procs*params%number_blocks, 3 ) )
    sfc_com_list = -1

    ! allocate space filling curve list, maximal number of elements = max number of blocks
    allocate( sfc_sorted_list( lgt_n, 2) )

    ! allocate buffer arrays here
    allocate( buffer_data( size(hvy_block,1), size(hvy_block,2), &
              size(hvy_block,3), size(hvy_block,4), size(hvy_block,5) ) )
    allocate( buffer_light( params%number_blocks ) )

    ! number of blocks
    N = params%number_blocks

    ! reset block count
    num_blocks = 0

    ! light data start line
    my_light_start = rank*params%number_blocks

    ! size of data array, use for readability
    data_size = size(hvy_block,1) * size(hvy_block,2) * size(hvy_block,3) * size(hvy_block,4)


!---------------------------------------------------------------------------------------------
! main body

    if (params%number_procs == 1) then
        ! on only one proc, no balancing is required
        return
    endif

    !---------------------------------------------------------------------------------
    ! First step: define how many blocks each mpirank should have.
    !---------------------------------------------------------------------------------
    call set_desired_num_blocks_per_rank(params, dist_list, opt_dist_list, lgt_n, hvy_n)
    call write_block_distribution( dist_list, params )
    ! at this point, we know how many blocks a mpirank has: "dist_list(myrank+1)"
    ! and how many it should have, if equally distributed: "opt_dist_list(myrank+1)"


    select case(distribution)
        case("equal")! simple uniformly distribution

            allocate( affinity(1:params%number_blocks) )
            allocate( friends( 1:number_procs, 1:number_procs ))

            ! at this point, we know how many blocks a mpirank has: "dist_list(myrank+1)"
            ! and how many it should have, if equally distributed: "opt_dist_list(myrank+1)"
            if (maxval(abs(dist_list-opt_dist_list))==0) then
                ! the distribution is fine, nothing to do.
                return
            endif

            ! determine matrix of number of neighbor relations between mpiranks
            call compute_friends_table(params, hvy_neighbor, friends, hvy_active, hvy_n)
!
!            !---------------------------------------------------------------------------------
!            ! second step: create plan for communication
!            !---------------------------------------------------------------------------------
!            ! column
!            !    1     sender proc
!            !    2     receiver proc
!            !    3     number of blocks to send
!
!            com_i = 1
!            ! loop over all procs
!            do id_send = 1, number_procs
!                send_deficit = opt_dist_list(id_send) - dist_list(id_send)
!                ! find a receiver, if this is a SENDER
!                do while (send_deficit < 0)
!                    ! from my list of friends, gather my currently best friend
!                    tmp = maxloc( friends(id_send,:) )
!                    id_recv = tmp(1)
!                    ! "remove" the friend we're looking at from the list
!                    friends(id_send, id_recv) = -2
!
!                    ! can this proc take some of my blocks?
!                    recv_deficit = opt_dist_list(id_recv) - dist_list(id_recv)
!                    ! proc can receive data
!                    if ( ( recv_deficit > 0 ) .and. ( id_recv /= id_send ) ) then
!                        ! calculate number of send/receive blocks
!                        com_N = minval( abs((/send_deficit, recv_deficit/)) )
!
!                        ! create com plan
!                        com_plan( com_i, 1 ) = id_send - 1
!                        com_plan( com_i, 2 ) = id_recv - 1
!                        com_plan( com_i, 3 ) = com_N
!
!                        ! change distribution list
!                        dist_list(id_send) = dist_list(id_send) - com_N
!                        dist_list(id_recv) = dist_list(id_recv) + com_N
!
!                        ! new com was created
!                        com_i = com_i + 1
!                    end if
!                    ! recompute senders deficit (how many blocks sender still wants to get rid of)
!                    send_deficit = opt_dist_list(id_send) - dist_list(id_send)
!                end do
!            end do
!            ! we counted one too far
!            com_i = com_i -1
!
!            !---------------------------------------------------------------------------------
!            ! third: actually send/receive data
!            !---------------------------------------------------------------------------------
!            l = 1
!            ! loop over all planed communications
!            do k = 1, com_i
!                !*************** SEND CASE
!                if ( com_plan(k, 1) == rank ) then
!                    ! yes, and I'll send "com_plan(k, 3)" blocks to proc "com_plan(k, 2)"
!                    !---------------------------------------------------------------------------------------
!                    ! affinity list, HEAVY DATA ARRAY
!                    call compute_affinity(params, my_block_list, hvy_neighbor, rank, com_plan(k, 2), affinity, hvy_active, hvy_n)
!
!                    com_N = 1
!                    do while ( com_N <= com_plan(k, 3) )
!                        ! fetch most desirable block
!                        tmp = maxloc(affinity)
!                        heavy_id = tmp(1)
!                        ! we now use this block and must be sure not to use it again
!                        affinity(heavy_id) = -99
!                        light_id = rank*params%number_blocks + heavy_id
!
!                        if (my_block_list(light_id, 1) /= -1) then! block is active
!                            ! fill buffer, heavy data
!                            buffer_data(:, :, :, com_N) = hvy_block(:, :, :, heavy_id)
!                            ! ... light data
!                            buffer_light(com_N) = light_id
!                            ! count com
!                            com_N = com_N + 1
!                            ! delete heavy data
!                            hvy_block(:, :, :, heavy_id) = 0.0_rk
!                            ! delete light data
!                            my_block_list( light_id, : ) = -1
!                        else
!                          write(*,*) "ERROR: load balancing error 576879"
!                          stop "unforeseen 576879"
!                        end if
!                    end do
!
!                    ! error case
!                    if ( com_N-1 /= com_plan(k, 3) ) then
!                        stop "ERROR: load balancing error 003857"
!                    end if
!
!                    ! send data
!                    call MPI_Send( buffer_light, (com_N-1), MPI_INTEGER4, com_plan(k, 2), tag, WABBIT_COMM, ierr)
!                    call MPI_Send( buffer_data, data_size*(com_N-1), MPI_REAL8, com_plan(k, 2), tag, WABBIT_COMM, ierr)
!
!                !*************** RECV CASE
!                elseif ( com_plan(k, 2) == rank ) then
!                    ! proc have to receive data
!                    ! receive data
!                    call MPI_Recv( buffer_light, com_plan(k, 3), MPI_INTEGER4, com_plan(k, 1), tag, WABBIT_COMM, status, ierr)
!                    call MPI_Recv( buffer_data, data_size*com_plan(k, 3), MPI_REAL8, com_plan(k, 1), tag, WABBIT_COMM, status, ierr)
!
!                    ! loop over all received blocks
!                    do l = 1,  com_plan(k, 3)
!
!                        ! find free "light id", work on reduced light data, so free id is heavy id
!                        call get_free_light_id( hvy_free_id, my_block_list( my_light_start+1 : my_light_start+params%number_blocks , 1 ), params%number_blocks )
!                        ! calculate light id
!                        lgt_free_id = my_light_start + hvy_free_id
!
!                        ! write light data
!                        my_block_list( lgt_free_id, :) = lgt_block( buffer_light(l), : )
!
!                        ! write heavy data
!                        hvy_block(:, :, :, hvy_free_id) = buffer_data(:, :, :, l)
!
!                        if (my_block_list(lgt_free_id, 1)<0 .or. my_block_list(lgt_free_id, 1)>3) then
!                          write(*,*) "For some reason, someone sent me an empty block (code: 7712345)"
!                          stop
!                        endif
!                    end do
!                end if
!            end do
!
!            ! synchronize light data
!            call synchronize_lgt_data( params, lgt_block )

            deallocate( friends, affinity )


        case("sfc_hilbert","sfc_z")
            !---------------------------------------------------------------------------------
            ! 1st: calculate space filling curve index for all blocks
            !---------------------------------------------------------------------------------
            do k = 1, lgt_n
                select case (distribution)
                case("sfc_z")
                    !-----------------------------------------------------------
                    ! Z - curve
                    !-----------------------------------------------------------
                    if (params%threeD_case) then
                        call treecode_to_sfc_id_3D( sfc_id, lgt_block( lgt_active(k), 1:params%max_treelevel ), params%max_treelevel )
                    else
                        call treecode_to_sfc_id_2D( sfc_id, lgt_block( lgt_active(k), 1:params%max_treelevel ), params%max_treelevel )
                    endif

                case("sfc_hilbert")
                    !-----------------------------------------------------------
                    ! Hilbert curve
                    !-----------------------------------------------------------
                    if (params%threeD_case) then
                        ! transfer treecode to hilbertcode
                        call treecode_to_hilbertcode_3D( lgt_block( lgt_active(k), 1:params%max_treelevel ), hilbertcode, params%max_treelevel)
                        ! calculate sfc position from hilbertcode
                        call treecode_to_sfc_id_3D( sfc_id, hilbertcode, params%max_treelevel )
                    else
                        ! transfer treecode to hilbertcode
                        call treecode_to_hilbertcode_2D( lgt_block( lgt_active(k), 1:params%max_treelevel ), hilbertcode, params%max_treelevel)
                        ! calculate sfc position from hilbertcode
                        call treecode_to_sfc_id_2D( sfc_id, hilbertcode, params%max_treelevel )
                    endif

                end select

                ! fill sfc list
                sfc_sorted_list(k, 1) = sfc_id
                sfc_sorted_list(k, 2) = lgt_active(k)
            end do

            ! sort sfc_list according to the first dimension, thus the position on
            ! the space filling curve (this was a bug, fixed: Thomas, 13/03/2018)
            if (lgt_n > 1) then
                call quicksort_ik(sfc_sorted_list, 1, lgt_n, 1, 2)
            end if

            !---------------------------------------------------------------------------------
            ! 2nd: plan communication
            !---------------------------------------------------------------------------------
            ! proc_dist_id: process responsible for current part of sfc
            ! proc_data_id: process who stores data of sfc element

            ! we start the loop on the root rank (0), then assign the first elements
            ! of the SFC, then to second rank, etc. (thus: proc_dist_id is a loop variable)
            proc_dist_id = 0

            ! communication counter. each communication (=send and receive) is stored
            ! in a long list
            com_i = 1

            ! loop over sfc_list
            do k = 1, lgt_n
                ! if the current owner of the SFC is supposed to have zero blocks
                ! then it does not really own this part of the SFC. So we look for the
                ! first rank which is supposed to hold at least one block, and declare it as owner
                ! of this part. NOTE: as we try to minimize communication during send/recv in
                ! load balancing, it may well be that the list of active mpiranks (ie those
                ! which have nonzero number of blocks) is non contiguous, i.e.
                ! opt_dist_list = 1 1 1 0 0 0 0 1 0 1
                ! can happen.
                do while ( opt_dist_list(proc_dist_id+1) == 0 )
                    proc_dist_id = proc_dist_id + 1
                end do

                ! find out on which mpirank lies the block that we're looking at
                call lgt_id_to_proc_rank( proc_data_id, sfc_sorted_list(k,2), params%number_blocks )

                ! does this block lie on the right mpirank, i.e., the current part of the
                ! SFC? if so, nothing needs to be done. otherwise, the following if is active
                if ( proc_dist_id /= proc_data_id ) then
                    ! as this block is one the wrong rank, it will be sent away from its
                    ! current owner (proc_data_id) to the owner of this part of the
                    ! SFC (proc_dist_id)

                    ! save this send+receive operation in the list of planned communications
                    ! column
                    !    1     sender proc
                    !    2     receiver proc
                    !    3     block light data id
                    sfc_com_list(com_i, 1) = proc_data_id           ! sender mpirank
                    sfc_com_list(com_i, 2) = proc_dist_id           ! receiver mpirank
                    sfc_com_list(com_i, 3) = sfc_sorted_list(k,2)   ! light id of block

                    ! increase communication index
                    com_i = com_i + 1
                end if

                ! The opt_dist_list defines how many blocks this rank should have, and
                ! we just treated one (which either already was on the mpirank or will be on
                ! it after communication), so remove one item from the opt_dist_list
                opt_dist_list( proc_dist_id+1 ) = opt_dist_list( proc_dist_id+1 ) - 1

                ! if there is no more blocks to be checked, increase mpirank counter by one
                if ( opt_dist_list( proc_dist_id+1 ) == 0 ) then
                    proc_dist_id = proc_dist_id + 1
                end if
            end do

            !---------------------------------------------------------------------------------
            ! 3rd: actual communication (send/recv)
            !---------------------------------------------------------------------------------
            ! loop over com list, and sen / recv data.
            ! note: delete com list elements after send/receive and if proc not
            ! responsible for com list entry.
            ! NOTE: we try to send as many blocks as possible at one time, not just one
            ! block at a time. The com_list already has a very nice structure, there is no
            ! need to sort them (Thomas, 13/03/2018), on the contrary, sorting ended oup with
            ! more communications
            do k = 1, com_i
                ! com list element is active
                if ( sfc_com_list(k, 1) /= -1 ) then

                    !-----------------------------------------------------------
                    ! I am the sender in this operation?
                    !-----------------------------------------------------------
                    if ( sfc_com_list(k, 1) == rank ) then

                        ! create send buffer, search list, send next l blocks to the same receiver
                        l = 0
                        do while ( (sfc_com_list(k+l, 1) == sfc_com_list(k, 1)) .and. (sfc_com_list(k+l, 2) == sfc_com_list(k, 2)) )
                            lgt_id = sfc_com_list(k+l, 3)
                            ! calculate heavy id from light id
                            call lgt_id_to_hvy_id( heavy_id, lgt_id, rank, params%number_blocks )

                            ! send buffer: fill buffer, heavy data
                            buffer_data(:, :, :, :, l+1 ) = hvy_block(:, :, :, :, heavy_id )
                            ! ... light data
                            buffer_light( l+1 ) = lgt_id

                            !..then delete what I just got rid of
                            hvy_block(:, :, :, :, heavy_id) = 0.0_rk
                            lgt_block(lgt_id, : ) = -1

                            ! go to next element
                            l = l + 1
                        end do

                        ! send data
                        call MPI_Send( buffer_light, l, MPI_INTEGER4, sfc_com_list(k, 2), tag, WABBIT_COMM, ierr)
                        call MPI_Send( buffer_data, data_size*l, MPI_REAL8, sfc_com_list(k, 2), tag, WABBIT_COMM, ierr)

                        ! delete all com list elements
                        sfc_com_list(k:k+l-1, :) = -1

                    !-----------------------------------------------------------
                    ! I am the receiver in this operation?
                    !-----------------------------------------------------------
                    elseif ( sfc_com_list(k, 2) == rank ) then

                        ! count received data sets, recv next l blocks from the sender
                        l = 1
                        do while ( (sfc_com_list(k+l, 1) == sfc_com_list(k, 1)) .and. (sfc_com_list(k+l, 2) == sfc_com_list(k, 2)) )

                            ! delete element
                            sfc_com_list(k+l, :) = -1

                            ! go to next element
                            l = l + 1
                        end do

                        ! receive data
                        call MPI_Recv( buffer_light, l, MPI_INTEGER4, sfc_com_list(k, 1), tag, WABBIT_COMM, status, ierr)
                        call MPI_Recv( buffer_data, data_size*l, MPI_REAL8, sfc_com_list(k, 1), tag, WABBIT_COMM, status, ierr)

                        ! delete first com list element after receiving data
                        sfc_com_list(k, :) = -1

                        ! save comm count
                        com_N = l

                        ! loop over all received blocks
                        do l = 1,  com_N
                            ! fetch a free light id slot on my rank
                            call get_free_local_light_id( params, rank, lgt_block, lgt_free_id )
                            call lgt_id_to_hvy_id( hvy_free_id, lgt_free_id, rank, N )

                            ! copy the data from the buffers
                            lgt_block( lgt_free_id, :) = lgt_block( buffer_light(l), : )
                            hvy_block(:, :, :, :, hvy_free_id) = buffer_data(:, :, :, :, l)

                        end do
                    else
                        ! I am not concerned by this operation, delete element
                        sfc_com_list(k, :) = -1
                    end if
                end if
            end do

            !---------------------------------------------------------------------------------
            ! 4th: synchronize light data
            !---------------------------------------------------------------------------------
            call synchronize_lgt_data( params, lgt_block, refinement_status_only=.false. )

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: block distribution scheme is unknown"
            write(*,*) distribution
            stop

    end select

    ! clean up
    deallocate( opt_dist_list )
    deallocate( dist_list )
    deallocate( com_plan )
    deallocate( sfc_com_list )
    deallocate( buffer_data )
    deallocate( buffer_light )
    deallocate( sfc_sorted_list )

    ! timing
    call toc( params, "balance_load", MPI_wtime()-t0 )
end subroutine balance_load
