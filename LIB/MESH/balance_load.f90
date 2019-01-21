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

subroutine balance_load( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n, hvy_tmp)

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
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank, proc_dist_id, proc_data_id
    ! number of processes
    integer(kind=ik)                    :: number_procs
    ! block distribution lists
    integer(kind=ik), allocatable, save :: opt_dist_list(:), dist_list(:), friends(:,:), affinity(:)
    ! loop variables
    integer(kind=ik)                    :: k, N, l, com_i, com_N, heavy_id, sfc_id, neq
    ! free light/heavy data id
    integer(kind=ik)                    :: lgt_free_id, hvy_free_id
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0, t1
    ! space filling curve list
    integer(kind=ik), allocatable, save :: sfc_com_list(:,:), sfc_sorted_list(:,:)
    ! hilbert code
    integer(kind=ik)                    :: hilbertcode(params%max_treelevel)

!---------------------------------------------------------------------------------------------
! variables initialization

    if (params%number_procs == 1) then
        ! on only one proc, no balancing is required
        return
    endif

    ! start time
    t0 = MPI_Wtime()

    ! MPI_parameters
    rank = params%rank
    number_procs = params%number_procs

    ! allocate block to proc lists
    if (.not.allocated(opt_dist_list)) allocate( opt_dist_list(1:number_procs))
    if (.not.allocated(dist_list)) allocate( dist_list(1:number_procs))
    ! allocate sfc com list, maximal number of communications is when every proc wants to send all of his blocks
    ! NOTE: it is not necessary or wise to reset this array (it is large!)
    if (.not.allocated(sfc_com_list)) allocate( sfc_com_list( number_procs*params%number_blocks, 3 ) )

    ! allocate space filling curve list, number of elements is the number of active blocks
    ! and for each block, we store the space-filling-curve-index and the lgt ID
    if (.not.allocated(sfc_sorted_list)) allocate( sfc_sorted_list( size(lgt_block,1), 2) )

    ! number of blocks
    N = params%number_blocks
    neq = params%n_eqn

!---------------------------------------------------------------------------------------------
! main body


    !---------------------------------------------------------------------------------
    ! First step: define how many blocks each mpirank should have.
    !---------------------------------------------------------------------------------
    call set_desired_num_blocks_per_rank(params, dist_list, opt_dist_list, lgt_n, hvy_n)
!    call write_block_distribution( params, dist_list, "block_dist.dat" )

    ! at this point, we know how many blocks a mpirank has: "dist_list(myrank+1)"
    ! and how many it should have, if equally distributed: "opt_dist_list(myrank+1)"


    select case(params%block_distribution)
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
            t1 = MPI_wtime()
            do k = 1, lgt_n
                select case (params%block_distribution)
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
            call toc( params, "balance_load (SFC+sort)", MPI_wtime()-t1 )

            !---------------------------------------------------------------------------------
            ! 2nd: plan communication (fill list of blocks to transfer)
            !---------------------------------------------------------------------------------
            t1 = MPI_wtime()
            ! proc_dist_id: process responsible for current part of sfc
            ! proc_data_id: process who stores data of sfc element

            ! we start the loop on the root rank (0), then assign the first elements
            ! of the SFC, then to second rank, etc. (thus: proc_dist_id is a loop variable)
            proc_dist_id = 0

            ! communication counter. each communication (=send and receive) is stored
            ! in a long list
            com_i = 0

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
                    com_i = com_i + 1
                    sfc_com_list(com_i, 1) = proc_data_id           ! sender mpirank
                    sfc_com_list(com_i, 2) = proc_dist_id           ! receiver mpirank
                    sfc_com_list(com_i, 3) = sfc_sorted_list(k,2)   ! light id of block
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
            call block_xfer( params, sfc_com_list, com_i, lgt_block, hvy_block, lgt_active, &
                 lgt_n, lgt_sortednumlist, hvy_tmp )
            call toc( params, "balance_load (comm)", MPI_wtime()-t1 )

        case default
            call abort(2009182147, "[balance_load.f90] ERROR: block distribution scheme is unknown")

    end select

    ! timing
    call toc( params, "balance_load (TOTAL)", MPI_wtime()-t0 )
end subroutine balance_load
