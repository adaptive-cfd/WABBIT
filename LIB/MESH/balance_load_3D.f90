! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: balance_load_3D.f90
! version: 0.5
! author: msr
!
! balance the load
!
! input:    - params, light and heavy data, neighbor data, lists of active blocks
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 09/02/17    - create
!
! ********************************************************************************************

subroutine balance_load_3D( params, lgt_block, hvy_block, lgt_active, lgt_n)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)

    ! list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n

    ! send/receive buffer, note: size is equal to block data array, because if a block want to send all his data
    real(kind=rk)                       :: buffer_data( size(hvy_block,1), size(hvy_block,2), size(hvy_block,3), size(hvy_block,4), size(hvy_block,5) )
    integer(kind=ik)                    :: buffer_light( params%number_blocks )

    ! light data list for working
    integer(kind=ik)                    :: my_block_list( size(lgt_block, 1), params%max_treelevel+2)
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

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! loop variables
    integer(kind=ik)                    :: k, N, num_blocks, l, com_i, com_N, heavy_id, sfc_id

    ! com plan
    integer(kind=ik), allocatable       :: com_plan(:,:)

    ! size of data array
    integer(kind=ik)                    :: data_size

    ! free light/heavy data id
    integer(kind=ik)                    :: free_light_id, free_heavy_id

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

    ! space filling curve list
    integer(kind=ik), allocatable       :: sfc_list(:), sfc_com_list(:,:)

    ! hilbert code
    !integer(kind=ik)                    :: hilbertcode(params%max_treelevel)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameter
    rank            = params%rank
    number_procs    = params%number_procs

    ! start time
    sub_t0 = MPI_Wtime()

    tag = 0

    distribution    = params%block_distribution

    ! allocate block to proc lists
    allocate( opt_dist_list(1:number_procs), dist_list(1:number_procs))

    allocate( affinity(1:params%number_blocks) )

    allocate( friends( 1:number_procs, 1:number_procs ))

    ! allocate com plan, maximal number of communications: every proc send every other proc something
    allocate( com_plan( number_procs*(number_procs-1), 3 ), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if
    com_plan = -1

    ! allocate sfc com list, maximal number of communications: every proc send all of his blocks
    allocate( sfc_com_list( number_procs*params%number_blocks, 3 ), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if
    sfc_com_list = -1

    ! allocate space filling curve list, maximal number of elements = max number of blocks
    allocate( sfc_list( 8**params%max_treelevel ), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    ! number of light data (since the light data is redunantly stored on all CPU,
    ! this number corresponds usually to the maximum number of blocks possible in
    ! the entire system)
    N = size(lgt_block, 1)

    ! reset block count
    num_blocks = 0

    ! light data start line
    my_light_start = rank*params%number_blocks

    ! set light data list for working, only light data coresponding to proc are not zero
    my_block_list = 0
    my_block_list( my_light_start+1: my_light_start+params%number_blocks, :) = lgt_block( my_light_start+1: my_light_start+params%number_blocks, :)

    ! size of data array, use for readability
    data_size = size(hvy_block,1) * size(hvy_block,2) * size(hvy_block,3) * size(hvy_block,4)

    ! reset send/receive buffer
    buffer_data = 9.0e9_rk

!---------------------------------------------------------------------------------------------
! main body

    select case(distribution)

        case("sfc_z")

            ! current block distribution
            call set_desired_num_blocks_per_rank(params, dist_list, opt_dist_list, lgt_active, lgt_n)
            ! write debug infos: current distribution list
            if ( params%debug ) then
                call write_block_distribution( dist_list )
            end if

            !---------------------------------------------------------------------------------
            ! first: calculate space filling curve
            !---------------------------------------------------------------------------------
            ! reset old lists
            sfc_list  = -1
            dist_list = 0

            ! loop over active blocks
            do k = 1, lgt_n
                ! calculate sfc position
                call treecode_to_3D_z_curve( sfc_id, lgt_block( lgt_active(k), 1:params%max_treelevel ), params%max_treelevel )
                ! fill sfc list with light data id
                sfc_list(sfc_id+1) = lgt_active(k)
            end do

            !---------------------------------------------------------------------------------
            ! second: distribute all blocks
            !---------------------------------------------------------------------------------
            ! equal distribution
            dist_list = lgt_n / number_procs
            do k = 1, mod(lgt_n, number_procs)
                dist_list( k ) = dist_list( k ) + 1
            end do

            !---------------------------------------------------------------------------------
            ! third: create com list
            !---------------------------------------------------------------------------------
            ! column
            !    1     sender proc
            !    2     receiver proc
            !    3     block light data id

            ! proc_dist_id: process responsible for current part of sfc
            ! proc_data_id: process who stores data of sfc element
            proc_dist_id = 0

            com_i = 1
            ! loop over sfc_list
            do k = 1, size(sfc_list,1)
                ! sfc element is active
                if ( sfc_list(k) /= -1 ) then

                    ! process with heavy data
                    call lgt_id_to_proc_rank( proc_data_id, sfc_list(k), params%number_blocks )

                    ! data has to send
                    if ( proc_dist_id /= proc_data_id ) then
                        ! create com plan
                        sfc_com_list(com_i, 1) = proc_data_id
                        sfc_com_list(com_i, 2) = proc_dist_id
                        sfc_com_list(com_i, 3) = sfc_list(k)
                        com_i = com_i + 1

                    else
                        ! nothing to do, block is allready on correct proc

                    end if

                    ! next scf element, so check proc id, switch if last block is distributed
                    dist_list( proc_dist_id+1 ) = dist_list( proc_dist_id+1 ) - 1
                    if ( dist_list( proc_dist_id+1 ) == 0 ) then
                        proc_dist_id = proc_dist_id + 1
                    end if

                end if
            end do

            ! stop load balancing, if nothing to do
            if ( com_i == 1 ) then
                ! the distribution is fine, nothing to do.
                return
            endif

            !---------------------------------------------------------------------------------
            ! fourth: communicate
            !---------------------------------------------------------------------------------
            ! loop over com list, create send buffer and send/receive data
            ! note: delete com list elements after send/receive and if proc not
            ! responsible for com list entry -> this means: no extra com plan needed
            do k = 1, com_i
                ! com list element is active
                if ( sfc_com_list(k, 1) /= -1 ) then

                    ! proc send data
                    if ( sfc_com_list(k, 1) == rank ) then

                        ! create send buffer, search list
                        l = 0
                        do while ( (sfc_com_list(k+l, 1) == sfc_com_list(k, 1)) .and. (sfc_com_list(k+l, 2) == sfc_com_list(k, 2)) )

                            ! calculate heavy id from light id
                            call lgt_id_to_hvy_id( heavy_id, sfc_com_list(k+l, 3), rank, params%number_blocks )

                            ! send buffer: fill buffer, heavy data
                            buffer_data(:, :, :, :, l+1 ) = hvy_block(:, :, :, :, heavy_id )
                            ! ... light data
                            buffer_light( l+1 ) = sfc_com_list(k+l, 3)

                            ! delete heavy data
                            hvy_block(:, :, :, :, heavy_id) = 0.0_rk
                            ! delete light data
                            my_block_list( sfc_com_list(k+l, 3) , : ) = -1

                            ! go to next element
                            l = l + 1

                        end do

                        ! send data
                        call MPI_Send( buffer_light, (l), MPI_INTEGER4, sfc_com_list(k, 2), tag, MPI_COMM_WORLD, ierr)
                        call MPI_Send( buffer_data, data_size*(l), MPI_REAL8, sfc_com_list(k, 2), tag, MPI_COMM_WORLD, ierr)

                        ! delete all com list elements
                        sfc_com_list(k:k+l-1, :) = -1

                    ! proc receive data
                    elseif ( sfc_com_list(k, 2) == rank ) then

                        ! count received data sets
                        l = 1
                        do while ( (sfc_com_list(k+l, 1) == sfc_com_list(k, 1)) .and. (sfc_com_list(k+l, 2) == sfc_com_list(k, 2)) )

                            ! delete element
                            sfc_com_list(k+l, :) = -1

                            ! go to next element
                            l = l + 1

                        end do

                        ! receive data
                        call MPI_Recv( buffer_light, (l), MPI_INTEGER4, sfc_com_list(k, 1), tag, MPI_COMM_WORLD, status, ierr)
                        call MPI_Recv( buffer_data, data_size*(l), MPI_REAL8, sfc_com_list(k, 1), tag, MPI_COMM_WORLD, status, ierr)

                        ! delete first com list element after receiving data
                        sfc_com_list(k, :) = -1

                        ! save comm count
                        com_N = l

                        ! loop over all received blocks
                        do l = 1,  com_N

                            ! find free "light id", work on reduced light data, so free id is heavy id
                            call get_free_light_id( free_heavy_id, my_block_list( my_light_start+1 : my_light_start+params%number_blocks , 1 ), params%number_blocks )
                            ! calculate light id
                            free_light_id = my_light_start + free_heavy_id

                            ! write light data
                            my_block_list( free_light_id, :) = lgt_block( buffer_light(l), : )

                            ! write heavy data
                            hvy_block(:, :, :, :, free_heavy_id) = buffer_data(:, :, :, :, l)

                            ! error case
                            if (my_block_list(free_light_id, 1)<0 .or. my_block_list(free_light_id, 1)>7) then
                              write(*,*) "For some reason, someone sent me an empty block (code: 7712345)"
                              stop
                            endif

                        end do


                    ! nothing to do
                    else
                        ! delete com list element
                        ! note: only to have a clean list
                        sfc_com_list(k, :) = -1

                    end if

                end if
            end do

            !---------------------------------------------------------------------------------
            ! sixth: synchronize light data
            !---------------------------------------------------------------------------------
            lgt_block = 0
            call MPI_Allreduce(my_block_list, lgt_block, size(lgt_block,1)*size(lgt_block,2), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: block distribution scheme is unknown"
            write(*,*) distribution
            stop

    end select

    ! clean up
    deallocate( friends, affinity )
    deallocate( opt_dist_list, stat=allocate_error )
    deallocate( dist_list, stat=allocate_error )
    deallocate( com_plan, stat=allocate_error )
    deallocate( sfc_list, stat=allocate_error )
    deallocate( sfc_com_list, stat=allocate_error )

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "balance_load" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "balance_load"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine balance_load_3D
