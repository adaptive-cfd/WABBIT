! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: balance_load.f90
! version: 0.4
! author: msr
!
! balance the load
!
! input:    - params, light and heavy data
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
! 16/11/2016  Avoid some communication by more carefully distributing the excess blocks
!
! = TODO ======================================================================================
! * generalzie output so that it can hanlde any number of mpirnaks to ascii file
! ********************************************************************************************

subroutine balance_load( params, block_list, block_data, neighbor_list )

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    ! heavy data array - neifghbor data
    integer(kind=ik), intent(inout)       :: neighbor_list(:)

    ! send/receive buffer, note: size is equal to block data array, because if a block want to send all his data
    real(kind=rk)                       :: buffer_data( size(block_data,1), size(block_data,2), size(block_data,3), size(block_data,4) )
    integer(kind=ik)                    :: buffer_light( params%number_blocks )

    ! light data list for working
    integer(kind=ik)                    :: my_block_list( size(block_list, 1), params%max_treelevel+2)
    ! light id start
    integer(kind=ik)                    :: my_light_start

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank, proc_id
    ! number of processes
    integer(kind=ik)                    :: number_procs
    ! MPI message tag
    integer(kind=ik)                    :: tag
    ! MPI status
    integer                             :: status(MPI_status_size)

    ! distribution type
    character(len=80)                   :: distribution

    ! block distribution lists
    integer(kind=ik), allocatable       :: opt_dist_list(:), dist_list(:), friends_loc(:,:), friends(:,:), affinity(:)

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! loop variables
    integer(kind=ik)                    :: k, N, num_blocks, l, com_i, com_N, excess_blocks, avg_blocks, &
                                           id_send, id_recv, send_deficit, recv_deficit, tmp(1), q, r, light_id, heavy_id, &
                                           direc_id

    ! com plan
    integer(kind=ik), allocatable       :: com_plan(:,:)

    ! size of data array
    integer(kind=ik)                    :: data_size

    ! free light/heavy data id
    integer(kind=ik)                    :: free_light_id, free_heavy_id

!---------------------------------------------------------------------------------------------
! interfaces
interface
  subroutine compute_affinity(params, block_list, block_data, neighbor_list, rank, rank_partner, affinity)
    use mpi
    use module_params
    implicit none

    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(inout)     :: block_list(:, :)
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    integer(kind=ik), intent(out)       :: neighbor_list(:), affinity(:)
    integer(kind=ik), intent(in) :: rank, rank_partner
  end subroutine
  subroutine set_desired_num_blocks_per_rank(params, block_list, dist_list, opt_dist_list)
    use module_params
    use mpi
    implicit none
    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(in)     :: block_list(:, :)
    integer(kind=ik), intent(out)    :: dist_list(:), opt_dist_list(:)
  end subroutine
  subroutine compute_friends_table(params, neighbor_list, friends)
    use module_params
    use mpi
    type (type_params), intent(in) :: params
    integer(kind=ik),intent(inout) :: friends(:,:)
    ! heavy data array - neifghbor data
    integer(kind=ik), intent(in) :: neighbor_list(:)
  end subroutine
end interface

!---------------------------------------------------------------------------------------------
! variables initialization

    tag = 0

    distribution    = params%block_distribution

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    ! determinate process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)

    ! allocate block to proc lists
    allocate( opt_dist_list(1:number_procs), dist_list(1:number_procs))
    allocate( affinity(1:params%number_blocks) )
    allocate( friends( 1:number_procs, 1:number_procs ))
    ! allocate com plan, maximal number of communications: every proc send every other proc something
    allocate( com_plan( number_procs*(number_procs-1), 3 ), stat=allocate_error )
    com_plan = -1

    ! number of light data (since the light data is redunantly stored on all CPU,
    ! this number corresponds usually to the maximum number of blocks possible in
    ! the entire system)
    N = size(block_list, 1)

    ! reset block count
    num_blocks = 0


    ! light data start line
    my_light_start = rank*params%number_blocks

    ! set light data list for working, only light data coresponding to proc are not zero
    my_block_list = 0
    my_block_list( my_light_start+1: my_light_start+params%number_blocks, :) = block_list( my_light_start+1: my_light_start+params%number_blocks, :)

    ! size of data array, use for readability
    data_size = size(block_data,1) * size(block_data,2) * size(block_data,3)

    ! reset send/receive buffer
    buffer_data = 9.0e9_rk

!---------------------------------------------------------------------------------------------
! main body

    select case(distribution)
    case("equal")! simple uniformly distribution
            !---------------------------------------------------------------------------------
            ! First step: define how many blocks each mpirank should have.
            !---------------------------------------------------------------------------------
            call set_desired_num_blocks_per_rank(params, block_list, dist_list, opt_dist_list)
            ! at this point, we know how many blocks a mpirank has: "dist_list(myrank)"
            ! and how many it should have, if equally distributed: "opt_dist_list(myrank)"

            ! determine matrix of number of neighbor relations between mpiranks
            call compute_friends_table(params, neighbor_list, friends)

            !---------------------------------------------------------------------------------
            ! second step: create plan for communication
            !---------------------------------------------------------------------------------
            ! column
            !    1     sender proc
            !    2     receiver proc
            !    3     number of blocks to send

            com_i = 1
            ! loop over all procs
            do id_send = 1, number_procs
                send_deficit = opt_dist_list(id_send) - dist_list(id_send)
                ! find a receiver, if this is a SENDER
                do while (send_deficit < 0)
                    ! from my list of friends, gather my currently best friend
                    tmp = maxloc( friends(id_send,:) )
                    id_recv = tmp(1)
                    ! "remove" the friend we're looking at from the list
                    friends(id_send, id_recv) = -2

                    ! can this proc take some of my blocks?
                    recv_deficit = opt_dist_list(id_recv) - dist_list(id_recv)
                    ! proc can receive data
                    if ( ( recv_deficit > 0 ) .and. ( id_recv /= id_send ) ) then
                        ! calculate number of send/receive blocks
                        com_N = minval( abs((/send_deficit, recv_deficit/)) )

                        ! create com plan
                        com_plan( com_i, 1 ) = id_send - 1
                        com_plan( com_i, 2 ) = id_recv - 1
                        com_plan( com_i, 3 ) = com_N

                        ! change distribution list
                        dist_list(id_send) = dist_list(id_send) - com_N
                        dist_list(id_recv) = dist_list(id_recv) + com_N

                        ! new com was created
                        com_i = com_i + 1
                    end if
                    ! recompute senders deficit (how many blocks sender still wants to get rid of)
                    send_deficit = opt_dist_list(id_send) - dist_list(id_send)
                end do
            end do
            ! we counted one too far
            com_i = com_i -1

            ! NOTE
            ! The following relations hold
            ! light_id(heavy_id) = rank*params%number_blocks + heavy_id
            ! heavy_id(light_id) = light_id - rank*params%number_blocks
            ! heavy_id(direc_id) = (direc_id-1)/16 + 1
            ! direc_id(heavy_id,dir) = (heavy_id-1)*16 + dir  (and dir=[1,16])

            !---------------------------------------------------------------------------------
            ! third: actually send/receive data
            !---------------------------------------------------------------------------------
            l = 1
            ! loop over all planed communications
            do k = 1, com_i
                !*************** SEND CASE
                if ( com_plan(k, 1) == rank ) then
                    ! yes, and I'll send "com_plan(k, 3)" blocks to proc "com_plan(k, 2)"
                    !---------------------------------------------------------------------------------------
                    ! affinity list, HEAVY DATA ARRAY
                    call compute_affinity(params, my_block_list, block_data, neighbor_list, rank, com_plan(k, 2), affinity)

                    com_N = 1
                    do while ( com_N <= com_plan(k, 3) )
                        ! fetch most desirable block
                        tmp = maxloc(affinity)
                        heavy_id = tmp(1)
                        ! we now use this block and must be sure not to use it again
                        affinity(heavy_id) = -99
                        light_id = rank*params%number_blocks + heavy_id

                        if (my_block_list(light_id, 1) /= -1) then! block is active
                            ! fill buffer, heavy data
                            buffer_data(:, :, :, com_N) = block_data(:, :, :, heavy_id)
                            ! ... light data
                            buffer_light(com_N) = light_id
                            ! count com
                            com_N = com_N + 1
                            ! delete heavy data
                            block_data(:, :, :, heavy_id) = 0.0_rk
                            ! delete light data
                            my_block_list( light_id, : ) = -1
                            ! write(*,'("rank=",i2," sends light=",i4," heayv=",i4," to ",i2," affinity=",i4)')&
                            ! rank,light_id,heavy_id,com_plan(k, 2), tmp(1)
                        else
                          write(*,*) my_block_list(light_id, 1), block_list(light_id, 1)
                          stop "unforeseen 576879"
                        end if
                    end do

                    ! error case
                    if ( com_N-1 /= com_plan(k, 3) ) then
                        stop "ERROR: load balancing error 003857"
                    end if

                    ! send data
                    call MPI_Send( buffer_light, (com_N-1), MPI_INTEGER4, com_plan(k, 2), tag, MPI_COMM_WORLD, ierr)
                    call MPI_Send( buffer_data, data_size*(com_N-1), MPI_REAL8, com_plan(k, 2), tag, MPI_COMM_WORLD, ierr)

                !*************** RECV CASE
                elseif ( com_plan(k, 2) == rank ) then
                    ! proc have to receive data
                    ! receive data
                    call MPI_Recv( buffer_light, com_plan(k, 3), MPI_INTEGER4, com_plan(k, 1), tag, MPI_COMM_WORLD, status, ierr)
                    call MPI_Recv( buffer_data, data_size*com_plan(k, 3), MPI_REAL8, com_plan(k, 1), tag, MPI_COMM_WORLD, status, ierr)

                    ! loop over all received blocks
                    do l = 1,  com_plan(k, 3)

                        ! find free light id
                        call get_free_light_id( free_heavy_id, my_block_list( my_light_start+1 : my_light_start+params%number_blocks , 1 ), params%number_blocks )
                        ! calculate light id
                        free_light_id = my_light_start + free_heavy_id

                        ! write light data
                        my_block_list( free_light_id, :) = block_list( buffer_light(l), : )

                        ! write heavy data
                        block_data(:, :, :, free_heavy_id) = buffer_data(:, :, :, l)

                    end do
                end if !end recv case
                ! synchronize light data
                block_list = 0
                call MPI_Allreduce(my_block_list, block_list, size(block_list,1)*size(block_list,2), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

            end do

            ! synchronize light data
            block_list = 0
            call MPI_Allreduce(my_block_list, block_list, size(block_list,1)*size(block_list,2), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

        case("sfc1")

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

end subroutine balance_load

!-------------------------------------------------------------------------------
! compute the affinity for all my blocks to a given rank.
! That means, if a block has many neighbor relations with the target rank, it has a high
! value in the array, otherwise a low value or -10 if the block is not active at all
!
! This list is at the core of heuristic load balancing: we know whom to send to, and using this list
! we decide which blocks will be sent.
!-------------------------------------------------------------------------------
subroutine compute_affinity(params, block_list, block_data, neighbor_list, rank, rank_partner, affinity)
  use mpi
  use module_params
  implicit none
  ! user defined parameter structure
  type (type_params), intent(in)      :: params
  ! light data array
  integer(kind=ik), intent(inout)     :: block_list(:, :)
  ! heavy data array - block data
  real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
  ! heavy data array - neifghbor data
  integer(kind=ik), intent(out)       :: neighbor_list(:), affinity(:)
  integer(kind=ik), intent(in)        :: rank, rank_partner
  integer :: heavy_id, light_id,direc_id, proc_id, q

  ! affinity list, HEAVY DATA ARRAY
  affinity = 0

  if(size(affinity)/=params%number_blocks) then
    write(*,*)"unforeseen"
    stop
  endif

  ! loop over all heavy data on my rank
  do heavy_id = 1,params%number_blocks
    ! get corresponding light_id
    light_id = rank*params%number_blocks + heavy_id
    if (block_list(light_id,1) == -1) then
      ! inactive block have extremely low affinity to be transferred
      affinity(heavy_id) = -1
    else
      ! loop over all directions and count how many neighbor relations I do have with the rank_partner
      do q = 1, 16 ! q is relative direction
        direc_id = (heavy_id-1)*16 + q ! where is this direction for this heavy_id in the neighbor list
        proc_id = (neighbor_list(direc_id) / params%number_blocks) + 1 ! one based proc_id
        if (proc_id-1 == rank_partner) then !receiver? ZERO BASED
          ! a shared border with the target rank is a high priority
          affinity(heavy_id) = affinity(heavy_id)+20
        elseif (proc_id-1 /= rank) then
          ! so I dont share this border with the target rank, but it is an mpi border
          ! when I'm out of good candidates, I should at least send blocks that are not surrounded only by my blocks
          affinity(heavy_id) = affinity(heavy_id)+1
        endif
      enddo
    end if
  end do

end subroutine compute_affinity




subroutine set_desired_num_blocks_per_rank(params, block_list, dist_list, opt_dist_list)
  use module_params
  use mpi
  implicit none
  ! user defined parameter structure
  type (type_params), intent(in)      :: params
  ! light data array
  integer(kind=ik), intent(in)     :: block_list(:, :)
  integer(kind=ik), intent(out)    :: dist_list(:), opt_dist_list(:)
  ! integer(kind=ik), intent(out)    :: currnt_BlocksPerRank(:), desired_BlocksPerRank(:)
  integer :: k, num_blocks, proc_id, avg_blocks, number_procs, rank, excess_blocks, ierr
  integer :: N_light

  ! determinate process rank
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  ! determinate process number
  call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)

  dist_list = 0
  opt_dist_list = 0
  N_light = size(block_list,1)
  ! count number of active blocks and current block distribution
  do k = 1, N_light
      ! block is active? we can tell by the first treecode entry. if it is -1
      ! then the block is not active
      if ( block_list(k, 1) /= -1 ) then
          ! count block for proc corresponing to light data
          ! the array dist_list has the length of the number of mpi procs and holds the number
          ! of blocks used on this process. From the light data ID "k", we can compute the mpirank
          ! of the process holding the block "proc_id"
          proc_id = (k-1) / params%number_blocks + 1 ! one based index
          dist_list( proc_id ) = dist_list( proc_id ) + 1
      end if
  end do
  ! count global number of blocks on all mpiranks
  num_blocks = sum(dist_list)

  ! optimal distribution of blocks per mpirank. The simple division of "num_blocks" by "number_procs" actually
  ! yields a double (since it is not guaranteed that all mpiranks hold the exact same number of blocks)
  ! using the integer division, decimal places are cut
  avg_blocks = num_blocks / number_procs
  opt_dist_list(:) = avg_blocks

  ! some blocks are missing due to the integer division
  excess_blocks = num_blocks - sum(opt_dist_list)
  ! distribute remaining blocks (the excess blocks, if we have some)
  do while ( excess_blocks>0 )
      ! first we try to be clever and increase the counter of "desired" blocks for
      ! procs that already have more blocks than they should (by one)
      do proc_id = 1, number_procs
        ! check if this proc_id has more blocks than it is supposed to and if so, we attribute it one of the excess blocks
        if ( dist_list(proc_id) > avg_blocks) then
            opt_dist_list(proc_id) = opt_dist_list(proc_id) + 1
            ! we got rid of one excess block
            excess_blocks = excess_blocks - 1
            ! no more blocks to distribute?
            if (excess_blocks==0) exit
        end if
      end do

      ! no more blocks to distribute?
      if (excess_blocks==0) exit

      ! second, it may be that this is not enough: there are still bocks to be
      ! distributed. so now we repeat the loop, but look for mpiranks that have
      ! enough blocks and give them one more.
      do proc_id = 1, number_procs
        if ( dist_list(proc_id) == avg_blocks) then
            opt_dist_list(proc_id) = opt_dist_list(proc_id) + 1
            excess_blocks = excess_blocks - 1
            if (excess_blocks==0) exit
        end if
      end do

      ! no more blocks to distribute?
      if (excess_blocks==0) exit

      ! third, it may still not be enough...so just pick some
      do proc_id = 1, number_procs
        if ( dist_list(proc_id) < avg_blocks) then
            opt_dist_list(proc_id) = opt_dist_list(proc_id) + 1
            excess_blocks = excess_blocks - 1
            if (excess_blocks==0) exit
        end if
      end do
  end do ! end of excess block distribution

  if (rank==0) then
    ! error checking. the sum of newly distributed blocks must of course be
    ! the same as the number we had before distribution
    if (sum(opt_dist_list)/=num_blocks .or. maxval(abs(opt_dist_list-avg_blocks))>1) then
        write(*,*) "something went wrong - during balancing, we lost or gained some blocks", excess_blocks
        write(*,*) "or we have more than +-1 block difference among them"
        write(*,*) opt_dist_list
        stop
    end if
    ! transfer information. for each time step dump to disk how many blocks
    ! we have to transfer
    open(14,file='load_balancing.t',status='unknown',position='append')
    write(14,'(10(i5,1x))') opt_dist_list-dist_list
    close(14)

    open(14,file='blocks_per_rank.t',status='unknown',position='append')
    write(14,'(10(i5,1x))') dist_list
    close(14)
  end if

end subroutine


! compute friends table. This is a mpisize x mpisize array of integers, and it counts
! how many neighbor relations the mpiranks have to each other mpirank. the matrix is
! (most of the time) symmetric, exceptions from symmetry can occur in the case of finer/coarser
! relations
! In the friends array, we thus just know who's my most important neighbor. The array can look like this: (mpisize=4)
! -1    3    2    4
!  3   -1   22   12
!  2   22   -1   27
!  4   12   27   -1
! Note we set (-1) on the diagonal (since an mpirank is likely its own best friend, they're egocentric in that regard)
subroutine compute_friends_table(params, neighbor_list, friends)
  use module_params
  use mpi
  implicit none
  type (type_params), intent(in) :: params
  integer(kind=ik),intent(inout) :: friends(:,:)
  integer(kind=ik),intent(in) :: neighbor_list(:)


  integer(kind=ik), allocatable :: friends_loc(:,:)
  integer(kind=ik) :: k, proc_id, number_procs, rank, ierr
  
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)

  ! We need the buffer for the friends array:
  allocate( friends_loc( 1:number_procs, 1:number_procs ))
  friends_loc = 0
  friends = 0


  ! loop over all entries in the neighbor array. note this is currently one-D and contains
  ! 16 entries per heavy data.
  do k = 1, 16*params%number_blocks
    ! is the block active?
    if (neighbor_list(k)/=-1) then
      ! one-based indexing
      proc_id = (neighbor_list(k) / params%number_blocks) + 1
      friends_loc(rank+1, proc_id) = friends_loc(rank+1, proc_id)+1
    endif
  enddo
  ! overwrite my own.. (rank is zero-based)
  friends_loc(rank+1,rank+1) = -1

  ! the friends must be known on all mpiranks
  call MPI_Allreduce(friends_loc, friends, number_procs**2, MPI_INTEGER,MPI_SUM, MPI_COMM_WORLD, ierr)

  deallocate(friends_loc)
end subroutine
