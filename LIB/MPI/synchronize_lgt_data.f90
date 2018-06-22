!> \file
!> \brief
!! Synchronize ligt data array among all mpiranks
!> \details
!-------------------------------------------------------------------------------
! the light data array looks like this (here for 2 mpiranks) (x=used .=unused)
!
! mpirank=* x..xx.xx...xx....................xx.x.x.x..x........................
!
! if both mpiranks do something with their light data, it changes:
!
! mpirank=0 x..xxxxxxxxxx.......................................................
! mpirank=1 .................................xx.x.x.............................
!
! not you *could* use a simple MPI_ALLREDUCE with lgt_N complexity (using MPI_MAX
! or MPI_SUM), but this can get very expensive.
! better do something smart.
!
!   NOTE: one idea implemented here is that the lgt_block data is actually only
!         small numbers (0..7 for octree, +-1 for refinement) so we can save a
!         lot of transfer volume by converting it to kind=1 integer. This should
!         one day be done throughout the code.
!
!   NOTE: since, by definition, the light data is OUTDATED at this point, you
!         cannot use active lists here since they are likewise outdated.
!
!-------------------------------------------------------------------------------
subroutine synchronize_lgt_data( params, lgt_block, refinement_status_only )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> INPUT: light data array, locally modified (=valid only on my section of the array)
    !> OUTPUT: light data array, synchronized on all procs (=valid everywhere)
    integer(kind=ik), intent(inout) :: lgt_block(:, :)
    !> Some operations change the treecodes etc and have to sync the entire array
    !> of light data. But sometimes, only the refinement status is altered: in that case
    !> we should communicate only that, of course, to save time
    logical, intent(in) :: refinement_status_only

    ! local variables
    integer(kind=ik) :: mpisize, mpirank, N, lgt_start, lgt_end, lgt_id, ierr, &
    buffer_size, lgt_num, buffer_start, k, R
    ! send/receive buffer for data synchronization
    integer(kind=1), allocatable, save      :: my_lgt_block_send_buffer(:,:), my_lgt_block_recv_buffer(:,:)
    ! kind=1 integer copy of light data, which will only hold my data
    integer(kind=1), allocatable, save      :: my_lgt_block(:, :)
    integer(kind=ik), allocatable, save     :: proc_lgt_num(:), my_proc_lgt_num(:)

    mpirank = params%rank
    mpisize = params%number_procs
    N = params%number_blocks
    R = params%max_treelevel+2

    if (.not.allocated(proc_lgt_num)) allocate( proc_lgt_num(1:mpisize) )
    if (.not.allocated(my_proc_lgt_num)) allocate( my_proc_lgt_num(1:mpisize) )
    if (.not.allocated(my_lgt_block)) allocate( my_lgt_block( size(lgt_block,1), size(lgt_block,2)) )
    if (.not.allocated(my_lgt_block_send_buffer)) allocate( my_lgt_block_send_buffer( size(lgt_block,1), size(lgt_block,2) ) )
    if (.not.allocated(my_lgt_block_recv_buffer)) allocate( my_lgt_block_recv_buffer( size(lgt_block,1), size(lgt_block,2) ) )

    ! ==========================================================================
    ! view on lgt_block on one mpirank:
    !
    ! mpirank=7 ........................xx.x.x.xx.x...........................
    !                                   ^         ^                          ^
    !                                   |         |                          |
    !                                lgt_start   lgt_end                lgt_end2
    !
    ! Thus I have to send
    !       lgt_num = lgt_end - lgt_start + 1
    ! data points. NOTE: this is not the same as
    !       N = lgt_end2 - lgt_start + 1
    ! since it includes possibile holes.
    ! ==========================================================================

    call proc_to_lgt_data_start_id( lgt_start, mpirank, N )
    ! first we copy the entire part to the new array (here lgt_end === lgt_end2 in the ascii art)
    ! even though we could also try to just copy the interior part
    lgt_end = lgt_start + N - 1
    ! copy and reduce precision to kind=1 integer
    my_lgt_block(lgt_start:lgt_end, :) = int( lgt_block(lgt_start:lgt_end, :), kind=1)

    ! fetch last used light id, that is the value of lgt_end in the above ascii art
    ! this will be the interval which we have to communicate via MPI
    do lgt_id = lgt_end, lgt_start, -1
        if (my_lgt_block(lgt_id,1) /= -1) exit
    enddo
    lgt_end = lgt_id
    lgt_num = lgt_end - lgt_start + 1


    ! Next, we figure out how much the union of all the different subsets (lgt_start:lgt_end)
    ! on all the mpiranks.
    my_proc_lgt_num(:) = 0
    my_proc_lgt_num(mpirank+1) = lgt_num
    ! synchronize array
    call MPI_Allreduce(my_proc_lgt_num, proc_lgt_num, mpisize, MPI_INTEGER4, MPI_SUM, WABBIT_COMM, ierr)
    ! this is the global buffer size
    buffer_size = sum(proc_lgt_num)


    ! ! now we can allocate send/receive buffer arrays
    ! allocate( my_lgt_block_send_buffer( 1:buffer_size, size(lgt_block,2) ) )
    ! allocate( my_lgt_block_recv_buffer( buffer_size, size(lgt_block,2) ) )

    ! ==========================================================================
    ! The data in my_lgt_block:
    ! mpirank=0 aaaa............................................................
    ! mpirank=1 ..............bbbb..............................................
    ! mpirank=2 ..........................cccccc................................
    ! mpirank=3 .....................................................dd.........
    !
    ! We remove most of the holes (in my_lgt_block_send_buffer)
    ! mpirank=0 aaaa............
    ! mpirank=1 ....bbbb........
    ! mpirank=2 ........cccccc..
    ! mpirank=3 ..............dd
    !
    ! And synchronize that (in my_lgt_block_recv_buffer)
    ! buffer = aaaabbbbccccccdd
    !
    ! NOTE: aaaa can contain some wholes, but is not expected to be very hollow
    !
    ! ==========================================================================

    ! fill the buffer
    buffer_start = sum(proc_lgt_num(1:mpirank))
    my_lgt_block_send_buffer(1:buffer_size,:) = 0
    my_lgt_block_send_buffer( buffer_start+1 : buffer_start+lgt_num, : ) = my_lgt_block( lgt_start:lgt_end, :)

    ! synchronize the buffer
    if (refinement_status_only) then
        call MPI_Allreduce(my_lgt_block_send_buffer(1:buffer_size,R), my_lgt_block_recv_buffer(1:buffer_size,R),&
        buffer_size, MPI_INTEGER1, MPI_SUM, WABBIT_COMM, ierr)
    else
        call MPI_Allreduce(my_lgt_block_send_buffer(1:buffer_size,:), my_lgt_block_recv_buffer(1:buffer_size,:), buffer_size*size(lgt_block,2), &
        MPI_INTEGER1, MPI_SUM, WABBIT_COMM, ierr)
    endif

    ! we need to delete the old lgt_block array to avoid any rotting corpses somewhere.
    ! it is a little tricky to see why this is the case, but we found it to be necessary.
    ! basically, if on one mpirank, the last elements get removed, then the synchronized part
    ! gets smaller and the last elements are not overwritten with the copy statement below
    if (refinement_status_only) then
        !> \todo We should avoid resetting the entire array, as it can be expensive.
        !> maybe one can figure out using the input array what must be deleted
        !lgt_block(:,size(lgt_block,2)) = 0
        ! not necessary, only refinement status has changed
    else
        !> \todo We should avoid resetting the entire array, as it can be expensive.
        !lgt_block(:,:) = -1
        ! reset only first column
        lgt_block(:,1) = -1
    endif


    ! unpack synchronized buffer into the light data array.
    ! loop over number of procs and reset lgt_block array
    do k = 1, mpisize
        ! proc k-1 has sent data
        if ( proc_lgt_num(k) /= 0 ) then
            ! write received light data, recover original int precision (kind=ik)
            if (refinement_status_only) then
                lgt_block( (k-1)*N+1 : (k-1)*N + proc_lgt_num(k), R ) =  &
                int( my_lgt_block_recv_buffer( sum(proc_lgt_num(1:k-1))+1 : sum(proc_lgt_num(1:k-1))+proc_lgt_num(k), R ), kind=ik)
            else
                lgt_block( (k-1)*N+1 : (k-1)*N + proc_lgt_num(k), : ) =  &
                int( my_lgt_block_recv_buffer( sum(proc_lgt_num(1:k-1))+1 : sum(proc_lgt_num(1:k-1))+proc_lgt_num(k), : ), kind=ik)
            endif
        end if
    end do
end subroutine
