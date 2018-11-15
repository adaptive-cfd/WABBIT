! block transfer between MPIranks
! The idea is that you provide a list of transfers (american english: xfer) which will be
! executed using non-blocking MPI operations here.
!
! Be sure to gather all possible xfers before calling this routine - do not call
! it several times for each block.
!
! All lists are outdated at the end of the routine. - call create_active_and_sorted_lists
!
! All xfers are completed when the code leaves this routine.
! Q: Wouldn't it be nicer to wait elsewhere, in order to do some stuff while we wait?
! A: well, the question is how much you can do and where you wait. I don't think the
! performance is better. note you have to wait somewhere! always!
!
! NOTE: We expect the xfer_list to be identical on all ranks
subroutine block_xfer( params, xfer_list, N_xfers, lgt_block, hvy_block, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_tmp )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> list of transfers (:,1)=sender (:,2)=recv (:,3)=lgt_id
    integer(kind=ik), intent(inout)     :: xfer_list(:,:)
    integer(kind=ik), intent(in)        :: N_xfers
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)

    integer(kind=ik) :: k, lgt_id, mpirank_recver, mpirank_sender, myrank
    integer(kind=ik) :: lgt_id_new, hvy_id_new, hvy_id, npoints, ierr, tag

    ! array of mpi requests, taken from stack. for extremely large N_xfers, that can cause stack problems
    integer(kind=ik) :: requests(1:N_xfers)
    integer(kind=ik) :: ireq

    logical :: not_enough_memory
    integer(kind=ik) :: counter, k_start

    ! if the list of xfers is empty, then we just return.
    if (N_xfers==0) return

    myrank = params%rank
    counter = 0
    k_start = 1
    ! size of one block, in points
    npoints = size(hvy_block,1)*size(hvy_block,2)*size(hvy_block,3)*size(hvy_block,4)

! the routine can resume from the next line on (like in the olden days, before screens were invented)
1   ireq = 0
    not_enough_memory = .false.

    do k = k_start, N_xfers
        ! we will transfer this block:
        lgt_id = xfer_list(k,3)

        ! from its current owner: "mpirank_sender"
        mpirank_sender = xfer_list(k,1)

        ! to its new owner:
        mpirank_recver = xfer_list(k,2)

        ! its new light id will be "lgt_id_new"
        call get_free_local_light_id( params, mpirank_recver, lgt_block, lgt_id_new, ignore_error=.true. )

        ! the idea is now if w do not have enough memory (ie no free block on target rank) we can
        ! wait until the current requests are finnished. Then we retry the loop.
        if (lgt_id_new == -1) then
            k_start = k
            not_enough_memory = .true.
            counter = counter + 1
            ! interrupt do loop
            exit
        endif

        ! use  unique tag for each message
        tag = lgt_id

        ! Am I the target rank who receives this block of data?
        if (myrank == mpirank_recver) then
            !-------------------------------------------------------------------
            ! RECV CASE
            !-------------------------------------------------------------------
            ! get hvy id where to store the data
            call lgt_id_to_hvy_id( hvy_id_new, lgt_id_new, myrank, params%number_blocks )

            ireq = ireq + 1

            ! open channel to receive one block.
            call MPI_irecv( hvy_block(:,:,:,:,hvy_id_new), npoints, MPI_DOUBLE_PRECISION, mpirank_sender, &
            tag, WABBIT_COMM, requests(ireq), ierr)

            if (ierr /= MPI_SUCCESS) call abort(1809181531, "[block_xfer.f90] MPI_irecv failed!")

        ! Am I the owner of this block, so will I have to send data?
    elseif (myrank == mpirank_sender) then
            !-------------------------------------------------------------------
            ! SEND CASE
            !-------------------------------------------------------------------
            ! what heavy ID (on its owner proc, which is me) does the block have?
            call lgt_id_to_hvy_id( hvy_id, lgt_id, mpirank_sender, params%number_blocks )

            ireq = ireq + 1

            ! send the block to the receiver. Note unfortunately we cannot delete it right away, since
            ! we have to wait for the MPI_REQUEST to be finnished.
            call MPI_isend( hvy_block(:,:,:,:,hvy_id), npoints, MPI_DOUBLE_PRECISION, mpirank_recver, tag, &
            WABBIT_COMM, requests(ireq), ierr)

            if (ierr /= MPI_SUCCESS) call abort(1809181532, "[block_xfer.f90] MPI_isend failed!")
        endif

        ! AVOID race condition:
        ! all cpu copy the light ID: this way, we ensure that our newly xfered block
        ! is not deleted again by other send/recv ops (it is hence excluded that lgt_id_new
        ! in the next loop is again the same value as it is right now)
        lgt_block( lgt_id_new, : ) = lgt_block( lgt_id, : )
        ! note you mut NOT delete the original block now. only after xfer is completed.
    enddo

    ! wait for all my xfers to be completed.
    if (ireq > 0) then
        call MPI_waitall( ireq, requests(1:ireq), MPI_STATUSES_IGNORE, ierr )
        if (ierr /= MPI_SUCCESS) call abort(1809181533, "[block_xfer.f90] MPI_waitall failed!")
    endif

    ! now all data is transfered on this rank, and we will no longer change it.
    ! so now, we can safely delete the original blocks
    do k = 1, N_xfers
        ! delete block. it has been moved previously, and now we delete the original
        lgt_block( xfer_list(k,3), 1 ) = -1
    enddo

    ! it may happen that we cannot execute a xfer because the target ranks has no more memory left
    ! but it s possible that this is only due to non-completed xfers! so now, if not_enough_memory=.true.
    ! we waited once for all xfers and the temporary blocks are freed. Now we can try again! We count
    ! to avoid infinite loops - in this case, we simply *really* ran out of memory.
    if (not_enough_memory) then
        if (counter < 10) then
            ! like in the golden age of computer programming:
            goto 1
        else
            call abort(1909181808, "[block_xfer.f90]: we tried but we simply do not have enough memory to complete transfer. &
            & you need to allocate more --memory")
        endif
    endif

end subroutine block_xfer
