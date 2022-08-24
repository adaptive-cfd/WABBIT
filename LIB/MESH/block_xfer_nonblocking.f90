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
subroutine block_xfer( params, xfer_list, N_xfers, hvy_block, msg )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> list of transfers (:,1)=sender (:,2)=recv (:,3)=lgt_id
    integer(kind=ik), intent(inout)     :: xfer_list(:,:)
    integer(kind=ik), intent(in)        :: N_xfers
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    character(len=*), intent(in), optional :: msg

    integer(kind=ik) :: k, lgt_id, mpirank_recver, mpirank_sender, myrank, i, Nxfer_done, Nxfer_total, Nxfer_notPossibleNow
    integer(kind=ik) :: lgt_id_new, hvy_id_new, hvy_id, npoints, npoints2
    integer :: ierr, tag
    logical :: xfer_started(1:N_xfers)
    logical :: source_block_deleted(1:N_xfers)

    ! array of mpi requests, taken from stack. for extremely large N_xfers, that can cause stack problems
    integer(kind=ik) :: requests(1:2*N_xfers)
    integer(kind=ik) :: ireq

    logical :: not_enough_memory
    integer(kind=ik) :: counter
    real(kind=rk) :: t0
    character(len=40) :: msg2

    ! if the list of xfers is empty, then we just return.
    if (N_xfers==0) return

    t0 = MPI_wtime()
    myrank = params%rank
    counter = 0
    xfer_started = .false.
    source_block_deleted = .false.
    ! size of one block, in points
    npoints = size(hvy_block,1)*size(hvy_block,2)*size(hvy_block,3)*size(hvy_block,4)

    Nxfer_done = 0
    Nxfer_total = N_xfers

    msg2 = ""
    if (present(msg)) msg2=msg

! the routine can resume from the next line on (like in the olden days, before screens were invented)
1   ireq = 0
    not_enough_memory = .false.

    Nxfer_notPossibleNow = 0

    do k = 1, N_xfers
        ! if this xfer has already been started, we are happy and go to the next one.
        if (xfer_started(k)) cycle

        ! we will transfer this block:
        lgt_id = xfer_list(k,3)

        ! from its current owner: "mpirank_sender"
        mpirank_sender = xfer_list(k,1)

        ! to its new owner:
        mpirank_recver = xfer_list(k,2)

        ! its new light id will be "lgt_id_new"
        call get_free_local_light_id( params, mpirank_recver, lgt_id_new, ignore_error=.true. )

        ! the idea is now if we do not have enough memory (ie no free block on target rank) we can
        ! wait until the current requests are finnished. Then we retry the loop.
        if (lgt_id_new == -1) then
            not_enough_memory = .true.
            ! this marks that we can NOT delete the source block after waiting for MPI xfer
            xfer_started(k) = .false.
            Nxfer_notPossibleNow = Nxfer_notPossibleNow + 1
            ! skip this xfer (it will be treated in the next iteration)
            cycle
        else
            ! this marks that we can delete the source block after waiting for MPI xfer
            xfer_started(k) = .true.
            Nxfer_done = Nxfer_done + 1
        endif



        ! Am I the target rank who receives this block of data?
        if (myrank == mpirank_recver) then
            !-------------------------------------------------------------------
            ! RECV CASE
            !-------------------------------------------------------------------
            ! get hvy id where to store the data
            call lgt2hvy( hvy_id_new, lgt_id_new, myrank, params%number_blocks )

            ireq = ireq + 1
            tag = 2*k ! use unique tag for each message

            ! open channel to receive one block.
            call MPI_irecv( hvy_block(:,:,:,:,hvy_id_new), npoints, MPI_DOUBLE_PRECISION, mpirank_sender, &
            tag, WABBIT_COMM, requests(ireq), ierr)

            if (ierr /= MPI_SUCCESS) call abort(1809181531, "[block_xfer.f90] "//trim(adjustl(msg2))//" MPI_irecv failed!")

        ! Am I the owner of this block, so will I have to send data?
        elseif (myrank == mpirank_sender) then
            !-------------------------------------------------------------------
            ! SEND CASE
            !-------------------------------------------------------------------
            ! what heavy ID (on its owner proc, which is me) does the block have?
            call lgt2hvy( hvy_id, lgt_id, mpirank_sender, params%number_blocks )

            ireq = ireq + 1
            tag = 2*k ! use unique tag for each message

            ! send the block to the receiver. Note unfortunately we cannot delete it right away, since
            ! we have to wait for the MPI_REQUEST to be finnished.
            call MPI_isend( hvy_block(:,:,:,:,hvy_id), npoints, MPI_DOUBLE_PRECISION, mpirank_recver, tag, &
            WABBIT_COMM, requests(ireq), ierr)

            if (ierr /= MPI_SUCCESS) call abort(1809181532, "[block_xfer.f90] "//trim(adjustl(msg2))//" MPI_isend failed!")

        endif

        ! AVOID race condition:
        ! all cpu copy the light ID: this way, we ensure that our newly xfered block
        ! is not deleted again by other send/recv ops (it is hence excluded that lgt_id_new
        ! in the next loop is again the same value as it is right now)
        lgt_block( lgt_id_new, : ) = lgt_block( lgt_id, : )
        ! note you mut NOT delete the original block now. only after xfer is completed.
    enddo

    ! wait for all my xfers to be completed. (ireq is different on each mpirank)
    ! both send/recv are waited for
    if (ireq > 0) then
        call MPI_waitall( ireq, requests(1:ireq), MPI_STATUSES_IGNORE, ierr )
        if (ierr /= MPI_SUCCESS) call abort(1809181533, "[block_xfer.f90] "//trim(adjustl(msg2))//" MPI_waitall failed!")
    endif


    ! now all data is transfered on this rank, and we will no longer change it.
    ! so now, we can safely delete the original blocks, if their xfer was started.
    ! Note that the waiting is different for each rank, and they do arrive here at different times,
    ! but they all eventually delete the source blocks only after they have been transfered.
    do i = 1, N_xfers
        ! delete block. it has been moved previously, and now we delete the original
        if (xfer_started(i) .and. .not. source_block_deleted(i)) then
            lgt_block( xfer_list(i,3), : ) = -1
            ! pay attention here: you must not delete the block twice. after the first delete, it is
            ! free and can be used for a newly xferd block: then it is no longer free and must not be
            ! deleted again
            source_block_deleted(i) = .true.
        endif
    enddo

    ! it may happen that we cannot execute a xfer because the target ranks has no more memory left
    ! but it s possible that this is only due to non-completed xfers! so now, if not_enough_memory=.true.
    ! we waited once for all xfers and the temporary blocks are freed. Now we can try again! We count
    ! to avoid infinite loops - in this case, we simply *really* ran out of memory.
    if (not_enough_memory) then
        if (counter < 200000) then
            counter = counter + 1
            if (myrank==0) then
                ! t0 is just an identifier for the call
                call append_t_file( "block_xfer.t", (/t0, real(counter,kind=rk), real(Nxfer_total,kind=rk), &
                real(Nxfer_done,kind=rk), real(Nxfer_notPossibleNow,kind=rk) /) )
            endif
            ! like in the golden age of computer programming:
            goto 1
        else
            call close_all_t_files()
            call abort(1909181808, "[block_xfer.f90]: "//trim(adjustl(msg2))//" not enough memory to complete transfer.")
        endif
    endif

end subroutine block_xfer
