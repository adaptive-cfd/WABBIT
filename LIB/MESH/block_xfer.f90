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

subroutine block_xfer( params, xfer_list, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !
    integer(kind=ik), intent(inout) :: xfer_list(:,:)
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

    integer(kind=ik) :: N_xfers, k, lgt_id, mpirank_newowner, mpirank_owner, myrank
    integer(kind=ik) :: lgt_id_new, hvy_id_new, hvy_id, npoints, ierr, tag

    myrank = params%rank
    ! size of one block, in points
    npoints = size(hvy_block,1)*size(hvy_block,2)*size(hvy_block,3)*size(hvy_block,4)

    N_xfers = size(xfer_list,1)

    xfer_list(1:N_xfers,4) = MPI_REQUEST_NULL

    do k = 1, N_xfers
        ! we will transfer this block:
        lgt_id = xfer_list(k,1)

        ! from its current owner: "mpirank_owner"
        call lgt_id_to_proc_rank( mpirank_owner, lgt_id, params%number_blocks )

        ! to its new owner:
        mpirank_newowner = xfer_list(k,2)

        ! its new light id will be "lgt_id_new"
        call get_free_local_light_id( params, mpirank_newowner, lgt_block, lgt_id_new )
        ! we also return the new light id in the xfer list to the caller. that is helpful
        ! for block merging and avoids looking again for their sisters.
        xfer_list(k, 3) = lgt_id_new

        tag = lgt_id

        ! Am I the target rank who receives all the data?
        if (myrank == mpirank_newowner) then
            !------------------------
            ! RECV CASE
            !------------------------
            ! get hvy id where to store the data
            call lgt_id_to_hvy_id( hvy_id_new, lgt_id_new, myrank, params%number_blocks )

            ! open channel to receive one block
            call MPI_irecv( hvy_block(:,:,:,:,hvy_id_new), npoints, MPI_REAL8, mpirank_owner, &
            tag, WABBIT_COMM, xfer_list(k, 4), ierr)

        elseif (myrank == mpirank_owner) then
            ! Am I the owner of this block, so will I have to send data?
            !------------------------
            ! SEND CASE
            !------------------------
            ! what heavy ID (on this proc) does the block have?
            call lgt_id_to_hvy_id( hvy_id, lgt_id, mpirank_owner, params%number_blocks )

            ! send the block to the receiver
            call MPI_isend( hvy_block(:,:,:,:,hvy_id), npoints, MPI_REAL8, mpirank_newowner, tag, &
            WABBIT_COMM,  xfer_list(k, 4), ierr)
        endif

        ! all cpu copy the light ID: this way, we ensure that our newly xfered block
        ! is not deleted again by other send/recv ops (it is hence excluded that lgt_id_new
        ! in the next loop is again the same value as it is right now)
        lgt_block( lgt_id_new, : ) = lgt_block( lgt_id, : )
        ! note you mut NOT delete the original block now. only after xfer is completed.
    enddo

    ! wait for all my xfers to be completed.
    call MPI_waitall( N_xfers, xfer_list(:,4), MPI_STATUSES_IGNORE, ierr )

    ! now all data is transfered on this rank, and we will no longer change it.
    ! so now, we can safely delete the original blocks
    do k = 1, N_xfers
        ! delete block. it has been moved previously, and now we delete the original
        lgt_block( xfer_list(k,1), 1 ) = -1
    enddo

end subroutine block_xfer
