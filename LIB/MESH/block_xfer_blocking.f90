! block transfer between MPIranks
! The idea is that you provide a list of transfers (american english: xfer) which will be
! executed using blocking MPI operations here.
!
! Be sure to gather all possible xfers before calling this routine - do not call
! it several times for each block.
!
! Q: There is a nonblocking version of this routine. Shouldn't it be better?
! A: YES! but it does not work on the BLueGene/Q. No idea why.
!
! NOTE: We expect the xfer_list to be identical on all ranks
subroutine block_xfer( params, xfer_list, N_xfers, lgt_block, hvy_block, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_work )
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
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)

    integer(kind=ik) :: k, l, lgt_id, neq, heavy_id, com_N, data_size, &
    hvy_free_id, ierr, lgt_free_id, myrank, tag, mpirank_recver, mpirank_sender
    integer(kind=ik), allocatable, save :: buffer_light( : )
    integer :: status(MPI_status_size)

    ! if the list of xfers is empty, then we just return.
    if (N_xfers==0) return

    if (.not.allocated(buffer_light)) allocate( buffer_light( params%number_blocks ) )

    myrank = params%rank
    neq = params%n_eqn

    ! size of data array, use for readability
    data_size = size(hvy_block,1) * size(hvy_block,2) * size(hvy_block,3) * size(hvy_block,4)


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
    do k = 1, N_xfers
        tag = 0
        ! check if the transfer is already taken care of: in this case, it is deleted by a -1

        mpirank_sender = xfer_list(k, 1)
        mpirank_recver = xfer_list(k, 2)

        if ( xfer_list(k, 1) /= -1 ) then

            !-----------------------------------------------------------
            ! I am the sender in this operation?
            !-----------------------------------------------------------
            if ( mpirank_sender == myrank ) then

                ! create send buffer, search list, send next l blocks to the same receiver
                ! (if they belong there)
                l = 0
                do while ( (xfer_list(k+l, 1) == mpirank_sender) .and. (xfer_list(k+l, 2) == mpirank_recver) )
                    lgt_id = xfer_list(k+l, 3)
                    ! calculate heavy id from light id
                    call lgt_id_to_hvy_id( heavy_id, lgt_id, myrank, params%number_blocks )

                    ! as the blocks which are consecutive in the list are not necessarily
                    ! consecutive in the hvy_block array, we first collect them in the bufffer
                    ! array. Then they are contiguous.
                    hvy_work(:, :, :, 1:neq, l+1 ) = hvy_block(:, :, :, 1:neq, heavy_id )
                    ! ... light data
                    buffer_light( l+1 ) = lgt_id

                    ! ... then delete the block I just got rid of
                    lgt_block(lgt_id, : ) = -1

                    ! go to next element
                    l = l + 1
                end do

                ! send data
                call MPI_Send( buffer_light(1:l), l, MPI_INTEGER4, mpirank_recver, tag, WABBIT_COMM, ierr)
                call MPI_Send( hvy_work(:,:,:,1:neq,1:l), data_size*l, MPI_REAL8, mpirank_recver, tag, WABBIT_COMM, ierr)

                ! delete all com list elements
                xfer_list(k:k+l-1, :) = -1

            !-----------------------------------------------------------
            ! I am the receiver in this operation?
            !-----------------------------------------------------------
            elseif ( mpirank_recver == myrank ) then

                ! count received data sets, recv next l blocks from the sender
                l = 1
                do while ( (xfer_list(k+l, 1) == mpirank_sender) .and. (xfer_list(k+l, 2) == mpirank_recver) )
                    ! delete element
                    xfer_list(k+l, :) = -1
                    l = l + 1
                end do

                ! receive data
                call MPI_Recv( buffer_light(1:l), l, MPI_INTEGER4, mpirank_sender, tag, WABBIT_COMM, status, ierr)
                call MPI_Recv( hvy_work(:,:,:,1:neq,1:l), data_size*l, MPI_REAL8, mpirank_sender, tag, WABBIT_COMM, status, ierr)

                ! delete first com list element after receiving data
                xfer_list(k, :) = -1

                ! save comm count
                com_N = l

                ! loop over all received blocks
                do l = 1,  com_N
                    ! fetch a free light id slot on my myrank
                    call get_free_local_light_id( params, myrank, lgt_block, lgt_free_id )
                    call lgt_id_to_hvy_id( hvy_free_id, lgt_free_id, myrank, params%number_blocks )

                    ! copy the data from the buffers
                    lgt_block( lgt_free_id, :) = lgt_block( buffer_light(l), : )
                    hvy_block(:, :, :, 1:neq, hvy_free_id) = hvy_work(:, :, :, 1:neq, l)

                end do
            else
                ! I am not concerned by this operation, delete element
                xfer_list(k, :) = -1
            end if
        end if
    end do


    call synchronize_lgt_data( params, lgt_block, refinement_status_only=.false. )

end subroutine block_xfer
