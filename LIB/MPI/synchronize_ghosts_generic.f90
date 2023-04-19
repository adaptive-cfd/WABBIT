subroutine sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, &
    syncSameLevelOnly1, lgt_BlocksToSync  )

    implicit none

    type (type_params), intent(in) :: params
    !> light data array
    integer(kind=ik), intent(in)   :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)   :: hvy_neighbor(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)   :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)   :: hvy_n
    logical, optional, intent(in)  :: syncSameLevelOnly1
    logical, optional, intent(in)  :: lgt_BlocksToSync(:)

    integer(kind=ik)   :: myrank, mpisize, ii0, ii1, Bs(1:3)
    integer(kind=ik)   :: N, k, neighborhood, level_diff, Nstages
    integer(kind=ik)   :: recver_lgtID, recver_rank, recver_hvyID
    integer(kind=ik)   :: sender_hvyID, sender_lgtID
    logical :: syncSameLevelOnly

    integer(kind=ik) :: ijk(2,3), isend, irecv
    integer(kind=ik) :: bounds_type, istage, inverse
    real(kind=rk) :: t0

    t0 = MPI_wtime()
    if (present(syncSameLevelOnly1)) then
        syncSameLevelOnly = syncSameLevelOnly1
    else
        syncSameLevelOnly = .false.
    endif

    if (.not. ghost_nodes_module_ready) then
        ! in order to keep the syntax clean, buffers are module-global and need to be
        ! allocated here.
        call init_ghost_nodes( params )
    endif

    ! if this mpirank has no active blocks, it has nothing to do here.
    if (hvy_n == 0) return

!
! New Idea (09 Apr 2023)
!
! During adapt_mesh, which is the most performance-hungry part of the algorithm with biorthogonal wavelets
! we should sync certain levels only. This we need also for the complete wavelet transform and denoising/CVS.
! Synching a level should mean we syn (J,J-1). We start from Jmax. This means: sync all neighborhoods that
! involve level J, also on blocks of level J-1. It does not mean a block on (J-1) is completely synced, just those
! neighborhoods.
!
! 2nd IDEA: stage-free ghost nodes.
! It turn out, if the coarser block sends not-correctly interpolatedy points, they can be corrected on the
! receiving fine block. Only a few are affected; the ones toward the interface. The ones further away
! are directly interpolated correctly. 19 apr 2023: This can be made work, but it is tedious and does not yield an
! immense speedup neither. On my local machine, its ~5%, but on large scale parallel sims, it may be more significant.
! Idea is described in inskape notes.

    Bs    = params%Bs
    N     = params%number_blocks
    myrank  = params%rank
    mpisize = params%number_procs
    ! default is two stages: first, copy&decimate, then interpolate
    Nstages = 2
    ! but if we sync only sameLevel neighbors, then we don't need the second stage
    if (syncSameLevelOnly) Nstages = 1

    ! call reset_ghost_nodes(  params, hvy_block, hvy_active, hvy_n )

    ! We require two stages: first, we fill all ghost nodes which are simple copy (including restriction),
    ! then in the second stage we can use interpolation and fill the remaining ones.
    do istage = 1, Nstages
        !***************************************************************************
        ! (i) stage initialization
        !***************************************************************************
        int_pos(:) = 1
        real_pos(:) = 0
        internalNeighbor_pos = 1

        ! compute, locally from the grid info, how much data I recv from and send to all
        ! other mpiranks. this gives us the start indices of each rank in the send/recv buffers. Note
        ! we do not count our internal nodes (.false. as last argument), as they are not put in the
        ! buffer at any time.
        call get_my_sendrecv_amount_with_ranks(params, lgt_block, hvy_neighbor, hvy_active, hvy_n, &
             Data_recvCounter, Data_sendCounter, MetaData_recvCounter, MetaData_sendCounter, istage, &
             count_internal=.false., ncomponents=size(hvy_block,4))


        ! reset iMetaData_sendBuffer, but only the parts that will actually be treated.
        do k = 1, params%number_procs
            ii0 = sum(MetaData_sendCounter(0:(k-1)-1)) + 1
            ii1 = ii0 + MetaData_sendCounter(k-1)
            iMetaData_sendBuffer(ii0:ii1) = -99

            ii0 = sum(MetaData_recvCounter(0:(k-1)-1)) + 1
            ii1 = ii0 + MetaData_recvCounter(k-1)
            iMetaData_recvBuffer(ii0:ii1) = -99
        enddo


        !***************************************************************************
        ! (ii) prepare data for sending
        !***************************************************************************
        do k = 1, hvy_n
            ! calculate light id
            sender_hvyID = hvy_active(k)
            call hvy2lgt( sender_lgtID, sender_hvyID, myrank, N )

            ! loop over all neighbors
            do neighborhood = 1, size(hvy_neighbor, 2)
                ! neighbor exists
                if ( hvy_neighbor( sender_hvyID, neighborhood ) /= -1 ) then
                    ! neighbor light data id
                    recver_lgtID = hvy_neighbor( sender_hvyID, neighborhood )
                    ! calculate neighbor rank
                    call lgt2proc( recver_rank, recver_lgtID, N )
                    ! neighbor heavy id
                    call lgt2hvy( recver_hvyID, recver_lgtID, recver_rank, N )
                    ! define level difference: sender - receiver, so +1 means sender on higher level
                    level_diff = lgt_block( sender_lgtID, params%Jmax + IDX_MESH_LVL ) - lgt_block( recver_lgtID, params%Jmax + IDX_MESH_LVL )

                    ! leveldiff = -1 : sender coarser than recver, interpolation on sender side
                    ! leveldiff =  0 : sender is same level as recver
                    ! leveldiff = +1 : sender is finer than recver, restriction is applied on sender side
                    if ( istage == 1 ) then
                        ! This is preparation for sending: hence in phase 1, level_diff=-1 skips and in phase 2 [0,+1 skip]
                        if (level_diff == -1) cycle
                        if (syncSameLevelOnly .and. level_diff/=0) cycle
                    else
                        ! in stage two leveldiff +1 and 0 are already done
                        if (level_diff ==  0) cycle
                        if (level_diff == +1) cycle
                    endif

                    if ( myrank == recver_rank ) then
                        ! internal relation (no communication) has its own buffer (to avoid senseless copying
                        ! from send to recv buffer)
                        internalNeighborSyncs( internalNeighbor_pos:internalNeighbor_pos+4-1 ) = (/sender_hvyID, recver_hvyID, neighborhood, level_diff/)
                        internalNeighbor_pos = internalNeighbor_pos + 4
                    else
                        ! external relation (MPI communication)
                        call send_prepare_external( params, recver_rank, hvy_block, sender_hvyID, recver_hvyID, neighborhood, level_diff )
                    end if ! (myrank==recver_rank)
                end if ! neighbor exists
            end do ! loop over all possible  neighbors
        end do ! loop over all heavy active

        !***************************************************************************
        ! (iii) transfer part (send/recv)
        !***************************************************************************
        call start_xfer_mpi( params, iMetaData_sendBuffer, rData_sendBuffer, iMetaData_recvBuffer, rData_recvBuffer, isend, irecv )


        !***************************************************************************
        ! (iv) Unpack received data in the ghost node layers
        !***************************************************************************
        ! process-internal ghost points (direct copy)
        call unpack_ghostlayers_internal( params, hvy_block )

        ! before unpacking the data we received from other ranks, we wait for the transfer
        ! to be completed
        call finalize_xfer_mpi(params, isend, irecv)

        ! process-external ghost points (copy from buffer)
        call unpack_ghostlayers_external( params, hvy_block )

    end do ! loop over stages 1,2

    call toc( "WRAPPER: sync ghosts", MPI_wtime()-t0 )

end subroutine sync_ghosts

! subroutine sync_ghosts_nostages( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
!
!     implicit none
!
!     type (type_params), intent(in) :: params
!     !> light data array
!     integer(kind=ik), intent(in)   :: lgt_block(:, :)
!     !> heavy data array - block data
!     real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)
!     !> heavy data array - neighbor data
!     integer(kind=ik), intent(in)   :: hvy_neighbor(:,:)
!     !> list of active blocks (heavy data)
!     integer(kind=ik), intent(in)   :: hvy_active(:)
!     !> number of active blocks (heavy data)
!     integer(kind=ik), intent(in)   :: hvy_n
!
!     integer(kind=ik)   :: myrank, mpisize, ii0, ii1, Bs(1:3)
!     integer(kind=ik)   :: N, k, neighborhood, level_diff
!     integer(kind=ik)   :: recver_lgtID, recver_rank, recver_hvyID
!     integer(kind=ik)   :: sender_hvyID, sender_lgtID
!
!     integer(kind=ik) :: ijk(2,3), isend, irecv
!     integer(kind=ik) :: bounds_type, inverse
!     real(kind=rk) :: t0
!
!     t0 = MPI_wtime()
!
!     if (.not. ghost_nodes_module_ready) then
!         ! in order to keep the syntax clean, buffers are module-global and need to be
!         ! allocated here.
!         call init_ghost_nodes( params )
!     endif
!
!     ! if this mpirank has no active blocks, it has nothing to do here.
!     if (hvy_n == 0) return
!
! !
! ! New Idea (09 Apr 2023)
! !
! ! During adapt_mesh, which is the most performance-hungry part of the algorithm with biorthogonal wavelets
! ! we should sync certain levels only. This we need also for the complete wavelet transform and denoising/CVS.
! ! Synching a level should mean we syn (J,J-1). We start from Jmax. This means: sync all neighborhoods that
! ! involve level J, also on blocks of level J-1. It does not mean a block on (J-1) is completely synced, just those
! ! neighborhoods.
! !
! ! 2nd IDEA: stage-free ghost nodes.
! ! It turn out, if the coarser block sends not-correctly interpolatedy points, they can be corrected on the
! ! receiving fine block. Only a few are affected; the ones toward the interface. The ones further away
! ! are directly interpolated correctly.
!
!     Bs    = params%Bs
!     N     = params%number_blocks
!     myrank  = params%rank
!     mpisize = params%number_procs
!
!     ! call reset_ghost_nodes(  params, hvy_block, hvy_active, hvy_n )
!
!         !***************************************************************************
!         ! (i) stage initialization
!         !***************************************************************************
!         int_pos(:) = 1
!         real_pos(:) = 0
!         internalNeighbor_pos = 1
!
!         ! compute, locally from the grid info, how much data I recv from and send to all
!         ! other mpiranks. this gives us the start indices of each rank in the send/recv buffers. Note
!         ! we do not count our internal nodes (.false. as last argument), as they are not put in the
!         ! buffer at any time.
!         call get_my_sendrecv_amount_with_ranks_nostages(params, lgt_block, hvy_neighbor, hvy_active, hvy_n, &
!              Data_recvCounter, Data_sendCounter, MetaData_recvCounter, MetaData_sendCounter, &
!              count_internal=.false., ncomponents=size(hvy_block,4))
!
!
!         ! reset iMetaData_sendBuffer, but only the parts that will actually be treated.
!         do k = 1, params%number_procs
!             ii0 = sum(MetaData_sendCounter(0:(k-1)-1)) + 1
!             ii1 = ii0 + MetaData_sendCounter(k-1)
!             iMetaData_sendBuffer(ii0:ii1) = -99
!
!             ii0 = sum(MetaData_recvCounter(0:(k-1)-1)) + 1
!             ii1 = ii0 + MetaData_recvCounter(k-1)
!             iMetaData_recvBuffer(ii0:ii1) = -99
!         enddo
!
!
!         !***************************************************************************
!         ! (ii) prepare data for sending
!         !***************************************************************************
!         do k = 1, hvy_n
!             ! calculate light id
!             sender_hvyID = hvy_active(k)
!             call hvy2lgt( sender_lgtID, sender_hvyID, myrank, N )
!
!             ! loop over all neighbors
!             do neighborhood = 1, size(hvy_neighbor, 2)
!                 ! neighbor exists
!                 if ( hvy_neighbor( sender_hvyID, neighborhood ) /= -1 ) then
!                     ! neighbor light data id
!                     recver_lgtID = hvy_neighbor( sender_hvyID, neighborhood )
!                     ! calculate neighbor rank
!                     call lgt2proc( recver_rank, recver_lgtID, N )
!                     ! neighbor heavy id
!                     call lgt2hvy( recver_hvyID, recver_lgtID, recver_rank, N )
!                     ! define level difference: sender - receiver, so +1 means sender on higher level
!                     ! leveldiff = -1 : sender coarser than recver, interpolation on sender side
!                     ! leveldiff =  0 : sender is same level as recver
!                     ! leveldiff = +1 : sender is finer than recver, restriction is applied on sender side
!                     level_diff = lgt_block( sender_lgtID, params%Jmax + IDX_MESH_LVL ) - lgt_block( recver_lgtID, params%Jmax + IDX_MESH_LVL )
!
!                     if ( myrank == recver_rank ) then
!                         ! internal relation (no communication) has its own buffer (to avoid senseless copying
!                         ! from send to recv buffer)
!                         internalNeighborSyncs( internalNeighbor_pos:internalNeighbor_pos+4-1 ) = (/sender_hvyID, recver_hvyID, neighborhood, level_diff/)
!                         internalNeighbor_pos = internalNeighbor_pos + 4
!                     else
!                         ! external relation (MPI communication)
!                         call send_prepare_external( params, recver_rank, hvy_block, sender_hvyID, recver_hvyID, neighborhood, level_diff )
!                     end if ! (myrank==recver_rank)
!                 end if ! neighbor exists
!             end do ! loop over all possible  neighbors
!         end do ! loop over all heavy active
!
!         !***************************************************************************
!         ! (iii) transfer part (send/recv)
!         !***************************************************************************
!         call start_xfer_mpi( params, iMetaData_sendBuffer, rData_sendBuffer, iMetaData_recvBuffer, rData_recvBuffer, isend, irecv )
!
!
!         !***************************************************************************
!         ! (iv) Unpack received data in the ghost node layers
!         !***************************************************************************
!         ! process-internal ghost points (direct copy)
!         call unpack_ghostlayers_internal( params, hvy_block )
!
!         ! before unpacking the data we received from other ranks, we wait for the transfer
!         ! to be completed
!         call finalize_xfer_mpi(params, isend, irecv)
!
!         ! process-external ghost points (copy from buffer)
!         call unpack_ghostlayers_external( params, hvy_block )
!
!         ! call fixInterpolatedPoints_postSync( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)
!
!     call toc( "WRAPPER: sync ghosts", MPI_wtime()-t0 )
!
! end subroutine sync_ghosts_nostages



subroutine send_prepare_external( params, recver_rank, hvy_block, sender_hvyID, &
    recver_hvyID, neighborhood, level_diff )
    implicit none

    type (type_params), intent(in) :: params
    integer(kind=ik), intent(in)   :: recver_rank ! zero-based
    integer(kind=ik), intent(in)   :: sender_hvyID, recver_hvyID
    integer(kind=ik), intent(in)   :: neighborhood
    integer(kind=ik), intent(in)   :: level_diff
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)

    ! merged information of level diff and an indicator that we have a historic finer sender
    integer(kind=ik)   :: buffer_size
    integer(kind=ik)   :: ijk1(2,3), nc, s(1:4)

    nc = size(hvy_block,4)

    if (size(res_pre_data,4) < size(hvy_block,4)) then
        s(1) = size( res_pre_data, 1 )
        s(2) = size( res_pre_data, 2 )
        s(3) = size( res_pre_data, 3 )
        s(4) = size( hvy_block, 4 )
        deallocate( res_pre_data )
        allocate( res_pre_data(s(1), s(2), s(3), s(4)) )
    endif


    ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
    ! They depend on the neighbor-relation, level difference and the bounds type.
    ! The last index is 1-sender 2-receiver 3-restricted/predicted.
    if ( level_diff == 0 ) then
        ! simply copy the ghost node layer (no interpolation or restriction here) to a line buffer, which
        ! we will send to our neighbor mpirank
        ijk1 = ijkGhosts(:,:, neighborhood, level_diff, SENDER)

        call GhostLayer2Line( params, line_buffer, buffer_size, &
        hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), :, sender_hvyID) )

    else
        ! up/downsample data first, then flatten to 1D buffer
        ijk1 = ijkGhosts(:,:, neighborhood, level_diff, SENDER)

        call restrict_predict_data( params, res_pre_data, ijk1, neighborhood, level_diff, hvy_block, sender_hvyID )

        ijk1 = ijkGhosts(:,:, neighborhood, level_diff, RESPRE)

        call GhostLayer2Line( params, line_buffer, buffer_size, &
        res_pre_data( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), 1:nc) )
    end if

    ! the chunk of data is added to the MPI buffers (preparation for sending)
    call AppendLineToBuffer( iMetaData_sendBuffer, rData_sendBuffer, buffer_size, recver_rank, line_buffer, &
    recver_hvyID, neighborhood, level_diff )

end subroutine send_prepare_external


subroutine unpack_ghostlayers_external( params, hvy_block, lgt_BlocksToSync )
    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    logical, optional, intent(in)  :: lgt_BlocksToSync(:)

    integer(kind=ik) :: recver_rank ! zero-based
    integer(kind=ik) :: l, recver_hvyID, neighborhood
    integer(kind=ik) :: level_diff, bounds_type, buffer_position, buffer_size, rank_destination
    integer(kind=ik) :: ijk1(2,3), i0, nc

    nc = size(hvy_block,4)

    do recver_rank = 0, params%number_procs-1 ! zero-based
        ! skip my own rank, skip ranks that did not send me anything
        if (recver_rank /= params%rank .and. Data_recvCounter(recver_rank) /= 0) then
            ! start index of this mpirank in the int_buffer
            l = sum(MetaData_recvCounter(0:recver_rank-1)) + 1

            do while ( iMetaData_recvBuffer(l) > -99 )
                ! unpack the description of the next data chunk
                recver_hvyID     = iMetaData_recvBuffer(l)
                neighborhood     = iMetaData_recvBuffer(l+1)
                level_diff       = iMetaData_recvBuffer(l+2)
                buffer_position  = iMetaData_recvBuffer(l+3)
                buffer_size      = iMetaData_recvBuffer(l+4)
                rank_destination = iMetaData_recvBuffer(l+5)

#ifdef DEV
                if (rank_destination /= params%rank) then
                    write(*,*) "rank=", params%rank, "dest=", rank_destination, "co", MetaData_recvCounter, &
                    "l=", l, "recver_rank=", recver_rank
                    call abort(7373872, "EXT this data seems to be not mine!")
                endif
#endif

                ! copy data to line buffer. we now need to extract this to the ghost nodes layer (2D/3D)
                i0 = sum(Data_recvCounter(0:recver_rank-1)) + buffer_position
                line_buffer(1:buffer_size) = rData_recvBuffer( i0 : i0+buffer_size-1 )


                ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
                ! They depend on the neighbor-relation, level difference and the bounds type.
                ! The last index is 1-sender 2-receiver 3-restricted/predicted.
                call Line2GhostLayer( params, line_buffer, ijkGhosts(:,:, neighborhood, level_diff, RECVER), hvy_block, recver_hvyID )

                ! increase buffer position marker
                l = l + 6
            end do
        end if
    end do

end subroutine unpack_ghostlayers_external



subroutine unpack_ghostlayers_internal( params, hvy_block, lgt_BlocksToSync )
    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    logical, optional, intent(in)  :: lgt_BlocksToSync(:)

    integer(kind=ik) :: l, recver_hvyID, neighborhood
    integer(kind=ik) :: sender_hvyID, level_diff
    integer(kind=ik) :: ijk1(2,3), ijk2(2,3), nc

    nc = size(hvy_block,4)

    do l = 1, internalNeighbor_pos-1, 4 ! note this is one too far (because the pointer is ready for the next element) hence -1
        ! required info:  sender_hvyID, recver_hvyID, neighborhood, level_diff
        sender_hvyID = internalNeighborSyncs(l)
        recver_hvyID = internalNeighborSyncs(l+1)
        neighborhood = internalNeighborSyncs(l+2)
        level_diff   = internalNeighborSyncs(l+3)

        if ( level_diff == 0 ) then
            ! simply copy from sender block to receiver block (NOTE: both are on the same MPIRANK)
            ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
            ! They depend on the neighbor-relation and level difference
            ! The last index is 1-sender 2-receiver 3-restricted/predicted.
            ijk1 = ijkGhosts(:,:, neighborhood, level_diff, RECVER)
            ijk2 = ijkGhosts(:,:, neighborhood, level_diff, SENDER)

            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), 1:nc, recver_hvyID ) = &
            hvy_block( ijk2(1,1):ijk2(2,1), ijk2(1,2):ijk2(2,2), ijk2(1,3):ijk2(2,3), 1:nc, sender_hvyID)
        else
            ! interpolation or restriction before inserting
            call restrict_predict_data( params, res_pre_data, ijkGhosts(1:2,1:3, neighborhood, level_diff, SENDER), &
            neighborhood, level_diff, hvy_block, sender_hvyID )

            ! copy interpolated / restricted data to ghost nodes layer
            ! NOTE: the indices of ghost nodes data chunks are stored globally in the ijkGhosts array (see module_MPI).
            ! They depend on the neighbor-relation, level difference and the bounds type.
            ! The last index is 1-sender 2-receiver 3-restricted/predicted.

            ijk1 = ijkGhosts(:, :, neighborhood, level_diff, RECVER)
            ijk2 = ijkGhosts(:, :, neighborhood, level_diff, RESPRE)

            hvy_block( ijk1(1,1):ijk1(2,1), ijk1(1,2):ijk1(2,2), ijk1(1,3):ijk1(2,3), 1:nc, recver_hvyID ) = &
            res_pre_data( ijk2(1,1):ijk2(2,1), ijk2(1,2):ijk2(2,2), ijk2(1,3):ijk2(2,3), 1:nc)
        end if
    end do

end subroutine unpack_ghostlayers_internal


subroutine GhostLayer2Line( params, line_buffer, buffer_counter, hvy_data )
    implicit none

    type (type_params), intent(in)   :: params
    !> data buffer
    real(kind=rk), intent(inout)     :: line_buffer(:)
    ! buffer size
    integer(kind=ik), intent(out)    :: buffer_counter
    !> heavy block data, all data fields
    real(kind=rk), intent(inout)     :: hvy_data(:, :, :, :)

    ! loop variable
    integer(kind=ik) :: i, j, k, dF
    ! reset buffer size
    buffer_counter = 0

    ! loop over all data fields
    do dF = 1, size(hvy_data,4)
        do k = 1, size(hvy_data, 3) ! third dimension, note: for 2D cases k is always 1
            do j = 1, size(hvy_data, 2)
                do i = 1, size(hvy_data, 1)
                    ! increase buffer size
                    buffer_counter = buffer_counter + 1
                    ! write data buffer
                    line_buffer(buffer_counter)   = hvy_data( i, j, k, dF )
                end do
            end do
        end do
    end do

end subroutine GhostLayer2Line

!############################################################################################################

subroutine Line2GhostLayer( params, line_buffer, data_bounds, hvy_block, hvy_id )
    implicit none

    type (type_params), intent(in)  :: params
    !> data buffer
    real(kind=rk), intent(inout)    :: line_buffer(:)
    !> data_bounds
    integer(kind=ik), intent(inout) :: data_bounds(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)    :: hvy_block(:, :, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)    :: hvy_id

    ! loop variable
    integer(kind=ik) :: i, j, k, dF, buffer_i

    buffer_i = 1
    ! loop over all data fields
    do dF = 1, size(hvy_block,4)
        do k = data_bounds(1,3), data_bounds(2,3) ! third dimension, note: for 2D cases k is always 1
            do j = data_bounds(1,2), data_bounds(2,2)
                do i = data_bounds(1,1), data_bounds(2,1)
                    ! write data buffer
                    hvy_block( i, j, k, dF, hvy_id ) = line_buffer( buffer_i )
                    buffer_i = buffer_i + 1
                end do
            end do
        end do
    end do

end subroutine Line2GhostLayer




subroutine AppendLineToBuffer( iMetaData_sendBuffer, rData_sendBuffer, buffer_size, recver_rank, line_buffer, &
    hvy_id, neighborhood, level_diff )

    implicit none

    !> send buffers, integer and real
    integer(kind=ik), intent(inout) :: iMetaData_sendBuffer(:)
    real(kind=rk), intent(inout)    :: rData_sendBuffer(:)
    ! data buffer size
    integer(kind=ik), intent(in)    :: buffer_size
    ! id integer
    integer(kind=ik), intent(in)    :: recver_rank ! zero-based
    ! restricted/predicted data buffer
    real(kind=rk), intent(inout)    :: line_buffer(:)
    ! data buffer intergers, receiver heavy id, neighborhood id, level difference
    integer(kind=ik), intent(in)    :: hvy_id, neighborhood, level_diff

    integer(kind=ik)                :: buffer_position, i0, l0

    buffer_position = real_pos(recver_rank+1) + 1

    i0 = sum(Data_sendCounter(0:recver_rank-1)) + buffer_position
    l0 = sum(MetaData_sendCounter(0:recver_rank-1)) + 0 ! note as int_pos is 1-based, we here use 0-based offset

    ! real data
    if (buffer_size>0) then
        rData_sendBuffer( i0:i0+buffer_size-1 ) = line_buffer(1:buffer_size)
    endif

    ! save position of NEXT patch in real buffer
    real_pos(recver_rank+1) = real_pos(recver_rank+1) + buffer_size

    ! save: neighbor id, neighborhood, level difference, buffer size
    iMetaData_sendBuffer( l0+int_pos(recver_rank+1)   ) = hvy_id
    iMetaData_sendBuffer( l0+int_pos(recver_rank+1)+1 ) = neighborhood
    iMetaData_sendBuffer( l0+int_pos(recver_rank+1)+2 ) = level_diff
    iMetaData_sendBuffer( l0+int_pos(recver_rank+1)+3 ) = buffer_position
    iMetaData_sendBuffer( l0+int_pos(recver_rank+1)+4 ) = buffer_size
    iMetaData_sendBuffer( l0+int_pos(recver_rank+1)+5 ) = recver_rank ! zero-based ! FIVE

    ! mark end of buffer with -99, will be overwritten by next element if it is not the last one
    iMetaData_sendBuffer( l0+int_pos(recver_rank+1)+6 ) = -99

    int_pos(recver_rank+1) = int_pos(recver_rank+1) + 6 ! FIVE

end subroutine AppendLineToBuffer


!############################################################################################################

subroutine start_xfer_mpi( params, iMetaData_sendBuffer, rData_sendBuffer, iMetaData_recvBuffer, rData_recvBuffer, isend, irecv )

    implicit none

    type (type_params), intent(in)  :: params
    integer(kind=ik), intent(inout) :: iMetaData_sendBuffer(:)
    integer(kind=ik), intent(inout) :: iMetaData_recvBuffer(:)
    real(kind=rk), intent(inout)    :: rData_sendBuffer(:)
    real(kind=rk), intent(inout)    :: rData_recvBuffer(:)
    integer(kind=ik), intent(out)   :: isend, irecv

    integer(kind=ik) :: length_realBuffer, mpirank_partner
    integer(kind=ik) :: length_intBuffer, rank, ierr, tag
    integer(kind=ik) :: k, i0, l0

    rank = params%rank
    isend = 0
    irecv = 0

    ! these two arrays are module-global.
    recv_request = MPI_REQUEST_NULL
    send_request = MPI_REQUEST_NULL

    do mpirank_partner = 0, params%number_procs-1 ! zero based

        if ( MetaData_recvCounter(mpirank_partner) > 1 .and. mpirank_partner/=rank ) then
            ! length of integer buffer
            length_intBuffer = MetaData_recvCounter(mpirank_partner)
            i0 = sum(MetaData_recvCounter(0:mpirank_partner-1)) + 1 ! note exclude k of course do not run 0:mpirank_partner

            ! increase communication counter
            irecv = irecv + 1

            ! receive data
            tag = mpirank_partner
            call MPI_Irecv( iMetaData_recvBuffer(i0:i0+length_intBuffer-1), length_intBuffer, MPI_INTEGER4, &
            mpirank_partner, tag, WABBIT_COMM, recv_request(irecv), ierr)
        endif

        if ( MetaData_sendCounter(mpirank_partner) > 1 .and. mpirank_partner/=rank  ) then
            ! length of integer buffer
            length_intBuffer = MetaData_sendCounter(mpirank_partner)
            i0 = sum(MetaData_sendCounter(0:mpirank_partner-1)) + 1 ! note exclude k of course do not run 0:mpirank_partner

            ! increase communication counter
            isend = isend + 1

            ! send data
            tag = rank
            call MPI_Isend( iMetaData_sendBuffer(i0:i0+length_intBuffer-1), length_intBuffer, MPI_INTEGER4, &
            mpirank_partner, tag, WABBIT_COMM, send_request(isend), ierr)
        end if

    end do

    ! Note: internal ghost nodes do not pass by the buffer, they only require the integer
    ! buffer above. Therefore, no data (a message of zero length) will be sent to myself
    do mpirank_partner = 0, params%number_procs-1 ! zero based

        if (Data_sendCounter(mpirank_partner) > 0) then
            ! increase communication counter
            isend = isend + 1

            ! the amount of data is pre-computed in get_my_sendrecv_amount_with_ranks
            ! hence we do know how much data we will receive
            length_realBuffer = Data_sendCounter(mpirank_partner)

            i0 = sum(Data_sendCounter(0:mpirank_partner-1)) + 1 ! note exclude k of course do not run 0:mpirank_partner

            ! send data
            tag = rank
            call MPI_Isend( rData_sendBuffer(i0:i0+length_realBuffer-1), length_realBuffer, MPI_REAL8, &
            mpirank_partner, tag, WABBIT_COMM, send_request(isend), ierr)

        end if

        if (Data_recvCounter(mpirank_partner) > 0) then
            ! increase communication counter
            irecv = irecv + 1

            ! the amount of data is pre-computed in get_my_sendrecv_amount_with_ranks
            ! hence we do know how much data we will receive
            length_realBuffer = Data_recvCounter(mpirank_partner)

            i0 = sum(Data_recvCounter(0:mpirank_partner-1)) + 1 ! note exclude k of course do not run 0:mpirank_partner

            ! receive data
            tag = mpirank_partner
            call MPI_Irecv( rData_recvBuffer(i0:i0+length_realBuffer-1), length_realBuffer, MPI_REAL8, &
            mpirank_partner, MPI_ANY_TAG, WABBIT_COMM, recv_request(irecv), ierr)
        end if
    end do
end subroutine start_xfer_mpi


subroutine finalize_xfer_mpi(params, isend, irecv)
    implicit none
    type (type_params), intent(in) :: params
    integer(kind=ik), intent(in) :: isend, irecv
    integer(kind=ik) :: ierr

    ! synchronize non-blocking communications
    if (isend>0) then
        call MPI_Waitall( isend, send_request(1:isend), MPI_STATUSES_IGNORE, ierr)
    end if
    if (irecv>0) then
        call MPI_Waitall( irecv, recv_request(1:irecv), MPI_STATUSES_IGNORE, ierr)
    end if
end subroutine


! returns two lists with numbers of points I send to all other procs and how much I
! receive from each proc. note: strictly locally computed, NO MPI comm involved here
subroutine get_my_sendrecv_amount_with_ranks(params, lgt_block, hvy_neighbor, hvy_active,&
     hvy_n, Data_recvCounter, Data_sendCounter, MetaData_recvCounter, MetaData_sendCounter, istage, &
     count_internal, ncomponents)

    implicit none

    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n, ncomponents, istage
    integer(kind=ik), intent(inout)     :: Data_recvCounter(0:), Data_sendCounter(0:)
    integer(kind=ik), intent(inout)     :: MetaData_recvCounter(0:), MetaData_sendCounter(0:)
    logical, intent(in)                 :: count_internal

    integer(kind=ik) :: k, sender_hvyID, sender_lgtID, myrank, N, neighborhood, recver_rank
    integer(kind=ik) :: ijk(2,3), inverse, ierr, recver_hvyID, recver_lgtID,level_diff, status, new_size

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    N = params%number_blocks

    Data_recvCounter(:) = 0
    Data_sendCounter(:) = 0
    MetaData_recvCounter(:) = 0
    MetaData_sendCounter(:) = 0

    do k = 1, hvy_n
        ! calculate light id
        sender_hvyID = hvy_active(k)
        call hvy2lgt( sender_lgtID, sender_hvyID, myrank, N )

        ! loop over all neighbors
        do neighborhood = 1, size(hvy_neighbor, 2)
            ! neighbor exists
            if ( hvy_neighbor( sender_hvyID, neighborhood ) /= -1 ) then
                ! neighbor light data id
                recver_lgtID = hvy_neighbor( sender_hvyID, neighborhood )
                ! calculate neighbor rank
                call lgt2proc( recver_rank, recver_lgtID, N )
                ! neighbor heavy id
                call lgt2hvy( recver_hvyID, recver_lgtID, recver_rank, N )

                ! define level difference: sender - receiver, so +1 means sender on higher level
                ! leveldiff = -1 : sender coarser than recver, interpolation on sender side
                ! leveldiff =  0 : sender is same level as recver
                ! leveldiff = +1 : sender is finer than recver, restriction is applied on sender side
                level_diff = lgt_block( sender_lgtID, params%Jmax + IDX_MESH_LVL ) - lgt_block( recver_lgtID, params%Jmax + IDX_MESH_LVL )


                if (recver_rank /= myrank .or. count_internal) then
                    ! it now depends on the stage if we have to sent this data
                    ! or not.
                    ! In stage 1, only level_diff = {+1, 0} is treated
                    ! In stage 2, only level_diff = -1

                    ! send counter. how much data will I send to other mpiranks?
                    if ((istage==1 .and. (level_diff==+1 .or. level_diff==0)) .or. (istage==2.and.level_diff==-1)) then
                        ! why is this RECVER and not sender? Well, complicated. The amount of data on the sender patch
                        ! is not the same as in the receiver patch, because we interpolate or downsample. We effectively
                        ! transfer only the data the recver wants - not the extra data.
                        ijk = ijkGhosts(:, :, neighborhood, level_diff, RECVER)

                        Data_sendCounter(recver_rank) = Data_sendCounter(recver_rank) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)
                    endif

                    ! recv counter. how much data will I recv from other mpiranks?
                    ! This is NOT the same number as before
                    if ((istage==1 .and. (-1*level_diff==+1 .or. -1*level_diff==0)) .or. (istage==2.and.-1*level_diff==-1)) then

                        inverse = inverse_neighbor(neighborhood, dim)

                        ijk = ijkGhosts(:, :, inverse, -1*level_diff, RECVER)

                        Data_recvCounter(recver_rank) = Data_recvCounter(recver_rank) + &
                        (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)
                    endif

                    ! counter for integer buffer: for each neighborhood, we send 6 integers as metadata
                    ! as this is a fixed number it does not depend on the type of neighborhood etc, so
                    ! technically one would need only one for send/recv
                    if ((istage==1 .and. (level_diff==+1 .or. level_diff==0)) .or. (istage==2.and.level_diff==-1)) then
                        MetaData_sendCounter(recver_rank) = MetaData_sendCounter(recver_rank) + 6 ! FIVE
                    endif

                    if ((istage==1 .and. (-1*level_diff==+1 .or. -1*level_diff==0)) .or. (istage==2.and.-1*level_diff==-1)) then
                        MetaData_recvCounter(recver_rank) = MetaData_recvCounter(recver_rank) + 6 ! FIVE
                    endif
                endif



            end if ! neighbor exists
        end do ! loop over all possible  neighbors
    end do ! loop over all heavy active



    ! NOTE: for the int buffer, we mosly start at some index l0 and then loop unitl
    ! we find a -99 indicating the end of the buffer. this could be avoided by using
    ! for instead of while loops in the main routines, but I do not have time now.
    !
    ! In the meantime, notice we extent the amount of data by one, to copy the last -99
    ! to the buffers
    MetaData_recvCounter(:) = MetaData_recvCounter(:) + 1
    MetaData_sendCounter(:) = MetaData_sendCounter(:) + 1


    ! NOTE ACTUAL SEND / RECV DATA IS NEQN
    Data_recvCounter(:) = Data_recvCounter(:) * ncomponents
    Data_sendCounter(:) = Data_sendCounter(:) * ncomponents

    ! NOTE: this feature is against wabbits memory policy: we try to allocate the
    ! whole memory of the machine on startup, then work with that. however, we have to
    ! reserve portions of that memory for the state vector, the RHS slots, etc, and the ghost nodes
    ! buffer. However, estimating those latter is difficult: it depends on the grid and the parallelization
    if (sum(Data_recvCounter) > size(rData_recvBuffer, 1)) then
        ! out-of-memory case: the preallocated buffer is not large enough.
        write(*,'("rank=",i4," OOM for ghost nodes and increases its buffer size to 125%")') myrank
        new_size = size(rData_recvBuffer,1)*125/100
        deallocate(rData_recvBuffer)
        allocate( rData_recvBuffer(1:new_size), stat=status )
        if (status /= 0) call abort(999992, "Buffer allocation failed. Not enough memory?")
    endif

    if (sum(Data_sendCounter) > size(rData_sendBuffer, 1)) then
        ! out-of-memory case: the preallocated buffer is not large enough.
        write(*,'("rank=",i4," OOM for ghost nodes and increases its buffer size to 125%")') myrank
        new_size = size(rData_sendBuffer,1)*125/100
        deallocate(rData_sendBuffer)
        allocate( rData_sendBuffer(1:new_size), stat=status )
        if (status /= 0) call abort(999993, "Buffer allocation failed. Not enough memory?")
    endif
end subroutine get_my_sendrecv_amount_with_ranks


! ! returns two lists with numbers of points I send to all other procs and how much I
! ! receive from each proc. note: strictly locally computed, NO MPI comm involved here
! subroutine get_my_sendrecv_amount_with_ranks_nostages(params, lgt_block, hvy_neighbor, hvy_active,&
!      hvy_n, Data_recvCounter, Data_sendCounter, MetaData_recvCounter, MetaData_sendCounter, &
!      count_internal, ncomponents)
!
!     implicit none
!
!     type (type_params), intent(in)      :: params
!     !> light data array
!     integer(kind=ik), intent(in)        :: lgt_block(:, :)
!     !> heavy data array - neighbor data
!     integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
!     !> list of active blocks (heavy data)
!     integer(kind=ik), intent(in)        :: hvy_active(:)
!     !> number of active blocks (heavy data)
!     integer(kind=ik), intent(in)        :: hvy_n, ncomponents
!     integer(kind=ik), intent(inout)     :: Data_recvCounter(0:), Data_sendCounter(0:)
!     integer(kind=ik), intent(inout)     :: MetaData_recvCounter(0:), MetaData_sendCounter(0:)
!     logical, intent(in)                 :: count_internal
!
!     integer(kind=ik) :: k, sender_hvyID, sender_lgtID, myrank, N, neighborhood, recver_rank
!     integer(kind=ik) :: ijk(2,3), inverse, ierr, recver_hvyID, recver_lgtID,level_diff, status, new_size
!
!     call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
!     N = params%number_blocks
!
!     Data_recvCounter(:) = 0
!     Data_sendCounter(:) = 0
!     MetaData_recvCounter(:) = 0
!     MetaData_sendCounter(:) = 0
!
!     do k = 1, hvy_n
!         ! calculate light id
!         sender_hvyID = hvy_active(k)
!         call hvy2lgt( sender_lgtID, sender_hvyID, myrank, N )
!
!         ! loop over all neighbors
!         do neighborhood = 1, size(hvy_neighbor, 2)
!             ! neighbor exists
!             if ( hvy_neighbor( sender_hvyID, neighborhood ) /= -1 ) then
!                 ! neighbor light data id
!                 recver_lgtID = hvy_neighbor( sender_hvyID, neighborhood )
!                 ! calculate neighbor rank
!                 call lgt2proc( recver_rank, recver_lgtID, N )
!                 ! neighbor heavy id
!                 call lgt2hvy( recver_hvyID, recver_lgtID, recver_rank, N )
!
!                 ! define level difference: sender - receiver, so +1 means sender on higher level
!                 ! leveldiff = -1 : sender coarser than recver, interpolation on sender side
!                 ! leveldiff =  0 : sender is same level as recver
!                 ! leveldiff = +1 : sender is finer than recver, restriction is applied on sender side
!                 level_diff = lgt_block( sender_lgtID, params%Jmax + IDX_MESH_LVL ) - lgt_block( recver_lgtID, params%Jmax + IDX_MESH_LVL )
!
!
!                 if (recver_rank /= myrank .or. count_internal) then
!                     ! it now depends on the stage if we have to sent this data
!                     ! or not.
!                     ! In stage 1, only level_diff = {+1, 0} is treated
!                     ! In stage 2, only level_diff = -1
!
!                     ! send counter. how much data will I send to other mpiranks?
!                     ! why is this RECVER and not sender? Well, complicated. The amount of data on the sender patch
!                     ! is not the same as in the receiver patch, because we interpolate or downsample. We effectively
!                     ! transfer only the data the recver wants - not the extra data.
!                     ijk = ijkGhosts(:, :, neighborhood, level_diff, RECVER)
!
!                     Data_sendCounter(recver_rank) = Data_sendCounter(recver_rank) + &
!                     (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)
!
!                     ! recv counter. how much data will I recv from other mpiranks?
!                     ! This is NOT the same number as before
!                     inverse = inverse_neighbor(neighborhood, dim)
!
!                     ijk = ijkGhosts(:, :, inverse, -1*level_diff, RECVER)
!
!                     Data_recvCounter(recver_rank) = Data_recvCounter(recver_rank) + &
!                     (ijk(2,1)-ijk(1,1)+1) * (ijk(2,2)-ijk(1,2)+1) * (ijk(2,3)-ijk(1,3)+1)
!
!                     ! counter for integer buffer: for each neighborhood, we send 6 integers as metadata
!                     ! as this is a fixed number it does not depend on the type of neighborhood etc, so
!                     ! technically one would need only one for send/recv
!                     MetaData_sendCounter(recver_rank) = MetaData_sendCounter(recver_rank) + 6 ! FIVE
!
!                     MetaData_recvCounter(recver_rank) = MetaData_recvCounter(recver_rank) + 6 ! FIVE
!                 endif
!
!
!
!             end if ! neighbor exists
!         end do ! loop over all possible  neighbors
!     end do ! loop over all heavy active
!
!
!
!     ! NOTE: for the int buffer, we mosly start at some index l0 and then loop unitl
!     ! we find a -99 indicating the end of the buffer. this could be avoided by using
!     ! for instead of while loops in the main routines, but I do not have time now.
!     !
!     ! In the meantime, notice we extent the amount of data by one, to copy the last -99
!     ! to the buffers
!     MetaData_recvCounter(:) = MetaData_recvCounter(:) + 1
!     MetaData_sendCounter(:) = MetaData_sendCounter(:) + 1
!
!
!     ! NOTE ACTUAL SEND / RECV DATA IS NEQN
!     Data_recvCounter(:) = Data_recvCounter(:) * ncomponents
!     Data_sendCounter(:) = Data_sendCounter(:) * ncomponents
!
!     ! NOTE: this feature is against wabbits memory policy: we try to allocate the
!     ! whole memory of the machine on startup, then work with that. however, we have to
!     ! reserver portions of that memory for the state vector, the RHS slots, etc, and the ghost nodes
!     ! buffer. However, estimating those latter is difficult: it depends on the grid and the parallelization
!     if (sum(Data_recvCounter) > size(rData_recvBuffer, 1)) then
!         ! out-of-memory case: the preallocated buffer is not large enough.
!         write(*,'("rank=",i4," OOM for ghost nodes and increases its buffer size to 125%")') myrank
!         new_size = size(rData_recvBuffer,1)*125/100
!         deallocate(rData_recvBuffer)
!         allocate( rData_recvBuffer(1:new_size), stat=status )
!         if (status /= 0) call abort(999992, "Buffer allocation failed. Not enough memory?")
!     endif
!
!     if (sum(Data_sendCounter) > size(rData_sendBuffer, 1)) then
!         ! out-of-memory case: the preallocated buffer is not large enough.
!         write(*,'("rank=",i4," OOM for ghost nodes and increases its buffer size to 125%")') myrank
!         new_size = size(rData_sendBuffer,1)*125/100
!         deallocate(rData_sendBuffer)
!         allocate( rData_sendBuffer(1:new_size), stat=status )
!         if (status /= 0) call abort(999993, "Buffer allocation failed. Not enough memory?")
!     endif
! end subroutine
