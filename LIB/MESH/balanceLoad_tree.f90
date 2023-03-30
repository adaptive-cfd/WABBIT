!> \image html balancing.svg "Load balancing" width=400
! ********************************************************************************************
subroutine balanceLoad_tree( params, hvy_block, tree_ID)
    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID
    !=============================================================================

    integer(kind=ik)  :: myrank, mpirank, mpirank_owner, ierr, &
                         k, N, l, com_i, com_N, heavy_id, sfc_id, neq, &
                         lgt_free_id, hvy_free_id, hilbertcode(params%Jmax), total_weight, mpisize
    integer(kind=ik), allocatable, save :: sfc_com_list(:,:), sfc_sorted_list(:,:), cpu_weight(:)
    real(kind=rk) :: t0, t1, avg_weight, remaining_weight, ideal_weight
    logical       :: is_predictable

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    ! on only one proc, no balancing is required
    if (params%number_procs == 1) return

    t0 = MPI_Wtime()
    myrank = params%rank
    mpisize = params%number_procs
    N = params%number_blocks
    neq = params%n_eqn

    ! allocate block to proc lists
    if (.not.allocated(cpu_weight)) allocate( cpu_weight(1:mpisize))
    ! allocate sfc com list, maximal number of communications is when every proc wants to send all of his blocks
    ! NOTE: it is not necessary or wise to reset this array (it is huge!)
    if (.not.allocated(sfc_com_list)) allocate( sfc_com_list( mpisize*params%number_blocks, 3 ) )
    ! allocate space filling curve list, number of elements is the number of active blocks
    ! and for each block, we store the space-filling-curve-index and the lgt ID
    if (.not.allocated(sfc_sorted_list)) allocate( sfc_sorted_list( size(lgt_block,1), 3) )

    ! for each block compute its weight. Note for technical reasons, weight is an integer
    ! but ne can assign a default weight of say 10 for ordinary blocks, which means we can
    ! in pratice also assign say 13/10 weight for a block. Total_weight is the sum of all blocks
    ! to be distributed: this number is divided my N_cpu to get the approx. local portion of the weight.
    call assignLoadbalancingWeight_tree(params, tree_ID, total_weight)


    !---------------------------------------------------------------------------------
    ! 1st: calculate space filling curve index for all blocks
    !---------------------------------------------------------------------------------
    ! COLLECTIVE OPERATION (performed redundantly on all ranks, no comm)
    t1 = MPI_wtime()
    select case (params%block_distribution)
    case("sfc_z")
        !-----------------------------------------------------------
        ! Z - curve
        !-----------------------------------------------------------
        if (params%dim == 3) then
            do k = 1, lgt_n(tree_ID)
                call treecode_to_sfc_id_3D( sfc_id, lgt_block( lgt_active(k, tree_ID), 1:params%Jmax ), params%Jmax )
                sfc_sorted_list(k, 1:3) = (/sfc_id, lgt_active(k, tree_ID), lgt_block(lgt_active(k, tree_ID), params%Jmax+IDX_REFINE_STS)/)
            end do
        else
            do k = 1, lgt_n(tree_ID)
                call treecode_to_sfc_id_2D( sfc_id, lgt_block( lgt_active(k, tree_ID), 1:params%Jmax ), params%Jmax )
                sfc_sorted_list(k, 1:3) = (/sfc_id, lgt_active(k, tree_ID), lgt_block(lgt_active(k, tree_ID), params%Jmax+IDX_REFINE_STS)/)
            end do
        endif
    case("sfc_hilbert")
        !-----------------------------------------------------------
        ! Hilbert curve
        !-----------------------------------------------------------
        if (params%dim == 3) then
            do k = 1, lgt_n(tree_ID)
                ! transfer treecode to hilbertcode
                call treecode_to_hilbertcode_3D( lgt_block( lgt_active(k, tree_ID), 1:params%Jmax ), hilbertcode, params%Jmax)
                ! calculate sfc position from hilbertcode
                call treecode_to_sfc_id_3D( sfc_id, hilbertcode, params%Jmax )
                sfc_sorted_list(k, 1:3) = (/sfc_id, lgt_active(k, tree_ID), lgt_block(lgt_active(k, tree_ID), params%Jmax+IDX_REFINE_STS)/)
            end do
        else
            do k = 1, lgt_n(tree_ID)
                ! transfer treecode to hilbertcode
                call treecode_to_hilbertcode_2D( lgt_block( lgt_active(k, tree_ID), 1:params%Jmax ), hilbertcode, params%Jmax)
                ! calculate sfc position from hilbertcode
                call treecode_to_sfc_id_2D( sfc_id, hilbertcode, params%Jmax )
                sfc_sorted_list(k, 1:3) = (/sfc_id, lgt_active(k, tree_ID), lgt_block(lgt_active(k, tree_ID), params%Jmax+IDX_REFINE_STS)/)
            end do
        endif

    case default
        call abort(210309, "The block_dist is unkown"//params%block_distribution)

    end select

    ! sort sfc_list according to the first dimension, thus the position on
    ! the space filling curve (this was a bug, fixed: Thomas, 13/03/2018)
    if (lgt_n(tree_ID) > 1) then
        call quicksort_ik(sfc_sorted_list, 1, lgt_n(tree_ID), 1, 3)
    end if
    call toc( "balanceLoad_tree (SFC+sort)", MPI_wtime()-t1 )


    ! We now know what the total weight to be distributed is, and each block has been assigned
    ! a unique index on the space filling curve. This list of blocks is sorted by the
    ! space filling curve index.
    ! -> It is now a matter of distributing (contiguous) chunks of blocks to the
    !    mpiranks, and doing so as well as possible: We want the maximum number of load
    !    on a mpirank to be as small as possible. The slowest CPU determines execution speed
    !    (if one idles that is not as bad as if one is very slow because then all other CPU idle)

    !---------------------------------------------------------------------------------
    ! 2nd: plan communication (fill list of blocks to transfer)
    !---------------------------------------------------------------------------------
    t1 = MPI_wtime()
    ! mpirank: process responsible for current part of sfc
    ! mpirank_owner: process who stores data of sfc element

    ! we start the loop on the root rank (0), then assign the first elements
    ! of the SFC, then to second rank, etc. (thus: mpirank is a loop variable)
    mpirank = 0

    ! communication counter. each communication (=send and receive) is stored
    ! in a long list
    com_i = 0
    cpu_weight = 0_ik

    ideal_weight = dble(total_weight) / dble(mpisize)

remaining_weight = real( sum(sfc_sorted_list(1:lgt_n(tree_ID), 3)), kind=rk)
avg_weight = remaining_weight / real( mpisize-mpirank, kind=rk)
! write(*,*) avg_weight

    ! prepare lists for transfering of blocks
    ! COLLECTIVE OPERATION (performed redundantly on all ranks, no comm)
    do k = 1, lgt_n(tree_ID)
        ! if (mpirank /= mpisize-1) then
        !     ! for the last rank, the weight is not recomputed (it gets the remaining
        !     ! blocks anyways). Otherwise, from the remaining blocks, recompute the desired
        !     ! average weight per mpirank:
        !     remaining_weight = real( sum(sfc_sorted_list(k:lgt_n(tree_ID), 3)), kind=rk)
        !     avg_weight = remaining_weight / real( mpisize-mpirank, kind=rk)
        ! endif


        ! if ((cpu_weight(mpirank+1) + sfc_sorted_list(k,3) < nint(avg_weight)) .or. (mpirank==mpisize-1)) then
        !     ! cpu can take more and is still not full:
        !     ! (or we are at the last mpirank, which gets the rest)
        !     cpu_weight(mpirank+1) = cpu_weight(mpirank+1) + sfc_sorted_list(k,3)
        ! else
        !     ! if adding one more block is effectively only one more weight unit (ie the weight of a
        !     ! standard block)
        !     if (cpu_weight(mpirank+1) + sfc_sorted_list(k,3) - nint(avg_weight) <= 20_ik) then
        !         cpu_weight(mpirank+1) = cpu_weight(mpirank+1) + sfc_sorted_list(k,3)
        !     endif
        !     ! Now, the mpirank is full -> the next blocks will be distributed to the next mpirank
        !     ! increase cpu counter, but only until we are at the last rank
        !     mpirank = min(mpirank+1, mpisize-1)
        ! endif

        if (cpu_weight(mpirank+1) + sfc_sorted_list(k,3) >= nint(avg_weight) + 10_ik) then
            ! Now, the mpirank is full -> the next blocks will be distributed to the next mpirank
            ! increase cpu counter, but only until we are at the last rank
            mpirank = min(mpirank+1, mpisize-1)

            remaining_weight = real( sum(sfc_sorted_list(k:lgt_n(tree_ID), 3)), kind=rk)
            avg_weight = remaining_weight / real( mpisize-mpirank, kind=rk)
            ! write(*,*) avg_weight
        endif

        cpu_weight(mpirank+1) = cpu_weight(mpirank+1) + sfc_sorted_list(k,3)

        ! find out on which mpirank lies the block that we're looking at
        call lgt2proc( mpirank_owner, sfc_sorted_list(k,2), params%number_blocks )

        ! does this block lie on the right mpirank, i.e., the current part of the
        ! SFC? if so, nothing needs to be done. otherwise, the following if is active
        if ( mpirank /= mpirank_owner ) then
            ! as this block is one the wrong rank, it will be sent away from its
            ! current owner (mpirank_owner) to the owner of this part of the
            ! SFC (mpirank)

            ! save this send+receive operation in the list of planned communications
            ! column
            !    1     sender mpirank
            !    2     receiver mpirank
            !    3     block light data id
            com_i = com_i + 1
            sfc_com_list(com_i, 1) = mpirank_owner        ! sender mpirank
            sfc_com_list(com_i, 2) = mpirank              ! receiver mpirank
            sfc_com_list(com_i, 3) = sfc_sorted_list(k,2) ! light id of block
        end if
    end do

    ! if (myrank==0) then
    !     write(*,'(300(i2,1x))') sfc_sorted_list(1:lgt_n(tree_ID),3)
    !     write(*,*) "ideal_weight", ideal_weight, lgt_n(tree_ID), total_weight
    !     write(*,'("cpu_weight=",16(i3,1x))') cpu_weight
    !     write(*,*) "xfers", com_i
    ! endif

    !---------------------------------------------------------------------------------
    ! 3rd: actual communication (send/recv)
    !---------------------------------------------------------------------------------
    call block_xfer( params, sfc_com_list, com_i, hvy_block, msg="balanceLoad_tree" )
    call toc( "balanceLoad_tree (comm)", MPI_wtime()-t1 )

    ! the block xfer changes the light data, and afterwards active lists are outdated.
    ! NOTE: an idea would be to also xfer the neighboring information (to save the updateNeighbors_tree
    ! call) but that is tricky: the neighbor list contains light ID of the neighbors, and those
    ! also change with the xfer.
    call updateMetadata_tree(params, tree_ID)

    ! timing
    call toc( "balanceLoad_tree (TOTAL)", MPI_wtime()-t0 )
end subroutine balanceLoad_tree


subroutine assignLoadbalancingWeight_tree(params, treeID, total_weight)
    implicit none
    type (type_params), intent(in)      :: params
    integer(kind=ik) :: total_weight
    integer(kind=ik) :: k, hvyID, lgtID, treeID, n, lgtID_neighbor, level_me, level_neighbor

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    !
    !
    ! NOTE Functionality currently disabled NOTE
    !
    ! we just assign the same weight for all blocks

    ! ! Problem: hvy_weight is a distributed dataset and not available on all procs.
    ! ! Solution: use the slot for "refinement status" temporarily for the wheight and also
    ! ! sync it with refinementStatusOnly
    ! do k = 1, hvy_n(treeID)
    !     hvyID = hvy_active(k, treeID)
    !     call hvy2lgt(lgtID, hvyID, params%rank, params%number_blocks)
    !
    !     ! default weight
    !     lgt_block( lgtID, params%Jmax+IDX_REFINE_STS ) = 10_ik
    !
    !     ! do n = 1, size(hvy_neighbor, 2)
    !     !     if (hvy_neighbor(hvyID, n) /= -1_ik) then
    !     !         ! neighbor light data id
    !     !         lgtID_neighbor = hvy_neighbor( hvyID, n )
    !     !         level_me       = lgt_block( lgtID, params%Jmax + IDX_MESH_LVL )
    !     !         level_neighbor = lgt_block( lgtID_neighbor, params%Jmax + IDX_MESH_LVL )
    !     !
    !     !         if (level_neighbor < level_me) then
    !     !             ! block affect by coarse extension
    !     !             lgt_block( lgtID, params%Jmax+IDX_REFINE_STS ) = 40_ik
    !     !             exit ! the neighbor loop
    !     !         endif
    !     !     endif
    !     ! enddo
    ! enddo
    !
    ! ! call synchronize_lgt_data(params, refinement_status_only=.true.)

    !!!!!!!!!!!!!!!!!!!
    do k = 1, lgt_n(treeID)
        lgtID = lgt_active(k, treeID)
        ! default weight
        lgt_block( lgtID, params%Jmax+IDX_REFINE_STS ) = 10_ik
    enddo
    !!!!!!!!!!!!!!!!!!!

    total_weight = 0_ik
    do k = 1, lgt_n(treeID)
        total_weight = total_weight + lgt_block( lgt_active(k, treeID), params%Jmax+IDX_REFINE_STS )
    enddo
end subroutine
