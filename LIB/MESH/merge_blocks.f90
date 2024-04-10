!> \brief Merge 4 or 8 blocks into a new one (restriction operator).
!
!! This routine merges either 4 or 8 blocks into one new block with half the resolution
!! It is supposed that all blocks are on the same mpirank, thus gather_blocks_on_proc has
!! to be called prior to merging.
!! Note ghost nodes are not copied (thus synching is required afterwards). Why you ask? Because the coarser
!! ghost node layer is physically larger than the fine one. Hence, it cannot be filled just from the data
!! that we have on entry in ths routine.
!! Note we keep the light data synchronized among CPUS, so that after moving, all CPU are up-to-date with their light data.
!! However, the active lists are outdated after this routine.
! ********************************************************************************************
subroutine merge_blocks( params, hvy_block, lgt_blocks_to_merge )
    implicit none

    type (type_params), intent(in)      :: params                         !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)       !> heavy data array - block data
    integer(kind=ik), intent(inout)     :: lgt_blocks_to_merge(:)         !> list of blocks to merge, can contain 4 or 8 blocks
    real(kind=rk), ALLOCATABLE, save    :: tmpblock(:,:,:,:)
    integer(kind=ik)                    :: N_merge                        ! number of blocks to be merged, can be 4 or 8
    ! what CPU is responsible for merging:
    integer(kind=ik)                    :: data_rank(8)
    integer(kind=ik)                    :: heavy_ids(8), tree_ID          ! list of block ids, proc ranks

    integer(kind=tsize)                 :: treecode
    integer(kind=ik) :: i1, i2, im, i, g, level, lgt_merge_id, Jmax, hvy_merge_id, N
    integer(kind=ik), dimension(3)      ::  icoars1, icoars2, icoarsm, Bs, ifine1, ifine2

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors and tree_N are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    ! number of blocks to be merged
    N_merge = size(lgt_blocks_to_merge,1)
    Bs = params%Bs
    g  = params%g
    Jmax = params%Jmax
    ! details of merged block
    level = lgt_block( lgt_blocks_to_merge(1), IDX_MESH_LVL )
    tree_ID = lgt_block( lgt_blocks_to_merge(1), IDX_TREE_ID )
    treecode = get_tc(lgt_block( lgt_blocks_to_merge(1), IDX_TC_1 : IDX_TC_2 ))

    ! Check which CPU holds the blocks. The CPU will also hold the merged, new block
    do i = 1, N_merge
        call lgt2proc( data_rank(i), lgt_blocks_to_merge(i), params%number_blocks )
    enddo

!-------------------------------------------------------------------------------
#ifdef DEV
    if ( N_merge /= 4 .and. N_merge /= 8) then
        call abort('You try to merge neither n=4 or 8 blocks...this cannot work.')
    endif

    ! Check if all blocks lie on the same rank
    if ( maxval(data_rank(1:N_merge)-data_rank(1)) /= 0 ) then
        call abort("You try to merge blocks on different ranks, but you must call gather_ranks before.")
    endif

    do i = 1, size(lgt_blocks_to_merge)
        treecode = get_tc(lgt_blocks_to_merge(i), IDX_TC_1 : IDX_TC_2)
        if (tc_get_digit_at_level_b( treecode, params%dim, level, params%Jmax) /= i-1) then
            call abort(647483," You try to merge blocks which do not belong together")
        endif
        do i1 = 1, level-1
            if (tc_get_digit_at_level_b( treecode, params%dim, level, params%Jmax) /=
                tc_get_digit_at_level_b( get_tc(lgt_blocks_to_merge(1), IDX_TC_1 : IDX_TC_2), params%dim, level, params%Jmax)) then
                call abort(647483," You try to merge blocks which do not belong together")
            endif
        enddo
    enddo
#endif
!-------------------------------------------------------------------------------



    ! allocate tmp buffer (used only for biorthogonal wavelets, not harten-style multiresolution)
    if ( .not. allocated(tmpblock) ) then
        allocate( tmpblock(size(hvy_block,1), size(hvy_block,2), size(hvy_block,3), size(hvy_block,4)) )
        tmpblock = 0.0_rk
    else
        ! it may happen that we use this routine for fields with different # components: in this case,
        ! tmp block has to change size.
        if (size(tmpblock,1)/=size(hvy_block,1) .or.size(tmpblock,2)/=size(hvy_block,2) &
        .or.size(tmpblock,4)/=size(hvy_block,4) .or.size(tmpblock,3)/=size(hvy_block,3) ) then
            deallocate(tmpblock)
            allocate( tmpblock(size(hvy_block,1), size(hvy_block,2), size(hvy_block,3), size(hvy_block,4)) )
        endif
    endif


    !-------------------------------------------------------------------------------
    ! merge the specified blocks into one new block. Merging is done in two steps,
    ! first for light data (which all CPUS do redundantly, so light data is kept synched)
    ! Then only the responsible rank will perform the heavy data merging.

    !-------------------------------------------------------------------------------
    ! a) light data (collective operation)
    ! fetch a free light ID for the merged blocks
    call get_free_local_light_id(params, data_rank(1), lgt_merge_id, message="merge_blocks")
    ! create light data entry for the new block
    lgt_block( lgt_merge_id, : ) = -1
    call set_tc(lgt_block( lgt_merge_id, IDX_TC_1:IDX_TC_2), tc_clear_until_level_b(treecode, &
        dim=params%dim, level=level-1, max_level=params%Jmax))
    lgt_block( lgt_merge_id, IDX_MESH_LVL ) = level-1
    lgt_block( lgt_merge_id, idx_refine_sts ) = 0
    lgt_block( lgt_merge_id, IDX_TREE_ID ) = tree_ID

    !-------------------------------------------------------------------------------
    ! b) heavy data merging (individual operation)
    if ( data_rank(1) == params%rank) then
        ! loop over the sisters and store the corresponding mpirank and heavy index
        ! in a list (still, we are dealing with all four sisters)
        do i = 1, N_merge
            call lgt2hvy( heavy_ids(i), lgt_blocks_to_merge(i), data_rank(1), params%number_blocks )
        enddo

        ! get heavy id of merge block
        call lgt2hvy( hvy_merge_id, lgt_merge_id, data_rank(1), params%number_blocks )

        ! detail is not computed (yet)
        hvy_details(:, hvy_merge_id) = -1.0_rk

        ! indices of subblocks on new, merged, coarser block
        do i = 1, params%dim
            icoars1(i) = g+1
            icoars2(i) = Bs(i)+g
            icoarsm(i) = Bs(i)/2 + g
        enddo

        ! indices on the finer blocks which we merge to a coarse one
        do i = 1,  params%dim
            ! [INKSCAPE]: Please note the fortran code always runs
            ! until X:2:Bs+g the important number is the
            ! starting index
            ! BS even
            ifine1(i) = g+1 ! start point of first index, low range
            ifine2(i) = g+1 ! start point of second index, high range
        enddo


        if (N_merge == 4) then
            ! ************ 2D case ***********************
            ! biorthogonal wavelets: apply a low-pass filter (called H) prior to decimation
            ! sister 0
            call restriction_prefilter_vct(params, hvy_block( :,:,:,:, heavy_ids(1) ), tmpblock) ! low-pass filtering
            hvy_block(icoars1(1):icoarsm(1), icoars1(2):icoarsm(2), :, :, hvy_merge_id) = tmpblock( ifine1(1):Bs(1)+g:2, ifine1(2):Bs(2)+g:2, :,:)

            ! sister 1
            call restriction_prefilter_vct(params, hvy_block( :,:,:,:, heavy_ids(2) ), tmpblock)
            hvy_block(icoars1(1):icoarsm(1), icoarsm(2)+1:icoars2(2), :, :, hvy_merge_id) = tmpblock( ifine1(1):Bs(1)+g:2, ifine2(2):Bs(2)+g:2, :,:)

            ! sister 2
            call restriction_prefilter_vct(params, hvy_block( :,:,:,:, heavy_ids(3) ), tmpblock)
            hvy_block(icoarsm(1)+1:icoars2(1), icoars1(2):icoarsm(2), :, :, hvy_merge_id) = tmpblock( ifine2(1):Bs(1)+g:2, ifine1(2):Bs(2)+g:2, :,:)

            ! sister 3
            call restriction_prefilter_vct(params, hvy_block( :,:,:,:, heavy_ids(4) ), tmpblock)
            hvy_block(icoarsm(1)+1:icoars2(1), icoarsm(2)+1:icoars2(2), :, :, hvy_merge_id) = tmpblock( ifine2(1):Bs(1)+g:2, ifine2(2):Bs(2)+g:2, :,:)

        elseif (N_merge == 8) then
            ! ************ 3D case ***********************
            ! sister 0
            call restriction_prefilter_vct(params, hvy_block(:,:,:,:,heavy_ids(1)), tmpblock)
            hvy_block(icoars1(1):icoarsm(1)  , icoars1(2):icoarsm(2)  , icoars1(3):icoarsm(3),  :, hvy_merge_id )  = tmpblock( ifine1(1):Bs(1)+g:2, ifine1(2):Bs(2)+g:2, ifine1(3):Bs(3)+g:2, :)

            ! sister 1
            call restriction_prefilter_vct(params, hvy_block(:,:,:,:,heavy_ids(2)), tmpblock)
            hvy_block(icoars1(1):icoarsm(1)  , icoarsm(2)+1:icoars2(2), icoars1(3):icoarsm(3),  :, hvy_merge_id )  = tmpblock( ifine1(1):Bs(1)+g:2, ifine2(2):Bs(2)+g:2, ifine1(3):Bs(3)+g:2, :)

            ! sister 2
            call restriction_prefilter_vct(params, hvy_block(:,:,:,:,heavy_ids(3)), tmpblock)
            hvy_block(icoarsm(1)+1:icoars2(1), icoars1(2):icoarsm(2)  , icoars1(3):icoarsm(3),  :, hvy_merge_id )  = tmpblock( ifine2(1):Bs(1)+g:2, ifine1(2):Bs(2)+g:2, ifine1(3):Bs(3)+g:2, :)

            ! sister 3
            call restriction_prefilter_vct(params, hvy_block(:,:,:,:,heavy_ids(4)), tmpblock)
            hvy_block(icoarsm(1)+1:icoars2(1), icoarsm(2)+1:icoars2(2), icoars1(3):icoarsm(3),  :, hvy_merge_id )  = tmpblock( ifine2(1):Bs(1)+g:2, ifine2(2):Bs(2)+g:2, ifine1(3):Bs(3)+g:2, :)

            ! sister 4
            call restriction_prefilter_vct(params, hvy_block(:,:,:,:,heavy_ids(5)), tmpblock)
            hvy_block(icoars1(1):icoarsm(1)  , icoars1(2):icoarsm(2)  , icoarsm(3)+1:icoars2(3), :, hvy_merge_id ) = tmpblock( ifine1(1):Bs(1)+g:2, ifine1(2):Bs(2)+g:2, ifine2(3):Bs(3)+g:2, :)

            ! sister 5
            call restriction_prefilter_vct(params, hvy_block(:,:,:,:,heavy_ids(6)), tmpblock)
            hvy_block(icoars1(1):icoarsm(1)  , icoarsm(2)+1:icoars2(2), icoarsm(3)+1:icoars2(3), :, hvy_merge_id ) = tmpblock( ifine1(1):Bs(1)+g:2, ifine2(2):Bs(2)+g:2, ifine2(3):Bs(3)+g:2, :)

            ! sister 6
            call restriction_prefilter_vct(params, hvy_block(:,:,:,:,heavy_ids(7)), tmpblock)
            hvy_block(icoarsm(1)+1:icoars2(1), icoars1(2):icoarsm(2)  , icoarsm(3)+1:icoars2(3), :, hvy_merge_id ) = tmpblock( ifine2(1):Bs(1)+g:2, ifine1(2):Bs(2)+g:2, ifine2(3):Bs(3)+g:2, :)

            ! sister 7
            call restriction_prefilter_vct(params, hvy_block(:,:,:,:,heavy_ids(8)), tmpblock)
            hvy_block(icoarsm(1)+1:icoars2(1), icoarsm(2)+1:icoars2(2), icoarsm(3)+1:icoars2(3), :, hvy_merge_id ) = tmpblock( ifine2(1):Bs(1)+g:2, ifine2(2):Bs(2)+g:2, ifine2(3):Bs(3)+g:2, :)
        endif
    endif

    ! merging is complete now, remove the original blocks from light data:
    do i = 1, N_merge
        lgt_block( lgt_blocks_to_merge(i), : ) = -1
    enddo
end subroutine merge_blocks
