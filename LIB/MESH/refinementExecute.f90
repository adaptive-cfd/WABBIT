!> \brief Refine block by interpolating it to the next finer level and creating
!!  four(2D) or eight(3D) new blocks, each carrying a part of the interpolated data.
!> Interpolated data is set in array, and last daughter overwrites motherblock
subroutine refineBlock(params, hvy_block, hvyID, tree_ID)
    implicit none

    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)             !> heavy data array - block data
    integer(kind=ik), intent(in)        :: hvyID                                !> Block to be reined
    integer(kind=ik), intent(in)        :: tree_ID

    integer(kind=ik)                    :: k, N, dF, k_daughter                 ! loop variables
    integer(kind=ik)                    :: rank, dim
    integer(kind=ik)                    :: g                                    ! grid parameter
    integer(kind=ik), dimension(3)      :: Bs
    real(kind=rk), allocatable, save    :: data_predict_fine(:,:,:,:)  ! data fields for interpolation
    integer(kind=ik)                    :: lgt_free_id, free_hvy_id, lgt_id, level, idx_d(2,3)
    integer(kind=tsize)                 :: treecode

    dim  = params%dim
    N    = params%number_blocks
    rank = params%rank
    Bs   = params%Bs
    g    = params%g

    ! Allocate data buffer for interpolation (data_predict_fine).
    ! The prediction below turns a block with size Bs+2*g (including ghost nodes) into a refined block
    ! of size 2*(Bs+2*g)-1.
    !
    ! NOTE: the predictor for the refinement acts on the extended blocks i.e. it
    ! includes the ghost nodes layer. Therefore, you MUST call sync_ghosts before this routine.
    ! The datafield for prediction is one level up, i.e. it contains Bs+g + (Bs+2g-1) points.
    !
    ! The refinement routine may be called for different trees, for example the flow tree and the mask tree
    ! and those trees may differ in their number of components. Ensure the buffer has the right size, if not, 
    ! dealloacte it (and subsequently re-alloacte it correctly)
    if (allocated(data_predict_fine)) then
        if (size(data_predict_fine, 4) < size(hvy_block,4)) deallocate(data_predict_fine)
    endif

    ! actual allocation of data buffer for interpolation
    if (dim == 2) then
        if (.not. allocated(data_predict_fine)) allocate( data_predict_fine( 2*(Bs(1)+2*g)-1, 2*(Bs(2)+2*g)-1, 1, size(hvy_block,4)) )
    else
        if (.not. allocated(data_predict_fine)) allocate( data_predict_fine( 2*(Bs(1)+2*g)-1, 2*(Bs(2)+2*g)-1, 2*(Bs(3)+2*g)-1, size(hvy_block,4)) )
    endif

    ! extract treecode and mesh level
    call hvy2lgt(lgt_id, hvyID, rank, N)
    treecode = get_tc(lgt_block( lgt_id, IDX_TC_1 : IDX_TC_2 ))
    level    = lgt_block( lgt_id, IDX_MESH_LVL )

    !------------------------------------------------------------------------------------------------------
    ! 1st Step: interpolate block data
    !------------------------------------------------------------------------------------------------------

    ! loop over all components
    do dF = 1, size(hvy_block,4)
        ! NOTE: the refinement interpolation acts on the entire block including ghost nodes.
        call prediction(hvy_block(:, :, :, dF, hvyID), data_predict_fine(:, :, :, df), params%order_predictor)
    end do

    !------------------------------------------------------------------------------------------------------
    ! 2nd Step: split new data and write into new blocks
    !------------------------------------------------------------------------------------------------------

    ! create all 4 or 8 new daughters
    do k_daughter = 0, 2**(params%dim) -1

        ! find a free light id on this rank
        if (k_daughter < 2**(params%dim) -1) then
            call get_free_local_light_id( params, rank, lgt_free_id, message="refinement_execute")
            call lgt2hvy( free_hvy_id, lgt_free_id, rank, N )
        else
            ! last daughter overwrites mother
            free_hvy_id = hvyID
            call hvy2lgt( lgt_free_id, free_hvy_id, rank, N )
        endif

        treecode = tc_set_digit_at_level_b(treecode, k_daughter, dim=params%dim, level=level+1, max_level=params%Jmax)

        ! init array - needed to change values if never adressed
        lgt_block( lgt_free_id, : ) = -1

        ! write new light data
        ! new treecode
        call set_tc(lgt_block( lgt_free_id, IDX_TC_1 : IDX_TC_2 ), treecode)
        ! new block is on (level + 1)
        lgt_block( lgt_free_id, IDX_MESH_LVL ) = level+1
        ! new blocks have refinement_status==0 (STAY)
        lgt_block( lgt_free_id, idx_refine_sts ) = REF_FRESHLY_REFINED
        ! the tree_ID is the same as the one of the mother block
        lgt_block( lgt_free_id, IDX_TREE_ID ) = tree_ID

        ! save interpolated data - select data from correct quadrant
        ! k=0 -> X_l, Y_l; k=1 -> X_l, Y_r; k=2 -> X_r, Y_l; k=3 -> X_r, Y_r;
        ! compute indices of new block, include ghost points
        idx_d(1,1) = 1*g +     1 + Bs(1) * modulo(k_daughter/2,2)
        idx_d(2,1) = 3*g + Bs(1) + Bs(1) * modulo(k_daughter/2,2)
        idx_d(1,2) = 1*g +     1 + Bs(2) * modulo(k_daughter,2)
        idx_d(2,2) = 3*g + Bs(1) + Bs(2) * modulo(k_daughter,2)
        if (dim == 3) then
            idx_d(1,3) = 1*g +     1 + Bs(3) * (k_daughter/4)
            idx_d(2,3) = 3*g + Bs(3) + Bs(3) * (k_daughter/4)
        endif

        if (dim == 2) then
            hvy_block( :, :, 1, 1:size(hvy_block,4), free_hvy_id ) = \
                data_predict_fine( idx_d(1,1):idx_d(2,1), idx_d(1,2):idx_d(2,2), 1, 1:size(hvy_block,4))
        else
            hvy_block( :, :, :, 1:size(hvy_block,4), free_hvy_id ) = \
                data_predict_fine( idx_d(1,1):idx_d(2,1), idx_d(1,2):idx_d(2,2), idx_d(1,3):idx_d(2,3), 1:size(hvy_block,4))
        endif
    end do

end subroutine refineBlock



!> \brief Refine block directly to spaghetti decomposition
!!    SC is copied to correct quadrant / octant
!!    WC is set to 0
!> All daughter blocks are created and at last mother is set as inactive
!  Module MPI could be used in order to do the copying more easily, as we have the indices there in ijkPatches
subroutine refineBlock2SpaghettiWD(params, hvy_block, hvyID, tree_ID)
    implicit none

    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)             !> heavy data array - block data
    integer(kind=ik), intent(in)        :: hvyID                                !> Block to be reined
    integer(kind=ik), intent(in)        :: tree_ID

    integer(kind=ik)                    :: k, N, dF, k_daughter                 ! loop variables
    integer(kind=ik)                    :: rank, dim, g, idx_daughter(2,3), idx_mother(2,3)
    integer(kind=ik), dimension(3)      :: Bs
    integer(kind=ik)                    :: lgt_free_id, free_hvy_id, lgt_id
    integer(kind=tsize)                 :: treecode
    integer(kind=ik)                    :: level, lgtID
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: tmp_wd  ! used for data juggling

    if (allocated(tmp_wd)) then
        if (size(tmp_wd, 4) < size(hvy_block, 4)) deallocate(tmp_wd)
    endif
    if (.not. allocated(tmp_wd)) allocate(tmp_wd(1:size(hvy_block, 1), 1:size(hvy_block, 2), 1:size(hvy_block, 3), 1:size(hvy_block, 4)) )


    dim = params%dim
    N = params%number_blocks
    rank = params%rank
    Bs = params%Bs
    g  = params%g
    call hvy2lgt(lgtID, hvyID, rank, N)

    ! extract treecode and mesh level
    treecode = get_tc(lgt_block( lgt_id, IDX_TC_1 : IDX_TC_2 ))
    level    = lgt_block( lgt_id, IDX_MESH_LVL )

    ! ------------------------------------------------------------------------------------------------------
    ! second: split new data and write into new blocks
    !--------------------------
    ! create all four daughters
    do k_daughter = 0,2**(params%dim) -1
        ! find a free light id on this rank
        if (k_daughter < 2**(params%dim) -1) then
            call get_free_local_light_id( params, rank, lgt_free_id, message="refinement_execute")
            call lgt2hvy( free_hvy_id, lgt_free_id, rank, N )
        endif

        treecode = tc_set_digit_at_level_b(treecode, k_daughter, dim=params%dim, level=level+1, max_level=params%Jmax)

        ! init array - needed to change values if never adressed
        lgt_block( lgt_free_id, : ) = -1

        ! write new light data
        ! new treecode
        call set_tc(lgt_block( lgt_free_id, IDX_TC_1 : IDX_TC_2 ), treecode)
        ! new level + 1
        lgt_block( lgt_free_id, IDX_MESH_LVL ) = level+1
        ! new blocks have refinement_status==0 (STAY)
        lgt_block( lgt_free_id, idx_refine_sts ) = 0
        ! the tree_ID is the same as the one of the mother block
        lgt_block( lgt_free_id, IDX_TREE_ID ) = tree_ID

        ! save interpolated data, loop over all datafields - select data from correct quadrant
        ! k=0 -> X_l, Y_l; k=1 -> X_l, Y_r; k=2 -> X_r, Y_l; k=3 -> X_r, Y_r;
        do dF = 1, size(hvy_block,4)
            ! Init values as 0 to set all WC to 0
            tmp_wd(:, :, :, df) = 0.0_rk

            ! compute bounds from which quadrant / octant to copy
            idx_mother = -1
            idx_daughter = -1
            call get_indices_of_ghost_patch(params%Bs, params%g, params%dim, -k_daughter-1, idx_daughter, gminus=g, gplus=g, lvl_diff=+1)
            call get_indices_of_ghost_patch(params%Bs, params%g, params%dim, -k_daughter-9, idx_mother, gminus=g, gplus=g, lvl_diff=-1)

            ! now copy values from this quadrant / octant to SC area in Mallat format
            if (dim == 2) then
                tmp_wd(idx_daughter(1,1):idx_daughter(2,1), idx_daughter(1,2):idx_daughter(2,2), 1, df) = &
                    hvy_block(idx_mother(1,1):idx_mother(2,1), idx_mother(1,2):idx_mother(2,2), 1, df, hvyID)
            else
                tmp_wd(idx_daughter(1,1):idx_daughter(2,1), idx_daughter(1,2):idx_daughter(2,2), &
                    idx_daughter(1,3):idx_daughter(2,3), df) = &
                    hvy_block(idx_mother(1,1):idx_mother(2,1), idx_mother(1,2):idx_mother(2,2), &
                    idx_mother(1,3):idx_mother(2,3), df, hvyID)
            endif

        end do
        ! transform to spaghetti form as this is what we usually work with
        call Mallat2Spaghetti_block(params, tmp_wd(:, :, :, 1:size(hvy_block,4)), hvy_block(:, :, :, 1:size(hvy_block,4), free_hvy_id))
    end do

    ! delete mother
    lgt_block(lgtID, :) = -1
    lgt_block(lgtID, IDX_REFINE_STS) = 0

end subroutine refineBlock2SpaghettiWD



!> \brief Refine mesh (3D version). All cpu loop over their heavy data and check if the refinement
!! flag +1 is set on the block. If so, we take this block, interpolate it to the next finer
!! level and create four(2D) or eight(3D) new blocks, each carrying a part of the interpolated data.
!! As all CPU first work individually, the light data array is synced afterwards.
!! Basically, this function is just a wrapper 
!
!> \note The interpolation (or prediction) operator here is applied to a block INCLUDING
!! any ghost nodes. You must sync first.
!
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays \n
! ********************************************************************************************
subroutine refinement_execute_tree( params, hvy_block, tree_ID, time )

    implicit none

    type (type_params), intent(in)      :: params                               !< user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)             !< heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID                              !< tree to work on
    real(kind=rk), intent(in), optional :: time                                 !< current simulation time, used for development output only

    integer(kind=ik)                    :: i_b, lgt_id, ierr
    integer(kind=ik)                    :: blocks_refined
    integer(kind=ik), allocatable, save :: blocks_refined_list(:)
    character(len=clong) :: format_string, string_prepare

    ! every proc loop over its active heavy data array
    blocks_refined = 0
    do i_b = 1, hvy_n(tree_ID)
        ! light data id
        call hvy2lgt( lgt_id, hvy_active(i_b, tree_ID), params%rank, params%number_blocks )

        ! block wants to refine
        if ( (lgt_block( lgt_id, idx_refine_sts) == +1) ) then
            blocks_refined = blocks_refined + 1
            call refineBlock(params, hvy_block, hvy_active(i_b, tree_ID), tree_ID)
        endif
    end do

    if (params%debug_refinement) then
        if (.not. allocated(blocks_refined_list)) allocate(blocks_refined_list(1:params%number_procs))

        ! debug output to see how many blocks have been refined, gather information on rank 0 and print to file
        call MPI_GATHER(blocks_refined, 1, MPI_INTEGER, blocks_refined_list, 1, MPI_INTEGER, 0, WABBIT_COMM, ierr)
        if (params%rank == 0) then
            open(unit=99, file=trim("debug_refinement.csv"), status="unknown", position="append")
            string_prepare = "-1.0E+00,"  ! set negative time, just to have csv with the same length in every row
            if (present(time)) write(string_prepare,'(es16.6,",")') time
            if (params%number_procs == 1) then
                write(99,'(A, i0)') trim(adjustl(string_prepare)), blocks_refined_list(1)
            else
                write(format_string, '("(A, i0,",i0,"("","",i0))")') params%number_procs - 1
                write(99,format_string) trim(adjustl(string_prepare)), blocks_refined_list(1), blocks_refined_list(2:params%number_procs)
            endif
            close(99)
        endif
    endif

    ! synchronize light data
    call synchronize_lgt_data( params, refinement_status_only=.false. )

end subroutine refinement_execute_tree



!> \brief Refine mesh. All cpu loop over their heavy data and check if the refinement
!! flag +1 is set on the block. If so, we either take this block, interpolate it to the next finer
!! level and create four(2D) or eight(3D) new blocks, each carrying a part of the interpolated data,
!! or copy them to Mallat decomposed blocks if they are on the input level
!! As all CPU first work individually, the light data array is synced afterwards.
!> This function is used when parts of the blocks are already waveled decomposed in adapt_tree,
!! so for block on the given level it refines directly to mallat wavelet decomposition (and then spaghetti),
!! for other blocks it does normal refinement
!
!> \note The interpolation (or prediction) operator here is applied to a block INCLUDING
!! any ghost nodes. You must sync first. This does not apply to blocks which will be in MallatWD though
!
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays \n
! ********************************************************************************************
subroutine refinementExecute_lvl2SpaghettiWD( params, hvy_block, tree_ID, level )

    implicit none

    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)             !> heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID                              !> Tree to refine
    integer(kind=ik), intent(in)        :: level                                !> Level, where blocks will be refined to MallatWD

    integer(kind=ik)                    :: k, N, dF, k_daughter                 ! loop variables
    integer(kind=ik)                    :: rank                                 ! process rank
    integer(kind=ik)                    :: g                                    ! grid parameter
    integer(kind=ik), dimension(3)      :: Bs
    integer(kind=ik)                    :: lgt_free_id, free_hvy_id, lgt_id 
    integer(kind=tsize)                 :: treecode
    integer(kind=ik)                    :: level_me, hvyID

    N = params%number_blocks
    rank = params%rank
    Bs = params%Bs
    g  = params%g

    ! every proc loop over its active heavy data array
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k, tree_ID)

        ! light data id
        call hvy2lgt( lgt_id, hvyID, rank, N )

        ! block wants to refine
        if ( (lgt_block( lgt_id, idx_refine_sts) == +1) ) then
            level_me = lgt_block( lgt_id, IDX_MESH_LVL)
            
            ! Daughters created as Mallat WDed blocks
            if (level_me == level) then
                call refineBlock2SpaghettiWD(params, hvy_block, hvyID, tree_ID)
            ! Daughters refined normally
            else
                call refineBlock(params, hvy_block, hvyID, tree_ID)
            endif

        end if
    end do

    ! synchronize light data
    call synchronize_lgt_data( params, refinement_status_only=.false. )

end subroutine refinementExecute_lvl2SpaghettiWD
