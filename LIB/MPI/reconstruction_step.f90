subroutine coarseExtensionUpdate_level( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, hvy_n, lgt_n, &
    inputDataSynced, level)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n, lgt_n
    ! if the input data are sync'ed we do not do it here: otherwise, call
    ! ghost nodes synchronization
    logical, intent(in) :: inputDataSynced
    integer(kind=ik), intent(in) :: level

! real(kind=rk), allocatable :: WCtmp(:,:,:,:,:,:) ! code used to verify that FWT after manip yields same coeffs
    integer(kind=ik) :: N, k, neighborhood, level_diff, hvyID, lgtID, hvyID_neighbor, lgtID_neighbor, level_me, level_neighbor
    integer(kind=ik) :: nx,ny,nz,nc, g, Bs(1:3), ii, Nreconl, Nreconr, nnn, p, ierr, g_this, g_spaghetti, Nwcl,Nwcr, idx(2,3)
!    integer(kind=ik) :: ix,iy,iz,ic,iwc ! code used to verify that FWT after manip yields same coeffs

    ! The WC array contains SC (scaling function coeffs) as well as all WC (wavelet coeffs)
    ! Note: the precise naming of SC/WC is not really important. we just apply
    ! the correct decomposition/reconstruction filters - thats it.
    !
    ! INDEX            2D     3D     LABEL
    ! -----            --    ---     ---------------------------------
    ! wc(:,:,:,:,1)    HH    HHH     sc scaling function coeffs
    ! wc(:,:,:,:,2)    HG    HGH     wcx wavelet coeffs
    ! wc(:,:,:,:,3)    GH    GHH     wcy wavelet coeffs
    ! wc(:,:,:,:,4)    GG    GGH     wcxy wavelet coeffs
    ! wc(:,:,:,:,5)          HHG     wcz wavelet coeffs
    ! wc(:,:,:,:,6)          HGG     wcxz wavelet coeffs
    ! wc(:,:,:,:,7)          GHG     wcyz wavelet coeffs
    ! wc(:,:,:,:,8)          GGG     wcxyz wavelet coeffs
    !
    real(kind=rk), allocatable, dimension(:,:,:,:,:), save :: wc
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: tmp_reconst
    real(kind=rk) :: t0
    logical, allocatable, save :: toBeManipulated(:)
    character(len=80) :: fname

    if (.not. params%isLiftedWavelet) return

    t0 = MPI_Wtime()
    nx = size(hvy_block, 1)
    ny = size(hvy_block, 2)
    nz = size(hvy_block, 3)
    nc = size(hvy_block, 4)
    g  = params%g
    Bs = params%bs
    Nreconl = params%Nreconl
    Nreconr = params%Nreconr 

    if ((Bs(1) < Nreconr-g).or.(Bs(2) < Nreconr-g)) then
        ! NOTE: Nreconr > Nreconl always.
        write(*,*) params%wavelet, "Bs=", Bs, "Bs_min=", Nreconr-g
        call abort(991234, "For the chosen wavelet, block size is too small!")
    endif

    ! it turns out, when the coefficients are spaghetti-ordered,
    ! we can sync one less point and still have enough coefficients.
    select case(params%wavelet)
    case("CDF44")
        g_spaghetti = 6
    case("CDF42")
        g_spaghetti = 4
    case("CDF22")
        g_spaghetti = 2
    case default
        ! full sync of all points
        g_spaghetti = g
    end select

    ! HUGE ! code used to verify that FWT after manip yields same coeffs
    ! allocate(WCtmp(1:nx, 1:ny, 1:nz, 1:nc, 1:8, 1:size(hvy_block,5))) ! code used to verify that FWT after manip yields same coeffs

    if (allocated(tmp_reconst)) then
        if (.not. areArraysSameSize(hvy_block(:,:,:,:,1), tmp_reconst) ) deallocate(tmp_reconst)
    endif
    if (allocated(wc)) then
        if (.not. areArraysSameSize(hvy_block(:,:,:,:,1), wc(:,:,:,:,1)) ) deallocate(wc)
    endif
    if (.not. allocated(wc)) allocate(wc(1:nx, 1:ny, 1:nz, 1:nc, 1:8) )
    if (.not. allocated(tmp_reconst)) allocate(tmp_reconst(1:nx, 1:ny, 1:nz, 1:nc) )

    !---------------------------------------------------------------------------
    ! create the list of blocks that will be affected by coarseExtension.
    ! It is a heavy array (distributed). This avoids checking again and again
    ! and improves code readability - performance-wise it does not matter
    if (allocated(toBeManipulated)) then
        if (size(toBeManipulated,1) < hvy_n) deallocate(toBeManipulated)
    endif
    if (.not. allocated(toBeManipulated)) allocate(toBeManipulated(1:hvy_n))

    ! default no block is to be manipulated
    toBeManipulated(1:hvy_n) = .false.

    ! check for all blocks if they have coarser neighbors - in this case they'll be manipulated
    ! nnn = 0
    do k = 1, hvy_n
        hvyID = hvy_active(k)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )

        do neighborhood = 1, size(hvy_neighbor, 2)
            ! neighbor exists ?
            if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                ! neighbor light data id
                lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                level_me       = lgt_block( lgtID, IDX_MESH_LVL )
                level_neighbor = lgt_block( lgtID_neighbor, IDX_MESH_LVL )

                ! we proceed level-wise
                if ((level_neighbor < level_me).and.(level_me==level)) then
                    toBeManipulated(k) = .true.
                    ! nnn = nnn + 1
                    ! its enough if one neighborhood is true
                    exit
                endif
            endif
        enddo
    end do
    call toc( "coarseExtension 1 (toBeManipulated list)", 1010, MPI_Wtime()-t0 )
    ! write(*,*) "rank", params%rank, "N", nnn, hvy_n, lgt_n, "level=", level
    !---------------------------------------------------------------------------


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 1st, we sync the ghost nodes, in order to apply the decomposition filters
    ! HD and GD to the data. It may be that this sync'ing is done before calling.
    if (.not. inputDataSynced) then
        t0 = MPI_Wtime()
        g_this = max(ubound(params%HD,1), ubound(params%GD,1))! ubound GD is the largest (GD is not symmetric but HD is)
        call sync_level_with_all_neighbours( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, &
            level, g_minus=g_this, g_plus=g_this)
        call toc( "coarseExtension (sync 1)", 1011, MPI_Wtime()-t0 )
    endif


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 2nd. We compute the FWT decomposition by applying the HD GD filters to the data on level "level".
    ! Data are filtered, decimated and then stored in the block in Spaghetti order (the counterpart
    ! to Mallat ordering). Coefficients SC and WC are computed only inside the block
    ! not in the ghost nodes (that is trivial: we cannot compute the filters in the ghost node
    ! layer). The first point (g+1),(g+1) is a scaling function coefficient, regardless of g
    ! odd or even. Bs must be even. We also make a copy of the sync'ed data n hvy_work - we use this
    ! to fix up the SC in the coarse extension case. 
    ! NOTE: we indeed need to tranform more blocks than we'll actually modify with coarseExt. Strictly speaking,
    ! FWT is required for blocks 
    ! (1) are affected by coarse Ext
    ! (2) blocks that are neighbors to blocks from (1)
    ! because we need to sync their WC. However, this logic is not implemented, and we FWT all blocks on the level.
    t0 = MPI_Wtime()
    do k = 1, hvy_n
        hvyID = hvy_active(k)

        ! We compute detail coefficients on the fly here, for all blocks
        ! on the level.
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
        level_me = lgt_block( lgtID, IDX_MESH_LVL )

        ! FWT required for a block that is on the level
        if (level_me == level) then
            ! hvy_work now is a copy with sync'ed ghost points.
            hvy_work(:,:,:,1:nc,hvyID) = hvy_block(:,:,:,1:nc,hvyID)

            ! data WC/SC now in Spaghetti order
            call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
        endif
    end do
    call toc( "coarseExtension 2 (FWT)", 1012, MPI_Wtime()-t0 )


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 3rd we sync the WC coefficients, but only on the same level. This is
    ! required as several fine blocks hit a coarse one, and this means we have to sync
    ! those in order to be able to reconstruct with modified WC/SC. Note it does not matter
    ! if all blocks are sync'ed: we'll not reconstruct on the coarse block anyways, and the
    ! ghost nodes on the fine bloc (WC/SC) are overwritten in the coarse-extension assumption anyways.
    t0 = MPI_Wtime()
    call sync_level_only( params, lgt_block, hvy_block, hvy_neighbor, &
    hvy_active, hvy_n, level, g_minus=g_spaghetti, g_plus=g_spaghetti)
    ! Note we tested it and syncSameLevelOnly1=.true. is indeed slightly faster (compared to full sync)
    call toc( "coarseExtension 3 (sync 2)", 1013, MPI_Wtime()-t0 )


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 3B) : Not all blocks will be modified by coarseExtension, but are still FWT'ed.
    ! We restore those blocks now (by copying back the backup and avoiding the expensive IWT).
    ! As we already have the FWT here, we might as well compute the largest WC (detail for coarseningIndicator)
    do k = 1, hvy_n
        hvyID = hvy_active(k)

        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
        level_me = lgt_block( lgtID, IDX_MESH_LVL )
        
        if (.not. toBeManipulated(k) .and. (level_me==level)) then
            ! this block has been FWT'ed but is not modified, howver, we can extract its details now
            call spaghetti2inflatedMallat_block(params, hvy_block(:,:,:,:,hvyID), wc)

            ! ! extract largest wavelet coeffcienct
            ! do p = 1, nc
            !     if (params%dim==3) then
            !         hvy_details(p, hvyID) = maxval( abs(wc(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, p, 2:8)) )
            !     else
            !         hvy_details(p, hvyID) = maxval( abs(wc(g+1:Bs(1)+g, g+1:Bs(2)+g, :, p, 2:4)) )
            !     endif
            ! enddo

            ! restore original data (this block is not modified, but currently transformed to wavelet space)
            hvy_block(:,:,:,1:nc,hvyID) = hvy_work(:,:,:,1:nc,hvyID)
        endif
    enddo   


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 4th. Loop over all blocks and check if they have *coarser* neighbors
    ! routine operates on a single tree (not a forest)
    t0 = MPI_Wtime()
    do k = 1, hvy_n
        hvyID = hvy_active(k)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )

        ! this block just remains the same as before..
        if (.not. toBeManipulated(k)) cycle
        
        ! ===> the following code is only executed for blocks that are manipulated by CE

        ! copy transformed data to "inflated mallat ordering":
        call spaghetti2inflatedMallat_block(params, hvy_block(:,:,:,:,hvyID), wc)

        ! loop over all relevant neighbors
        do neighborhood = 1, size(hvy_neighbor, 2)
            ! neighbor exists ?
            if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                ! neighbor light data id
                lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                level_me       = lgt_block( lgtID, IDX_MESH_LVL )
                level_neighbor = lgt_block( lgtID_neighbor, IDX_MESH_LVL )

                if (level_neighbor < level_me) then
                    ! manipulation of coeffs
                    call coarseExtensionManipulateWC_block(params, wc, neighborhood)
                    call coarseExtensionManipulateSC_block(params, wc, hvy_work(:,:,:,:,hvyID), neighborhood)
                elseif (level_neighbor > level_me) then
                    ! it is actually possible for a block to have both finer and coarser neighbors. If its
                    ! coarser, (level_neighbor<level_me), then coarseExtension is applied to this block.
                    ! However for small Bs (or large wavelets, i.e. large Nreconr Nreconl), the coarseExt
                    ! reconstruction step also takes into account the ghost nodes layer on the adjacent 
                    ! side. If this adjacent block is then finer (level_neighbor > level_me), those WC
                    ! should be zero and this is what we assure here.
                    call coarseExtensionManipulateWC_block(params, wc, neighborhood, params%g, params%g)
                    ! note if the neighboring blocks are on the same level, the WC are sync'ed between the 
                    ! blocks, and we must NOT delete them.
                    
                    ! What about the SC?
                    ! See the extensive comment in WaveDecomposition_dim1 on that. 
                    ! "Magically", the SC are correct.
                endif
            endif
        enddo

        ! WCtmp(:,:,:,:,1:8,hvyID) = wc ! code used to verify that FWT after manip yields same coeffs

        ! ! evaluate detail for blocks that were affected by coarseExtension (blocks on
        ! ! the level J which are not affected by it are computed above)
        ! do p = 1, nc
        !     if (params%dim==3) then
        !         hvy_details(p, hvyID) = maxval( abs(wc(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, p, 2:8)) )
        !     else
        !         hvy_details(p, hvyID) = maxval( abs(wc(g+1:Bs(1)+g, g+1:Bs(2)+g, :, p, 2:4)) )
        !     endif
        ! enddo

        ! copy back original data, then fix the coarse extension parts by
        ! copying tmp_reconst (the inverse of the manipulated SC/WC) into the patches
        ! relevant for coarse extension
        hvy_block(:,:,:,1:nc,hvyID) = hvy_work(:,:,:,1:nc,hvyID)

        ! reconstruct from the manipulated coefficients
        call inflatedMallat2spaghetti_block(params, wc, tmp_reconst)
        call waveletReconstruction_block(params, tmp_reconst)

        ! reconstruction part. We manipulated the data and reconstructed them on the entire block with modified coeffs.
        ! Now, we copy those reconstructed data back to the original block - this is
        ! the actual coarseExtension.
        ! JB: In theory we synched the SC and WC of the same level so we could overwrite the whole domain
        !     however, as with blocks with both finer and coarser neighbours the finer neighbours do not synch WC and SC
        !     the WR reconstructs false values at the border there
        !     For CVS this problem should solve itself as every block can find same-level neighbours for all directions
        do neighborhood = 1, size(hvy_neighbor,2)
            if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                ! neighbor light data id
                lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                level_me       = lgt_block( lgtID, IDX_MESH_LVL )
                level_neighbor = lgt_block( lgtID_neighbor, IDX_MESH_LVL )

                if (level_neighbor < level_me) then
                    ! coarse extension case (neighbor is coarser)
                    idx(:, :) = 1
                    call get_indices_of_modify_patch(params, neighborhood, idx, (/ nx, ny, nz/), (/Nreconl, Nreconl, Nreconl/), (/Nreconr, Nreconr, Nreconr/))

                    hvy_block(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc, hvyID) = &
                        tmp_reconst(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc)
                endif
            endif
        enddo
    end do
    call toc( "coarseExtension 4 (manipulation loop)", 1014, MPI_Wtime()-t0 )

    ! unfortunately, the above loop affects the load balancing. in the sync_ghosts
    ! step, CPUS will be in sync again, but since they arrive at different times at this line of
    ! code, idling occurs -> bad for performance.

    ! final sync step. The reason is: we have kept some SC and zeroed some WC, and
    ! reconstructed the signal. The ghost nodes are outdated now: the reconstruction
    ! step altered the signal.
    t0 = MPI_Wtime()
    g_this = max(ubound(params%HD,1),ubound(params%GD,1))
    call sync_level_with_all_neighbours( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, level, g_minus=g_this, g_plus=g_this)
    call toc( "coarseExtension 5 (sync 3)", 1015, MPI_Wtime()-t0 )


    ! code used to verify that FWT after manip yields same coeffs
    ! do k = 1, hvy_n
    !     if (.not. toBeManipulated(k)) cycle
    !
    !     hvyID = hvy_active(k)
    !     call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
    !
    !     hvy_work(:,:,:,1:nc,hvyID) = hvy_block(:,:,:,1:nc,hvyID)
    !     call waveletDecomposition_block(params, hvy_work(:,:,:,:,hvyID))
    !
    !     call spaghetti2inflatedMallat_block(params, hvy_work(:,:,:,:,hvyID), wc)
    !
    !     wc = wc - WCtmp(:,:,:,:,1:8,hvyID)
    !
    ! write(fname,'(i1,"_",i0.6,".diff")') params%rank, hvyID
    ! write(*,*) fname
    !
    !     open(14, file=fname, status='replace')
    !     do ix = g+1, Bs(1)+g!nx
    !     do iy = g+1, Bs(1)+g!ny
    !     do iz = g+1, Bs(1)+g!nz
    !     do iwc = 1, 8
    !         write(14,'(es15.8))') wc(ix,iy,iz,1,iwc)
    !     enddo
    !     enddo
    !     enddo
    !     enddo
    !     close(14)
    !     abort(197)
    ! enddo
    !
    ! deallocate(WCtmp)
end subroutine



!> \brief Modify the SC and WC of a wavelet decomposed blocks at fine/coarse interfaces
!> This routine assumes that the input is already wavelet decomposed in spaghetti form
!> It will work on all blocks or which have the REF_TMP_UNTREATED flag in refinement status for leaf-wise operation
subroutine coarse_extension_modify_tree(params, lgt_block, hvy_data, hvy_tmp, hvy_neighbor, hvy_active, hvy_n, lgt_n, tree_ID, sc_skip_ghosts)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)             !< light data array
    real(kind=rk), intent(inout)        :: hvy_data(:, :, :, :, :)     !< heavy data array, WDed in Spaghetti form
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)      !< heavy work data array - block data.
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)           !< heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_active(:)               !< list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n, lgt_n                !< number of active blocks (heavy and light data)
    integer(kind=ik), intent(in)        :: tree_ID                     !< Tree to be investigated
    logical, optional, intent(in)       :: sc_skip_ghosts              !< for second CE modify we can skip ghosts points

    integer(kind=ik)                    :: iteration, k, neighborhood, lgtID, hvyID, tree_me
    integer(kind=ik)                    :: nx,ny,nz,nc, level_me, level_neighbor, lgtID_neighbor
    logical                             :: toBeManipulated, scSkipGhosts
    real(kind=rk), allocatable, dimension(:,:,:,:,:), save :: wc

    nx = size(hvy_data, 1)
    ny = size(hvy_data, 2)
    nz = size(hvy_data, 3)
    nc = size(hvy_data, 4)

    scSkipGhosts = .false.
    if (present(sc_skip_ghosts)) scSkipGhosts = sc_skip_ghosts

    if (allocated(wc)) then
        if (size(wc, 4) < nc) deallocate(wc)
    endif
    if (.not. allocated(wc)) allocate(wc(1:nx, 1:ny, 1:nz, 1:nc, 1:8) )

    do k = 1, hvy_n
        ! Set toBeManipulated - This is one of those sneaky errors I searched 10hours for - JB
        toBeManipulated = .false.

        hvyID = hvy_active(k)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
        level_me       = lgt_block( lgtID, IDX_MESH_LVL )
        tree_me        = lgt_block( lgtID, IDX_TREE_ID )

        ! check if this block is to be modified
        do neighborhood = 1, size(hvy_neighbor, 2)
            ! neighbor exists ?
            if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                ! neighbor light data id
                lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                level_neighbor = lgt_block( lgtID_neighbor, IDX_MESH_LVL )

                ! we proceed level-wise
                if ((level_neighbor < level_me) .and. (tree_me==tree_ID)) then
                    toBeManipulated = .true.
                    ! nnn = nnn + 1
                    ! its enough if one neighborhood is true
                    exit
                endif
            endif
        enddo

        if (toBeManipulated) then
            ! transform to inflated mallat
            call spaghetti2inflatedMallat_block(params, hvy_data(:,:,:,:,hvyID), wc)

            ! loop over all relevant neighbors
            do neighborhood = 1, size(hvy_neighbor, 2)
                ! neighbor exists ?
                if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                    ! neighbor light data id
                    lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                    level_neighbor = lgt_block( lgtID_neighbor, IDX_MESH_LVL )


                    if (level_neighbor < level_me) then
                        ! manipulation of coeffs
                        call coarseExtensionManipulateWC_block(params, wc, neighborhood)
                        call coarseExtensionManipulateSC_block(params, wc, hvy_tmp(:,:,:,:,hvyID), neighborhood, scSkipGhosts)
                    elseif (level_neighbor > level_me) then
                        ! it is actually possible for a block to have both finer and coarser neighbors. If its
                        ! coarser, (level_neighbor<level_me), then coarseExtension is applied to this block.
                        ! However for small Bs (or large wavelets, i.e. large Nreconr Nreconl), the coarseExt
                        ! reconstruction step also takes into account the ghost nodes layer on the adjacent 
                        ! side. If this adjacent block is then finer (level_neighbor > level_me), those WC
                        ! should be zero and this is what we assure here, as they have not been correctly computed
                        call coarseExtensionManipulateWC_block(params, wc, neighborhood, params%g, params%g)
                        ! note if the neighboring blocks are on the same level, the WC are sync'ed between the 
                        ! blocks, and we must NOT delete them.
                        
                        ! What about the SC?
                        ! See the extensive comment in WaveDecomposition_dim1 on that. 
                        ! "Magically", the SC are correct.
                    endif
                endif
            enddo

            ! transform wc from inflated mallat back to spaghetti
            call inflatedMallat2spaghetti_block(params, wc, hvy_data(:,:,:,:,hvyID))
        endif
    enddo
end subroutine



!> \brief Apply CE for all blocks in a tree. This copies back the old values and then overwrites those on affected patches
subroutine coarse_extension_reconstruct_tree(params, lgt_block, hvy_data, hvy_tmp, hvy_neighbor, hvy_active, hvy_n, lgt_n)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)             !< light data array
    real(kind=rk), intent(inout)        :: hvy_data(:, :, :, :, :)     !< heavy data array
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)      !< heavy work data array - block data.
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)           !< heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_active(:)               !< list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n, lgt_n                !< number of active blocks (heavy and light data)

    integer(kind=ik)                    :: iteration, k, neighborhood, lgtID, hvyID, Nreconl, Nreconr
    integer(kind=ik)                    :: nx,ny,nz,nc, level_me, level_neighbor, lgtID_neighbor, idx(2,3)
    logical                             :: toBeManipulated
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: tmp_reconst

    integer(kind=ik) :: iy

    nx = size(hvy_data, 1)
    ny = size(hvy_data, 2)
    nz = size(hvy_data, 3)
    nc = size(hvy_data, 4)

    Nreconl = params%Nreconl
    Nreconr = params%Nreconr 

    if (allocated(tmp_reconst)) then
        if (size(tmp_reconst, 4) < nc) deallocate(tmp_reconst)
    endif
    if (.not. allocated(tmp_reconst)) allocate(tmp_reconst(1:nx, 1:ny, 1:nz, 1:nc) )

    do k = 1, hvy_n
        ! Set toBeManipulated - This is one of those sneaky errors I searched 10hours for - JB
        toBeManipulated = .false.

        hvyID = hvy_active(k)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
        level_me       = lgt_block( lgtID, IDX_MESH_LVL )

        ! check if this block is to be modified
        do neighborhood = 1, size(hvy_neighbor, 2)
            ! neighbor exists ?
            if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                ! neighbor light data id
                lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                level_neighbor = lgt_block( lgtID_neighbor, IDX_MESH_LVL )

                ! check if this patch is a CE patch
                if ((level_neighbor < level_me)) then
                    toBeManipulated = .true.
                    exit
                endif
            endif
        enddo

        if (toBeManipulated) then
        
            ! reconstruct from the manipulated coefficients
            tmp_reconst = hvy_data(:,:,:,1:nc,hvyID)
            call waveletReconstruction_block(params, tmp_reconst)

            ! if (lgt_block(lgtID, IDX_TC_2) == 4032 .and. lgt_block(lgtID, IDX_MESH_LVL)==5) then
            !     write(*, '("113 lvl ", i0)') level
            !     do iy = 1, 19
            !         write(*, '(19(es8.1))') tmp_reconst(1:19, iy, 1, 1)
            !     enddo
            ! endif

            hvy_data(:,:,:,1:nc,hvyID) = hvy_tmp(:,:,:,1:nc,hvyID)

            ! reconstruction part. We manipulated the data and reconstructed them on the entire block with modified coeffs.
            ! Now, we copy those reconstructed data back to the original block - this is
            ! the actual coarseExtension.
            ! JB: In theory we synched the SC and WC of the same level so we could overwrite the whole domain
            !     however, as with blocks with both finer and coarser neighbours the finer neighbours do not synch WC and SC
            !     the WR reconstructs false values at the border there
            !     For CVS this problem should solve itself as every block can find same-level neighbours for all directions
            do neighborhood = 1, size(hvy_neighbor,2)
                if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                    ! neighbor light data id
                    lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                    level_neighbor = lgt_block( lgtID_neighbor, IDX_MESH_LVL )

                    if (level_neighbor < level_me) then
                        ! coarse extension case (neighbor is coarser)
                        idx(:, :) = 1
                        call get_indices_of_modify_patch(params, neighborhood, idx, (/ nx, ny, nz/), (/Nreconl, Nreconl, Nreconl/), (/Nreconr, Nreconr, Nreconr/))

                        hvy_data(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc, hvyID) = &
                            tmp_reconst(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), 1:nc)
                    endif
                endif
            enddo

            ! if (lgt_block(lgtID, IDX_TC_2) == 4032 .and. lgt_block(lgtID, IDX_MESH_LVL)==5) then
            !     write(*, '("113 data lvl ", i0)') level
            !     do iy = 1, 19
            !         write(*, '(19(es8.1))') hvy_data(1:19, iy, 1, 1, hvyID)
            !     enddo
            ! endif

        else  ! block is not modified, rewrite old values
            hvy_data(:,:,:,1:nc,hvyID) = hvy_tmp(:,:,:,1:nc,hvyID)
        endif
    enddo
end subroutine