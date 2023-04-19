subroutine coarseExtensionUpdate_tree( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, hvy_n, lgt_n, &
    inputDataSynced )
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

    integer(kind=ik) :: N, k, neighborhood, level_diff, hvyID, lgtID, hvyID_neighbor, lgtID_neighbor, level_me, level_neighbor, Nwcl
    integer(kind=ik) :: nx,ny,nz,nc, g, Bs(1:3), Nwcr, ii, Nscl, Nscr, Nreconl, Nreconr, nnn
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

    if ((params%wavelet=="CDF40").or.(params%wavelet=="CDF20")) return

    nx = size(hvy_block, 1)
    ny = size(hvy_block, 2)
    nz = size(hvy_block, 3)
    nc = size(hvy_block, 4)
    g  = params%g
    Bs = params%bs
    Nscl    = params%Nscl
    Nscr    = params%Nscr
    Nwcl    = params%Nwcl
    Nwcr    = params%Nwcr
    Nreconl = params%Nreconl
    Nreconr = params%Nreconr

    ! if (.not. areArraysSameSize(sc, hvy_block(:,:,:,:)))

    if (allocated(tmp_reconst)) then
        if (.not. areArraysSameSize(hvy_block(:,:,:,:,1), tmp_reconst) ) deallocate(tmp_reconst)
    endif
    if (allocated(wc)) then
        if (.not. areArraysSameSize(hvy_block(:,:,:,:,1), wc(:,:,:,:,1)) ) deallocate(wc)
    endif
    if (.not. allocated(wc)) allocate(wc(1:nx, 1:ny, 1:nz, 1:nc, 1:8) )
    if (.not. allocated(tmp_reconst)) allocate(tmp_reconst(1:nx, 1:ny, 1:nz, 1:nc) )

    !---------------------------------------------------------------------------
    t0 = MPI_Wtime()

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
                level_me       = lgt_block( lgtID, params%Jmax + IDX_MESH_LVL )
                level_neighbor = lgt_block( lgtID_neighbor, params%Jmax + IDX_MESH_LVL )

                if (level_neighbor < level_me) then
                    toBeManipulated(k) = .true.
                    ! its enough if one is true
                    ! nnn = nnn + 1
                    exit
                endif
            endif
        enddo
    end do
    call toc( "coarseExtension (toBeManipulated list)", MPI_Wtime()-t0 )
    !---------------------------------------------------------------------------
    ! write(*,*) "rank", params%rank, "Nblocksforreon", nnn, hvy_n, lgt_n


    ! First, we sync the ghost nodes, in order to apply the decomposition filters
    ! HD and GD to the data. It may be that this is done before calling.
    if (.not. inputDataSynced) then
        t0 = MPI_Wtime()
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
        call toc( "coarseExtension (sync 1)", MPI_Wtime()-t0 )
    endif


    ! 2nd. We compute the decomposition by applying the HD GD filters to the data.
    ! Data are decimated and then stored in the block in Spaghetti order (the counterpart
    ! to Mallat ordering). Coefficients SC and WC are computed only inside the block
    ! not in the ghost nodes (that is trivial: we cannot compute the filters in the ghost node
    ! layer). The first point (g+1),(g+1) is a scaling function coefficient, regardless of g
    ! odd or even. Bs must be even. We also make a copy of the sync'ed data n hvy_work - we use this
    ! to fix up the SC in the coarse extension case.
    t0 = MPI_Wtime()
    do k = 1, hvy_n
        hvyID = hvy_active(k)

        if (toBeManipulated(k)) then
            ! hvy_work now is a copy with sync'ed ghost points.
            hvy_work(:,:,:,1:nc,hvyID) = hvy_block(:,:,:,1:nc,hvyID)
            ! data WC/SC now in Spaghetti order
            call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
        endif
    end do
    call toc( "coarseExtension (FWT)", MPI_Wtime()-t0 )


    ! 3rd we sync the decompose coefficients, but only on the same level. This is
    ! required as several fine blocks hit a coarse one, and this means we have to sync
    ! those in order to be able to reconstruct with modified WC/SC. Note it does not matter
    ! if all blocks are sync'ed: we'll not reconstruct on the coarse block anyways, and the
    ! ghost nodes on the fine bloc (WC/SC) are overwritten in the coarse-extension assumption anyways.
    t0 = MPI_Wtime()
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, &
    hvy_active, hvy_n, syncSameLevelOnly1=.true. )
    ! Note we tested it and syncSameLevelOnly1=.true. is indeed slightly faster (compared to full sync) 
    call toc( "coarseExtension (sync 2)", MPI_Wtime()-t0 )


    ! routine operates on a single tree (not a forest)
    ! 4th. Loop over all blocks and check if they have *coarser* neighbors. Neighborhood 1..4 are
    ! always equidistant.
    t0 = MPI_Wtime()
    do k = 1, hvy_n
        ! this block just remains the same as before..
        if (.not. toBeManipulated(k)) cycle

        ! ===> the following code is only executed for blocks that are manipulated by CE
        hvyID = hvy_active(k)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )

        ! copy transform data in "inflated mallat ordering":
        call spaghetti2inflatedMallat_block(params, hvy_block(:,:,:,:,hvyID), wc)

        ! loop over all relevant neighbors
        do neighborhood = 1, size(hvy_neighbor, 2)
            ! neighbor exists ?
            if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                ! neighbor light data id
                lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                level_me       = lgt_block( lgtID, params%Jmax + IDX_MESH_LVL )
                level_neighbor = lgt_block( lgtID_neighbor, params%Jmax + IDX_MESH_LVL )

                if (level_neighbor < level_me) then
                    ! manipulation of coeffs
                    call coarseExtensionManipulateWC_block(params, wc, neighborhood)
                    call coarseExtensionManipulateSC_block(params, wc, hvy_work(:,:,:,:,hvyID), neighborhood)
                endif
            endif
        enddo

        ! copy back original data, then fix the coarse extension parts by
        ! copying tmp_reconst (the inverse of the manipulated SC/WC) into the patches
        ! relevant for coarse extension
        hvy_block(:,:,:,1:nc,hvyID) = hvy_work(:,:,:,1:nc,hvyID)

        ! reconstruct from the manipulated coefficients
        call mallat2spaghetti_block(params, wc, tmp_reconst)
        call waveletReconstruction_block(params, tmp_reconst)

        ! reconstruction part. We manipulated the data and reconstructed it in some regions.
        ! Now, we copy those reconstructed data back to the original block - this is
        ! the actual coarseExtension.
        do neighborhood = 1, size(hvy_neighbor,2)
            if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                ! neighbor light data id
                lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                level_me       = lgt_block( lgtID, params%Jmax + IDX_MESH_LVL )
                level_neighbor = lgt_block( lgtID_neighbor, params%Jmax + IDX_MESH_LVL )

                if (level_neighbor < level_me) then
                    ! coarse extension case (neighbor is coarser)
                    if (params%dim == 2) then
                        select case (neighborhood)
                        case (9:10)
                            ! -x
                            hvy_block(1:Nreconl, :, :, 1:nc, hvyID) = tmp_reconst(1:Nreconl,:,:,1:nc)
                        case (11:12)
                            ! +x
                            hvy_block(nx-Nreconr:nx, :, :, 1:nc, hvyID) = tmp_reconst(nx-Nreconr:nx,:,:,1:nc)
                        case (13:14)
                            ! +y
                            hvy_block(:, ny-Nreconr:ny, :, 1:nc, hvyID) = tmp_reconst(:, ny-Nreconr:ny,:,1:nc)
                        case (15:16)
                            ! -y
                            hvy_block(:, 1:Nreconl, :, 1:nc, hvyID) = tmp_reconst(:, 1:Nreconl,:,1:nc)
                        case (5)
                            hvy_block(1:Nreconl, ny-Nreconr:ny, :, 1:nc, hvyID) = tmp_reconst(1:Nreconl, ny-Nreconr:ny,:,1:nc)
                        case (6)
                            hvy_block(1:Nreconl, 1:Nreconl, :, 1:nc, hvyID) = tmp_reconst(1:Nreconl, 1:Nreconl,:,1:nc)
                        case (7)
                            ! top right corner
                            hvy_block(nx-Nreconr:nx, ny-Nreconr:ny, :, 1:nc, hvyID) = tmp_reconst(nx-Nreconr:nx, ny-Nreconr:ny,:,1:nc)
                        case (8)
                            hvy_block(nx-Nreconr:nx, 1:Nreconl, :, 1:nc, hvyID) = tmp_reconst(nx-Nreconr:nx, 1:Nreconl,:,1:nc)
                        end select
                    else
                        select case(neighborhood)
                            ! ---faces---
                        case (35:38)
                            ! +x
                            hvy_block(nx-Nscr:nx, :, :, 1:nc, hvyID) = tmp_reconst(nx-Nscr:nx, :, :, 1:nc)
                        case (43:46)
                            ! -x
                            hvy_block(1:Nscl, :, :, 1:nc, hvyID) = tmp_reconst(1:Nscl, :, :, 1:nc)
                        case (39:42)
                            ! +y
                            hvy_block(:, ny-Nscr:ny, :, 1:nc, hvyID) = tmp_reconst(:, ny-Nscr:ny, :, 1:nc)
                        case (31:34)
                            ! -y
                            hvy_block(:, 1:Nscl, :, 1:nc, hvyID) = tmp_reconst(:, 1:Nscl, :, 1:nc)
                        case (27:30)
                            ! +z
                            hvy_block(:, :, nz-Nscr:nz, 1:nc, hvyID) = tmp_reconst(:, :, nz-Nscr:nz, 1:nc)
                        case (47:50)
                            ! -z
                            hvy_block(:, :, 1:Nscl, 1:nc, hvyID) = tmp_reconst(:, :, 1:Nscl, 1:nc)
                            ! --- corners ---
                        case (26)
                            hvy_block(1:Nscl, 1:Nscl, 1:Nscl, 1:nc, hvyID) = tmp_reconst(1:Nscl, 1:Nscl, 1:Nscl, 1:nc)
                        case (23)
                            hvy_block(nx-Nscr:nx, 1:Nscl, 1:Nscl, 1:nc, hvyID) = tmp_reconst(nx-Nscr:nx, 1:Nscl, 1:Nscl, 1:nc)
                        case (22)
                            hvy_block(1:Nscl, 1:Nscl, nz-Nscr:nz, 1:nc, hvyID) = tmp_reconst(1:Nscl, 1:Nscl, nz-Nscr:nz, 1:nc)
                        case (19)
                            hvy_block(nx-Nscr:nx, 1:Nscl, nz-Nscr:nz, 1:nc, hvyID) = tmp_reconst(nx-Nscr:nx, 1:Nscl, nz-Nscr:nz, 1:nc)
                        case (25)
                            hvy_block(1:Nscl, ny-Nscr:ny, 1:Nscl, 1:nc, hvyID) = tmp_reconst(1:Nscl, ny-Nscr:ny, 1:Nscl, 1:nc)
                        case (21)
                            hvy_block(1:Nscl, ny-Nscr:ny, nz-Nscr:nz, 1:nc, hvyID) = tmp_reconst(1:Nscl, ny-Nscr:ny, nz-Nscr:nz, 1:nc)
                        case (20)
                            hvy_block(nx-Nscr:nx, ny-Nscr:ny, nz-Nscr:nz, 1:nc, hvyID) = tmp_reconst(nx-Nscr:nx, ny-Nscr:ny, nz-Nscr:nz, 1:nc)
                        case (24)
                            hvy_block(nx-Nscr:nx, ny-Nscr:ny, 1:Nscl, 1:nc, hvyID) = tmp_reconst(nx-Nscr:nx, ny-Nscr:ny, 1:Nscl, 1:nc)
                            ! ---(partial) edges---
                        case (51:52)
                            hvy_block(:, 1:Nscl, nz-Nscr:nz, 1:nc, hvyID) = tmp_reconst(:, 1:Nscl, nz-Nscr:nz, 1:nc)
                        case (53:54)
                            hvy_block(nx-Nscr:nx, :, nz-Nscr:nz, 1:nc, hvyID) = tmp_reconst(nx-Nscr:nx, :, nz-Nscr:nz, 1:nc)
                        case (55:56)
                            hvy_block(:, ny-Nscr:ny, nz-Nscr:nz, 1:nc, hvyID) = tmp_reconst(:, ny-Nscr:ny, nz-Nscr:nz, 1:nc)
                        case (57:58)
                            hvy_block(1:Nscl, :, nz-Nscr:nz, 1:nc, hvyID) = tmp_reconst(1:Nscl, :, nz-Nscr:nz, 1:nc)
                        case (59:60)
                            hvy_block(:, 1:Nscl, 1:Nscl, 1:nc, hvyID) = tmp_reconst(:, 1:Nscl, 1:Nscl, 1:nc)
                        case (61:62)
                            hvy_block(nx-Nscr:nx, :, 1:Nscl, 1:nc, hvyID) = tmp_reconst(nx-Nscr:nx, :, 1:Nscl, 1:nc)
                        case (63:64)
                            hvy_block(:, ny-Nscr:ny, 1:Nscl, 1:nc, hvyID) = tmp_reconst(:, ny-Nscr:ny, 1:Nscl, 1:nc)
                        case (65:66)
                            hvy_block(1:Nscl, :, 1:Nscl, 1:nc, hvyID) = tmp_reconst(1:Nscl, :, 1:Nscl, 1:nc)
                        case (67:68)
                            hvy_block(nx-Nscr:nx, 1:Nscl, :, 1:nc, hvyID) = tmp_reconst(nx-Nscr:nx, 1:Nscl, :, 1:nc)
                        case (69:70)
                            hvy_block(1:Nscl, 1:Nscl, :, 1:nc, hvyID) = tmp_reconst(1:Nscl, 1:Nscl, :, 1:nc)
                        case (71:72)
                            hvy_block(nx-Nscr:nx, ny-Nscr:ny, :, 1:nc, hvyID) = tmp_reconst(nx-Nscr:nx, ny-Nscr:ny, :, 1:nc)
                        case (73:74)
                            hvy_block(1:Nscl, ny-Nscr:ny, :, 1:nc, hvyID) = tmp_reconst(1:Nscl, ny-Nscr:ny, :, 1:nc)
                        end select
                    endif
                endif
            endif
        enddo
    end do
    call toc( "coarseExtension (manipulation loop)", MPI_Wtime()-t0 )
    ! unfortunately, the above loop affects the load balancing. in the sync_ghosts
    ! step, CPUS will be in sync again, but since they arrive at different times at this line of
    ! code, idling occurs -> bad for performance.

    ! final sync step. The reason is: we have kept some SC and zeroed some WC, and
    ! reconstructed the signal. The ghost nodes are outdated now: the reconstruction
    ! step altered the signal.
    t0 = MPI_Wtime()
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
    call toc( "coarseExtension (sync 3)", MPI_Wtime()-t0 )
end subroutine
