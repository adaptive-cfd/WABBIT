subroutine substitution_step( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, hvy_n )
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
    integer(kind=ik), intent(in)        :: hvy_n

    integer(kind=ik) :: N, k, neighborhood, level_diff, hvyID, lgtID, hvyID_neighbor, lgtID_neighbor, level_me, level_neighbor, Nwcl
    integer(kind=ik) :: nx,ny,nz,nc, g, Bs(1:3), Nwcr, ii, Nscl, Nscr, Nreconl, Nreconr
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: sc, wcx, wcy, wcxy, tmp_reconst
    logical :: manipulated

    if ((params%wavelet=="CDF40").or.(params%wavelet=="CDF20")) return

    nx = size(hvy_block, 1)
    ny = size(hvy_block, 2)
    nz = size(hvy_block, 3)
    nc = size(hvy_block, 4)
    g = params%g
    Bs = params%bs
    Nscl    = params%Nscl
    Nscr    = params%Nscr
    Nwcl    = params%Nwcl
    Nwcr    = params%Nwcr
    Nreconl = params%Nreconl
    Nreconr = params%Nreconr

    if (.not. allocated(sc  )) allocate(  sc(1:nx, 1:ny, 1:nz, 1:nc) )
    if (.not. allocated(wcx )) allocate( wcx(1:nx, 1:ny, 1:nz, 1:nc) )
    if (.not. allocated(wcy )) allocate( wcy(1:nx, 1:ny, 1:nz, 1:nc) )
    if (.not. allocated(wcxy)) allocate(wcxy(1:nx, 1:ny, 1:nz, 1:nc) )
    if (.not. allocated(tmp_reconst)) allocate(tmp_reconst(1:nx, 1:ny, 1:nz, 1:nc) )

    ! maybe this call is redundant - FIXME
    ! First, we sync the ghost nodes, in order to apply the decomposition filters
    ! HD and GD to the data. It may be that this is done before calling.
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )


    ! 2nd. We compute the decomposition by applying the HD GD filters to the data.
    ! Data are decimated and then stored in the block in Spaghetti order (the counterpart
    ! to Mallat ordering). Coefficients SC and WC are computed only inside the block
    ! not in the ghost nodes (that is trivial: we cannot compute the filters in the ghost node
    ! layer). The first point (g+1),(g+1) is a scaling function coefficient, regardless of g
    ! odd or even. Bs must be even. We also make a copy of the sync'ed data n hvy_work - we use this
    ! to fix up the SC in the coarse extension case.
    do k = 1, hvy_n
        hvyID = hvy_active(k)
        ! hvy_work now is a copy with sync'ed ghost points.
        hvy_work(:,:,:,1:nc,hvyID) = hvy_block(:,:,:,1:nc,hvyID)

        ! data WC/SC now in Spaghetti order
        call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
    end do


    ! 3rd we sync the decompose coefficients, but only on the same level. This is
    ! required as several fine blocks hit a coarse one, and this means we have to sync
    ! those in order to be able to reconstruct with modified WC/SC. Note it does not matter
    ! if all blocks are sync'ed: we'll not reconstruct on the coarse block anyways, and the
    ! ghost nodes on the fine bloc (WC/SC) are overwritten in the coarse-extension assumption anyways.
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, &
    hvy_active, hvy_n, syncSameLevelOnly1=.true. )

    ! routine operates on a single tree (not a forest)
    ! 4th. Loop over all blocks and check if they have *coarser* neighbors. Neighborhood 1..4 are
    ! always equidistant.
    do k = 1, hvy_n
        hvyID = hvy_active(k)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )

        call spaghetti2mallat_block(params, hvy_block(:,:,:,:,hvyID), sc, wcx, wcy, wcxy)
        manipulated = .false.

        ! loop over all relevant neighbors
        do neighborhood = 5, 16
            ! neighbor exists ?
            if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                ! neighbor light data id
                lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                level_me       = lgt_block( lgtID, params%Jmax + IDX_MESH_LVL )
                level_neighbor = lgt_block( lgtID_neighbor, params%Jmax + IDX_MESH_LVL )

                if (level_neighbor < level_me) then
                    manipulated = .true.
                    ! manipulation of coeffs
                    call coarseExtensionManipulateWC_block(params, wcx, wcy, wcxy, neighborhood)
                    call coarseExtensionManipulateSC_block(params, sc, hvy_work(:,:,:,:,hvyID), neighborhood)
                endif
            endif
        enddo

        ! copy back original data, then fix the coarse extension parts by
        ! copying tmp_reconst (the inverse of the manipulated SC/WC) into the patches
        ! relevant for coarse extension
        hvy_block(:,:,:,1:nc,hvyID) = hvy_work(:,:,:,1:nc,hvyID)

        if (manipulated) then
            ! ensures that 3/4 of the numbers are zero - required for reconstruction
            ! note when copying Spaghetti to Mallat, this is automatically done, but
            ! when manipulating coefficients, it may happen that we set nonzero values
            ! where a zero should be. Here: only SC (WC are set to zero anyways)
            call setRequiredZerosWCSC_block(params, sc)
            ! wavelet reconstruction - we do not call the routine in module_interplation
            ! because this requires us to copy data back to Spaghetti-ordering (which is
            ! unnecessary here, even though it would not hurt)
            call blockFilterCustom_vct( params, sc  , sc  , "HR", "HR", "--" )
            call blockFilterCustom_vct( params, wcx , wcx , "HR", "GR", "--" )
            call blockFilterCustom_vct( params, wcy , wcy , "GR", "HR", "--" )
            call blockFilterCustom_vct( params, wcxy, wcxy, "GR", "GR", "--" )
            tmp_reconst = sc + wcx + wcy + wcxy

            ! reconstruction part. We manipulated the data and reconstructed it in some regions.
            ! Now, we copy those reconstructed data back to the original block - this is
            ! the actual coarseExtension.
            do neighborhood = 5, 16
                if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                    ! neighbor light data id
                    lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                    level_me       = lgt_block( lgtID, params%Jmax + IDX_MESH_LVL )
                    level_neighbor = lgt_block( lgtID_neighbor, params%Jmax + IDX_MESH_LVL )

                    if (level_neighbor < level_me) then
                        ! coarse extension case (neighbor is coarser)
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
                    endif
                endif
            enddo
        end if
    end do

    ! this is probably an optional sync step
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

end subroutine
