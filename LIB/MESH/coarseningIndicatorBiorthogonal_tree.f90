subroutine coarseningIndicatorBiorthogonal_tree( time, params, level_this, hvy_block, hvy_tmp, tree_ID)

    use module_indicators

    implicit none
    real(kind=rk), intent(in)           :: time
    type (type_params), intent(in)      :: params                         !> user defined parameter structure
    integer(kind=ik), intent(in)        :: level_this                     !> current level to look at (in the case of biorthogonal wavelets, not in "harten-multiresolution")
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)       !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)         !> heavy work data array - block data.
    integer(kind=ik), intent(in)        :: tree_ID

    ! local variables
    integer(kind=ik) :: k, Jmax, neq, lgtID, g, mpierr, hvyID, p, N_thresholding_components, &
    tags, ierr, level, neighborhood, level_me, level_neighbor, lgtID_neighbor
    integer(kind=ik), dimension(3) :: Bs
    ! local block spacing and origin
    real(kind=rk) :: dx(1:3), x0(1:3), R
    real(kind=rk), allocatable :: norm(:), wc_mag(:)
    integer(kind=ik) :: nx,ny,nz,nc, Nwcl, Nwcr, Nscl, Nscr, Nreconl, Nreconr
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: sc, wcx, wcy, wcxy

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    ! in the default case we threshold all statevector components
    Jmax = params%Jmax
    neq = params%n_eqn
    Bs = params%Bs
    g = params%g

    !> reset refinement status to "stay" on all blocks
    do k = 1, lgt_n(tree_ID)
        lgtID = lgt_active(k, tree_ID)
        lgt_block( lgtID, Jmax + IDX_REFINE_STS ) = 0
    enddo


    !---------------------------------------------------------------------------
    !> Compute normalization for eps, if desired.
    !---------------------------------------------------------------------------
    !! versions <30.05.2018 used fixed eps for all qtys of the state vector, but that is not very smart
    !! as each qty can have different mangitudes. If the switch eps_normalized is on, we compute here
    !! the vector of normalization factors for each qty that adaptivity will be based on (state vector
    !! or vorticity). The nor is specified in params%eps_norm, default is Linfty.

    !! default norm (e.g. for compressible navier-stokes) is 1 so in this case eps
    !! is an absolute value.
    allocate(norm(1:neq))
    allocate(wc_mag(1:neq))
    norm = 1.0_rk

    ! if we coarsen randomly or everywhere, well, why compute the norm ?
    if ( params%eps_normalized ) then
        ! Apply thresholding directly to the statevector (hvy_block), not to derived quantities
        call component_wise_tree_norm(params, hvy_block, tree_ID, params%eps_norm, norm)

        ! avoid division by zero (corresponds to using an absolute eps if the norm is very small)
        do p = 1, neq
            if (norm(p) <= 1.0e-9_rk) norm(p) = 1.0_rk
        enddo

        ! during dev is it useful to know what the normalization is, if that is active
        call append_t_file('eps_norm.t', (/time, norm, params%eps/))
    endif

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

    ! maybe this call is redundant - FIXME
    ! First, we sync the ghost nodes, in order to apply the decomposition filters
    ! HD and GD to the data. It may be that this is done before calling.
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )


    ! 2nd. We compute the decomposition by applying the HD GD filters to the data.
    ! Data are decimated and then stored in the block in Spaghetti order (the counterpart
    ! to Mallat ordering). Coefficients SC and WC are computed only inside the block
    ! not in the ghost nodes (that is trivial: we cannot compute the filters in the ghost node
    ! layer). The first point (g+1),(g+1) is a scaling function coefficient, regardless of g
    ! odd or even. Bs must be even. We also make a copy of the sync'ed data n hvy_work - we use this
    ! to fix up the SC in the coarse extension case.
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        ! data WC/SC now in Spaghetti order
        call waveletDecomposition_block(params, hvy_tmp(:,:,:,1:nc,hvyID))
    end do


    ! 3rd we sync the decompose coefficients, but only on the same level. This is
    ! required as several fine blocks hit a coarse one, and this means we have to sync
    ! those in order to be able to reconstruct with modified WC/SC. Note it does not matter
    ! if all blocks are sync'ed: we'll not reconstruct on the coarse block anyways, and the
    ! ghost nodes on the fine bloc (WC/SC) are overwritten in the coarse-extension assumption anyways.
    call sync_ghosts( params, lgt_block, hvy_tmp, hvy_neighbor, &
    hvy_active(:,tree_ID), hvy_n(tree_ID), syncSameLevelOnly1=.true. )

    ! routine operates on a single tree (not a forest)
    ! 4th. Loop over all blocks and check if they have *coarser* neighbors. Neighborhood 1..4 are
    ! always equidistant.
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k, tree_ID)
        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )

        level_me = lgt_block( lgtID, params%Jmax + IDX_MESH_LVL )

        !!!!!!!!!!!!!!!!!
        if (level_me /= level_this) cycle
        !!!!!!!!!!!!!!!!!

        call spaghetti2mallat_block(params, hvy_tmp(:,:,:,1:nc,hvyID), sc, wcx, wcy, wcxy)

        ! loop over all relevant neighbors
        do neighborhood = 5, 16
            ! neighbor exists ?
            if ( hvy_neighbor(hvyID, neighborhood) /= -1 ) then
                ! neighbor light data id
                lgtID_neighbor = hvy_neighbor( hvyID, neighborhood )
                level_neighbor = lgt_block( lgtID_neighbor, params%Jmax + IDX_MESH_LVL )

                if (level_neighbor < level_me) then
                    ! manipulation of coeffs
                    call coarseExtensionManipulateWC_block(params, wcx, wcy, wcxy, neighborhood)
                endif
            endif
        enddo

        do p = 1, nc
            ! if all details are smaller than C_eps, we can coarsen.
            ! check interior WC only
            wc_mag(p) = maxval(abs(wcx(g+1:Bs(1)+g,g+1:Bs(2)+g,:,p)))
            wc_mag(p) = max( wc_mag(p), maxval(abs(wcy(g+1:Bs(1)+g,g+1:Bs(2)+g,:,p))) )
            wc_mag(p) = max( wc_mag(p), maxval(abs(wcxy(g+1:Bs(1)+g,g+1:Bs(2)+g,:,p))) )

            wc_mag(p) = wc_mag(p) / norm(p)
        enddo

        if (maxval(wc_mag)< params%eps) then
            lgt_block(lgtID, params%Jmax+IDX_REFINE_STS) = -1_ik
        endif

    enddo

    ! force blocks on maximum refinement level to coarsen, if parameter is set
    if (params%force_maxlevel_dealiasing) then
        do k = 1, lgt_n(tree_ID)
            lgtID = lgt_active(k, tree_ID)
            if (lgt_block(lgtID, Jmax + IDX_MESH_LVL) == params%Jmax) then
                ! force blocks on maxlevel to coarsen
                lgt_block(lgtID, Jmax + IDX_REFINE_STS) = -1
            endif
        enddo
    endif

    deallocate(norm)


    !> after modifying all refinement flags, we need to synchronize light data
    call synchronize_lgt_data( params,  refinement_status_only=.true. )

end subroutine
