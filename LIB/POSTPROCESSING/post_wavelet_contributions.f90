subroutine post_wavelet_contributions(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData
    use module_wavelets

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: fname_in, fname_out
    real(kind=rk)                      :: time
    integer(kind=ik)                   :: iteration, k, lgt_ID, tc_length, g
    integer(kind=ik), dimension(3)     :: Bs

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvy_ID, ix, iy, iz, even_odd, idim, hvy_tmp_neqn
    integer(kind=ik)                   :: level_me, target_level, Jmax_active, Jmin_active
    integer(kind=ik)                   :: sub0_levels, sub_level, max_sub0_levels
    integer(kind=ik), external         :: count_div2
    integer(kind=ik)                   :: pos_underscore
    logical                            :: norms_only
    logical                            :: use_box
    logical                            :: logical_has_extent, logical_has_end
    real(kind=rk)                      :: norms(1:7)
    real(kind=rk), dimension(3)        :: box_origin, box_extent, box_end
    real(kind=rk), dimension(3)        :: box_end_cmd

    real(kind=rk), dimension(3)        :: domain

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename)
    call get_command_argument(2, fname_in)

    ! does the user need help?
    if (fname_in=='--help' .or. fname_in=='--h' .or. fname_in=='-h') then
        if (params%rank==0) then
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " Wabbit postprocessing: compute wavelet contributions per level"
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " This tool decomposes the input field using wavelets, then"
            write(*,'(A)') " reconstructs the field keeping only wavelet coefficients"
            write(*,'(A)') " from each level separately. Outputs one file per level."
            write(*,'(A)') ""
            write(*,'(A)') " ./wabbit-post --wavelet-contributions ux_000.h5 --wavelet=CDF40"
            write(*,'(A)') ""
            write(*,'(A)') " Output: ux-WC0_000.h5, ux-WC1_000.h5, ..., ux-WCJ_000.h5"
            write(*,'(A)') " With sub0 enabled: ux-WC0s1_000.h5, ux-WC0s2_000.h5, ..."
            write(*,'(A)') ""
            write(*,'(A)') " To only compute norms per level (no files saved):"
            write(*,'(A)') " ./wabbit-post --wavelet-contributions ux_000.h5 --wavelet=CDF40 --norms-only"
            write(*,'(A)') ""
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " --wavelet=CDF44   - Wavelet type (default: CDF40)"
            write(*,'(A)') " --Jmin=N           - Minimum contribution level to process (default: 0)"
            write(*,'(A)') " --norms-only       - Only compute and print norms, do not save files"
            write(*,'(A)') " --sub0-levels=N    - Additional WC-only periodic sub0 contributions"
            write(*,'(A)') "                     (requires active Jmin=0, max by repeated /2 of Bs)"
            write(*,'(A)') " --box-origin-x/y/z - Box origin for local contributions"
            write(*,'(A)') " --box-extent-x/y/z - Box extent for local contributions"
            write(*,'(A)') " --box-end-x/y/z    - Box end point for local contributions"
            write(*,'(A)') "                     (outside coefficients are set to zero, slice-wise)"
            write(*,'(A)') "                     Per dimension: use either --box-extent OR --box-end"
            write(*,'(A)') "                     For each dimension, specify all 3 (x,y,z) consistently"
            write(*,'(A)') "-----------------------------------------------------------"
        end if
        return
    endif

    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF40" )
    call get_cmd_arg( "--Jmin", params%Jmin, default=0_ik )
    call get_cmd_arg_bool( "--norms-only", norms_only, default=.false. )
    call get_cmd_arg( "--sub0-levels", sub0_levels, default=0_ik )
    call get_cmd_arg( "--box-origin-x", box_origin(1), default=0.0_rk )
    call get_cmd_arg( "--box-origin-y", box_origin(2), default=0.0_rk )
    call get_cmd_arg( "--box-origin-z", box_origin(3), default=0.0_rk )
    call get_cmd_arg( "--box-extent-x", box_extent(1), default=-1.0_rk )
    call get_cmd_arg( "--box-extent-y", box_extent(2), default=-1.0_rk )
    call get_cmd_arg( "--box-extent-z", box_extent(3), default=-1.0_rk )
    call get_cmd_arg( "--box-end-x", box_end_cmd(1), default=-1.0_rk )
    call get_cmd_arg( "--box-end-y", box_end_cmd(2), default=-1.0_rk )
    call get_cmd_arg( "--box-end-z", box_end_cmd(3), default=-1.0_rk )

    ! initialize wavelet transform
    ! also, set number of ghost nodes params%G to minimal value for this wavelet
    call setup_wavelet(params, params%g)

    params%forest_size = 1
    params%n_eqn = 1

    call check_file_exists(trim(fname_in))

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(fname_in, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)
    
    ! Check box parameters: for each dimension, use either extent or end (not both)
    ! but allow mixing across dimensions
    do idim = 1, params%dim
        logical_has_extent = box_extent(idim) > 0.0_rk
        logical_has_end = box_end_cmd(idim) > 0.0_rk
        
        if (logical_has_extent .and. logical_has_end) then
            write(*,'(A,I0,A)') "ERROR: Cannot specify both --box-extent-", idim, &
                " and --box-end-", idim
            call abort(260107_ik, "Cannot use both --box-extent and --box-end for same dimension")
        endif
        
        ! If end point is specified, validate it's greater than origin
        if (logical_has_end .and. box_end_cmd(idim) <= box_origin(idim)) then
            write(*,'(A,I0,A,F16.8,A,F16.8)') "ERROR: box-end-", idim, " (", box_end_cmd(idim), &
                ") must be > box-origin-", idim, " (", box_origin(idim), ")"
            call abort(260109_ik, "Box end point must be greater than box origin")
        endif
        
        ! Convert extent to end point if extent was used
        if (logical_has_extent) then
            box_end(idim) = box_origin(idim) + box_extent(idim)
        else if (logical_has_end) then
            box_end(idim) = box_end_cmd(idim)
        endif
    enddo
    
    ! Check if all dimension has box parameters set
    use_box = all((box_extent(1:params%dim) > 0.0_rk) .or. (box_end_cmd(1:params%dim) > 0.0_rk))

    params%Jmax = tc_length
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs
    allocate(params%butcher_tableau(1,1))

    params%useCoarseExtension = params%isLiftedWavelet
    params%useSecurityZone = params%isLiftedWavelet

    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component = 1

    Bs = params%Bs
    g  = params%g

    ! for full wavelet operations we need 8/7 for 3D or 4/3 for 2D as blocks
    params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) * 2.0_rk**params%dim / (2.0_rk**params%dim - 1.0_rk) * 1.2_rk)+4

    ! for norms-only, we can skip reloading, as we do not need to prune
    params%n_eqn = 1
    if (norms_only) then
        hvy_tmp_neqn = 2
    else
        hvy_tmp_neqn = 1
    endif

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, neqn_hvy_tmp=hvy_tmp_neqn)

    ! read input data
    call readHDF5vct_tree( (/ fname_in /), params, hvy_block, tree_ID)

    call updateMetadata_tree(params, tree_ID, search_overlapping=.true.)

    ! perform wavelet decomposition
    if (params%rank == 0) then
        write(*,'(A)') "Performing wavelet decomposition..."
    endif
    call wavelet_decompose_full_tree(params, hvy_block, tree_ID, hvy_tmp)

    ! determine active level range
    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)

    if (params%Jmin < 0) then
        call abort(260100_ik, "--Jmin must be >= 0")
    endif
    if (params%Jmin > Jmax_active) then
        call abort(260104_ik, "--Jmin exceeds max active level in input")
    endif

    if (params%rank == 0) then
        write(*,'(A,I0,A,I0)') "Active levels: ", params%Jmin, " to ", Jmax_active
        write(*,'(A)') "Computing contributions for each level ..."
    endif

    if (sub0_levels < 0) then
        call abort(260101_ik, "--sub0-levels must be >= 0")
    endif
    if (sub0_levels > 0) then
        if (params%Jmin /= 0) then
            call abort(260105_ik, "--sub0-levels requires --Jmin=0")
        endif
        if (Jmin_active /= 0) then
            call abort(260102_ik, "sub0 requested but active Jmin is not 0")
        endif

        max_sub0_levels = huge(1_ik)
        max_sub0_levels = min(max_sub0_levels, count_div2(Bs(1)))
        max_sub0_levels = min(max_sub0_levels, count_div2(Bs(2)))
        if (params%dim == 3) max_sub0_levels = min(max_sub0_levels, count_div2(Bs(3)))

        if (sub0_levels > max_sub0_levels) then
            write(*,'(A,I0,A,I0)') "Requested --sub0-levels=", sub0_levels, " but max allowed is ", max_sub0_levels
            call abort(260103_ik, "sub0-levels exceeds divide-by-2 limit")
        endif
    endif

    ! even_odd pattern for identifying WC vs SC
    even_odd = mod(params%g + 1, 2)

    ! append ending before last underscore or if it doesnt exist before last .
    pos_underscore = index(fname_in, "_", back=.true.)
    if (pos_underscore < 0) pos_underscore = index(fname_in, '.', back=.true.)

    ! loop over all active levels
    do target_level = params%Jmin, Jmax_active
        if (params%rank == 0) then
            write(*,'(A,I0,A,I0)') "Processing level ", target_level, " / ", Jmax_active
        endif

        ! fresh start for each target contribution
        call reload_and_decompose(reload_from_file=.not. norms_only .or. target_level == params%Jmin)

        ! zero out wavelet coefficients on all levels except target_level
        ! currently we intentionally keep WC only (SC are discarded for all levels)
        ! future extension could add an additional SC-only output branch
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
            level_me = lgt_block(lgt_ID, IDX_MESH_LVL)

            if (level_me /= target_level) then
                ! zero out all values on this block, no matter if SC or WC
                hvy_block(:,:,:,1,hvy_ID) = 0.0_rk
            ! else if (target_level > Jmin_active .and. target_level > params%Jmin) then
            else
                ! on the target level, keep WC but zero SC
                do iz = merge(1, g+1, params%dim==2), merge(1, Bs(3)+g, params%dim==2)
                    do iy = g+1, Bs(2)+g
                        do ix = g+1, Bs(1)+g
                            ! check if this is a scaling coefficient point
                            if (mod(ix,2)==even_odd .and. mod(iy,2)==even_odd .and. &
                                mod(iz,2)==merge(1,even_odd,params%dim==2)) then
                                ! this is a SC point - zero it out
                                hvy_block(ix,iy,iz,1,hvy_ID) = 0.0_rk
                            endif
                        enddo
                    enddo
                enddo
            endif
        enddo

        if (target_level == 0 .and. sub0_levels > 0) then
            ! First, save the regular L0 contribution (WC0), then add sub-level distributions.
            call finalize_and_output(target_level, -1_ik)

            do sub_level = 1, sub0_levels
                ! fresh start for each sub-level contribution
                call reload_and_decompose(reload_from_file=.not. norms_only)

                ! keep only level-0 blocks before sub0 processing
                do k = 1, hvy_n(tree_ID)
                    hvy_ID = hvy_active(k, tree_ID)
                    call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
                    level_me = lgt_block(lgt_ID, IDX_MESH_LVL)

                    if (level_me /= 0) then
                        hvy_block(:,:,:,1,hvy_ID) = 0.0_rk
                    endif
                enddo

                ! hierarchical sub-level extraction:
                ! SC0 -> (SC-1,WC-1) -> ... -> (SC-sub_level,WC-sub_level)
                ! while existing coarser WC are kept untouched by sub-level decomposition.
                do k = 1, hvy_n(tree_ID)
                    hvy_ID = hvy_active(k, tree_ID)
                    call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
                    level_me = lgt_block(lgt_ID, IDX_MESH_LVL)

                    if (level_me == 0) then
                        ! input is L0-decomposed, we further decompose it until we get the desired sub-level decomposition
                        call sub0_wavelet_decompose_n(params, hvy_block(:,:,:,1:1,hvy_ID), 0_ik, sub_level, hvy_tmp(:,:,:,1:1,hvy_ID))

                        ! keep only WC-sub_level:
                        ! remove WC0..WC(sub_level-1), then remove SC-sub_level
                        do iy = 0, sub_level-1
                            call sub0_delete_coefficients(params, hvy_block(:,:,:,1:1,hvy_ID), iy, "WC")
                        enddo
                        call sub0_delete_coefficients(params, hvy_block(:,:,:,1:1,hvy_ID), sub_level, "SC")

                        ! reconstruction from sub-level to level-0 decomposition is deferred to finalize_and_output,
                        ! so it can happen after optional local box filtering
                    endif
                enddo

                call finalize_and_output(0_ik, sub_level)

            enddo

            cycle
        endif

        call finalize_and_output(target_level, -1_ik)
    enddo

    if (params%rank == 0) then
        write(*,'(A)') "Done! Created contributions for all levels."
    endif

contains

    subroutine finalize_and_output(level_id, sublevel_id)
        implicit none

        integer(kind=ik), intent(in) :: level_id
        integer(kind=ik), intent(in) :: sublevel_id
        integer(kind=ik) :: n_nonzero_wavelets, mpierr, g(3)

        g(:) = params%g
        if (params%dim == 2) g(3) = 0

        if (use_box) call apply_local_box_filter()

        ! ensure reconstruction helper data matches currently filtered coefficients
        ! also count non-zero wavelet coefficients for the current contribution (for output, not used in reconstruction)
        n_nonzero_wavelets = 0_ik
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            hvy_tmp(:,:,:,1,hvy_ID) = hvy_block(:,:,:,1,hvy_ID)

            n_nonzero_wavelets = n_nonzero_wavelets + count(abs(hvy_block(g(1)+1:g(1)+params%bs(1),g(2)+1:g(2)+params%bs(2),g(3)+1:g(3)+params%bs(3),:,hvy_ID)) > 1.0e-12_rk)
        enddo
        call MPI_ALLREDUCE(MPI_IN_PLACE, n_nonzero_wavelets, 1, MPI_INTEGER4, MPI_SUM, WABBIT_COMM, mpierr)

        !!!
        ! Reconstruction starts here. We have to reconstruct from sub0 to 0 (if necessary) and then do normal reconstruction.
        !!!

        ! For sub0 outputs, reconstruct from isolated WC-sub_level to level-0 decomposition
        if (sublevel_id > 0_ik) then
            do k = 1, hvy_n(tree_ID)
                hvy_ID = hvy_active(k, tree_ID)
                call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
                level_me = lgt_block(lgt_ID, IDX_MESH_LVL)

                if (level_me == 0) then
                    call sub0_wavelet_reconstruct_n(params, hvy_block(:,:,:,1:1,hvy_ID), sublevel_id, 0_ik, hvy_tmp(:,:,:,1:1,hvy_ID))
                endif
            enddo
        endif

        call wavelet_reconstruct_full_tree(params, hvy_block, hvy_tmp, tree_ID)

        call componentWiseNorm_tree(params, hvy_block, tree_ID, "L1",     norms(1:1), threshold_state_vector=.false., norm_case="leaf")
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2",     norms(2:2), threshold_state_vector=.false., norm_case="leaf")
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "L3",     norms(3:3), threshold_state_vector=.false., norm_case="leaf")
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "L4",     norms(4:4), threshold_state_vector=.false., norm_case="leaf")
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "Linfty", norms(5:5), threshold_state_vector=.false., norm_case="leaf")
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "H1",     norms(6:6), threshold_state_vector=.false., norm_case="leaf")
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "Mean",     norms(7:7), threshold_state_vector=.false., norm_case="leaf")

        if (params%rank == 0) then
            if (sublevel_id < 0) then
                write(*,'(A,I0,A,I0)') "  Level ", level_id, " / ", Jmax_active
            else
                write(*,'(A,I0,A,I0)') "  Level 0 sub-level ", sublevel_id, " / ", sub0_levels
            endif
            write(*,'(A,es16.8)') "    L1 norm:     ", norms(1)
            write(*,'(A,es16.8)') "    L2 norm:     ", norms(2)
            write(*,'(A,es16.8)') "    L3 norm:     ", norms(3)
            write(*,'(A,es16.8)') "    L4 norm:     ", norms(4)
            write(*,'(A,es16.8)') "    LInfty norm: ", norms(5)
            write(*,'(A,es16.8)') "    H1 norm:     ", norms(6)
            write(*,'(A,es16.8)') "    Mean:     ", norms(7)
            write(*,'(A,I0)')     "    Non-zero wavelets: ", n_nonzero_wavelets
        endif
        if (.not. norms_only) then
            if (sublevel_id < 0) then
                write(fname_out, '(A, A, I0, A)') fname_in(1:pos_underscore-1), "-WC", level_id, fname_in(pos_underscore:LEN_TRIM(fname_in))
            else
                write(fname_out, '(A, A, I0, A)') fname_in(1:pos_underscore-1), "-WC0s", sublevel_id, fname_in(pos_underscore:LEN_TRIM(fname_in))
            endif

            call prune_fulltree2leafs(params, tree_ID)
            call saveHDF5_tree(fname_out, time, iteration, 1, params, hvy_block, tree_ID )
            if (params%rank == 0) then
                write(*,'(A,A)') "  Saved: ", trim(fname_out)
            endif
        endif

    end subroutine

    subroutine reload_and_decompose(reload_from_file)
        implicit none
        logical, intent(in) :: reload_from_file

        integer(kind=ik) :: k_block, hvy_ID_block

        if (reload_from_file) then
            call delete_tree(params, tree_ID)
            call createActiveSortedLists_forest(params)
            call readHDF5vct_tree( (/ fname_in /), params, hvy_block, tree_ID)
            call updateMetadata_tree(params, tree_ID, search_overlapping=.true.)
            call wavelet_decompose_full_tree(params, hvy_block, tree_ID, hvy_tmp)

            ! let's save the wavelet decomposed values in hvy_tmp second entry
            if (size(hvy_tmp,4) > 1) then
                do k_block = 1, hvy_n(tree_ID)
                    hvy_ID_block = hvy_active(k_block, tree_ID)
                    hvy_tmp(:,:,:,2,hvy_ID_block) = hvy_block(:,:,:,1,hvy_ID_block)
                enddo
            endif
        else
            if (size(hvy_tmp,4) > 1) then
                do k_block = 1, hvy_n(tree_ID)
                    hvy_ID_block = hvy_active(k_block, tree_ID)
                    hvy_block(:,:,:,1,hvy_ID_block) = hvy_tmp(:,:,:,2,hvy_ID_block)
                enddo
            else
                call abort(260511, "Something is odd. I do not have what you are searching for.")
            endif
        endif
    end subroutine

    subroutine apply_local_box_filter()
        implicit none

        integer(kind=ik) :: k_loc, hvy_id_loc, lgt_id_loc
        integer(kind=ik) :: ix_loc, iy_loc, iz_loc
        real(kind=rk) :: x0_loc(3), dx_loc(3)
        real(kind=rk) :: x, y, z
        real(kind=rk) :: xmin, xmax, ymin, ymax, zmin, zmax
        logical :: outside_block

        do k_loc = 1, hvy_n(tree_ID)
            hvy_id_loc = hvy_active(k_loc, tree_ID)
            call hvy2lgt(lgt_id_loc, hvy_id_loc, params%rank, params%number_blocks)
            call get_block_spacing_origin(params, lgt_id_loc, x0_loc, dx_loc)

            xmin = x0_loc(1)
            xmax = x0_loc(1) + dble(Bs(1)-1) * dx_loc(1)
            ymin = x0_loc(2)
            ymax = x0_loc(2) + dble(Bs(2)-1) * dx_loc(2)
            zmin = x0_loc(3)
            zmax = x0_loc(3)
            if (params%dim == 3) zmax = x0_loc(3) + dble(Bs(3)-1) * dx_loc(3)

            outside_block = (xmax < box_origin(1) .or. xmin > box_end(1) .or. &
                             ymax < box_origin(2) .or. ymin > box_end(2))
            if (params%dim == 3) outside_block = outside_block .or. (zmax < box_origin(3) .or. zmin > box_end(3))

            if (outside_block) then
                hvy_block(:,:,:, 1, hvy_id_loc) = 0.0_rk
                cycle
            endif

            do iz_loc = merge(1, g+1, params%dim==2), merge(1, Bs(3)+g, params%dim==2)
                z = 0
                if (params%dim == 3) z = dble(iz_loc-(g+1)) * dx_loc(3) + x0_loc(3)

                if (params%dim == 3 .and. (z < box_origin(3) .or. z > box_end(3))) then
                    hvy_block(:, :, iz_loc, 1, hvy_id_loc) = 0.0_rk
                    cycle
                endif

                do iy_loc = g+1, Bs(2)+g
                    y = dble(iy_loc-(g+1)) * dx_loc(2) + x0_loc(2)

                    if (y < box_origin(2) .or. y > box_end(2)) then
                        hvy_block(:, iy_loc, iz_loc, 1, hvy_id_loc) = 0.0_rk
                        cycle
                    endif

                    do ix_loc = g+1, Bs(1)+g
                        x = dble(ix_loc-(g+1)) * dx_loc(1) + x0_loc(1)
                        if (x < box_origin(1) .or. x > box_end(1)) then
                            hvy_block(ix_loc, iy_loc, iz_loc, 1, hvy_id_loc) = 0.0_rk
                        endif
                    enddo
                enddo
            enddo
        enddo
    end subroutine

end subroutine


integer(kind=ik) function count_div2(n)
    use module_globals
    implicit none
    integer(kind=ik), intent(in) :: n
    integer(kind=ik) :: m

    count_div2 = 0
    m = n
    do while (m > 0 .and. modulo(m, 2_ik) == 0)
        count_div2 = count_div2 + 1
        m = m / 2
    enddo
end function
