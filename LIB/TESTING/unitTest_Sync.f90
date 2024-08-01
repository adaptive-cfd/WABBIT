subroutine unitTest_fill_linearly(params, hvy_block, hvy_work, hvy_tmp, tree_ID)
    implicit none
    type (type_params), intent(inout)       :: params                     !> user defined parameter structure
    real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(out)              :: hvy_tmp(:, :, :, :, :)
    !> heavy work array: used for RHS evaluation in multistep methods (like RK4: u0, k1, k2 etc)
    real(kind=rk), intent(out)              :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik), intent(in)            :: tree_ID

    integer(kind=ik)                        :: k, lgt_id, hvy_id, Jm_a
    real(kind=rk)                           :: ddx(1:3), xx0(1:3)
    integer(kind=ik)                        :: ix, iy, iz, g, Bs(3), g_min
    real(kind=rk)                           :: x, y, z

    Bs = params%Bs
    g  = params%g
    Jm_a = maxActiveLevel_tree(tree_ID)

    ! Set initial condition being a function in x,y,z: f(x,y,z) = 4*Bs(1)*x + 8*Bs(2)*y + 16*Bs(3)*z
    ! only set it at interioor of the blocks
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        ! compute block spacing and origin from treecode
        call get_block_spacing_origin( params, lgt_id, xx0, ddx )

        ! calculate f(x,y,z) for first datafield
        do iz = 1, size(hvy_block, 3)
            if (params%dim == 2) then
                z = 0
            else
                z = real(iz-(g+1), kind=rk) * ddx(3) + xx0(3)
            endif
            do iy = 1, size(hvy_block, 2)
                y = real(iy-(g+1), kind=rk) * ddx(2) + xx0(2)
                do ix = 1, size(hvy_block, 1)
                    x = real(ix-(g+1), kind=rk) * ddx(1) + xx0(1)
                    hvy_block(ix, iy, iz, 1, hvy_id) = 2**(Jm_a)*Bs(1)*x + 2**(Jm_a+1)*Bs(2)*y + 2**(Jm_a+2)*Bs(3)*z
                enddo
            enddo
        enddo

        ! now the entire block (incl ghost nodes) holds the exact solution: make a
        ! copy of the block for later comparison, but use work arrays usually used for RK4 substages
        ! so no additional memory is used.
        hvy_work(:,:,:,1,hvy_id,1) = hvy_block(:,:,:,1,hvy_id)

        ! in some rare cases the fact that we have now in fact filled the ghost
        ! nodes correctly introduces a bias: it means, we can now interpolate, even
        ! without the ghost nodes on the coarse block filled. hence, reset the
        ! ghost nodes of the testing blocks:
        !-- x-direction
        hvy_block(1:g, :, :, 1, hvy_id)           = -1.0
        hvy_block(Bs(1)+g+1:Bs(1)+2*g, :, :, 1, hvy_id) = -1.0
        !-- y-direction
        hvy_block(:, 1:g, :, 1, hvy_id)           = -1.0
        hvy_block(:, Bs(2)+g+1:Bs(2)+2*g, :, 1, hvy_id) = -1.0
        !-- z-direction
        if ( params%dim == 3 ) then
            hvy_block(:, :, 1:g, 1, hvy_id)           = -1.0
            hvy_block(:, :, Bs(3)+g+1:Bs(3)+2*g, 1, hvy_id) = -1.0
        end if
    end do
end subroutine

subroutine unitTest_Sync( params, hvy_block, hvy_work, hvy_tmp, tree_ID, abort_on_fail, verbose)

    implicit none
    type (type_params), intent(inout)  :: params                     !> user defined parameter structure
    real(kind=rk),  intent(inout)      :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(out)         :: hvy_tmp(:, :, :, :, :)
    !> heavy work array: used for RHS evaluation in multistep methods (like RK4: u0, k1, k2 etc)
    real(kind=rk), intent(out)         :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    logical, intent(in)                :: abort_on_fail
    logical, optional, intent(in)      :: verbose

    integer(kind=ik)           :: k, l, lgt_id, hvy_id, fail_crit, fail_normal, ierr
    integer(kind=ik)           :: rank, Bs(3)
    real(kind=rk)              :: ddx(1:3), xx0(1:3)
    integer(kind=ik)           :: g, ix, iy, iz, g_depth, p_f, p_t, p_f_old, g_min
    integer(kind=tsize)        :: treecode
    real(kind=rk)              :: x, y, z
    character(len=80)          :: file_dump
    logical                    :: apply_verbose

    apply_verbose = .false.
    if (present(verbose)) apply_verbose = verbose

    rank = params%rank

    if (rank == 0) then
        write(*,'(20("_/¯\"))')
        write(*,'("UNIT TEST: Beginning rudimentary ghost nodes test")')
        write(*,'("UNIT TEST: It tests if all ghost patches are correctly set up")')
    end if

    Bs = params%Bs
    g  = params%g
    ! distance of g_min*ddx at periodicity boundaries is ignored to exclude step-effects in comparison
    g_min = max(abs(ubound(params%HD, dim=1)), abs(ubound(params%HR, dim=1))+1)
    fail_crit = 0
    fail_normal = 0

    if (rank == 0) then
        write(*,'("UNIT TEST: testing Bs=",i4," x ",i4," x ",i4," blocks-per-mpirank=",i5)')  Bs(1),Bs(2),Bs(3), params%number_blocks
    end if

    !---------------------------------------------------------------------------
    ! Test syncing: adaptive with one of the test grids, createGrid_simple_adaptive is enough actually
    !---------------------------------------------------------------------------
    call createGrid_simple_adaptive(params, hvy_block, tree_ID)

    ! now test syncing
    do l = params%g, params%g_rhs, -1
        ! reset values
        call unitTest_fill_linearly(params, hvy_block, hvy_work, hvy_tmp, tree_ID)

        g_depth = l
        ! first condition: SC are copied for this length anyways
        ! second condition: somehow values on coarse blocks at coarse-fine interfaces where neighbouring two fine blocks meet need this criterion
        !    it is probably not linked to HR directly but rather for each X in CDFXY necessary but thats what I found
        if (g_depth > abs(lbound(params%HD, dim=1))/2 .and. g_depth > abs(lbound(params%HR, dim=1))) then
            if (rank == 0) then
                write(*, '("UNIT TEST: Testing g_sync=", i0 ,", using sync with full restriction filter")') g_depth
            endif
            call sync_ghosts_tree( params, hvy_block, tree_ID, g_minus=g_depth, g_plus=g_depth)
        else
            if (rank == 0) then
                write(*, '("UNIT TEST: Testing g_sync=", i0 ,", using sync with copy for restriction filter")') g_depth
            endif
            call sync_ghosts_RHS_tree( params, hvy_block, tree_ID, g_minus=g_depth, g_plus=g_depth)
        endif

        ! test blocks, only test points that are in the middle part of the domain 0.25 < xyz < 0.75 to avoid periodicity problems
        p_f = 0
        p_t = 0
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
            treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))

            ! compute block spacing and origin from treecode
            call get_block_spacing_origin( params, lgt_id, xx0, ddx )

            p_f_old = p_f

            ! check data
            do iz = 1, size(hvy_block, 3)
                if (params%dim == 2) then
                    z = 0.0
                else
                    z = real(iz-(g+1), kind=rk) * ddx(3) + xx0(3)
                endif
                do iy = 1, size(hvy_block, 2)
                    y = real(iy-(g+1), kind=rk) * ddx(2) + xx0(2)
                    do ix = 1, size(hvy_block, 1)
                        x = real(ix-(g+1), kind=rk) * ddx(1) + xx0(1)
                        ! buffer zone around domain borders are ignored as there are jumps in function values which interpolation might not be able to handle
                        if (x < g_min*ddx(1) .or. x >= 1-g_min*ddx(1) .or. y < g_min*ddx(2) .or. y >= 1-g_min*ddx(2) &
                            .or. ((z < g_min*ddx(3) .or. z >= 1-g_min*ddx(3)) .and. params%dim == 3)) then
                            hvy_tmp(ix, iy, iz, 1, hvy_id) = -99.0
                        elseif (((iz <= g-g_depth .or. iz > Bs(3)+g+g_depth) .and. params%dim==3) &
                         .or. iy <= g-g_depth .or. iy > Bs(2)+g+g_depth &
                         .or. ix <= g-g_depth .or. ix > Bs(1)+g+g_depth) then
                            ! check outside, these should stay as -1 (unset)
                            if (abs(hvy_block(ix, iy, iz, 1, hvy_id) - (-1.0)) < 1e-8) then
                                hvy_tmp(ix, iy, iz, 1, hvy_id) = 0.0
                            else
                                p_f = p_f + 1
                                hvy_tmp(ix, iy, iz, 1, hvy_id) = -10.0
                            endif
                        else
                            ! check inside - these should be the correct ghost point values, be genereous to the wavelets for interpolating
                            if (abs(hvy_block(ix, iy, iz, 1, hvy_id) - hvy_work(ix, iy, iz, 1, hvy_id, 1)) > 0.5) then
                                p_f = p_f + 1
                            endif
                            hvy_tmp(ix, iy, iz, 1, hvy_id) = (hvy_block(ix, iy, iz, 1, hvy_id) - hvy_work(ix, iy, iz, 1, hvy_id, 1)) * 100.0
                        endif
                    enddo
                enddo
            enddo

            if (p_f - p_f_old > 0) then
                write(*, '("UNIT TEST:    Difference in block ", i0, " with lgt_id= ", i0)') treecode, lgt_id
            endif
        enddo

        p_f_old = p_f
        call MPI_Allreduce(p_f_old, p_f, 1, MPI_INTEGER4, MPI_SUM, WABBIT_COMM, ierr)


        if (p_f > 0) then
            if (rank == 0) then
                write(*, '(A, i2)') "UNIT TEST:    ERROR! Did not sync all ghost nodes correctly for g_sync= ", g_depth
            endif
            ! those are used in simulations often and should be deemed as a critical failure if they do not work
            if (g_depth == params%g .or. g_depth == params%g_RHS) then
                fail_crit = fail_crit + 1
            else
                fail_normal = fail_normal + 1
            endif
        else
            if (rank == 0) then
                write(*,'("UNIT TEST:    HURRAYYY! Syncing works as expected for g_sync= ", i2)') g_depth
            endif
        endif

        ! dump block values for testing - this is 2D!
        if (apply_verbose) then
            if (any(l == (/ 1, 2, 3/))) then  ! change values if only debugging one case
                do k = 1, hvy_n(tree_ID)
                    hvy_ID = hvy_active(k, tree_ID)
                    call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
                    treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))

                    write(file_dump, '(A, i0, A, i0, A)') "block_dumped_g=", g_depth, "_tc=", treecode, ".t"
                    ! print block in a fancy way to file, digits is set to be size of largest 2D number +1
                    call dump_block_fancy(hvy_block(:, :, :, 1:1, hvy_id), file_dump, Bs, g, to_int=.true., digits=4)
                    ! call dump_block_fancy(hvy_block(:, :, :, 1:1, hvy_id), file_dump, Bs, g, to_int=.false., digits=4)

                    open(unit=32, file=file_dump, status='unknown', position='append')
                    write(32, '(A)') ""  ! new line
                    write(32, '(A, es8.2, A)') "1 corresponds to a deviation of 1 between correct value and synched one. Value deemed wrong at >50. Values with -99 were not considered"
                    close(32)

                    call dump_block_fancy(hvy_tmp(:, :, :, 1:1, hvy_id), file_dump, Bs, g, to_int=.true., digits=4, append=.true.)
                enddo
            endif

            ! write information so that we can debug neighbor relations
            if (l == params%g) then
                ! write neighbor_blocks file
                if (params%rank == 0) then
                    open(16,file='hvy_neighbors.dat',status='replace')
                    write(16,'(A)') "% lgt_id + 157*neighbors"
                    do k=1,hvy_n(tree_ID)
                        hvy_ID = hvy_active(k, tree_ID)
                        call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)

                        write(16,'(169(i0, 1x))') lgt_id, hvy_neighbor(k, :)
                    enddo
                    close(16)
                endif

                call saveHDF5_tree("gridsync_0001.h5", 0.0_rk, 1, 1, params, hvy_block, tree_ID, no_sync=.true.)
            endif
        endif
    enddo

    ! delete the grid we created for this subroutine
    call reset_tree(params, .true., tree_ID=tree_ID)

    ! print output
    if (fail_crit == 0 .and. fail_normal == 0) then
        if (rank == 0) then
            write(*,'("UNIT TEST: Numbers have been traded in succession. The customers seem pleased with their share.")')
        endif
    else
        if (rank == 0) then
            write(*,'("UNIT TEST: Numbers have been traded, but ", i0, " critical and ", i0, " normal customers are furious and reject their share.")') fail_crit, fail_normal
            if (fail_crit > 0 .and. abort_on_fail) then
                call abort(20240721, "We cannot continue our business with unhappy critical tests. Let's have a break for now and check what's going on!")
            endif
        endif
    endif

end subroutine
