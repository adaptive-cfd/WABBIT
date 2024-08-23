subroutine unit_test_waveletDecomposition_invertibility( params, hvy_block, hvy_work, hvy_tmp, tree_ID, verbose)

    implicit none
    type (type_params), intent(inout)       :: params                     !> user defined parameter structure
    real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(out)              :: hvy_tmp(:, :, :, :, :)
    !> heavy work array: used for RHS evaluation in multistep methods (like RK4: u0, k1, k2 etc)
    real(kind=rk), intent(out)              :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik), intent(in)            :: tree_ID
    logical, optional, intent(in)           :: verbose

    integer(kind=ik)                        :: k, hvy_id, lgt_id
    integer(kind=ik)                        :: g, ix, iy, iz, nc, ic, ii
    integer(kind=ik), dimension(3)          :: Bs
    real(kind=rk), allocatable :: norm_1(:), norm_ref(:), norm_2(:)
    integer(kind=tsize)        :: treecode
    character(len=80)                       :: file_dump
    logical                                 :: apply_verbose

    apply_verbose = .false.
    if (present(verbose)) apply_verbose = verbose

    if (params%rank == 0) then
        write(*, '("")')  ! newline
        write(*,'(20("_/¯\"))')
        write(*,'("UNIT TEST: Testing if adapt(adapt(U)) = adapt(U), performed on a non-equidistant grid.")')
        write(*,'("UNIT TEST: It checks if the implementation of wavelet decomposition and reconstruction in adapt_tree is correct.")')
        write(*,'("UNIT TEST: Due to CVS or coarse extension the first adaption alter values, but later application should not alter them.")')
    end if

    Bs = params%Bs
    g  = params%g
    nc = params%n_eqn

    allocate(norm_1(1:params%n_eqn))
    allocate(norm_2(1:params%n_eqn))
    allocate(norm_ref(1:params%n_eqn))

    !----------------------------------------------------------------------------
    ! Construct a random grid for testing
    !----------------------------------------------------------------------------
    ! this parameter controls roughly how dense the random grid is, i.e., in % of the
    ! complete memory.
    if (params%Jmax/=params%Jmin .or. params%Jmax>1) then
        params%max_grid_density = 0.10_rk
        ! perform at min 3 iterations of random refinement/coarsening
        call createRandomGrid_tree( params, hvy_block, hvy_tmp, level_init=params%Jmin, verbosity=.true., iterations=min(3, params%Jmax-params%Jmin), tree_ID=tree_ID )
    else
        if (params%rank == 0) then
            write(*, '("UNIT TEST: With these grid settings we need to do this test with an equidistand grid on level ", i0, ".")') params%Jmin
        endif
        call createEquidistantGrid_tree( params, hvy_block, params%Jmin, .true., tree_ID )
    endif

    !----------------------------------------------------------------------------
    ! create just some data...
    !----------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)
        call random_data(hvy_block(:,:,:,:,hvy_id))
    end do

    call sync_ghosts_tree( params, hvy_block, tree_ID )
    call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm_1)

    ! We adapt two times, the first time coarse extension and CVS alters the data, but the second time it should do the same
    ! adaption (first time)
    call adapt_tree( 0.0_rk, params, hvy_block, tree_ID, params%coarsening_indicator, hvy_tmp)
    ! this norm is the reference norm, as we compare to this one
    call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm_ref)

    ! print out some debugging infos to files
    if (apply_verbose) then
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k,tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
            write(file_dump, '(A, i0, A, i0, A)') "block_dumped_tc=", treecode, ".t"
            call dump_block_fancy(hvy_block(:, :, :, 1:1, hvy_id), file_dump, Bs, g)
        enddo
    endif

    ! adaption (second time)
    call adapt_tree( 0.0_rk, params, hvy_block, tree_ID, params%coarsening_indicator, hvy_tmp)
    ! we compare norms - this is not the most elegant way, but it'll do the trick for now.
    call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm_2)

    ! print out some debugging infos to files
    if (apply_verbose) then
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k,tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
            write(file_dump, '(A, i0, A, i0, A)') "block_dumped_tc=", treecode, ".t"
            open(unit=32, file=file_dump, status='unknown', position='append')
            write(32, '(A)') ""  ! new line
            write(32, '(A)') "After adapt-adapt:"
            write(32, '(A)') ""  ! new line
            close(32)

            call dump_block_fancy(hvy_tmp(:, :, :, 1:1, hvy_id), file_dump, Bs, g, append=.true.)
        enddo
    endif

    norm_1 = abs(norm_1 / norm_ref - 1.0_rk)
    norm_2 = abs(norm_2 / norm_ref - 1.0_rk)

    if (params%rank==0) write(*,'(A, es15.8)') "UNIT TEST: Relative L2 error between original data and after first adaption : ", norm_1(1)
    if (params%rank==0) write(*,'(A, es15.8)') "UNIT TEST: Relative L2 error between first and second adaption (should be 0): ", norm_2(1)

    if (norm_2(1)>1.0e-14_rk) then
        call abort(230306608, "Error! adapt(adapt(U)) /= adapt(U)! Call the police! Danger!!" )
    else
        if (params%rank==0) then
            write(*,'(20("_/¯\"))')
            write(*,'(A)') "           ( ("
            write(*,'(A)') "            ) )"
            write(*,'(A)') "          ........      How lovely that the invertibility test succeeded!"
            write(*,'(A)') "          |      |]       You've earned yourself a refreshing beverage."
            write(*,'(A)') "          \      /"
            write(*,'(A)') "           `----'"
        endif
    endif

    if (params%rank == 0) then
        write(*,'(20("_/¯\"))')
    end if

    ! delete the grid we created for this subroutine
    call reset_tree(params, .true., tree_ID=tree_ID)
end subroutine
