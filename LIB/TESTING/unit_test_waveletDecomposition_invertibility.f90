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

    integer(kind=ik)                        :: k, hvy_id, lgt_id, i_adapt, it_random, l_init
    integer(kind=ik)                        :: g, ix, iy, iz, nc, ic, ii, block_dump_max, Bs(1:3), io(1:3)
    real(kind=rk), allocatable :: norm_1(:), norm_ref(:), norm_2(:)
    real(kind=rk)                           :: x0(1:3), dx(1:3)
    integer(kind=tsize)        :: treecode
    character(len=80)                       :: file_dump
    logical                                 :: apply_verbose, problem, grid_is_equidistant

    apply_verbose = .false.
    if (present(verbose)) apply_verbose = verbose
    block_dump_max = 10

    if (params%rank == 0) then
        write(*, '("")')  ! newline
        write(*,'(20("_/¯\"))')
        write(*,'("UNIT TEST: Testing if adapt(adapt(U)) = adapt(U), performed on a non-equidistant grid.")')
        write(*,'("UNIT TEST: It checks if the implementation of wavelet decomposition and reconstruction in adapt_tree is correct.")')
        write(*,'("UNIT TEST: Due to CVS or coarse extension the values after the first adaption should be changed for an adaptive grid, but later applications should give the same results.")')
        write(*,'("UNIT TEST: This test is also used for verifying the minimum blocksize needed for the adapt_tree loop.")')
    end if

    Bs = params%Bs
    g  = params%g
    nc = params%n_eqn
    ! for odd block sizes, we have an overlap of the points from the center line and want to ignore those
    io = 0
    do k = 1,params%dim
        if (modulo(Bs(k),2) == 1) io(k) = 1
    enddo

    allocate(norm_1(1:params%n_eqn))
    allocate(norm_2(1:params%n_eqn))
    allocate(norm_ref(1:params%n_eqn))

    !----------------------------------------------------------------------------
    ! Construct a random grid for testing
    !----------------------------------------------------------------------------
    ! this parameter controls roughly how dense the random grid is, i.e., in % of the
    ! complete memory.
    if (params%Jmax/=params%Jmin .and. params%Jmax>1) then
        ! perform random refinement/coarsening
        it_random = max(4, params%Jmax-params%Jmin)
        l_init = max(min(3, params%Jmax), params%Jmin)  ! init on level 3 but adhere to Jmin Jmax restrictions
        call createRandomGrid_tree( params, hvy_block, hvy_tmp, level_init=l_init, verbosity=.true., iterations=it_random, tree_ID=tree_ID )
    else
        if (params%rank == 0) then
            write(*, '("UNIT TEST: With these grid settings we need to do this test with an equidistand grid on level ", i0, ".")') params%Jmin
        endif
        call createEquidistantGrid_tree( params, hvy_block, params%Jmin, .true., tree_ID )
    endif

    ! call createGrid_simple_adaptive(params, hvy_block, tree_ID)

    ! sometimes by chance or on purpose the resulting grid is equidistant, then every adaption should not alter the data
    grid_is_equidistant = maxActiveLevel_tree(tree_ID) == minActiveLevel_tree(tree_ID)

    !----------------------------------------------------------------------------
    ! create just some data...
    !----------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)
        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)

        call get_block_spacing_origin( params, lgt_id, x0, dx )

        if (all(io == 0)) then
            call random_data(hvy_block(:,:,:,:,hvy_id))
        else
            ! redundant grid, so we have to make sure that redundant points receive the same values
            call random_data_unique( hvy_block(:,:,:,:,hvy_id), x0, dx, (/params%g,params%g,params%g/), params%domain_size, params%Jmax, params%BS)
        endif
    end do

    call sync_ghosts_tree( params, hvy_block, tree_ID )

    ! print out some debugging infos to files
    if (apply_verbose) then
        ! if there are not too many blocks we print the values directly
        if (lgt_n(tree_ID) <= block_dump_max) then
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k,tree_ID)
                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
                write(file_dump, '(A, i0, A, i0, A)') "block_dumped_it0_tc=", treecode, ".t"
                open(unit=32, file=file_dump, status='unknown', position='append')
                write(32, '(A)') ""  ! new line
                write(32, '(A)') "Beginning:"
                write(32, '(A)') ""  ! new line
                close(32)
                call dump_block_fancy(hvy_block(:, :, :, 1:1, hvy_id), file_dump, Bs, g)
            enddo
        endif
        ! call saveHDF5_tree("WR-invertibility1_0000.h5", 0.0_rk, 0, 1, params, hvy_block, tree_ID, no_sync=.true., save_ghosts=.false.)
    endif

    ! We adapt two times (or more), the first time coarse extension and CVS alters the data, but afterwards it should do the same
    problem = .false.  ! init
    do i_adapt = 1,2

        ! this norm is the reference norm, as we compare to this one
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm_ref)

        call adapt_tree( dble(i_adapt), params, hvy_block, tree_ID, params%coarsening_indicator, hvy_tmp)
        ! this norm is the reference norm, as we compare to this one
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm_1)
        do k = 1, nc
            if (norm_ref(k) > 1.0e-12) then
                norm_1(k) = abs(norm_1(k) / norm_ref(k) - 1.0_rk)
            endif
        enddo

        ! prepare for test results, any test after the second should pass
        if (i_adapt > 1 .or. grid_is_equidistant) then
            problem = norm_1(1)>1.0e-14_rk .or. problem
        endif

        ! print out some debugging infos to files
        if (apply_verbose) then
            ! if there are not too many blocks we print the values directly
            if (lgt_n(tree_ID) <= block_dump_max) then
                do k = 1, hvy_n(tree_ID)
                    hvy_id = hvy_active(k,tree_ID)
                    call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                    treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
                    write(file_dump, '(A, i0, A,i0, A, i0, A)') "block_dumped_it", i_adapt,"_tc=", treecode, ".t"
                    open(unit=32, file=file_dump, status='unknown', position='append')
                    write(32, '(A)') ""  ! new line
                    write(32, '(A)') "After adapt:"
                    write(32, '(A)') ""  ! new line
                    close(32)
                    call dump_block_fancy(hvy_block(:, :, :, 1:1, hvy_id), file_dump, Bs, g, append=.true.)
                enddo
            endif
            write(file_dump, '(A, i0, A)') "WR-invertibility", i_adapt, "_0000.h5"
            call saveHDF5_tree(file_dump, 0.0_rk, 0, 1, params, hvy_block, tree_ID, no_sync=.true., save_ghosts=.false.)
        endif

        if (params%rank==0) then
            write(*,'(A, i0, A, es15.8)') "UNIT TEST: Relative L2 error between before and after adaption no ", i_adapt, " : ", norm_1(1)
        endif
    enddo

    if (apply_verbose) then
        call summarize_profiling( WABBIT_COMM )
    endif

    ! report to terminal / log and abort if critical
    if (params%rank==0) then
        if (problem) then
            call abort(230306608, "Error! adapt(adapt(U)) /= adapt(U)! Call the police! Danger!!" )
        endif
            
        write(*,'(20("_/¯\"))')
        write(*,'(A)') "           ( ("
        write(*,'(A)') "            ) )"
        write(*,'(A)') "          ........      How lovely that the invertibility test succeeded!"
        write(*,'(A)') "          |      |]       You've earned yourself a refreshing beverage."
        write(*,'(A)') "          \      /"
        write(*,'(A)') "           `----'"
        write(*,'(20("_/¯\"))')
    end if

    ! delete the grid we created for this subroutine
    call reset_tree(params, .true., tree_ID=tree_ID)
end subroutine
