subroutine unit_test_refineCoarsen( params, hvy_block, hvy_work, hvy_tmp, tree_ID, verbose)

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
    real(kind=rk), allocatable :: norm(:), norm_ref(:)
    integer(kind=tsize)        :: treecode
    character(len=80)                       :: file_dump
    logical                                 :: apply_verbose, error_OOM

    apply_verbose = .false.
    if (present(verbose)) apply_verbose = verbose

    if (params%rank == 0) then
        write(*, '("")')  ! newline
        write(*,'(20("_/¯\"))')
        write(*,'("UNIT TEST: Testing if Coarsen(Refine(U)) = U, performed on an equidistant grid.")')
        write(*,'("UNIT TEST: It checks if the implementation of Interpolation, Refinement and Block Merging are correct.")')
    end if

    Bs = params%Bs
    g  = params%g
    nc = params%n_eqn

    allocate(norm(1:params%n_eqn))
    allocate(norm_ref(1:params%n_eqn))

    if (params%Jmax<2) then
        if (params%rank==0) write(*,'(A)') "UNIT TEST: Test cannot be performed for Jmax<2, skipping it."
        return
    endif

    if (params%Jmax==params%Jmin) then
        if (params%rank==0) write(*,'(A)') "UNIT TEST: Test cannot be performed for equidistant grids, skipping it."
        return
    endif

    !----------------------------------------------------------------------------
    ! create an equidistant grid on level J=1 (and not Jmin, because that may well be 0)
    !----------------------------------------------------------------------------
    call createEquidistantGrid_tree( params, hvy_block, params%Jmin, .true., tree_ID )

    !----------------------------------------------------------------------------
    ! create just some data...
    !----------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)
        call random_data(hvy_block(:,:,:,:,hvy_id))
    end do

    call sync_ghosts_tree( params, hvy_block, tree_ID )

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

    call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm_ref)

    ! refine
    call refine_tree( params, hvy_block, hvy_tmp, "everywhere", tree_ID, error_OOM )

    if (error_OOM) call abort(2512118,"Refinement failed, out of memory. Try with more memory.")

    call sync_ghosts_tree( params, hvy_block, tree_ID )

    ! coarsening (back to the original level)
    call adapt_tree( 0.0_rk, params, hvy_block, tree_ID, "everywhere", hvy_tmp)

    ! we compare norms - this is not the most elegant way, but it'll do the trick for now.
    call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm)

    ! print out some debugging infos to files
    if (apply_verbose) then
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k,tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
            write(file_dump, '(A, i0, A, i0, A)') "block_dumped_tc=", treecode, ".t"
            open(unit=32, file=file_dump, status='unknown', position='append')
            write(32, '(A)') ""  ! new line
            write(32, '(A)') "After refine-coarsen:"
            write(32, '(A)') ""  ! new line
            close(32)

            call dump_block_fancy(hvy_tmp(:, :, :, 1:1, hvy_id), file_dump, Bs, g, append=.true.)
        enddo
    endif

    norm = abs(norm / norm_ref - 1.0_rk)

    if (params%rank==0) write(*,'(A, es15.8)') "UNIT TEST: Relative L2 error in Coarsen(Refine(u)) is: ", norm(1)

    if (norm(1)>1.0e-14_rk) then
        call abort(230306608, "Error! Coarsen(Refine(U)) /= U ! Call the police! Danger!!" )
    else
        if (params%rank==0) then
            write(*,'(20("_/¯\"))')
            write(*,'(A)') "           ( ("
            write(*,'(A)') "            ) )"
            write(*,'(A)') "          ........      How lovely that the refine-coarsen test succeeded!"
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
