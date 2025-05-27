subroutine unit_test_waveletDecomposition( params, hvy_block, hvy_work, hvy_tmp, tree_ID)

    implicit none
    type (type_params), intent(inout)       :: params                     !> user defined parameter structure
    real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(out)              :: hvy_tmp(:, :, :, :, :)
    !> heavy work array: used for RHS evaluation in multistep methods (like RK4: u0, k1, k2 etc)
    real(kind=rk), intent(out)              :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik), intent(in)            :: tree_ID

    integer(kind=ik)                        :: k, hvy_id, lgt_id
    integer(kind=ik)                        :: g, ix, iy, iz, nc, ic, ii
    integer(kind=ik), dimension(3)          :: Bs
    real(kind=rk), allocatable :: norm(:), norm_ref(:), wc(:,:,:,:,:)
    character(len=cshort)                   :: debug_name

    if (params%rank == 0) then
        write(*, '("")')  ! newline
        write(*,'(20("_/¯\"))')
        write(*,'("UNIT TEST: testing if IWT(FWT(U)) = U, performed on an equidistant grid.")')
        write(*,'("UNIT TEST: It checks if the filter banks HD,GD,HR,GR are correct.")')
        write(*,'("UNIT TEST: Wavelet=",A," g=", i2)') trim(adjustl(params%wavelet)), params%g
    end if

    Bs = params%Bs
    g = params%g
    nc = params%n_eqn

    allocate(norm(1:params%n_eqn))
    allocate(norm_ref(1:params%n_eqn))

    !----------------------------------------------------------------------------
    ! create an equidistant grid on level J=1 (and not Jmin, because that may well be 0)
    !----------------------------------------------------------------------------
    call createEquidistantGrid_tree( params, hvy_block, min(1, params%Jmax), .true., tree_ID )

    !----------------------------------------------------------------------------
    ! create just some data...
    !----------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)
        call random_data(hvy_block(:,:,:,:,hvy_id))
    end do

    call sync_ghosts_tree( params, hvy_block, tree_ID )

    call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm_ref)

    !----------------------------------------------------------------------------
    ! FWT
    !----------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        ! store original values
        hvy_tmp(:,:,:,1:nc,hvy_id) = hvy_block(:,:,:,1:nc,hvy_id)

        ! ! debug before WD
        ! write(debug_name, '(A,i0,A)') 'block_00_TC', lgt_block(lgt_id, IDX_TC_2), '.dat'
        ! call dump_block_fancy(hvy_block(:,:,:,1:nc,hvy_id), debug_name, params%Bs, params%g, digits=2)

        ! Wavelet decomposition
        call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvy_id))
    end do

    call sync_ghosts_tree( params, hvy_block, tree_ID )

    !---------------------------------------------------------------------------
    ! testing of conversion routines on FWT transformed data
    !---------------------------------------------------------------------------
    allocate(wc(1:size(hvy_block,1), 1:size(hvy_block,2), 1:size(hvy_block,3), 1:size(hvy_block,4), 1:8 ))
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)
        ! Test Sp -> InfMall -> Mall -> InfMall -> Sp
        call spaghetti2inflatedMallat_block(params, hvy_block(:,:,:,:,hvy_id), wc)
        call inflatedMallat2Mallat_block(params, wc, hvy_block(:,:,:,:,hvy_id))
        call Mallat2inflatedMallat_block(params, hvy_block(:,:,:,:,hvy_id), wc)
        call inflatedMallat2spaghetti_block(params, wc, hvy_block(:,:,:,:,hvy_id))

        ! Test Sp -> Mall -> Sp
        call spaghetti2Mallat_block(params, hvy_block(:,:,:,:,hvy_id), wc(:,:,:,:,1))
        call Mallat2spaghetti_block(params, wc(:,:,:,:,1), hvy_block(:,:,:,:,hvy_id))
    end do
    deallocate(wc)


    !----------------------------------------------------------------------------
    ! IWT
    !----------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        ! ! debug before WR
        ! write(debug_name, '(A,i0,A)') 'block_WD_TC', lgt_block(lgt_id, IDX_TC_2), '.dat'
        ! call dump_block_fancy(hvy_block(:,:,:,1:nc,hvy_id), debug_name, params%Bs, params%g, digits=2)

        ! Wavelet reconstruction
        call waveletReconstruction_block(params, hvy_block(:,:,:,:,hvy_id))

        ! ! debug after WR
        ! write(debug_name, '(A,i0,A)') 'block_WR_TC', lgt_block(lgt_id, IDX_TC_2), '.dat'
        ! call dump_block_fancy(hvy_block(:,:,:,1:nc,hvy_id), debug_name, params%Bs, params%g, digits=2)

        ! error IWT(FWT(u)) - u
        hvy_block(:,:,:,1:nc,hvy_id) = hvy_block(:,:,:,1:nc,hvy_id) - hvy_tmp(:,:,:,1:nc,hvy_id)
    end do

    ! compute norm of error
    call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm)

    do k = 1, nc
        if (norm_ref(k) > 1.0e-12) then
            norm(k) = norm(k) / norm_ref(k)
        endif
    enddo

    if (params%rank==0) write(*,'(A, es15.8)') "UNIT TEST: Relative L2 error in IWT(FWT(u)) is: ", norm(1)

    if (norm(1)>1.0e-14_rk) then
        ! do k = 1, hvy_n(tree_ID)
        !     hvy_id = hvy_active(k,tree_ID)
        !     call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        !     write(debug_name, '(A,i0,A)') 'block_TC', lgt_block(lgt_id, IDX_TC_2), '.dat'
        !     call dump_block_fancy(hvy_block(:,:,:,1:nc,hvy_id), debug_name, params%Bs, params%g, digits=2)
        ! end do
        call abort(230306608, "Error! IWT(FWT(U)) /= U! Call the police! Danger!!" )
    else
        if (params%rank==0) then
            write(*,'(20("_/¯\"))')
            write(*,'(A)') "           ( ("
            write(*,'(A)') "            ) )"
            write(*,'(A)') "          ........   How lovely that the wavelet decomposition test succeeded!"
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
