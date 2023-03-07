subroutine unitTest_waveletDecomposition( params, hvy_block, hvy_work, hvy_tmp, tree_ID)

    implicit none
    type (type_params), intent(inout)       :: params                     !> user defined parameter structure
    real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(out)              :: hvy_tmp(:, :, :, :, :)
    !> heavy work array: used for RHS evaluation in multistep methods (like RK4: u0, k1, k2 etc)
    real(kind=rk), intent(out)              :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik), intent(in)            :: tree_ID

    integer(kind=ik)                        :: k, hvyID
    integer(kind=ik)                        :: g, ix, iy, iz, nc
    integer(kind=ik), dimension(3)          :: Bs
    real(kind=rk), allocatable :: norm(:)

    if (params%rank == 0) then
        write(*,'(80("~"))')
        write(*,'("UNIT TEST: testing if IWT(FWT(u)) = Id")')
        write(*,'("This test is performed on an equidistant grid.")')
        write(*,'("It checks if the filter banks HD,GD,HR,GR are correct.")')
        write(*,'("Wavelet=",A," g=", i2)') trim(adjustl(params%wavelet)), params%g
    end if

    nc = 1

    allocate(norm(1:params%n_eqn))

    !----------------------------------------------------------------------------
    ! create an equidistant grid on level J=1 (and not Jmin, because that may well be 0)
    !----------------------------------------------------------------------------
    call createEquidistantGrid_tree( params, 1, .true., tree_ID )

    !----------------------------------------------------------------------------
    ! create just some data...
    !----------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        do iy = g+1, Bs(2)+g
            do ix = g+1, Bs(1)+g
                hvy_block(ix,iy,:,1:nc,hvyID) = rand_nbr()
            end do
        end do
    end do

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, &
    hvy_active(:,tree_ID), hvy_n(tree_ID) )

    !----------------------------------------------------------------------------
    ! FWT
    !----------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        hvy_tmp(:,:,:,1:nc,hvyID) = hvy_block(:,:,:,1:nc,hvyID)
        call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvyID))
    end do

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, &
    hvy_active(:,tree_ID), hvy_n(tree_ID) )

    !----------------------------------------------------------------------------
    ! IWT
    !----------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        call waveletReconstruction_block(params, hvy_block(:,:,:,:,hvyID))
        ! error IWT(FWT(u)) - u
        hvy_block(:,:,:,1:nc,hvyID) = hvy_block(:,:,:,1:nc,hvyID) - hvy_tmp(:,:,:,1:nc,hvyID)
    end do

    ! compute norm of error
    call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm)

    if (params%rank==0) write(*,*) "Error in IWT(FWT(u)) is: ", norm(1)

    if (norm(1)>1.0e-14_rk) then
        call abort(230306608, "Error in IWT(FWT(U)) is too large! Call the police. Dial 17 in France." )
    else
        if (params%rank==0) write(*,*) "How lovely that this test suceeded! You've earned yourself a refreshing beverage."
    endif

    if (params%rank == 0) then
        write(*,'(80("~"))')
    end if
end subroutine
