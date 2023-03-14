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
    integer(kind=ik)                        :: g, ix, iy, iz, nc, ic
    integer(kind=ik), dimension(3)          :: Bs
    real(kind=rk), allocatable :: norm(:), norm_ref(:)

    if (params%rank == 0) then
        write(*,'(80("~"))')
        write(*,'("UNIT TEST: testing if IWT(FWT(u)) = Id")')
        write(*,'("This test is performed on an equidistant grid.")')
        write(*,'("It checks if the filter banks HD,GD,HR,GR are correct.")')
        write(*,'("Wavelet=",A," g=", i2)') trim(adjustl(params%wavelet)), params%g
    end if

    Bs = params%Bs
    g = params%g
    nc = params%n_eqn

    allocate(norm(1:params%n_eqn))
    allocate(norm_ref(1:params%n_eqn))

    !----------------------------------------------------------------------------
    ! create an equidistant grid on level J=1 (and not Jmin, because that may well be 0)
    !----------------------------------------------------------------------------
    call createEquidistantGrid_tree( params, 1, .true., tree_ID )

    !----------------------------------------------------------------------------
    ! create just some data...
    !----------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvyID = hvy_active(k,tree_ID)
        do ic = 1, nc
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    hvy_block(ix,iy,:,ic,hvyID) = rand_nbr()
                end do
            end do
        end do
    end do

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, &
    hvy_active(:,tree_ID), hvy_n(tree_ID) )

    call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm_ref)

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

    norm = norm / norm_ref

    if (params%rank==0) write(*,*) "Relative L2 error in IWT(FWT(u)) is: ", norm

    if (norm(1)>1.0e-14_rk) then
        call abort(230306608, "Error in IWT(FWT(U)) is too large! Call the police. Dial 17 in France." )
    else
        if (params%rank==0) write(*,*) "How lovely that this test suceeded! You've earned yourself a refreshing beverage."
    endif

    if (params%rank == 0) then
        write(*,'(80("~"))')
    end if

    ! delete the grid we created for this subroutine
    call reset_tree(params, .true., tree_ID=tree_ID)
end subroutine
