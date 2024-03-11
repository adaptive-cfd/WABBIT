subroutine unitTest_refineCoarsen( params, hvy_block, hvy_work, hvy_tmp, tree_ID)

    implicit none
    type (type_params), intent(inout)       :: params                     !> user defined parameter structure
    real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(out)              :: hvy_tmp(:, :, :, :, :)
    !> heavy work array: used for RHS evaluation in multistep methods (like RK4: u0, k1, k2 etc)
    real(kind=rk), intent(out)              :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik), intent(in)            :: tree_ID

    integer(kind=ik)                        :: k, hvyID
    integer(kind=ik)                        :: g, ix, iy, iz, nc, ic, ii
    integer(kind=ik), dimension(3)          :: Bs
    real(kind=rk), allocatable :: norm(:), norm_ref(:), wc(:,:,:,:,:)

    if (params%rank == 0) then
        write(*,'(80("~"))')
        write(*,'("UNIT TEST: testing if Coarsen(Refine(u)) = Id")')
        write(*,'("This test is performed on an equidistant grid.")')
        write(*,'("It checks if the implementation of Interpolation, Refinement and Block Merging are correct.")')
    end if

    Bs = params%Bs
    g  = params%g
    nc = params%n_eqn

    allocate(norm(1:params%n_eqn))
    allocate(norm_ref(1:params%n_eqn))

    if (params%Jmax<2) then
        if (params%rank==0) write(*,*) "Test cannot be performed because of level restrictions: params%jmax=", params%jmax
        return
    endif

    if (params%Jmax==params%Jmin) then
        if (params%rank==0) write(*,*) "Test cannot be performed because of level restrictions: params%jmax=", params%jmax
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
        hvyID = hvy_active(k,tree_ID)
        if (params%dim == 3) then
            do ic = 1, nc
                do iz = g+1, Bs(3)+g
                    do iy = g+1, Bs(2)+g
                        do ix = g+1, Bs(1)+g
                            hvy_block(ix,iy,iz,ic,hvyID) = rand_nbr()
                        end do
                    end do
                end do
            end do
        else
            do ic = 1, nc
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        hvy_block(ix,iy,:,ic,hvyID) = rand_nbr()
                    end do
                end do
            end do
        endif
    end do

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, &
    hvy_active(:,tree_ID), hvy_n(tree_ID) )

    call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm_ref)

    ! refine
    call refine_tree( params, hvy_block, hvy_tmp, "everywhere", tree_ID  )

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, &
    hvy_active(:,tree_ID), hvy_n(tree_ID) )

    ! coarsening (back to the original level)
    call adapt_tree( 0.0_rk, params, hvy_block, tree_ID, "everywhere", hvy_tmp)

    ! we compare norms - this is not the most elegant way, but it'll do the trick for now.
    call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2", norm)

    norm = abs(norm / norm_ref - 1.0_rk)

    if (params%rank==0) write(*,*) "Relative L2 error in Coarsen(Refine(u)) is: ", norm

    if (norm(1)>1.0e-14_rk) then
        call abort(230306608, "Error in IWT(FWT(U)) is too large! Call the police! Danger!!" )
    else
        if (params%rank==0) then
            write(*,'(A)') "           ( ("
            write(*,'(A)') "            ) )"
            write(*,'(A)') "          ........"
            write(*,'(A)') "          |      |]"
            write(*,'(A)') "          \      /    How lovely that this test suceeded! You've earned yourself a refreshing beverage."
            write(*,'(A)') "           `----'"
        endif
    endif

    if (params%rank == 0) then
        write(*,'(80("~"))')
    end if

    ! delete the grid we created for this subroutine
    call reset_tree(params, .true., tree_ID=tree_ID)
end subroutine
