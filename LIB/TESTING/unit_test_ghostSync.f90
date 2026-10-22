subroutine unit_test_ghostSync( params, hvy_block, hvy_work, hvy_tmp, tree_ID, abort_on_fail, verbose)

    implicit none
    type (type_params), intent(inout) :: params                     !> user defined parameter structure
    real(kind=rk),  intent(inout)     :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(out)        :: hvy_tmp(:, :, :, :, :)
    !> heavy work array: used for RHS evaluation in multistep methods (like RK4: u0, k1, k2 etc)
    real(kind=rk), intent(out)        :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik), intent(in)      :: tree_ID
    logical, intent(in)               :: abort_on_fail
    logical, optional, intent(in)     :: verbose

    integer(kind=ik)                  :: k, it_random, lgt_id, hvy_id, l_init
    integer(kind=ik)                  :: rank, number_procs
    real(kind=rk)                     :: ddx(1:3), xx0(1:3)
    integer(kind=ik)                  :: g, number_blocks, ix, iy, iz, JmaxA, JminA
    integer(kind=ik), dimension(3)    :: Bs
    real(kind=rk)                     :: Lx, Ly, Lz, x, y, z
    integer(kind=ik)                  :: d
    integer(kind=ik)                  :: ifrequ

    ! error variable
    real(kind=rk)                     :: error2(1:32), error1(1:32), error_L2, error_Linfty, norm_L2, norm_Linfty, t0
    integer(kind=ik)                  :: ierr, ii
    logical                           :: test, apply_verbose

    apply_verbose = .false.
    if (present(verbose)) apply_verbose = verbose

    rank = params%rank
    Lx = params%domain_size(1)
    Ly = params%domain_size(2)
    Lz = params%domain_size(3)

    d = params%dim


    if (rank == 0) then
        write(*,'("")')  ! newline
        write(*,'(20("_/¯\"))')
        write(*,'("UNIT TEST: Beginning ghost nodes test")')
        write(*,'("UNIT TEST: It tests correctness of ghost patches and prediction order")')
    end if

    Bs = params%Bs
    g  = params%g
    number_procs  = params%number_procs
    number_blocks = params%number_blocks

    if (rank == 0) then
        write(*,'("UNIT TEST: testing Bs=",i3," x ",i3," x ",i3," blocks-per-mpirank=",i0)')  Bs(1),Bs(2),Bs(3), params%number_blocks
    end if

    !---------------------------------------------------------------------------
    ! Step 1: Construct a random grid for testing. Note we keep this grid
    ! and perform the same test for differnet frequencies (= resolutions) only on
    ! this one grid.
    !---------------------------------------------------------------------------
    ! this parameter controls roughly how dense the random grid is, i.e., in % of the
    ! complete memory.
    ! perform at min 3 iterations of random refinement/coarsening
    it_random = max(4, params%Jmax-params%Jmin)
    l_init = max(min(3, params%Jmax), params%Jmin)  ! init on level 3 but adhere to Jmin Jmax restrictions
    call createRandomGrid_tree( params, hvy_block, hvy_tmp, level_init=l_init, verbosity=.true., iterations=it_random, tree_ID=tree_ID )

    JmaxA = maxActiveLevel_tree(tree_ID)
    JminA = minActiveLevel_tree(tree_ID)
    if (JmaxA == JminA) then
        if (params%rank==0) write(*,'(A)') "UNIT TEST: By chance, generated an equidistant mesh: skipping ghost nodes test"
        ! delete the grid we created for this subroutine
        call reset_tree(params, .true., tree_ID=tree_ID)
        return
    endif

    !---------------------------------------------------------------------------
    ! Step 2: Actual testing of ghost node routines
    !---------------------------------------------------------------------------
    ! the entire test procedure is repeated for a bunch of frequencies, which is
    ! equivalent to using different block sizes, but way easier to program.
    ! These frequencies start at 1 and go to 2**(JmaxA - JminA)

    ! loop over frequencies
    do ifrequ = 1 , JmaxA - JminA + 1
        !-----------------------------------------------------------------------
        ! Fill the above constructed grid with the exact solution values
        !-----------------------------------------------------------------------
        ! loop over all active blocks
        do k = 1, hvy_n(tree_ID)
            ! hvy_id of the block we're looking at
            hvy_id = hvy_active(k, tree_ID)

            ! light id of this block
            call hvy2lgt( lgt_id, hvy_id, rank, params%number_blocks )

            ! compute block spacing and origin from treecode
            call get_block_spacing_origin( params, lgt_id, xx0, ddx )

            ! calculate f(x,y,z) for first datafield
            if ( params%dim == 3 ) then
                ! 3D:
                do iz = 1, Bs(3)+2*g
                    z = real(iz-(g+1), kind=rk) * ddx(3) + xx0(3)
                    do iy = 1, Bs(2)+2*g
                        y = real(iy-(g+1), kind=rk) * ddx(2) + xx0(2)
                        do ix = 1, Bs(1)+2*g
                            x = real(ix-(g+1), kind=rk) * ddx(1) + xx0(1)

                            ! use cos functions because theyre symmetric (symmetry BC)
                            hvy_block(ix, iy, iz, 1, hvy_id) &
                            = cos(2.0_rk**dble(ifrequ - 1)*x/Lx * 2.0_rk*pi) &
                            * cos(2.0_rk**dble(ifrequ - 1)*y/Ly * 2.0_rk*pi) &
                            * cos(2.0_rk**dble(ifrequ - 1)*z/Lz * 2.0_rk*pi)
                        enddo
                    enddo
                enddo
            else
                ! 2D:
                do iy = 1, Bs(2)+2*g
                    y = real(iy-(g+1), kind=rk) * ddx(2) + xx0(2)
                    do ix = 1, Bs(1)+2*g
                        x = real(ix-(g+1), kind=rk) * ddx(1) + xx0(1)

                        ! use cos functions because theyre symmetric (symmetry BC)
                        hvy_block(ix, iy, 1, 1, hvy_id) &
                        = cos(2.0_rk**dble(ifrequ - 1)*x/Lx * 2.0_rk*pi) &
                        * cos(2.0_rk**dble(ifrequ - 1)*y/Ly * 2.0_rk*pi)
                    enddo
                enddo
            end if

            ! now the entire block (incl ghost nodes) holds the exact solution: make a
            ! copy of the block for later comparison, but use work arrays usually used for RK4 substages
            ! so no additional memory is used.
            hvy_work(:,:,:,1,hvy_id,1) = hvy_block(:,:,:,1,hvy_id)

            ! in some rare cases the fact that we have now in fact filled the ghost
            ! nodes correctly introduces a bias: it means, we can now interpolate, even
            ! without the ghost nodes on the coarse block filled. hence, reset the
            ! ghost nodes of the testing blocks:
            !-- x-direction
            hvy_block(1:g, :, :, 1, hvy_id)           = 9.0e9_rk
            hvy_block(Bs(1)+g+1:Bs(1)+2*g, :, :, 1, hvy_id) = 9.0e9_rk
            !-- y-direction
            hvy_block(:, 1:g, :, 1, hvy_id)           = 9.0e9_rk
            hvy_block(:, Bs(2)+g+1:Bs(2)+2*g, :, 1, hvy_id) = 9.0e9_rk
            !-- z-direction
            if ( params%dim == 3 ) then
                hvy_block(:, :, 1:g, 1, hvy_id)           = 9.0e9_rk
                hvy_block(:, :, Bs(3)+g+1:Bs(3)+2*g, 1, hvy_id) = 9.0e9_rk
            end if
        end do


        !-----------------------------------------------------------------------
        ! synchronize ghost nodes (this is what we test here)
        !-----------------------------------------------------------------------
        call sync_ghosts_tree( params, hvy_block, tree_ID )

        ! verbose, only save for first iteration
        if (apply_verbose .and. ifrequ == 1) then
            call saveHDF5_tree("grid_0000.h5", 0.0_rk, 1, 1, params, hvy_block, tree_ID, no_sync=.true., save_ghosts=.true.)
        endif

        !-----------------------------------------------------------------------
        ! compute error (normalized, global, 2-norm)
        !-----------------------------------------------------------------------
        ! reset error
        error_L2 = 0.0_rk
        error_Linfty = 0.0_rk
        norm_L2 = 0.0_rk
        norm_Linfty = 0.0_rk

        ! loop over all active blocks and compute their error
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            error_L2     = error_L2 + sum( (hvy_block(:,:,:,1,hvy_id)-hvy_work(:,:,:,1,hvy_id,1))**2 )
            error_Linfty = max( error_Linfty, maxval(abs(hvy_block(:,:,:,1,hvy_id)-hvy_work(:,:,:,1,hvy_id,1))) )

            norm_L2     = norm_L2 + sum( (hvy_work(:,:,:,1,hvy_id,1))**2 )
            norm_Linfty = max( error_Linfty, maxval(abs(hvy_work(:,:,:,1,hvy_id,1))) )
        end do


        call MPI_Allreduce(error_L2, error1(ifrequ), 1, MPI_REAL8, MPI_SUM, WABBIT_COMM, ierr)
        error1(ifrequ) = sqrt(error1(ifrequ))

        call MPI_Allreduce(MPI_IN_PLACE, norm_L2, 1, MPI_REAL8, MPI_SUM, WABBIT_COMM, ierr)
        norm_L2 = sqrt(norm_L2)

        error1(ifrequ) = error1(ifrequ) / norm_L2

        call MPI_Allreduce(error_Linfty, error2(ifrequ), 1, MPI_REAL8, MPI_MAX, WABBIT_COMM, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, norm_Linfty, 1, MPI_REAL8, MPI_MAX, WABBIT_COMM, ierr)

        error2(ifrequ) = error2(ifrequ) / norm_L2

        ! output
        if (rank==0) then
            write(*,'("UNIT TEST: ghost nodes sync error_L2= ",es11.4,", error_L∞= ",es11.4,", frequ=",f5.1)')  &
            error1(ifrequ), error2(ifrequ), 2.0_rk**dble(ifrequ - 1)
        end if
    end do

    if (rank==0) then
        write(*,'(20("_/¯\"))')
        write(*,'(A)') "        .-."
        write(*,'(A)') "       ( o )"
        write(*,'(A)') "    /\_.' '._/\"
        write(*,'(A)') "    |         |"
        write(*,'(A)') "     \       /"
        write(*,'(A)') "      \    /`"
        write(*,'(A)') "    (__)  /"
        write(*,'(A)') "    `.__.'"
        write(*,'("UNIT TEST:", " L2 convergence order: ",32(f7.3,2x))') sqrt(error1(2:JmaxA-JminA+1) / error1(1:JmaxA-JminA))
        write(*,'("UNIT TEST:", " L2 mean conv.  order: ",f7.3)') sum(sqrt(error1(2:JmaxA-JminA+1) / error1(1:JmaxA-JminA))) / dble(JmaxA-JminA)
        write(*,'("UNIT TEST:", " L∞ convergence order: ",32(f7.3,2x))') sqrt(error2(2:JmaxA-JminA+1) / error2(1:JmaxA-JminA))
        write(*,'("UNIT TEST:", " L∞ mean conv.  order: ",f7.3)') sum(sqrt(error2(2:JmaxA-JminA+1) / error2(1:JmaxA-JminA))) / dble(JmaxA-JminA)
        write(*,'(20("¯\_/"))')

    endif

    if (abort_on_fail) then
        select case(params%order_predictor)
        case("multiresolution_2nd")
            if ((sum(sqrt(error1(2:JmaxA-JminA+1) / error1(1:JmaxA-JminA))) / dble(JmaxA-JminA)) < 1.50_rk) then
                call abort(70820231, "2nd order convergence not satisfied")
            endif
        case("multiresolution_4th")
            if ((sum(sqrt(error1(2:JmaxA-JminA+1) / error1(1:JmaxA-JminA))) / dble(JmaxA-JminA)) < 3.50_rk) then
                call abort(70820231, "4th order convergence not satisfied")
            endif
        case("multiresolution_6th")
            if ((sum(sqrt(error1(2:JmaxA-JminA+1) / error1(1:JmaxA-JminA))) / dble(JmaxA-JminA)) < 5.50_rk) then
                call abort(70820231, "6th order convergence not satisfied")
            endif
        case("multiresolution_8th")
            if ((sum(sqrt(error1(2:JmaxA-JminA+1) / error1(1:JmaxA-JminA))) / dble(JmaxA-JminA)) < 7.50_rk) then
                call abort(70820231, "8th order convergence not satisfied")
            endif
        case("multiresolution_10th")
            if ((sum(sqrt(error1(2:JmaxA-JminA+1) / error1(1:JmaxA-JminA))) / dble(JmaxA-JminA)) < 9.50_rk) then
                call abort(70820231, "10th order convergence not satisfied")
            endif
        case("multiresolution_12th")
            if ((sum(sqrt(error1(2:JmaxA-JminA+1) / error1(1:JmaxA-JminA))) / dble(JmaxA-JminA)) < 11.50_rk) then
                call abort(70820231, "12th order convergence not satisfied")
            endif
        case default
            call abort(251103, "I do not know this predictor order: " // params%order_predictor)
        end select
    endif

    ! delete the grid we created for this subroutine
    call reset_tree(params, .true., tree_ID=tree_ID)
end subroutine
