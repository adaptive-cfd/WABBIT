!> \file
!> \brief  unit test for ghost nodes synchronization
!> \version 0.5
!> \author msr, engels
!> \note input only params struct to this subroutine
!!       create new light/heavy data arrays here and deallocate them after this function
!! \details
!! \date 21/01/17 - create
!! \date 03/04/17 - major rewrite: no local memory allocation, convergence test is performed \n
!! \date 05/04/17 - use the renewed refine_mesh with random indicator
!
! ********************************************************************************************

subroutine unit_test_ghost_nodes_synchronization( params, lgt_block, hvy_block, hvy_work, hvy_tmp, &
    hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist, hvy_n, lgt_n )


    implicit none
    !> user defined parameter structure
    type (type_params), intent(inout)       :: params
    !> light data array
    integer(kind=ik),  intent(inout)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(out)              :: hvy_tmp(:, :, :, :, :)
    !> heavy work array: used for RHS evaluation in multistep methods (like RK4: u0, k1, k2 etc)
    real(kind=rk), intent(out)              :: hvy_work(:, :, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik),  intent(inout)        :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: lgt_active(:)
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: hvy_active(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)      :: lgt_sortednumlist(:,:)

    ! number of active blocks (heavy data)
    integer(kind=ik), intent(inout)         :: hvy_n
    ! number of active blocks (light data)
    integer(kind=ik), intent(inout)         :: lgt_n
    ! loop variables
    integer(kind=ik)                        :: k, l, lgt_id, hvy_id
    ! process rank
    integer(kind=ik)                        :: rank, number_procs
    ! spacing
    real(kind=rk)                           :: ddx(1:3), xx0(1:3)

    ! grid parameter
    integer(kind=ik)                        :: g, number_blocks, ix, iy, iz
    integer(kind=ik), dimension(3)          :: Bs
    real(kind=rk)                           :: Lx, Ly, Lz, x, y, z
    ! data dimensionality
    integer(kind=ik)                        :: d,  max_neighbors
    ! frequency of sin functions for testing:
    real(kind=rk)                           :: frequ(1:6)
    integer(kind=ik)                        :: ifrequ

    ! error variable
    real(kind=rk)                           :: error2(1:6), error1(1:6), error_L2, error_Linfty, norm_L2, norm_Linfty
    ! MPI error variable
    integer(kind=ik)                        :: ierr
    logical::test

    rank = params%rank
    Lx = params%domain_size(1)
    Ly = params%domain_size(2)
    Lz = params%domain_size(3)

    d = params%dim
    if ( params%dim == 3 ) then
        max_neighbors = 74
    else
        max_neighbors = 12
    endif


    if (rank == 0) then
        write(*,'(80("_"))')
        write(*,'("UNIT TEST: Beginning ghost nodes test")')
    end if

    Bs = params%Bs
    g  = params%n_ghosts
    number_procs  = params%number_procs
    number_blocks = params%number_blocks

    if (rank == 0) then
        write(*,'("UNIT TEST: testing Bs=",i4," x ",i4," x ",i4," blocks-per-mpirank=",i5)')  Bs(1),Bs(2),Bs(3), params%number_blocks
    end if

    !---------------------------------------------------------------------------
    ! Step 1: Construct a random grid for testing. Note we keep this grid
    ! and perform the same test for differnet frequencies (= resolutions) only on
    ! this one grid.
    !---------------------------------------------------------------------------
    ! this parameter controls roughly how dense the random grid is, i.e., in % of the
    ! complete memory.
    params%max_grid_density = 0.10_rk
    ! perform 5 iterations of random refinement/coarsening
    l = 5
    call create_random_grid( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n, 2, .true., l, tree_ID_flow )

    if (params%rank == 0) then
        write(*,'(80("-"))')
        write(*,'("UNIT TEST: performed ",i2," randomized refinement and coarsening steps")') l
        write(*,'(" done creating a random grid N_blocks=",i5, " Jmax=", i2)') lgt_n, max_active_level( lgt_block, lgt_active, lgt_n )
        write(*,'(" ready for testing.")')
    endif


    if (max_active_level( lgt_block, lgt_active, lgt_n ) == min_active_level( lgt_block, lgt_active, lgt_n )) then
        if (params%rank==0) write(*,*) "By chance, generated an equidistant mesh: skipping ghost nodes test"
        return
    endif

    !---------------------------------------------------------------------------
    ! Step 2: Actual testing of ghost node routines
    !---------------------------------------------------------------------------
    ! the entire test procedure is repeated for a bunch of frequencies, which is
    ! equivalent to using different block sizes, but way easier to program.
    ! These frequencies are tested:
    frequ=(/1.0_rk , 2.0_rk, 4.0_rk, 8.0_rk, 16.0_rk, 32.0_rk/)

    ! loop over frequencies
    do ifrequ = 1 , size(frequ)
        !-----------------------------------------------------------------------
        ! Fill the above constructed grid with the exact solution values
        !-----------------------------------------------------------------------
        ! loop over all active blocks
        do k = 1, hvy_n
            ! hvy_id of the block we're looking at
            hvy_id = hvy_active(k)
            ! light id of this block
            call hvy_id_to_lgt_id( lgt_id, hvy_id, rank, params%number_blocks )
            ! compute block spacing and origin from treecode
            call get_block_spacing_origin( params, lgt_id, lgt_block, xx0, ddx )

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
                            = cos(frequ(ifrequ)*x/Lx * 2.0_rk*pi) &
                            * cos(frequ(ifrequ)*y/Ly * 2.0_rk*pi) &
                            * cos(frequ(ifrequ)*z/Lz * 2.0_rk*pi)
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
                        = cos(frequ(ifrequ)*x/Lx * 2.0_rk*pi) &
                        * cos(frequ(ifrequ)*y/Ly * 2.0_rk*pi)
                    enddo
                enddo
            end if

            ! now the entire block (incl ghost nodes) holds the exact solution: make a
            ! copy of the block for later comparison, but use work arrays usually used for RK4 substages
            ! so no additional memory is used.
            hvy_work(:,:,:,1,hvy_id,1) = hvy_block(:,:,:,1,hvy_id)
        end do


        !-----------------------------------------------------------------------
        ! synchronize ghost nodes (this is what we test here)
        !-----------------------------------------------------------------------
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

        !-----------------------------------------------------------------------
        ! compute error (normalized, global, 2-norm)
        !-----------------------------------------------------------------------
        ! reset error
        error_L2 = 0.0_rk
        error_Linfty = 0.0_rk
        norm_L2 = 0.0_rk
        norm_Linfty = 0.0_rk

        ! loop over all active blocks and compute their error
        do k = 1, hvy_n
            error_L2     = error_L2 + sum((hvy_block(:,:,:,1,hvy_active(k))-hvy_work(:,:,:,1,hvy_active(k),1))**2)
            error_Linfty = max( error_Linfty, maxval(abs(hvy_block(:,:,:,1,hvy_active(k))-hvy_work(:,:,:,1,hvy_active(k),1))) )

            norm_L2     = norm_L2 + sum( (hvy_work(:,:,:,1,hvy_active(k),1))**2 )
            norm_Linfty = max( error_Linfty, maxval(abs(hvy_work(:,:,:,1,hvy_active(k),1))) )
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
            write(*,'(" done - ghost nodes synchronization error_L2 = ",es16.8," error_Linfty=",es16.8," frequ=",g12.4)')  &
            error1(ifrequ), error2(ifrequ), frequ(ifrequ)
        end if
    end do

    if (rank==0) then
        write(*,'(" done - L2 convergence order was ",6(g12.4,1x))')  sqrt(error1(2:6) / error1(1:5))
        write(*,'(" done - L2 mean convergence order was ",g12.4)')  sum(sqrt(error1(2:6) / error1(1:5))) / 5.0_rk
        write(*,'(" done - Linfty convergence order was ",6(g12.4,1x))')  sqrt(error2(2:6) / error2(1:5))
        write(*,'(" done - Linfty mean convergence order was ",g12.4)')  sum(sqrt(error2(2:6) / error2(1:5))) / 5.0_rk
    endif

end subroutine unit_test_ghost_nodes_synchronization
