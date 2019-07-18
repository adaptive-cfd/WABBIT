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
    ! coordinates vectors
    real(kind=rk), allocatable              :: coord_x(:), coord_y(:), coord_z(:)
    ! spacing
    real(kind=rk)                           :: ddx(1:3), xx0(1:3)

    ! grid parameter
    integer(kind=ik)                        :: g, number_blocks
    integer(kind=ik), dimension(3)          :: Bs
    real(kind=rk)                           :: Lx, Ly, Lz
    ! data dimensionality
    integer(kind=ik)                        :: d, dF, max_neighbors
    ! frequency of sin functions for testing:
    real(kind=rk)                           :: frequ(1:6)
    integer(kind=ik)                        :: ifrequ

    ! error variable
    real(kind=rk)                           :: error(1:6), my_error, norm, my_norm
    ! MPI error variable
    integer(kind=ik)                        :: ierr
    logical::test

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank = params%rank

    ! grid parameter
    Lx = params%domain_size(1)
    Ly = params%domain_size(2)
    Lz = params%domain_size(3)

    d = params%dim
    ! set data dimension
    if ( params%dim == 3 ) then
        max_neighbors = 74
    else
        max_neighbors = 12
    endif

!---------------------------------------------------------------------------------------------
! main body

    if (rank == 0) then
        write(*,'(80("_"))')
        write(*,'("UNIT TEST: Beginning ghost nodes test")')
    end if

    Bs = params%Bs
    g  = params%n_ghosts
    dF = params%n_eqn
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
    l = 5
    call create_random_grid( params, lgt_block, hvy_block, hvy_work, hvy_tmp, hvy_neighbor, lgt_active, &
        lgt_n, lgt_sortednumlist, hvy_active, hvy_n, 2, .true., l, tree_ID_flow )

    if (params%rank == 0) then
        write(*,'(80("-"))')
        write(*,'("UNIT TEST: performed ",i2," randomized refinement and coarsening steps")') l
        write(*,'(" done creating a random grid N_blocks=",i5, " Jmax=", i2)') lgt_n, max_active_level( lgt_block, lgt_active, lgt_n )
        write(*,'(" ready for testing.")')
    endif

    !---------------------------------------------------------------------------
    ! Step 2: Actual testing of ghost node routines
    !---------------------------------------------------------------------------
    ! the entire test procedure is repeated for a bunch of frequencies, which is
    ! equivalent to using different block sizes, but way easier to program.
    ! These frequencies are tested:
    frequ=(/1.0_rk , 2.0_rk, 4.0_rk, 8.0_rk, 16.0_rk, 32.0_rk/)
    allocate( coord_x( Bs(1) + 2*g ), coord_y( Bs(2) + 2*g ), coord_z( Bs(3) + 2*g ) )

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

          ! fill coordinate arrays, of course including ghost nodes
          do l = 1, Bs(1)+2*g
            coord_x(l) = real(l-(g+1), kind=rk) * ddx(1) + xx0(1)
          enddo
          do l = 1, Bs(2)+2*g
            coord_y(l) = real(l-(g+1), kind=rk) * ddx(2) + xx0(2)
          enddo
          do l = 1, Bs(3)+2*g
            coord_z(l) = real(l-(g+1), kind=rk) * ddx(3) + xx0(3)
          enddo

          ! calculate f(x,y,z) for first datafield
          if ( params%dim == 3 ) then
            ! 3D:
            call f_xyz_3D( coord_x, coord_y, coord_z, hvy_block(:, :, :, 1, hvy_id), Bs, g, Lx, Ly, Lz, frequ(ifrequ) )
          else
            ! 2D:
            call f_xy_2D( coord_x, coord_y, hvy_block(:, :, 1, 1, hvy_id), Bs, g, Lx, Ly, frequ(ifrequ)  )
          end if

        end do

        ! now the entire grid (incl ghost nodes) holds the exact solution: make a
        ! copy of the grid for later comparison, but use work arrays usually used for RK4 substages
        ! so no additional memory is used.
        hvy_work(:,:,:,1,:,1) = hvy_block(:,:,:,1,:)

        !-----------------------------------------------------------------------
        ! synchronize ghost nodes (this is what we test here)
        !-----------------------------------------------------------------------
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

        !-----------------------------------------------------------------------
        ! compute error (normalized, global, 2-norm)
        !-----------------------------------------------------------------------
        ! reset error
        my_error = 0.0_rk
        my_norm = 0.0_rk

        ! loop over all active blocks and compute their error
        do k = 1, hvy_n
          my_error = my_error + sqrt(sum((hvy_block(:,:,:,1,hvy_active(k))-hvy_work(:,:,:,1,hvy_active(k),1))**2 ))
          my_norm = my_norm  + sqrt(sum((hvy_work(:,:,:,1,hvy_active(k),1))**2 ))
        end do

        ! synchronize errors
        call MPI_Allreduce(my_error, error(ifrequ), 1, MPI_REAL8, MPI_SUM, WABBIT_COMM, ierr)
        call MPI_Allreduce(my_norm, norm, 1, MPI_REAL8, MPI_SUM, WABBIT_COMM, ierr)

        error(ifrequ) = error(ifrequ) / norm

        ! output
        if (rank==0) then
            write(*,'(" done - ghost nodes synchronization error = ",es16.8," frequ=",g12.4)')  error(ifrequ), frequ(ifrequ)
        end if
      end do

    if (rank==0) then
      write(*,'(" done - convergence order was ",6(g12.4,1x))')  sqrt(error(2:6) / error(1:5))
      write(*,'(" done - mean convergence order was ",g12.4)')  sum(sqrt(error(2:6) / error(1:5))) / 5.0_rk
    endif

    !---------------------------------------------------------------------------------------------
    ! last: clean up
    deallocate(coord_x, coord_y, coord_z)

end subroutine unit_test_ghost_nodes_synchronization
