! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: unit_test_ghost_nodes_synchronization.f90
! version: 0.5
! author: msr
!
! unit test for ghost nodes synchronization
! note: input only params struct to this subroutine
!       create new light/heavy data arrays here and deallocate them after this function
!
! input:    - params
! output:   -
!
! = log ======================================================================================
!
! 21/01/17 - create
!
! ********************************************************************************************

subroutine unit_test_ghost_nodes_synchronization( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    ! user defined parameter structure
    type (type_params), intent(inout)     :: params
    ! light data array
    integer(kind=ik),  intent(inout)      :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk),  intent(inout)         :: hvy_block(:, :, :, :, :)
    ! heavy work array  )
    real(kind=rk),  intent(inout)         :: hvy_work (:, :, :, :, :)
    ! neighbor array (heavy data)
    integer(kind=ik),  intent(inout)      :: hvy_neighbor(:,:)
    ! list of active blocks (light data)
    integer(kind=ik),  intent(inout)      :: lgt_active(:)
    ! list of active blocks (light data)
    integer(kind=ik),  intent(inout)      :: hvy_active(:)

    ! number of active blocks (heavy data)
    integer(kind=ik)                        :: hvy_n
    ! number of active blocks (light data)
    integer(kind=ik)                        :: lgt_n

    ! loop variables
    integer(kind=ik)                        :: k, l, lgt_id, hvy_id

    ! process rank
    integer(kind=ik)                        :: rank, number_procs

    ! coordinates vectors
    real(kind=rk), allocatable              :: coord_x(:), coord_y(:), coord_z(:)
    ! spacing
    real(kind=rk)                           :: ddx(1:3), xx0(1:3)

    ! grid parameter
    integer(kind=ik)                        :: Bs, g, Ds, number_blocks
    real(kind=rk)                           :: Lx, Ly, Lz
    ! data dimensionality
    integer(kind=ik)                        :: d, dF
    ! frequency of sin functions for testing:
    real(kind=rk)                           :: frequ(1:6)
    integer(kind=ik)                        :: ifrequ

    ! error variable
    real(kind=rk)                           :: error(1:6), my_error, norm, my_norm
    ! MPI error variable
    integer(kind=ik)                        :: ierr, allocate_error

    ! chance for block refinement, random number
    real(kind=rk)                           :: ref_chance, r

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank = params%rank

    ! grid parameter
    Lx = params%Lx
    Ly = params%Ly
    Lz = params%Lz

    ! set data dimension
    if ( params%threeD_case ) then
        d = 3
    else
        d = 2
    endif

!---------------------------------------------------------------------------------------------
! main body

    if (rank == 0) then
        write(*,'(80("_"))')
        write(*,'("UNIT TEST: Beginning ghost nodes test")')
    end if

    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes
    Ds = params%number_domain_nodes
    dF = params%number_data_fields
    number_procs  = params%number_procs
    number_blocks = params%number_blocks

    if (rank == 0) then
      write(*,'("UNIT TEST: testing Bs=",i4," blocks-per-mpirank=",i5)')  Bs, params%number_blocks
    end if

    !---------------------------------------------------------------------------
    ! Step 1: Construct a random grid for testing. Note we keep this grid
    ! and perform the same test for differnet frequencies (= resolutions) only on
    ! this one grid.
    !---------------------------------------------------------------------------
    ! allocate coord arrays
    allocate(coord_x( Bs + 2*g ),  stat=allocate_error )
    if ( allocate_error /= 0 ) call error_msg("ERROR: memory allocation fails")

    allocate(coord_y( Bs + 2*g ),  stat=allocate_error )
    if ( allocate_error /= 0 ) call error_msg("ERROR: memory allocation fails")

    allocate(coord_z( Bs + 2*g ),  stat=allocate_error )
    if ( allocate_error /= 0 ) call error_msg("ERROR: memory allocation fails")

    ! set all blocks to free (since if we call inicond twice, all blocks are used in the second call)
    lgt_block = -1
    lgt_active = -1; lgt_N = 0
    hvy_active = -1; hvy_N = 0

    ! init data array with zeros => after this: blocks have correct coordinate vectors
    call inicond_zeros( params, lgt_block, hvy_block )

    ! update lists of active blocks (light and heavy data)
    call create_hvy_active_list( lgt_block, hvy_active, hvy_n )
    call create_lgt_active_list( lgt_block, lgt_active, lgt_n )


    !---------------------------------------------------------------------------------------------
    ! second: refine some blocks (random)

    ! refinement chance
    ref_chance = 0.25_rk

    ! set random seed
    call init_random_seed()

    ! set refinement status
    do k = 1, lgt_n
        ! random number
        call random_number(r)
        ! set refinement status
        if ( r <= ref_chance ) then
            lgt_block( lgt_active(k), params%max_treelevel+2 ) = 1
        end if
    end do
    ! refine blocks
    ! check if block has reached maximal level
    call respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n )

    ! interpolate the new mesh
    if ( params%threeD_case ) then
        ! 3D:
        call refine_mesh_3D( params, lgt_block, hvy_block, hvy_active, hvy_n )
    else
        ! 2D:
        call refine_mesh_2D( params, lgt_block, hvy_block(:,:,1,:,:), hvy_active, hvy_n )
    end if

    ! update lists of active blocks (light and heavy data)
    call create_lgt_active_list( lgt_block, lgt_active, lgt_n )
    call create_hvy_active_list( lgt_block, hvy_active, hvy_n )

    ! at this point now, the grid used for testing is ready.

    !---------------------------------------------------------------------------
    ! Step 2: Actual testing of ghost node routines
    !---------------------------------------------------------------------------
    ! the entire test procedure is repeated for a bunch of frequencies, which is
    ! equivalent to using different block sizes, but way easier to program.
    ! These frequencies are tested:
    frequ=(/1.0_rk, 2.0_rk, 4.0_rk, 8.0_rk, 16.0_rk, 32.0_rk/)

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
          do l = 1, Bs+2*g
            coord_x(l) = real(l-(g+1), kind=rk) * ddx(1) + xx0(1)
            coord_y(l) = real(l-(g+1), kind=rk) * ddx(2) + xx0(2)
            coord_z(l) = real(l-(g+1), kind=rk) * ddx(3) + xx0(3)
          enddo

          ! calculate f(x,y,z) for first datafield
          if ( params%threeD_case ) then
            ! 3D:
            call f_xyz_3D( coord_x, coord_y, coord_z, hvy_block(:, :, :, 2, hvy_active(k)), Bs, g, Lx, Ly, Lz, frequ(ifrequ) )
          else
            ! 2D:
            call f_xy_2D( coord_x, coord_y, hvy_block(:, :, 1, 2, hvy_active(k)), Bs, g, Lx, Ly, frequ(ifrequ)  )
          end if

        end do

        ! now the entire grid (incl ghost nodes) holds the exact solution: make a
        ! copy of the grid for later comparison, but use work arrays ususally used for RK4 substages
        ! so no additional memory is used.
        hvy_work(:,:,:,1,:) = hvy_block(:,:,:,2,:)

        !-----------------------------------------------------------------------
        ! synchronize ghost nodes (this is what we test here)
        !-----------------------------------------------------------------------
        ! update neighbors
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

        ! synchronize ghost nodes
        call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

        !-----------------------------------------------------------------------
        ! compute error (normalized, global, 2-norm)
        !-----------------------------------------------------------------------
        ! reset error
        my_error = 0.0_rk
        my_norm = 0.0_rk

        ! loop over all active blocks and compute their error
        do k = 1, hvy_n
          my_error = my_error + sqrt(sum((hvy_block(:,:,:,2,hvy_active(k))-hvy_work(:,:,:,1,hvy_active(k)))**2 ))
          my_norm = my_norm  + sqrt(sum((hvy_work(:,:,:,1,hvy_active(k)))**2 ))
        end do

        ! synchronize errors
        call MPI_Allreduce(my_error, error(ifrequ), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(my_norm, norm, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        error(ifrequ) = error(ifrequ) / norm

        ! output
        if (rank==0) then
            write(*,'("UNIT TEST DONE: ghost nodes synchronization error = ",es16.8," frequ=",g12.4)')  error(ifrequ), frequ(ifrequ)
        end if
      end do

    if (rank==0) then
      write(*,'("UNIT TEST DONE: convergence order was ",6(g12.4,1x))')  sqrt(error(2:6) / error(1:5))
      write(*,'("UNIT TEST DONE: mean convergence order was ",g12.4)')  sum(sqrt(error(2:6) / error(1:5))) / 5.0_rk
    endif

    !---------------------------------------------------------------------------------------------
    ! last: clean up
    deallocate(coord_x, coord_y, coord_z)
end subroutine unit_test_ghost_nodes_synchronization
