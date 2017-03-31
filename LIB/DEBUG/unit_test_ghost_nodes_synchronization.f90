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

subroutine unit_test_ghost_nodes_synchronization( params )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)          :: params

    ! local copy of params struct
    ! user defined parameter structure
    type (type_params)                      :: params_loc

    ! light data array
    integer(kind=ik), allocatable           :: lgt_block_loc(:, :)
    ! heavy data array - block data
    real(kind=rk), allocatable              :: hvy_block_loc(:, :, :, :, :), hvy_block_loc_exact(:, :, :, :, :)

    ! list of active blocks (heavy data)
    integer(kind=ik), allocatable           :: hvy_active(:)
    ! number of active blocks (heavy data)
    integer(kind=ik)                        :: hvy_n

    ! list of active blocks (light data)
    integer(kind=ik), allocatable           :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik)                        :: lgt_n

    ! neighbor array
    integer(kind=ik), allocatable           :: hvy_neighbor_loc(:,:)

    ! allocation error variable
    integer(kind=ik)                        :: allocate_error

    ! loop variables
    integer(kind=ik)                        :: k, l, lgt_id, hvy_id

    ! process rank
    integer(kind=ik)                        :: rank

    ! coordinates vectors
    real(kind=rk), allocatable              :: coord_x(:), coord_y(:), coord_z(:)
    ! spacing
    real(kind=rk)                           :: ddx(1:3), xx0(1:3)

    ! grid parameter
    integer(kind=ik)                        :: Bs, g, Ds
    real(kind=rk)                           :: Lx, Ly, Lz
    ! data dimensionality
    integer(kind=ik)                        :: d
    ! frequency of sin functions for testing:
    real(kind=rk)                           :: frequ(1:6)
    integer(kind=ik)                        :: ifrequ

    ! error variable
    real(kind=rk)                           :: error(1:6), my_error, norm, my_norm
    ! MPI error variable
    integer(kind=ik)                        :: ierr

    ! chance for block refinement, random number
    real(kind=rk)                           :: ref_chance, r

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization


    ! local copy of params struct
    params_loc = params
    ! this ugly little flag supresses some verbosity
    params_loc%unit_test = .true.

    ! set MPI parameters
    rank = params_loc%rank

    ! grid parameter
    Lx = params_loc%Lx
    Ly = params_loc%Ly
    Lz = params_loc%Lz

    ! set data dimension
    if ( params_loc%threeD_case ) then
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

    ! choose bocksize and domainsize here:
    if ( params_loc%threeD_case ) then
        params_loc%number_block_nodes = 17
        params_loc%number_domain_nodes = 65
    else
        params_loc%number_block_nodes = 33
        params_loc%number_domain_nodes = 513
    endif

    Bs = params_loc%number_block_nodes
    g  = params_loc%number_ghost_nodes
    Ds = params_loc%number_domain_nodes

    ! determine the required number of blocks, given the current block
    ! size and the desired "full resolution" size "Ds"
    params_loc%number_blocks = ( ((Ds-1) / (Bs-1))**d ) * 2**d

    if (rank == 0) then
      write(*,'("UNIT TEST: testing Bs=",i4," blocks-per-mpirank=",i5)')  Bs, params_loc%number_blocks
    end if

    !---------------------------------------------------------------------------
    ! Step 1: Construct a random grid for testing. Note we keep this grid
    ! and perform the same test for differnet frequencies (= resolutions) only on
    ! this one grid.
    !---------------------------------------------------------------------------


    ! first: initializing new block data
    ! allocate block_list
    call allocate_block_list( params_loc, lgt_block_loc )
    ! allocate heavy data
    call allocate_block_data( params_loc, hvy_block_loc )
    call allocate_block_data( params_loc, hvy_block_loc_exact )

    ! active heavy block list
    allocate( hvy_active( params_loc%number_blocks ), stat=allocate_error )
    if ( allocate_error /= 0 ) call error_msg("ERROR: memory allocation fails")

    ! create active light block list
    allocate( lgt_active( size(lgt_block_loc, 1) ), stat=allocate_error )
    if ( allocate_error /= 0 ) call error_msg("ERROR: memory allocation fails")


    ! set all blocks to free (since if we call inicond twice, all blocks are used in the second call)
    lgt_active = -1; lgt_N = 0
    hvy_active = -1; hvy_N = 0

    ! init data array with zeros => after this: blocks have correct coordinate vectors
    call inicond_zeros( params_loc, lgt_block_loc, hvy_block_loc )

    ! update lists of active blocks (light and heavy data)
    call create_hvy_active_list( lgt_block_loc, hvy_active, hvy_n )
    call create_lgt_active_list( lgt_block_loc, lgt_active, lgt_n )

    ! init neighbor data array
    ! 2D: maximal 16 neighbors per block
    ! 3D: maximal 74 neighbors per block
    if ( params%threeD_case ) then
        ! 3D:
        allocate( hvy_neighbor_loc( params_loc%number_blocks, 74 ), stat=allocate_error )
        if ( allocate_error /= 0 ) call error_msg("ERROR: memory allocation fails")
    else
        ! 2D:
        allocate( hvy_neighbor_loc( params_loc%number_blocks, 16 ), stat=allocate_error )
        if ( allocate_error /= 0 ) call error_msg("ERROR: memory allocation fails")
    end if

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
            lgt_block_loc( lgt_active(k), params_loc%max_treelevel+2 ) = 1
        end if
    end do
    ! refine blocks
    ! check if block has reached maximal level
    call respect_min_max_treelevel( params_loc, lgt_block_loc, lgt_active, lgt_n )

    ! interpolate the new mesh
    if ( params%threeD_case ) then
        ! 3D:
        call refine_mesh_3D( params_loc, lgt_block_loc, hvy_block_loc, hvy_active, hvy_n )
    else
        ! 2D:
        call refine_mesh_2D( params_loc, lgt_block_loc, hvy_block_loc(:,:,1,:,:), hvy_active, hvy_n )
    end if

    ! update lists of active blocks (light and heavy data)
    call create_lgt_active_list( lgt_block_loc, lgt_active, lgt_n )
    call create_hvy_active_list( lgt_block_loc, hvy_active, hvy_n )

    !---------------------------------------------------------------------------------------------
    ! third: fill blocks with data

    ! allocate coord arrays
    allocate(coord_x( Bs + 2*g ),  stat=allocate_error )
    !if ( allocate_error /= 0 ) call error_msg("ERROR: memory allocation fails")
    if ( allocate_error /= 0 ) then
        write(*,*) "ERROR: memory allocation fails"
        write(*,*) "Code stops now."
        stop
    end if

    allocate(coord_y( Bs + 2*g ),  stat=allocate_error )
    !if ( allocate_error /= 0 ) call error_msg("ERROR: memory allocation fails")
    if ( allocate_error /= 0 ) then
        write(*,*) "ERROR: memory allocation fails"
        write(*,*) "Code stops now."
        stop
    end if

    allocate(coord_z( Bs + 2*g ),  stat=allocate_error )
    !if ( allocate_error /= 0 ) call error_msg("ERROR: memory allocation fails")
    if ( allocate_error /= 0 ) then
        write(*,*) "ERROR: memory allocation fails"
        write(*,*) "Code stops now."
        stop
    end if

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
          call hvy_id_to_lgt_id( lgt_id, hvy_id, rank, params_loc%number_blocks )
          ! compute block spacing and origin from treecode
          call get_block_spacing_origin( params_loc, lgt_id, lgt_block_loc, xx0, ddx )

          ! fill coordinate arrays, of course including ghost nodes
          do l = 1, Bs+2*g
            coord_x(l) = real(l-(g+1), kind=rk) * ddx(1) + xx0(1)
            coord_y(l) = real(l-(g+1), kind=rk) * ddx(2) + xx0(2)
            coord_z(l) = real(l-(g+1), kind=rk) * ddx(3) + xx0(3)
          enddo

          ! calculate f(x,y,z) for first datafield
          if ( params_loc%threeD_case ) then
            ! 3D:
            call f_xyz_3D( coord_x, coord_y, coord_z, hvy_block_loc(:, :, :, 2, hvy_active(k)), Bs, g, Lx, Ly, Lz, frequ(ifrequ) )
          else
            ! 2D:
            call f_xy_2D( coord_x, coord_y, hvy_block_loc(:, :, 1, 2, hvy_active(k)), Bs, g, Lx, Ly, frequ(ifrequ)  )
          end if

        end do

        ! now the entire grid (inlc ghost nodes) holds the exact solution: make a
        ! copy of the grid for later comparison (to be optimized as it consumes memory)
        hvy_block_loc_exact = hvy_block_loc

        !-----------------------------------------------------------------------
        ! synchronize ghost nodes (this is what we test here)
        !-----------------------------------------------------------------------
        ! update neighbors
        call update_neighbors( params_loc, lgt_block_loc, hvy_neighbor_loc, lgt_active, lgt_n, hvy_active, hvy_n )

        ! synchronize ghost nodes
        call synchronize_ghosts( params_loc, lgt_block_loc, hvy_block_loc, hvy_neighbor_loc, hvy_active, hvy_n )

        !-----------------------------------------------------------------------
        ! compute error (normalized, global, 2-norm)
        !-----------------------------------------------------------------------
        ! reset error
        my_error = 0.0_rk
        my_norm = 0.0_rk

        ! loop over all active blocks
        do k = 1, hvy_n
          my_error = my_error + sqrt(sum((hvy_block_loc(:,:,:,2,hvy_active(k))-hvy_block_loc_exact(:,:,:,2,hvy_active(k)))**2 ))
          my_norm = my_norm  + sqrt(sum((hvy_block_loc_exact(:,:,:,2,hvy_active(k)))**2 ))
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
    deallocate(coord_x, coord_y, coord_z, hvy_neighbor_loc)
    deallocate(lgt_block_loc, lgt_active, stat=allocate_error )
    deallocate(hvy_block_loc, stat=allocate_error )
    deallocate(hvy_block_loc_exact, stat=allocate_error )
    deallocate( hvy_active, stat=allocate_error )
end subroutine unit_test_ghost_nodes_synchronization
