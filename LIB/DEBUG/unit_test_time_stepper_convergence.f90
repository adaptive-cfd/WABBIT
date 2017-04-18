! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: unit_test_time_stepper_convergence.f90
! version: 0.5
! author: msr
!
! unit test for ghost nodes synchronization
! note: input only params struct to this subroutine
!       create new light/heavy data arrays here and deallocate them after this function
!
! input:    - params, empty ligth and heavy data arrays
! output:   -
!
! = log ======================================================================================
!
! 18/04/17 - create
!
! ********************************************************************************************

subroutine unit_test_time_stepper_convergence( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    ! user defined parameter structure
    type (type_params), intent(inout)       :: params
    ! light data array
    integer(kind=ik),  intent(inout)        :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)
    ! heavy work array  )
    real(kind=rk),  intent(inout)           :: hvy_work (:, :, :, :, :)
    ! neighbor array (heavy data)
    integer(kind=ik),  intent(inout)        :: hvy_neighbor(:,:)
    ! list of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: lgt_active(:)
    ! list of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: hvy_active(:)

    ! local user defined parameter structure - use to change settings
    type (type_params)                      :: params_loc

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
    integer(kind=ik)                        :: ierr

! allocation error variabel
    integer(kind=ik)                                :: allocate_error

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! copy params data to local struct
    params_loc = params

    ! set MPI parameters
    rank = params%rank

!    ! grid parameter
!    Lx = params%Lx
!    Ly = params%Ly
!    Lz = params%Lz
!
!    ! set data dimension
!    if ( params%threeD_case ) then
!        d = 3
!    else
!        d = 2
!    endif

!---------------------------------------------------------------------------------------------
! main body

    if (rank == 0) then
        write(*,'(80("_"))')
        write(*,'("UNIT TEST: Beginning time stepper convergence test")')
    end if

    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes
!    Ds = params%number_domain_nodes
!    dF = params%number_data_fields
!    number_procs  = params%number_procs
!    number_blocks = params%number_blocks

    if (rank == 0) then
      write(*,'("UNIT TEST: testing Bs=",i4," blocks-per-mpirank=",i5)')  Bs, params%number_blocks
    end if

    !---------------------------------------------------------------------------
    ! Step 1: Construct the test grid. Note: grid size is fixed on finest mesh
    ! level for testing different time steps
    !---------------------------------------------------------------------------
    ! allocate coord arrays
    allocate( coord_x( Bs + 2*g ), coord_y( Bs + 2*g ), coord_z( Bs + 2*g ) )

    ! set all blocks to free (since if we call inicond twice, all blocks are used in the second call)
    lgt_block = -1
    lgt_active = -1; lgt_N = 0
    hvy_active = -1; hvy_N = 0

    ! setup the finest grid level with some data (we don't care what data, we'll erase it)
    ! Note that active lists + neighbor relations are updated inside this routine as well, as
    ! the grid is modified
    call create_equidistant_base_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, params%max_treelevel )


!    !---------------------------------------------------------------------------------------------
!    ! second: refine some blocks (random), coarsen some blocks (random)
!    do l = 1, 5
!      call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, "random" )
!      call adapt_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, "random" )
!    end do
!
!    if (params%rank == 0) then
!        write(*,'(80("-"))')
!        write(*,'("UNIT TEST: performed ",i2," randomized refinement and coarsening steps")') l
!        write(*,'(" done creating a random grid N_blocks=",i5, " Jmax=", i2)') lgt_n, maxval(lgt_block(:,params%max_treelevel+1))
!        write(*,'(" ready for testing.")')
!    endif
!
    !---------------------------------------------------------------------------
    ! Step 2: time stepper testing
    !---------------------------------------------------------------------------
    ! i)   set physics to convection-diffusion physics
    ! ii)  calculate minimal dt for given cfl number
    ! iii) fill grid with exact (periodic) solution
    ! iv)  do one time step (repeat with smaller steps)

    ! set local physics to convection diffusion with proper parameters
    if ( params_loc%threeD_case ) then
        ! 3D:
        params_loc%physics_type         = '3D_convection_diffusion'
        ! only one datafield
        params_loc%number_data_fields   = 1
        ! choose diffusion coefficient and convection velocity
        if (allocated(params_loc%physics%u0)) then
            params_loc%physics%nu = 0.0_rk
            params_loc%physics%u0 = (/1.0_rk, 0.5_rk, 0.25_rk/)
        else
            allocate(params%physics%u0(3), params%physics%nu(1))
            params_loc%physics%nu = 0.0_rk
            params_loc%physics%u0 = (/1.0_rk, 0.5_rk, 0.25_rk/)
        end if
    else
        ! 2D:
        params_loc%physics_type         = '2D_convection_diffusion'
        ! only one datafield
        params_loc%number_data_fields   = 1
        ! choose diffusion coefficient and convection velocity
        if (allocated(params_loc%physics%u0)) then
            params_loc%physics%nu = 0.0_rk
            params_loc%physics%u0 = (/1.0_rk, 0.5_rk/)
        else
            !allocate(params%physics%u0(2), source=0.0_rk, stat=allocate_error) !params%physics%nu(1))
            allocate(params%physics%u0(2), stat=allocate_error)
            print*, allocate_error
            print*, params_loc%physics%u0
            !params_loc%physics%nu = 0.0_rk
            params_loc%physics%u0 = (/1.0_rk, 0.5_rk/)
        end if
    end if


print*, params_loc%physics%u0

    ! set dt parameter
    !params_loc%dt_type = cfl
    !params_loc%cfl = 0.5

    ! calculate dt
    !call calc_dt

    !params_loc%dt_type = fixed
    !params_loc%cfl = dt








!    frequ=(/1.0_rk, 2.0_rk, 4.0_rk, 8.0_rk, 16.0_rk, 32.0_rk/)
!
!    ! loop over frequencies
!    do ifrequ = 1 , size(frequ)
!
!        !-----------------------------------------------------------------------
!        ! Fill the above constructed grid with the exact solution values
!        !-----------------------------------------------------------------------
!        ! loop over all active blocks
!        do k = 1, hvy_n
!          ! hvy_id of the block we're looking at
!          hvy_id = hvy_active(k)
!          ! light id of this block
!          call hvy_id_to_lgt_id( lgt_id, hvy_id, rank, params%number_blocks )
!          ! compute block spacing and origin from treecode
!          call get_block_spacing_origin( params, lgt_id, lgt_block, xx0, ddx )
!
!          ! fill coordinate arrays, of course including ghost nodes
!          do l = 1, Bs+2*g
!            coord_x(l) = real(l-(g+1), kind=rk) * ddx(1) + xx0(1)
!            coord_y(l) = real(l-(g+1), kind=rk) * ddx(2) + xx0(2)
!            coord_z(l) = real(l-(g+1), kind=rk) * ddx(3) + xx0(3)
!          enddo
!
!          ! calculate f(x,y,z) for first datafield
!          if ( params%threeD_case ) then
!            ! 3D:
!            call f_xyz_3D( coord_x, coord_y, coord_z, hvy_block(:, :, :, 2, hvy_active(k)), Bs, g, Lx, Ly, Lz, frequ(ifrequ) )
!          else
!            ! 2D:
!            call f_xy_2D( coord_x, coord_y, hvy_block(:, :, 1, 2, hvy_active(k)), Bs, g, Lx, Ly, frequ(ifrequ)  )
!          end if
!
!        end do
!
!        ! now the entire grid (incl ghost nodes) holds the exact solution: make a
!        ! copy of the grid for later comparison, but use work arrays ususally used for RK4 substages
!        ! so no additional memory is used.
!        hvy_work(:,:,:,1,:) = hvy_block(:,:,:,2,:)
!
!        !-----------------------------------------------------------------------
!        ! synchronize ghost nodes (this is what we test here)
!        !-----------------------------------------------------------------------
!        call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
!
!        !-----------------------------------------------------------------------
!        ! compute error (normalized, global, 2-norm)
!        !-----------------------------------------------------------------------
!        ! reset error
!        my_error = 0.0_rk
!        my_norm = 0.0_rk
!
!        ! loop over all active blocks and compute their error
!        do k = 1, hvy_n
!          my_error = my_error + sqrt(sum((hvy_block(:,:,:,2,hvy_active(k))-hvy_work(:,:,:,1,hvy_active(k)))**2 ))
!          my_norm = my_norm  + sqrt(sum((hvy_work(:,:,:,1,hvy_active(k)))**2 ))
!        end do
!
!        ! synchronize errors
!        call MPI_Allreduce(my_error, error(ifrequ), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
!        call MPI_Allreduce(my_norm, norm, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
!
!        error(ifrequ) = error(ifrequ) / norm
!
!        ! output
!        if (rank==0) then
!            write(*,'(" done - ghost nodes synchronization error = ",es16.8," frequ=",g12.4)')  error(ifrequ), frequ(ifrequ)
!        end if
!      end do
!
!    if (rank==0) then
!      write(*,'(" done - convergence order was ",6(g12.4,1x))')  sqrt(error(2:6) / error(1:5))
!      write(*,'(" done - mean convergence order was ",g12.4)')  sum(sqrt(error(2:6) / error(1:5))) / 5.0_rk
!    endif
!
!    !---------------------------------------------------------------------------------------------
!    ! last: clean up
!    deallocate(coord_x, coord_y, coord_z)

end subroutine unit_test_time_stepper_convergence
