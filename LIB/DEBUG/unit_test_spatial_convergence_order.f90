! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: unit_test_spatial_convergence_order.f90
! version: 0.5
! author: msr
!
! unit test for time stepper
! note: need additional memory to save results for smallest time step
!
! input:    - params, empty light and heavy data arrays
! output:   -
!
! = log ======================================================================================
!
! 18/04/17 - create
!
! ********************************************************************************************

subroutine unit_test_spatial_convergence_order( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active )

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
    integer(kind=ik)                        :: k, l, lgt_id, hvy_id, hvy_id_old, lgt_id_old

    ! process rank
    integer(kind=ik)                        :: rank

    ! coordinates vectors
    real(kind=rk), allocatable              :: coord_x(:), coord_y(:), coord_z(:)
    ! spacing
    real(kind=rk)                           :: ddx(1:3), xx0(1:3)

    ! grid parameter
    integer(kind=ik)                        :: Bs, g
    real(kind=rk)                           :: Lx, Ly, Lz

    ! loop variable
    integer(kind=ik)                        :: idx

    ! error variable
    real(kind=rk), allocatable              :: error(:)
    real(kind=rk)                           :: my_error, norm, my_norm
    ! MPI error variable
    integer(kind=ik)                        :: ierr

    ! dummy time variable
    real(kind=rk)                           :: time

    ! array to save heavy data on finest mesh level
    real(kind=rk), allocatable              :: hvy_old (:, :, :, :, : )
    ! old light data array on finest mesh level
    integer(kind=ik)                        :: lgt_old( size(lgt_block,1), size(lgt_block,2))
    ! old list of active blocks (light data)
    integer(kind=ik)                        :: lgt_active_old( size(lgt_active) )
    ! old number of active blocks (light data)
    integer(kind=ik)                        :: lgt_n_old

    ! dummy variable for block existence
    logical                                 :: exists

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! copy params data to local struct
    params_loc = params

    ! set MPI parameters
    rank = params%rank

    ! grid parameter
    Lx = params_loc%Lx
    Ly = params_loc%Ly
    Lz = params_loc%Lz

    ! allocate error array
    allocate( error(params_loc%max_treelevel) )

!---------------------------------------------------------------------------------------------
! main body

    if (rank == 0) then
        write(*,'(80("_"))')
        write(*,'("UNIT TEST: Beginning spatial convergence test")')
    end if

    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    if (rank == 0) then
      write(*,'("UNIT TEST: testing Bs=",i4," blocks-per-mpirank=",i5)')  Bs, params_loc%number_blocks
    end if

    !---------------------------------------------------------------------------
    ! Step 1: calculate solution on finest grid and store result
    !---------------------------------------------------------------------------
    ! allocate coord arrays
    allocate( coord_x( Bs + 2*g ), coord_y( Bs + 2*g ), coord_z( Bs + 2*g ) )

    ! set all blocks to free (since if we call inicond twice, all blocks are used in the second call)
    lgt_block = -1
    lgt_active = -1; lgt_N = 0
    hvy_active = -1; hvy_N = 0

    ! setup the coarsest grid level with some data (we don't care what data, we'll erase it)
    ! Note that active lists + neighbor relations are updated inside this routine as well, as
    ! the grid is modified
    call create_equidistant_base_mesh( params_loc, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, 1, .true. )
    ! refine the mesh to max grid-level
    do l = 1, params_loc%max_treelevel
        call refine_mesh( params_loc, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, "everywhere" )
    end do

    ! reset time step to fixed small time step
    params_loc%time_step_method = 'fixed'
    params_loc%dt               = 1e-6_rk

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
            allocate(params_loc%physics%u0(3))
            allocate(params_loc%physics%nu(1))
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
            allocate(params_loc%physics%u0(2))
            allocate(params_loc%physics%nu(1))
            params_loc%physics%nu = 0.0_rk
            params_loc%physics%u0 = (/1.0_rk, 0.5_rk/)
        end if
    end if

    ! allocate heavy data array to save result for smallest dt
    if ( params_loc%threeD_case ) then
        ! 3D:
        allocate( hvy_old(Bs+2*g, Bs+2*g, Bs+2*g, 1, params_loc%number_blocks ) )
    else
        ! 2D:
        allocate( hvy_old(Bs+2*g, Bs+2*g, 1, 1, params_loc%number_blocks ) )
    end if

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
        call get_block_spacing_origin( params_loc, lgt_id, lgt_block, xx0, ddx )

        ! fill coordinate arrays, of course including ghost nodes
        do l = 1, Bs+2*g
            coord_x(l) = real(l-(g+1), kind=rk) * ddx(1) + xx0(1)
            coord_y(l) = real(l-(g+1), kind=rk) * ddx(2) + xx0(2)
            coord_z(l) = real(l-(g+1), kind=rk) * ddx(3) + xx0(3)
        enddo

        ! calculate f(x,y,z) for first datafield
        if ( params_loc%threeD_case ) then
            ! 3D:
            call f_xyz_3D( coord_x, coord_y, coord_z, hvy_block(:, :, :, 2, hvy_active(k)), Bs, g, Lx, Ly, Lz, 1.0_rk )
        else
            ! 2D:
            call f_xy_2D( coord_x, coord_y, hvy_block(:, :, 1, 2, hvy_active(k)), Bs, g, Lx, Ly, 1.0_rk  )
        end if

    end do

    !-----------------------------------------------------------------------
    ! do n time steps
    !-----------------------------------------------------------------------
    time = 0.0_rk
    do l = 1, 10
        call time_step_RK4( time, params_loc, lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, hvy_n )
    end do

    !-----------------------------------------------------------------------
    ! save result on finest grid
    !-----------------------------------------------------------------------
    ! save result, loop over all active blocks
    do k = 1, hvy_n
        hvy_old(:, :, :, 1, :)  = hvy_block(:, :, :, 2, :)
        lgt_old                 = lgt_block
        lgt_active_old          = lgt_active
        lgt_n_old               = lgt_n
    end do

    error(1) = 0.0_rk

    !-----------------------------------------------------------------------
    ! reset grid
    !-----------------------------------------------------------------------
    call reset_grid( params_loc, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, .false. )

    !---------------------------------------------------------------------------
    ! Step 2: spatial test - loop over all possible mesh levels
    !---------------------------------------------------------------------------
    ! i)   create grid on coarsest level
    ! ii)  refine to test grid-level - write exact solution
    ! iii) do n time steps
    ! iv)  move grid to finest level and calculate error
    ! v)   reset grid

    if (rank==0) then
        write(*,'(" running test ... may take some time ... ")')
    end if

    do idx = 2, params_loc%max_treelevel

        ! setup the coarsest grid level with some data (we don't care what data, we'll erase it)
        ! Note that active lists + neighbor relations are updated inside this routine as well, as
        ! the grid is modified
        call create_equidistant_base_mesh( params_loc, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, 1, .false. )
        ! refine the mesh to test grid-level
        do l = idx, params_loc%max_treelevel-1
            call refine_mesh( params_loc, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, "everywhere" )
        end do

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
            call get_block_spacing_origin( params_loc, lgt_id, lgt_block, xx0, ddx )

            ! fill coordinate arrays, of course including ghost nodes
            do l = 1, Bs+2*g
                coord_x(l) = real(l-(g+1), kind=rk) * ddx(1) + xx0(1)
                coord_y(l) = real(l-(g+1), kind=rk) * ddx(2) + xx0(2)
                coord_z(l) = real(l-(g+1), kind=rk) * ddx(3) + xx0(3)
            enddo

            ! calculate f(x,y,z) for first datafield
            if ( params_loc%threeD_case ) then
                ! 3D:
                call f_xyz_3D( coord_x, coord_y, coord_z, hvy_block(:, :, :, 2, hvy_active(k)), Bs, g, Lx, Ly, Lz, 1.0_rk )
            else
                ! 2D:
                call f_xy_2D( coord_x, coord_y, hvy_block(:, :, 1, 2, hvy_active(k)), Bs, g, Lx, Ly, 1.0_rk  )
            end if

        end do

        !-----------------------------------------------------------------------
        ! do n time steps
        !-----------------------------------------------------------------------
        time = 0.0_rk
        do l = 1, 10
            call time_step_RK4( time, params_loc, lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, hvy_n )
        end do

        !-----------------------------------------------------------------------
        ! move grid to finest level
        !-----------------------------------------------------------------------
        do l = 1, params_loc%max_treelevel
            call refine_mesh( params_loc, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, "everywhere" )
        end do

        !-----------------------------------------------------------------------
        ! calculate error
        ! compute error (normalized, global, 2-norm)
        !-----------------------------------------------------------------------
        ! reset error
        my_error = 0.0_rk
        my_norm = 0.0_rk

        ! loop over all active blocks and compute their error
        if ( params_loc%threeD_case ) then
            ! 3D:
            do k = 1, hvy_n

                ! light id for block k
                call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, params_loc%number_blocks )
                ! look for hvy_id in old data
                ! note: treecode is from current light data (lgt_block), heavy id is from old light data (lgt_old)
                call does_block_exist(lgt_block(lgt_id, 1:params_loc%max_treelevel), lgt_old, params_loc%max_treelevel, exists, lgt_id_old, lgt_active_old, lgt_n_old)

                ! old heavy id
                call lgt_id_to_hvy_id( hvy_id_old, lgt_id_old, rank, params_loc%number_blocks )

                my_error = my_error + sqrt( sum( ( hvy_block(g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 2, hvy_active(k)) - hvy_old(g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 1, k) )**2 ) )
                my_norm = my_norm  + sqrt(sum(( hvy_old(g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 1, k) )**2 ))

            end do
        else
            ! 2D:
            do k = 1, hvy_n

                ! light id for block k
                call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, params_loc%number_blocks )
                ! look for hvy_id in old data
                ! note: treecode is from current light data (lgt_block), heavy id is from old light data (lgt_old)
                call does_block_exist(lgt_block(lgt_id, 1:params_loc%max_treelevel), lgt_old, params_loc%max_treelevel, exists, lgt_id_old, lgt_active_old, lgt_n_old)

                ! old heavy id
                call lgt_id_to_hvy_id( hvy_id_old, lgt_id_old, rank, params_loc%number_blocks )

                ! calc error
                my_error = my_error + sqrt( sum( ( hvy_block(g+1:Bs+g, g+1:Bs+g, 1, 2, hvy_active(k)) - hvy_old(g+1:Bs+g, g+1:Bs+g, 1, 1, hvy_id_old) )**2 ) )
                my_norm = my_norm  + sqrt(sum(( hvy_old(g+1:Bs+g, g+1:Bs+g, 1, 1, hvy_id_old) )**2 ))

            end do
        end if

        ! synchronize errors
        call MPI_Allreduce(my_error, error(idx), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(my_norm, norm, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        error(idx) = error(idx) / norm

        ! output
        if (rank==0) then
            write(*,'(" done - spatial convergence test with dx = ",g12.4, " error = ", g12.4)')  ddx(1), error(idx)
        end if

        !-----------------------------------------------------------------------
        ! reset grid
        !-----------------------------------------------------------------------
        if (idx < params_loc%max_treelevel ) then
            call reset_grid( params_loc, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, .false. )
        end if

    end do

    if (rank==0) then
        write(*,'(" done - convergence order was ",6(g12.4,1x))')  sqrt(error(3:size(error)) / error(2:(size(error)-1)))
        write(*,'(" done - mean convergence order was ",g12.4)')  sum(sqrt(error(3:size(error)) / error(2:(size(error)-1)))) / real((size(error)-3)+1,kind=rk)
    endif

    !---------------------------------------------------------------------------------------------
    ! last: clean up
    deallocate(coord_x, coord_y, coord_z, hvy_old)

end subroutine unit_test_spatial_convergence_order
