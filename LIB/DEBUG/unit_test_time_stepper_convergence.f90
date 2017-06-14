!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name unit_test_time_stepper_convergence.f90
!> \version 0.5
!> \author msr
!
!> \brief unit test for time stepper
!> \note need additional memory to save results for smallest time step \n
!!
!!
!! input:    - params, empty light and heavy data arrays \n
!! output:   -                                           \n
!!
!!
!! = log ======================================================================================
!! \n
!! 18/04/17 - create
!
! ********************************************************************************************

subroutine unit_test_time_stepper_convergence( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist, com_lists, com_matrix )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> user defined parameter structure
    type (type_params), intent(inout)       :: params
    !> light data array
    integer(kind=ik),  intent(inout)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)
    !> heavy work array  )
    real(kind=rk),  intent(inout)           :: hvy_work (:, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik),  intent(inout)        :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: lgt_active(:)
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: hvy_active(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)      :: lgt_sortednumlist(:,:)

    ! communication lists:
    integer(kind=ik), intent(inout)     :: com_lists(:, :, :, :)
    ! communications matrix:
    integer(kind=ik), intent(inout)     :: com_matrix(:,:,:)

    ! local user defined parameter structure - use to change settings
    type (type_params)                      :: params_loc

    ! number of active blocks (heavy data)
    integer(kind=ik)                        :: hvy_n
    ! number of active blocks (light data)
    integer(kind=ik)                        :: lgt_n

    ! loop variables
    integer(kind=ik)                        :: k, l, lgt_id, hvy_id

    ! process rank
    integer(kind=ik)                        :: rank

    ! coordinates vectors
    real(kind=rk), allocatable              :: coord_x(:), coord_y(:), coord_z(:)
    ! spacing
    real(kind=rk)                           :: ddx(1:3), xx0(1:3)

    ! grid parameter
    integer(kind=ik)                        :: Bs, g
    real(kind=rk)                           :: Lx, Ly, Lz
    ! data dimensionality
    integer(kind=ik)                        :: d
    ! frequency of sin functions for testing:
    integer(kind=ik)                        :: num_dt(1:6)
    integer(kind=ik)                        :: idt

    ! error variable
    real(kind=rk)                           :: error(1:6), my_error, norm, my_norm
    ! MPI error variable
    integer(kind=ik)                        :: ierr

    ! time step, dx, time
    real(kind=rk)                           :: dt, dx, my_dx, time

    ! array to save heavy data on smallest time step
    real(kind=rk), allocatable              :: hvy_old (:, :, :, :, : )

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! copy params data to local struct
    params_loc = params

    ! set MPI parameters
    rank = params%rank

    ! grid parameter
    Lx = 1.0_rk
    Ly = 1.0_rk
    Lz = 1.0_rk

    ! set domain size
    params_loc%Lx = Lx
    params_loc%Ly = Ly
    params_loc%Lz = Lz

    ! set data dimension
    if ( params_loc%threeD_case ) then
        d = 3
    else
        d = 2
    end if

    ! reset dx
    my_dx = 9.0e9_rk
    dx    = 9.0e9_rk
    ! reset dt
    dt    = 9.0e9_rk

!---------------------------------------------------------------------------------------------
! main body

    if (rank == 0) then
        write(*,'(80("_"))')
        write(*,'("UNIT TEST: Beginning time stepper convergence test")')
    end if

    Bs = params_loc%number_block_nodes
    g  = params_loc%number_ghost_nodes

    if (rank == 0) then
      write(*,'("UNIT TEST: testing Bs=",i4," blocks-per-mpirank=",i5)')  Bs, params_loc%number_blocks
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
    call create_equidistant_base_mesh( params_loc, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, params_loc%min_treelevel, .true. )

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
        if (.not. allocated(params_loc%physics%u0)) then
            allocate(params_loc%physics%u0(2))
        end if
        if (.not. allocated(params_loc%physics%nu)) then
           allocate(params_loc%physics%nu(1))
        end if
        params_loc%physics%nu = 0.0_rk
        params_loc%physics%u0 = (/1.0_rk, 0.5_rk/)
    end if

    ! set dt parameter
    params_loc%time_step_method = 'CFL_cond'
    params_loc%CFL              = 0.5_rk

    ! FIXME: you could also look over light data, as ddx is available only from that. no mpi
    do k = 1, hvy_n
        ! light id of this block
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params_loc%rank, params_loc%number_blocks )
        ! compute blocks' spacing from treecode
        call get_block_spacing_origin( params_loc, lgt_id, lgt_block, xx0, ddx )
        ! find smallest dx of active blocks
        my_dx = min(my_dx, minval(ddx(1:d)) )

        ! HACK repair first datafield, as we're about to remove it
        hvy_block(:,:,:,1,hvy_active(k)) = 0.0_rk
        hvy_block(1,2,:,1,hvy_active(k)) = ddx(1)
        hvy_block(2,2,:,1,hvy_active(k)) = ddx(2)
        hvy_block(3,2,:,1,hvy_active(k)) = ddx(3)

    end do

    ! find globally smallest dx
    call MPI_Allreduce(my_dx, dx, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

    ! calculate dt
    call calculate_time_step( params_loc, dx, dt )

    ! set calculated dt as fixed dt
    params_loc%time_step_method = 'fixed'
    params_loc%dt = dt
    ! set time max
    params_loc%time_max = 1.0_rk

    ! set number array (numbers of dt)
    num_dt = (/ 32, 16, 8, 4, 2, 1 /)

    ! allocate heavy data array to save result for smallest dt
    if ( params_loc%threeD_case ) then
        ! 3D:
        allocate( hvy_old(Bs+2*g, Bs+2*g, Bs+2*g, 1, hvy_n ) )
    else
        ! 2D:
        allocate( hvy_old(Bs+2*g, Bs+2*g, 1, 1, hvy_n ) )
    end if

    if (rank==0) then
        write(*,'(" running test ... may take some time ... ")')
    end if

    ! loop over dts
    do idt = 1 , size(num_dt)

        ! set dt and time
        params_loc%dt  = dt/real(num_dt(idt),kind=rk)
        time           = 0.0_rk

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
        ! time stepper
        !-----------------------------------------------------------------------
        do l = 1, num_dt(idt)*10
            call time_stepper( time, params_loc, lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix )
        end do

        if ( idt == 1 ) then
            !-----------------------------------------------------------------------
            ! save result for smallest dt
            !-----------------------------------------------------------------------
            ! save result, loop over all active blocks
            do k = 1, hvy_n
                hvy_old(:, :, :, 1, k) = hvy_block(:, :, :, 2, hvy_active(k))
            end do

            error(idt) = 0.0_rk

        else
            !-----------------------------------------------------------------------
            ! compute error (normalized, global, 2-norm)
            !-----------------------------------------------------------------------
            ! reset error
            my_error = 0.0_rk
            my_norm = 0.0_rk

            ! loop over all active blocks and compute their error
            if ( params_loc%threeD_case ) then
                ! 3D:
                do k = 1, hvy_n
                    my_error = my_error + sqrt( sum( ( hvy_block(g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 2, hvy_active(k)) - hvy_old(g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 1, k) )**2 ) )
                    my_norm = my_norm  + sqrt(sum(( hvy_old(g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 1, k) )**2 ))
                end do
            else
                ! 2D:
                do k = 1, hvy_n
                    my_error = my_error + sqrt( sum( ( hvy_block(g+1:Bs+g, g+1:Bs+g, 1, 2, hvy_active(k)) - hvy_old(g+1:Bs+g, g+1:Bs+g, 1, 1, k) )**2 ) )
                    my_norm = my_norm  + sqrt(sum(( hvy_old(g+1:Bs+g, g+1:Bs+g, 1, 1, k) )**2 ))
                    !my_error = my_error + sqrt( sum( ( hvy_block(:, :, 1, 2, hvy_active(k)) - hvy_old(:, :, 1, 1, k) )**2 ) )
                    !my_norm = my_norm  + sqrt(sum(( hvy_old(:, :, 1, 1, k) )**2 ))
                end do
            end if

            ! synchronize errors
            call MPI_Allreduce(my_error, error(idt), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
            call MPI_Allreduce(my_norm, norm, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

            error(idt) = error(idt) / norm

        end if

        ! output
        if (rank==0) then
            write(*,'(" done - time stepper with dt = ",g12.4, " error = ", g12.4)')  dt/real(num_dt(idt),kind=rk), error(idt)
        end if

    end do

    if (rank==0) then
        write(*,'(" done - convergence order was ",6(g12.4,1x))')  sqrt(error(3:6) / error(2:5))
        write(*,'(" done - mean convergence order was ",g12.4)')  sum(sqrt(error(3:6) / error(2:5))) / 4.0_rk
    endif

    !---------------------------------------------------------------------------------------------
    ! last: clean up
    deallocate(coord_x, coord_y, coord_z, hvy_old)

end subroutine unit_test_time_stepper_convergence
