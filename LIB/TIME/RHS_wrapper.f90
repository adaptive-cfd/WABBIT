!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_wrapper.f90
!> \version 0.5
!> \author sm
!
!> \brief wrapper for RHS call in time step function
!
!>
!! calls RHS depending on physics
!!
!! input:
!!           - time variable
!!           - time step dt
!!           - params
!!           - heavy data an lgt_block
!!           - coefficients for Runge Kutta
!!           - loop variable from time stepper
!!
!! output:
!!           - hvy_work
!!
!! butcher table, e.g.
!!
!! |   |    |    |   |
!! |---|----|----|---|
!! | 0 | 0  | 0  |  0|
!! |c2 | a21| 0  |  0|
!! |c3 | a31| a32|  0|
!! | 0 | b1 | b2 | b3|
!!
!!
!! = log ======================================================================================
!! \n
!! 23/05/17 - create
!
!**********************************************************************************************

subroutine RHS_wrapper(time, dt, params, hvy_work, rk_coeff, j, lgt_block, hvy_active, hvy_n, hvy_block)

!----------------------------------------------------------------------------------------------
! modules

!----------------------------------------------------------------------------------------------
! variables

   implicit none

    !> time variable
    real(kind=rk), intent(in)           :: time
    !> dt
    real(kind=rk), intent(in)           :: dt

    !> user defined parameter structure, hvy_active
    type (type_params), intent(in)      :: params
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
    !> heavy data array - block data
    real(kind=rk), intent(in)           :: hvy_block(:, :, :, :, :)

    !> loop variable inside time step routine
    integer(kind=ik), intent(in)        :: j
    !> coefficient for time + coeff*dt
    real(kind=rk), intent(in)           :: rk_coeff

    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    !> global integral
    real(kind=rk), dimension(3)         :: volume_int

    !> spacing and origin of a block
    real(kind=rk), dimension(3)         :: dx, x0
    ! loop variables
    integer(kind=ik)                    :: k, dF, N_dF, lgt_id
    ! grid parameter, error variable
    integer(kind=ik)                    :: Bs, g


!---------------------------------------------------------------------------------------------
! variables initialization

    ! number of datafields
    N_dF  = params%number_data_fields

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

!---------------------------------------------------------------------------------------------
! main body

    ! RHS depends on physics
    select case(params%physics_type)

    case("ACM-new")
      !-------------------------------------------------------------------------
      ! 1st stage: init_stage.
      !-------------------------------------------------------------------------
      ! performs initializations in the RHS module, such as resetting integrals
      select case(params%physics_type)
      case ("ACM-new")
       ! this call is not done for all blocks, but only once, globally.
       call RHS_ACM( time+rk_coeff*dt, hvy_work(:,:,:,1:N_dF,hvy_active(1)), g, &
       x0, dx, hvy_work(:,:,:,1:N_dF,hvy_active(1)), "init_stage" )

      case default
       call abort(2152000, "physics_type is unknown"//params%physics_type)
      end select

      !-------------------------------------------------------------------------
      ! 2nd stage: integral_stage.
      !-------------------------------------------------------------------------
      ! For some RHS, the eqn depend not only on local, block based qtys, such as
      ! the state vector, but also on the entire grid, for example to compute a
      ! global forcing term (e.g. in FSI the forces on bodies). As the physics
      ! modules cannot see the grid, (they only see blocks), in order to encapsulate
      ! them nicer, two RHS stages have to be defined: integral / local stage.
      do k = 1, hvy_n
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        !---------- different physics modules ----------
        select case(params%physics_type)
        case ("ACM-new")
          ! input state vector: hvy_block, output RHS vector: hvy_work
          call RHS_ACM( time+rk_coeff*dt, hvy_block(:,:,:,1:N_dF, hvy_active(k)), g, &
               x0, dx, hvy_work(:,:,:,j*N_dF+1:(j+1)*N_dF, hvy_active(k)), "integral_stage" )

        case default
          call abort(2152000, "physics_type is unknown"//params%physics_type)
        end select
      enddo


      !-------------------------------------------------------------------------
      ! 3rd stage: post integral stage.
      !-------------------------------------------------------------------------
      ! in rhs module, used ror example for MPI_REDUCES
      select case(params%physics_type)
      case ("ACM-new")
        ! this call is not done for all blocks, but only once, globally.
        call RHS_ACM( time+rk_coeff*dt, hvy_work(:,:,:,1:N_dF,hvy_active(1)), g, &
        x0, dx, hvy_work(:,:,:,1:N_dF,hvy_active(1)), "post_stage" )

      case default
        call abort(2152000, "physics_type is unknown"//params%physics_type)
      end select


      !-------------------------------------------------------------------------
      ! 3rd stage: local evaluation of RHS on all blocks
      !-------------------------------------------------------------------------
      ! the second stage then is what you would usually do: evaluate local differential
      ! operators etc.
      do k = 1, hvy_n
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        !---------- different physics modules ----------
        select case(params%physics_type)
        case ("ACM-new")
          ! input state vector: hvy_block, output RHS vector: hvy_work
          call RHS_ACM( time+rk_coeff*dt, hvy_block(:,:,:,1:N_dF, hvy_active(k)), g, &
               x0, dx, hvy_work(:,:,:,j*N_dF+1:(j+1)*N_dF, hvy_active(k)), "local_stage" )

        case default
          call abort(2152000, "physics_type is unknown"//params%physics_type)
        end select
      enddo




        case('2D_convection_diffusion')
            ! loop over all data fields
            do dF = 1, N_dF
                ! loop over all active heavy data blocks
                do k = 1, hvy_n
                    ! copy ghost nodes to hvy_work
                    hvy_work( :, :, :, (dF-1)*5+j+1, hvy_active(k) ) = hvy_block(:, :, :, dF, hvy_active(k) )

                    ! convert given hvy_id to lgt_id for block spacing routine
                    call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

                    ! get block spacing for RHS
                    call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                    ! RHS (compute k-coefficients)
                    call RHS_2D_convection_diffusion( hvy_work( :, :, 1, (dF-1)*5+j+1, hvy_active(k) ), &
                                      dx(1), dx(2), g, Bs, &
                                      params%physics%u0( (dF-1)*2 + 1 ), params%physics%u0( (dF-1)*2 + 2 ), &
                                      params%physics%nu(dF), params%order_discretization  )
                 end do
            end do

       case('2D_navier_stokes')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n
                ! copy ghost nodes to hvy_work
                hvy_work( :, :, :, j*N_dF+1:(j+1)*N_dF, hvy_active(k) ) = hvy_block(:, :, :, 1:N_dF, hvy_active(k) )

                ! convert given hvy_id to lgt_id for block spacing routine
                call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

                ! get block spacing for RHS
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                ! RHS (compute k-coefficients)
                call RHS_2D_navier_stokes( params%physics_ns, g, Bs, &
                                  dx(1), dx(2), N_dF, &
                                  hvy_work( :, :, 1, j*N_dF+1:(j+1)*N_dF, hvy_active(k) ))
            end do

       case('3D_convection_diffusion')
           ! loop over all data fields
           do dF = 1, N_dF
               ! loop over all active heavy data blocks
               do k = 1, hvy_n
                   ! copy ghost nodes to hvy_work
                   hvy_work( :, :, :, (dF-1)*5+j+1, hvy_active(k) ) = hvy_block(:, :, :, dF, hvy_active(k) )

                   ! convert given hvy_id to lgt_id for block spacing routine
                   call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

                   ! get block spacing for RHS
                   call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                   ! RHS (compute k-coefficients)
                   call RHS_3D_convection_diffusion( hvy_work( :, :, 1, (dF-1)*5+j+1, hvy_active(k) ), &
                                     dx(1), dx(2), dx(3), g, Bs, &
                                     params%physics%u0( (dF-1)*2 + 1 ), params%physics%u0( (dF-1)*2 + 2 ), &
                                     params%physics%u0( (dF-1)*2 + 3 ), &
                                     params%physics%nu(dF), params%order_discretization  )
                end do
           end do

       case('3D_navier_stokes')
           ! loop over all active heavy data blocks
           do k = 1, hvy_n
               ! copy ghost nodes to hvy_work
               hvy_work( :, :, :, j*N_dF+1:(j+1)*N_dF, hvy_active(k) ) = hvy_block(:, :, :, 1:N_dF, hvy_active(k) )

               ! convert given hvy_id to lgt_id for block spacing routine
               call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

               ! get block spacing for RHS
               call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

               ! RHS (compute k-coefficients)
               call RHS_3D_navier_stokes( params%physics_ns, g, Bs, &
                              dx(1), dx(2), dx(3), N_dF, &
                              hvy_work( :, :, :, j*N_dF+1:(j+1)*N_dF, hvy_active(k) ))
           end do

       case('2D_advection')
           ! loop over all data fields
            do dF = 1, N_dF
                ! loop over all active heavy data blocks
                do k = 1, hvy_n
                    ! copy ghost nodes to hvy_work
                    hvy_work( :, :, :, (dF-1)*5+j+1, hvy_active(k) ) = hvy_block(:, :, :, dF, hvy_active(k) )

                    ! convert given hvy_id to lgt_id for block spacing routine
                    call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

                    ! get block spacing and origin for RHS
                    call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                    ! RHS (compute k-coefficients)
                    ! k_j = RHS((t+dt*c_j, data_field(t) + sum(a_jl*k_l)) (time-dependent rhs)
                    call RHS_2D_advection( hvy_work( :, :, 1, (dF-1)*5+j+1, hvy_active(k) ), &
                                       x0(1:2), dx(1:2), g, Bs, &
                                       time + rk_coeff*dt, &
                                       params%order_discretization  )
                 end do
            end do

        case default
            call abort(1717,"ERROR: physics type is unknown"//params%physics_type)

    end select


end subroutine RHS_wrapper
