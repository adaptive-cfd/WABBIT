
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \file
!> \callgraph
!> \brief wrapper for RHS call in time step function, computes RHS in work array
!! (inplace)
!> \version 0.5
!> \author sm
!! \date 23/05/17 - create
!!
!
!>\details
!! calls RHS depending on physics
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
!
!**********************************************************************************************

subroutine RHS_wrapper(time, params, hvy_state, hvy_rhs, hvy_grid, lgt_block, hvy_active, hvy_n, first_substep)

!----------------------------------------------------------------------------------------------
! modules

!----------------------------------------------------------------------------------------------
! variables

   implicit none

    !> time variable
    real(kind=rk), intent(in)           :: time
    !> user defined parameter structure, hvy_active
    type (type_params), intent(in)      :: params
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_rhs(:, :, :, :, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_state(:, :, :, :, :)
    !> hvy_grid are qtys that depend on grid and not explicitly on time
    real(kind=rk), intent(inout)        :: hvy_grid(:, :, :, :, :)
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n
    !> some operations might be done only in the first RK substep, hence we pass
    !! this flag to check if this is the first call at the current time level.
    logical, optional, intent(in)       :: first_substep

    !> global integral
    real(kind=rk), dimension(3)         :: volume_int

    !> spacing and origin of a block
    real(kind=rk), dimension(3)         :: dx, x0
    ! loop variables
    integer(kind=ik)                    :: k, dF, neqn, lgt_id
    ! grid parameter, error variable
    integer(kind=ik)                    :: g
    integer(kind=ik), dimension(3)      :: Bs

    logical :: first_substep2

    integer(kind=2)                    :: surface(3)=0

!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs    = params%Bs
    g     = params%n_ghosts

    ! the first_substep flag is optional and its default is "false"
    first_substep2 = .false.
    if (present(first_substep)) first_substep2=first_substep

!---------------------------------------------------------------------------------------------
! main body


    !-------------------------------------------------------------------------
    ! 1st stage: init_stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! performs initializations in the RHS module, such as resetting integrals
    call RHS_meta( params%physics_type, time, hvy_state(:,:,:,:,hvy_active(1)), g, x0, dx, &
        hvy_rhs(:,:,:,:,hvy_active(1)), hvy_grid(:,:,:,:,hvy_active(1)), "init_stage", &
        first_substep=first_substep2 )

    !-------------------------------------------------------------------------
    ! 2nd stage: integral_stage. (called for all blocks)
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

      call RHS_meta( params%physics_type, time, hvy_state(:,:,:,:, hvy_active(k)), g, x0, dx,&
          hvy_rhs(:,:,:,:,hvy_active(k)), hvy_grid(:,:,:,:,hvy_active(k)), &
          "integral_stage", first_substep=first_substep2 )
    enddo


    !-------------------------------------------------------------------------
    ! 3rd stage: post integral stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! in rhs module, used ror example for MPI_REDUCES
    call RHS_meta( params%physics_type, time, hvy_state(:,:,:,:, hvy_active(1)), g, x0, dx, &
        hvy_rhs(:,:,:,:,hvy_active(1)), hvy_grid(:,:,:,:,hvy_active(1)), &
        "post_stage", first_substep=first_substep2 )


    !-------------------------------------------------------------------------
    ! 4th stage: local evaluation of RHS on all blocks (called for all blocks)
    !-------------------------------------------------------------------------
    ! the second stage then is what you would usually do: evaluate local differential
    ! operators etc.
    do k = 1, hvy_n
      ! convert given hvy_id to lgt_id for block spacing routine
      call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
      ! get block spacing for RHS
      call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

      if ( .not. All(params%periodic_BC) ) then
        ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
        call get_adjacent_boundary_surface_normal(params, lgt_id, lgt_block, params%max_treelevel, surface)
      endif
      ! if (surface(1).ne. 0 .or. surface(2).ne.0) then
      !   write(*,*) "surface normal",lgt_block(lgt_id,1:params%max_treelevel)
      ! endif

      call RHS_meta( params%physics_type, time, hvy_state(:,:,:,:, hvy_active(k)), g, &
           x0, dx, hvy_rhs(:,:,:,:, hvy_active(k)), hvy_grid(:,:,:,:, hvy_active(k)), "local_stage", &
           boundary_flag=surface, first_substep=first_substep2 )
    enddo

end subroutine RHS_wrapper
