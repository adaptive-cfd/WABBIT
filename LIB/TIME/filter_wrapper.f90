
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \file
!> \callgraph
!> \brief wrapper for filter
!> \version 0.5
!> \author Pkrah
!! \date 30/04/18 - create
!!
!
!**********************************************************************************************

subroutine filter_wrapper(time, params, hvy_state, hvy_work, lgt_block, hvy_active, hvy_n)

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
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_state(:, :, :, :, :)
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
    integer(kind=ik)                    :: k, dF, neqn, lgt_id
    ! grid parameter, error variable
    integer(kind=ik)                    :: Bs, g


!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs    = params%Bs
    g     = params%nr_ghosts


    do k = 1, hvy_n
      ! convert given hvy_id to lgt_id for block spacing routine
      call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
      ! get block spacing for RHS
      call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

      call filter_meta(params%physics_type, time, hvy_state(:,:,:,:, hvy_active(k)), g, x0, dx,&
          hvy_work(:,:,:,:,hvy_active(k)))
    enddo

    !update state vector

end subroutine filter_wrapper
