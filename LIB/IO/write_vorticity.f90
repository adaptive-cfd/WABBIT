!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name write_vorticity.f90
!> \version 0.5
!> \author sm
!
!> \brief compute vorticity for time step t (for saving it on disk)
!
!>
!! input:
!!           - parameter array
!!           - light data array
!!           - heavy data array
!!
!! output:
!!           - 
!!
!!
!! = log ======================================================================================
!! \n
!! 24/07/17 - create
!
! ********************************************************************************************
subroutine write_vorticity( hvy_work, hvy_block, lgt_block, hvy_active, hvy_n, params, Bs, g, time, iteration, lgt_active, lgt_n)

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> physics parameter structure
    type (type_params), intent(in)                 :: params
    !> actual block data
    real(kind=rk), intent(in)                      :: hvy_block(:, :, :, :)
    !> hvy_work array to safe vorticity
    real(kind=rk), intent(inout)                   :: hvy_work(:, :, :, :, :)
    !> time
    real(kind=rk), intent(in)                      :: time

    integer(kind=ik), intent(in)                   :: Bs, g

    !> 
    integer(kind=ik), intent(in)                   :: hvy_n, lgt_n, iteration
    integer(kind=ik), intent(in)                   :: lgt_block(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)                   :: hvy_active(:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)                :: lgt_active(:)

    !> origin and spacing of the block
    real(kind=rk), dimension(3)                    :: dx, x0
    !> local datafields
    real(kind=rk), dimension(Bs+2*g, Bs+2*g)       :: u, v, vorticity
    ! loop variables
    integer(kind=ik)                               :: k, lgt_id
    ! file name
    character(len=80)                              :: fname

!---------------------------------------------------------------------------------------------
! variables initialization

    vorticity = 0.0_rk
    hvy_work  = 0.0_rk

!---------------------------------------------------------------------------------------------
! main body

    do k=1, hvy_n
       call hvy_id_to_lgt_id(lgt_id, hvy_active(k), params%rank, params%number_blocks)
       call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

       u = hvy_block(:, :, 1, hvy_active(k))
       v = hvy_block(:, :, 2, hvy_active(k))

      ! call compute_vorticity(params, u, v, dx, vorticity, Bs, g)

       hvy_work(:, :, 1, 1, hvy_active(k)) = vorticity(:,:)

   end do

   write( fname,'(a, "_", i12.12, ".h5")') 'vort', nint(time * 1.0e6_rk)

   call write_field(fname, time, iteration, 1, params, lgt_block, hvy_work(:,:,:,:,:), lgt_active, lgt_n, hvy_n)
   

end subroutine write_vorticity
