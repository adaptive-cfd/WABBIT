!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name write_mask.f90
!> \version 0.5
!> \author sm
!
!> \brief write mask at time step t to disk
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
!! 24/10/17 - create
!
! ********************************************************************************************
subroutine write_mask( hvy_work, lgt_block, hvy_active, hvy_n, params, time, iteration, lgt_active, lgt_n)

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> physics parameter structure
    type (type_params), intent(in)                 :: params
    !> hvy_work array to safe mask
    real(kind=rk), intent(in)                      :: hvy_work(:, :, :, :, :)
    !> time
    real(kind=rk), intent(in)                      :: time
    !> number of active blocks (heavy and light data)
    integer(kind=ik), intent(in)                   :: hvy_n, lgt_n, iteration
    !> light data array
    integer(kind=ik), intent(in)                   :: lgt_block(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)                   :: hvy_active(:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)                   :: lgt_active(:)

    ! origin and spacing of the block
    real(kind=rk), dimension(3)                    :: dx, x0
    ! loop variables
    integer(kind=ik)                               :: k, lgt_id
    ! file name
    character(len=80)                              :: fname
    !grid parameter
    integer(kind=ik)                               :: Bs, g

!---------------------------------------------------------------------------------------------
! variables initialization

    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes
!---------------------------------------------------------------------------------------------
! main body

    select case (params%physics_type)

        case('2D_acm')
 
            do k=1, hvy_n

                call hvy_id_to_lgt_id(lgt_id, hvy_active(k), params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
                ! safe mask in first field heavy work array
                call create_mask_2D(params, hvy_work(:,:,1,1,hvy_active(k)), x0, dx, Bs, g)
                
            end do

            write( fname,'(a, "_", i12.12, ".h5")') 'mask', nint(time * 1.0e6_rk)
            ! write field 1 of hvy_work array (mask) to disk
            call write_field(fname, time, iteration, 1, params, lgt_block, hvy_work, lgt_active, lgt_n, hvy_n)

        case('3D_acm')

            do k=1, hvy_n

                call hvy_id_to_lgt_id(lgt_id, hvy_active(k), params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
                ! safe mask in first field heavy work array
                call create_mask_3D(params, hvy_work(:,:,:,1,hvy_active(k)), x0, dx, Bs, g)
                
            end do

    end select

end subroutine write_mask
