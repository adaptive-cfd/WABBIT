!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name write_vorticity.f90
!> \version 0.5
!> \author sm
!
!> \brief compute vorticity for time step t (for writing it to disk)
!
!> \note routine uses hvy_work array for computation of vorticity, thus cannot be used during RK stages!
!!
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
subroutine write_vorticity( hvy_work, hvy_block, lgt_block, hvy_active, hvy_n, params, time, iteration, lgt_active, lgt_n)

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> physics parameter structure
    type (type_params), intent(in)                 :: params
    !> actual block data
    real(kind=rk), intent(in)                      :: hvy_block(:, :, :, :, :)
    !> hvy_work array to safe u,v(,w) and vorticity
    real(kind=rk), intent(inout)                   :: hvy_work(:, :, :, :, :)
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

    !> origin and spacing of the block
    real(kind=rk), dimension(3)                    :: dx, x0
    ! loop variables
    integer(kind=ik)                               :: k, lgt_id
    ! file name
    character(len=80)                              :: fname

    ! field numbers for navier stokes physics
    integer(kind=ik)                               :: df, pF, rhoF, UxF, UyF, UzF

!---------------------------------------------------------------------------------------------
! variables initialization


    pF   = 0
    rhoF = 0
    UxF  = 0
    UyF  = 0
    UzF  = 0

    if ( params%physics_type=='2D_navier_stokes' .or. params%physics_type=='3D_navier_stokes' ) then
        ! find fields
        do dF = 1, params%number_data_fields
            if ( params%physics_ns%names(dF) == "p" ) pF = dF
            if ( params%physics_ns%names(dF) == "rho" ) rhoF = dF
            if ( params%physics_ns%names(dF) == "Ux" ) UxF = dF
            if ( params%physics_ns%names(dF) == "Uy" ) UyF = dF
            if ( params%physics_ns%names(dF) == "Uz" ) UzF = dF
        end do
    end if

!---------------------------------------------------------------------------------------------
! main body

    ! calculate vorticity only for proper physic
    select case (params%physics_type)

        case('2D_navier_stokes')
            do k=1, hvy_n

                call hvy_id_to_lgt_id(lgt_id, hvy_active(k), params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
                ! store u,v in hvy_work array
                hvy_work(:,:,1,1,hvy_active(k)) = hvy_block(:,:,1,UxF,hvy_active(k))/hvy_block(:,:,1,rhoF,hvy_active(k))**2  ! u
                hvy_work(:,:,1,2,hvy_active(k)) = hvy_block(:,:,1,UyF,hvy_active(k))/hvy_block(:,:,1,rhoF,hvy_active(k))**2  ! v

                ! call compute_vorticity(params, hvy_work(:,:,:,1,hvy_active(k)), hvy_work(:,:,:,2,hvy_active(k)), hvy_work(:,:,:,3,hvy_active(k)), dx, hvy_work(:,:,:,4:6,hvy_active(k)))

            end do

            write( fname,'(a, "_", i12.12, ".h5")') 'vor', nint(time * 1.0e6_rk)
            ! write field 4 of hvy_work array (vorticity) to disk
            call write_field(fname, time, iteration, 4, params, lgt_block, hvy_work(:,:,:,:,:), lgt_active, lgt_n, hvy_n)

        case('3D_navier_stokes')

             do k=1, hvy_n

                call hvy_id_to_lgt_id(lgt_id, hvy_active(k), params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
                ! store u,v,w in hvy_work array
                hvy_work(:,:,:,1,hvy_active(k)) = hvy_block(:,:,:,UxF,hvy_active(k))/hvy_block(:,:,:,rhoF,hvy_active(k))**2  ! u
                hvy_work(:,:,:,2,hvy_active(k)) = hvy_block(:,:,:,UyF,hvy_active(k))/hvy_block(:,:,:,rhoF,hvy_active(k))**2  ! v
                hvy_work(:,:,:,3,hvy_active(k)) = hvy_block(:,:,:,UyF,hvy_active(k))/hvy_block(:,:,:,rhoF,hvy_active(k))**2  ! w

                ! call compute_vorticity(params, hvy_work(:,:,:,1,hvy_active(k)), hvy_work(:,:,:,2,hvy_active(k)), hvy_work(:,:,:,3,hvy_active(k)), dx, hvy_work(:,:,:,4:6,hvy_active(k)))

            end do

    end select

end subroutine write_vorticity
