!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name init_data.f90
!> \version 0.5
!> \author msr
!
!
!
! ********************************************************************************************
subroutine set_inicond_blocks(params, lgt_block, hvy_block, hvy_active, hvy_n, inicond)

  !---------------------------------------------------------------------------------------------
  ! variables

  implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)    :: params
    !> light data array
    integer(kind=ik), intent(inout)      :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)         :: hvy_block(:, :, :, :, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)      :: hvy_active(:)
    !> number of heavy and light active blocks
    integer(kind=ik), intent(inout)      :: hvy_n
    !> what function to use
    character(len=*), intent(in)         :: inicond
    ! loop variable
    integer(kind=ik)                     :: k, dF, pF, rhoF, UxF, UyF, UzF
    integer(kind=ik)                     :: hvy_id, lgt_id
    ! origin and spacing of blocks
    real(kind=rk)                        :: x0(1:3), dx(1:3)

    ! p0 value
    !> \todo get from ini file, rework gauss blob setup - see shear layer setup as template
    real(kind=rk)                        :: p0, rho0

  !---------------------------------------------------------------------------------------------
  ! interfaces

  !---------------------------------------------------------------------------------------------
  ! variables initialization

    pF   = 0
    rhoF = 0
    UxF  = 0
    UyF  = 0
    UzF  = 0

    rho0 = 1.0_rk
    p0   = 1.0e5_rk

  !---------------------------------------------------------------------------------------------
  ! main body

    !---------------------------------------------------------------------------
    ! on the grid, evaluate the initial condition
    !---------------------------------------------------------------------------
    ! loop over my active heavy data
    do k = 1, hvy_n
        ! hvy_id of the block we're looking at
        hvy_id = hvy_active(k)
        ! light id of this block
        call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! compute block spacing and origin from treecode
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        !> \note subroutine sets initial condition on first datafield
        ! set the initial condition on this block
        call initial_condition_on_block_wrapper( params, hvy_block(:,:,:,:,hvy_id), x0, dx, inicond )

!######****** HACK: this entire if-clause has to go!!!!! it makes no sense
        ! more than on datafield
        if ( params%number_data_fields > 1 ) then

            ! datafield writing depends on inicond
            select case( params%initial_cond )

                case ("gauss-blob","gauss_blob")
                    ! all other datafields get data from first field
                    do dF = 2, params%number_data_fields
                        hvy_block( :, :, :, dF, hvy_id) = hvy_block( :, :, :, 1, hvy_id)
                    end do

                case ("ns_pressure_blob")
                    ! pressure field gets gauss blob, velocity fields set to zero, density field to 1
                    ! find fields
                    do dF = 1, params%number_data_fields
                        if ( params%physics_ns%names(dF) == "p" ) pF = dF
                        if ( params%physics_ns%names(dF) == "rho" ) rhoF = dF
                        if ( params%physics_ns%names(dF) == "Ux" ) UxF = dF
                        if ( params%physics_ns%names(dF) == "Uy" ) UyF = dF
                        if ( params%physics_ns%names(dF) == "Uz" ) UzF = dF
                    end do

                    ! set fields
                    ! set gauss blob to p, note: shear layer is in field one
                    hvy_block( :, :, :, pF, hvy_id) = p0 + 1000.0_rk * hvy_block( :, :, :, 1, hvy_id)
                    ! set rho
                    hvy_block( :, :, :, rhoF, hvy_id) = rho0
                    ! set Ux
                    hvy_block( :, :, :, UxF, hvy_id) = 0.0_rk
                    ! set Uy
                    hvy_block( :, :, :, UyF, hvy_id) = 0.0_rk

                    if (params%threeD_case) then
                        ! set Uz to zero
                        hvy_block( :, :, :, UzF, hvy_id) = 0.0_rk
                    endif

                case ("shear_layer")
                    ! everything is done in initial_condition_on_block_wrapper subroutine

            end select

        else
            ! nothing to do
        end if

    enddo

end subroutine set_inicond_blocks
