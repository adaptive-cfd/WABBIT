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
subroutine set_inicond_all_blocks(params, lgt_block, hvy_block, hvy_active, hvy_n, inicond)

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
  integer(kind=ik)                     :: k, dF, pF, rhoF
  integer(kind=ik)                     :: hvy_id, lgt_id
  ! origin and spacing of blocks
  real(kind=rk)                        :: x0(1:3), dx(1:3)

  ! p0 value \todo get from ini file
  real(kind=rk)                        :: p0

  !---------------------------------------------------------------------------------------------
  ! interfaces

  !---------------------------------------------------------------------------------------------
  ! variables initialization

    pF = 0
    rhoF = 0

    p0 = 1.0e5_rk

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

        ! note: subroutine sets initial condition on first datafield
        ! set the initial condition on this block
        call initial_condition_on_block_wrapper( params, hvy_block(:,:,:,:,hvy_id), x0, dx, inicond )

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
                    ! find p field
                    do dF = 1, params%number_data_fields
                        if ( params%physics_ns%names(dF) == "p" ) pF = dF
                        if ( params%physics_ns%names(dF) == "rho" ) rhoF = dF
                    end do

                    do dF = 2, params%number_data_fields
                        if ( dF == pF ) then
                            ! set gauus blob to p, add p0
                            hvy_block( :, :, :, dF, hvy_id) = p0 + 1000.0_rk * hvy_block( :, :, :, 1, hvy_id)
                        elseif ( dF == rhoF ) then
                            ! ini rho field
                            hvy_block( :, :, :, dF, hvy_id) = 1.0_rk
                        else
                            ! velocity fields
                            hvy_block( :, :, :, dF, hvy_id) = 0.0_rk
                        end if
                    end do

                    ! set first datafield, note: still contains gauss blob
                    if ( pF == 1 ) then
                        ! first field is p field, add p0
                        hvy_block( :, :, :, 1, hvy_id) = 1000.0_rk * hvy_block( :, :, :, 1, hvy_id) + p0
                    elseif ( rhoF == 1 ) then
                        ! first field is rho, set to one
                        hvy_block( :, :, :, 1, hvy_id) = 1.0_rk
                    else
                        ! other field
                        hvy_block( :, :, :, 1, hvy_id) = 0.0_rk
                    end if

            end select

        else
            ! nothing to do
        end if

    enddo

end subroutine set_inicond_all_blocks
