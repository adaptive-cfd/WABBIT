!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name get_block_max_velocity_norm.f90
!> \version 0.5
!> \author msr
!
!> \brief calculate max velocity norm for dt calculation. \n
!
!> 
!! input:    
!!           - params
!!           - list of active hvy data
!!           - hvy data
!!           - datafield number of velocity fields
!!
!! output:   
!!           - norm of velocity
!!
!! = log ======================================================================================
!! \n
!! 15/08/17 - create
!
subroutine get_block_max_velocity_norm( params, hvy_block, hvy_active, hvy_n, dF_numbers, norm_u )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> heavy data array - block data
    real(kind=rk), intent(in)           :: hvy_block(:, :, :, :, :)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)

    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    !> the block in question
    integer(kind=ik), intent(in)        :: dF_numbers(3)

    !> output
    real(kind=rk), intent(out)          :: norm_u

    ! loop variables
    integer(kind=ik)                    :: k, g, Bs, ix, iy, iz

!---------------------------------------------------------------------------------------------
! variables initialization

    norm_u = 0.0_rk

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all active heavy data
    do k = 1, hvy_n

        if ( params%threeD_case ) then
            ! 3D
            ! loop over all grid point
            do ix = g+1, Bs+g
                do iy = g+1, Bs+g
                    do iz = g+1, Bs+g
                        norm_u = max(norm_u, &
                                     norm2( (/hvy_block( ix, iy, iz, dF_numbers(1), hvy_active(k) ), &
                                              hvy_block( ix, iy, iz, dF_numbers(2), hvy_active(k) ), &
                                              hvy_block( ix, iy, iz, dF_numbers(3), hvy_active(k) ) /) ) &
                                    )
                    end do
                end do
            end do

        else
            ! 2D
            ! loop over all grid point
            do ix = g+1, Bs+g
                do iy = g+1, Bs+g
                    norm_u = max(norm_u, &
                                 norm2( (/hvy_block( ix, iy, 1, dF_numbers(1), hvy_active(k) ), &
                                          hvy_block( ix, iy, 1, dF_numbers(2), hvy_active(k) ) /) ) &
                                )
                end do
            end do

        end if

    end do

end subroutine get_block_max_velocity_norm
