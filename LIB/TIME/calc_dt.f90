! ********************************
! WABBIT
! --------------------------------
!
! calculate timestep
!
! name: calc_dt.f90
! date: 30.09.2016
! author: msr
! version: 0.2
!
! ********************************

subroutine calc_dt(dt)

    use module_params
    use module_blocks

    implicit none

    real(kind=rk), intent(inout) 	:: dt

    real(kind=rk)                   :: dx
    integer                         :: k, N

    N = blocks_params%number_max_blocks

    ! minimal spacing
    dx = 9.0e9_rk

    ! loop over all blocks
    do k = 1, N
        if (blocks(k)%active) then
            dx = min(dx, blocks(k)%dx )
        end if
    end do

    ! time step
    dt = params%CFL * dx / norm2(params%u0)

end subroutine calc_dt
