! ********************************
! 2D AMR prototype
! --------------------------------
!
! calculate timestep
!
! name: calc_dt.f90
! date: 04.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine calc_dt(dt)

    use module_params
    use module_blocks

    implicit none

    real(kind=rk), intent(inout) 	:: dt

    real(kind=rk)                   :: dx
    integer                         :: k, N, block_num

    N = size(blocks_params%active_list, dim=1)

    ! minimal spacing
    dx = params%time_max
    do k = 1, N
        block_num = blocks_params%active_list(k)
        dx = min(dx, blocks(block_num)%dx )
    end do

    ! time step
    dt = params%CFL * dx / norm2(params%u0)
end subroutine calc_dt
