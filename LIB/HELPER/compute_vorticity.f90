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
subroutine compute_vorticity(params, u, v, dx, vorticity, Bs, g)

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> physics parameter structure
    type (type_params), intent(in)             :: params
    !> origin and spacing of the block
    real(kind=rk), dimension(3), intent(in)    :: dx
    !> local datafields
    real(kind=rk), dimension(Bs+2*g,Bs+2*g), intent(in)  :: u, v
    !> vorticity
    real(kind=rk), dimension(Bs+2*g,Bs+2*g), intent(out) :: vorticity
    
    !> grid parameters
    integer(kind=ik), intent(in)               :: Bs, g
    !> derivatives
    real(kind=rk)                              :: u_dy, v_dx, dx_inv, dy_inv
    ! loop variables
    integer(kind=ik)                           :: ix, iy
    ! coefficients for Tam&Webb
    real(kind=rk)                              :: a(-3:3)
!---------------------------------------------------------------------------------------------
! variables initialization

    vorticity = 0.0_rk

    ! Tam & Webb, 4th order optimized (for first derivative)
    a = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, &
         0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)

    dx_inv = 1.0_rk / (2.0_rk*dx(1))
    dy_inv = 1.0_rk / (2.0_rk*dx(2))

!    Bs = params%number_block_nodes
!    g  = params%number_ghost_nodes
!---------------------------------------------------------------------------------------------
! main body

    if (params%order_discretization == "FD_2nd_central" ) then
        do ix = g+1, Bs+g
            do iy = g+1, Bs+g
                u_dy = (u(ix,iy+1)-u(ix,iy-1))*dy_inv
                v_dx = (v(ix+1,iy)-v(ix-1,iy))*dx_inv
                vorticity(ix,iy) = v_dx - u_dy
             end do
        end do
    else if (params%order_discretization == "FD_4th_central_optimized") then
        do ix = g+1, Bs+g
            do iy = g+1, Bs+g
                u_dy = (a(-3)*u(ix,iy-3) + a(-2)*u(ix,iy-2) + a(-1)*u(ix,iy-1) + a(0)*u(ix,iy)&
              +  a(+1)*u(ix,iy+1) + a(+2)*u(ix,iy+2) + a(+3)*u(ix,iy+3))*dy_inv
                v_dx = (a(-3)*v(ix-3,iy) + a(-2)*v(ix-2,iy) + a(-1)*v(ix-1,iy) + a(0)*v(ix,iy)&
              +  a(+1)*v(ix+1,iy) + a(+2)*v(ix+2,iy) + a(+3)*v(ix+3,iy))*dx_inv
                vorticity(ix,iy) = v_dx - u_dy
            end do
        end do
    else
        write(*,*) "ERROR: discretization method in params%order_discretization is unknown"
        write(*,*) params%order_discretization
        stop
    end if

end subroutine compute_vorticity
