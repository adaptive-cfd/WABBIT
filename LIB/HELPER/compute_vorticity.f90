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
subroutine compute_vorticity(params, u, v, w, dx, vorticity)

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> physics parameter structure
    type (type_params), intent(in)                 :: params
    !> origin and spacing of the block
    real(kind=rk), dimension(3), intent(in)        :: dx
    !> local datafields
    real(kind=rk), dimension(:,:,:), intent(in)    :: u, v, w
    !> vorticity
    real(kind=rk), dimension(:,:,:,:), intent(out) :: vorticity
    
    !> grid parameters
    integer(kind=ik)                               :: Bs, g
    !> derivatives
    real(kind=rk)                                  :: u_dy, u_dz, v_dx, v_dz, w_dx, w_dy
    !> inverse of dx, dy, dz
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv
    ! loop variables
    integer(kind=ik)                               :: ix, iy, iz
    ! coefficients for Tam&Webb
    real(kind=rk)                                  :: a(-3:3)
!---------------------------------------------------------------------------------------------
! variables initialization

    vorticity = 0.0_rk

    ! Tam & Webb, 4th order optimized (for first derivative)
    a = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, &
         0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)

    dx_inv = 1.0_rk / (2.0_rk*dx(1))
    dy_inv = 1.0_rk / (2.0_rk*dx(2))

    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes
!---------------------------------------------------------------------------------------------
! main body

    if (params%threeD_case) then
        if (params%order_discretization == "FD_2nd_central" ) then
            do ix = g+1, Bs+g
                do iy = g+1, Bs+g
                    do iz = g+1, Bs+g
                        u_dy = (u(ix,iy+1,iz)-u(ix,iy-1,iz))*dy_inv
                        u_dz = (u(ix,iy,iz+1)-u(ix,iy,iz-1))*dz_inv
                        v_dx = (v(ix+1,iy,iz)-v(ix-1,iy,iz))*dx_inv
                        v_dz = (v(ix,iy,iz+1)-v(ix,iy,iz-1))*dz_inv
                        w_dx = (w(ix+1,iy,iz)-w(ix-1,iy,iz))*dx_inv
                        w_dy = (w(ix,iy+1,iz)-w(ix,iy-1,iz))*dy_inv
                    
                        vorticity(ix,iy,iz,1) = w_dy - v_dz
                        vorticity(ix,iy,iz,2) = u_dz - w_dx
                        vorticity(ix,iy,iz,3) = v_dx - u_dy
                    end do
                 end do
            end do
        else if (params%order_discretization == "FD_4th_central_optimized") then
            do ix = g+1, Bs+g
                do iy = g+1, Bs+g
                    do iz = g+1, Bs+g
                        u_dy = (a(-3)*u(ix,iy-3,iz) + a(-2)*u(ix,iy-2,iz) + a(-1)*u(ix,iy-1,iz) + a(0)*u(ix,iy,iz)&
                       +  a(+1)*u(ix,iy+1,iz) + a(+2)*u(ix,iy+2,iz) + a(+3)*u(ix,iy+3,iz))*dy_inv
                        u_dz = (a(-3)*u(ix,iy,iz-3) + a(-2)*u(ix,iy,iz-2) + a(-1)*u(ix,iy,iz-1) + a(0)*u(ix,iy,iz)&
                      +  a(+1)*u(ix,iy,iz+1) + a(+2)*u(ix,iy,iz+2) + a(+3)*u(ix,iy,iz+3))*dz_inv
                        v_dx = (a(-3)*v(ix-3,iy,iz) + a(-2)*v(ix-2,iy,iz) + a(-1)*v(ix-1,iy,iz) + a(0)*v(ix,iy,iz)&
                      +  a(+1)*v(ix+1,iy,iz) + a(+2)*v(ix+2,iy,iz) + a(+3)*v(ix+3,iy,iz))*dx_inv
                        v_dz = (a(-3)*v(ix,iy,iz-3) + a(-2)*v(ix,iy,iz-2) + a(-1)*v(ix,iy,iz-1) + a(0)*v(ix,iy,iz)&
                      +  a(+1)*v(ix,iy,iz+1) + a(+2)*v(ix,iy,iz+2) + a(+3)*v(ix,iy,iz+3))*dz_inv
                        w_dx = (a(-3)*w(ix-3,iy,iz) + a(-2)*w(ix-2,iy,iz) + a(-1)*w(ix-1,iy,iz) + a(0)*w(ix,iy,iz)&
                      +  a(+1)*w(ix+1,iy,iz) + a(+2)*w(ix+2,iy,iz) + a(+3)*w(ix+3,iy,iz))*dx_inv
                        w_dy = (a(-3)*w(ix,iy-3,iz) + a(-2)*w(ix,iy-2,iz) + a(-1)*w(ix,iy-1,iz) + a(0)*w(ix,iy,iz)&
                     +  a(+1)*w(ix,iy+1,iz) + a(+2)*w(ix,iy+2,iz) + a(+3)*w(ix,iy+3,iz))*dy_inv
                       
                        vorticity(ix,iy,iz,1) = w_dy - v_dz
                        vorticity(ix,iy,iz,2) = u_dz - w_dx
                        vorticity(ix,iy,iz,3) = v_dx - u_dy
                    end do
                end do
            end do
        else
            write(*,*) "ERROR: discretization method in params%order_discretization is unknown"
            write(*,*) params%order_discretization
            stop
        end if
    else
        if (params%order_discretization == "FD_2nd_central" ) then
            do ix = g+1, Bs+g
                do iy = g+1, Bs+g
                    u_dy = (u(ix,iy+1,1)-u(ix,iy-1,1))*dy_inv
                    v_dx = (v(ix+1,iy,1)-v(ix-1,iy,1))*dx_inv
                    vorticity(ix,iy,1,1) = v_dx - u_dy
                 end do
            end do
        else if (params%order_discretization == "FD_4th_central_optimized") then
            do ix = g+1, Bs+g
                do iy = g+1, Bs+g
                    u_dy = (a(-3)*u(ix,iy-3,1) + a(-2)*u(ix,iy-2,1) + a(-1)*u(ix,iy-1,1) + a(0)*u(ix,iy,1)&
                  +  a(+1)*u(ix,iy+1,1) + a(+2)*u(ix,iy+2,1) + a(+3)*u(ix,iy+3,1))*dy_inv
                    v_dx = (a(-3)*v(ix-3,iy,1) + a(-2)*v(ix-2,iy,1) + a(-1)*v(ix-1,iy,1) + a(0)*v(ix,iy,1)&
                  +  a(+1)*v(ix+1,iy,1) + a(+2)*v(ix+2,iy,1) + a(+3)*v(ix+3,iy,1))*dx_inv
                    vorticity(ix,iy,1,1) = v_dx - u_dy
                end do
            end do
        else
            write(*,*) "ERROR: discretization method in params%order_discretization is unknown"
            write(*,*) params%order_discretization
            stop
        end if
    end if

end subroutine compute_vorticity
