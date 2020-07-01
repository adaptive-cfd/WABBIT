!> \file
! WABBIT
!> \name gradient.f90
!> \version 0.1
!> \author dk
!
!> \brief computation of the gradient of a field
!
! = log ======================================================================================
!-----------------------------------------------------------------------------------------------------
subroutine gradient(fld, dx, Bs, g, discretization, grad_fld)

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> origin and spacing of the block
    real(kind=rk), dimension(3), intent(in)        :: dx
    !> local datafields
    real(kind=rk), dimension(:,:,:), intent(in)    :: fld
    !> vorticity
    real(kind=rk), dimension(:,:,:,:), intent(out) :: grad_fld
    character(len=*), intent(in)                   :: discretization
    !> grid parameters
    integer(kind=ik), intent(in)                   :: g
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    !> derivatives
    real(kind=rk)                                  :: fld_dx, fld_dy, fld_dz
    !> inverse of dx, dy, dz
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv
    ! loop variables
    integer(kind=ik)                               :: ix, iy, iz
    ! coefficients for Tam&Webb
    real(kind=rk)                                  :: a(-3:3)
!---------------------------------------------------------------------------------------------
! variables initialization

    fld_dx = 0.0_rk
    fld_dy = 0.0_rk
    fld_dz = 0.0_rk

    ! Tam & Webb, 4th order optimized (for first derivative)
    a = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, &
        0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)

!---------------------------------------------------------------------------------------------
! main body

    if (size(fld,3)>2) then
        !-----------------------------------------------------------------------
        ! 3D case 3D case 3D case 3D case
        !-----------------------------------------------------------------------
        dz_inv = 1.0_rk / dx(3)
        if (discretization == "FD_2nd_central" ) then
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    do iz = g+1, Bs(3)+g
                        fld_dx = (fld(ix+1,iy,iz)-fld(ix-1,iy,iz))*dx_inv*0.5_rk
                        fld_dy = (fld(ix,iy+1,iz)-fld(ix,iy-1,iz))*dy_inv*0.5_rk
                        fld_dz = (fld(ix,iy,iz+1)-fld(ix,iy,iz-1))*dz_inv*0.5_rk

                        grad_fld(ix,iy,iz,1) = fld_dx
                        grad_fld(ix,iy,iz,2) = fld_dy
                        grad_fld(ix,iy,iz,3) = fld_dz
                    end do
                end do
            end do
        else if (discretization == "FD_4th_central_optimized") then
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    do iz = g+1, Bs(3)+g
                        fld_dx = (a(-3)*fld(ix-3,iy,iz) + a(-2)*fld(ix-2,iy,iz) &
                        + a(-1)*fld(ix-1,iy,iz) + a(0)*fld(ix,iy,iz)&
                        + a(+1)*fld(ix+1,iy,iz) + a(+2)*fld(ix+2,iy,iz) &
                        + a(+3)*fld(ix+3,iy,iz))*dx_inv

                        fld_dy = (a(-3)*fld(ix,iy-3,iz) + a(-2)*fld(ix,iy-2,iz) &
                        + a(-1)*fld(ix,iy-1,iz) + a(0)*fld(ix,iy,iz)&
                        + a(+1)*fld(ix,iy+1,iz) + a(+2)*fld(ix,iy+2,iz) &
                        + a(+3)*fld(ix,iy+3,iz))*dy_inv
                        
                        fld_dz = (a(-3)*fld(ix,iy,iz-3) + a(-2)*fld(ix,iy,iz-2) &
                        + a(-1)*fld(ix,iy,iz-1) + a(0)*fld(ix,iy,iz) &
                        + a(+1)*fld(ix,iy,iz+1) + a(+2)*fld(ix,iy,iz+2) &
                        + a(+3)*fld(ix,iy,iz+3))*dz_inv

                        grad_fld(ix,iy,iz,1) = fld_dx
                        grad_fld(ix,iy,iz,2) = fld_dy
                        grad_fld(ix,iy,iz,3) = fld_dz
                    end do
                end do
            end do
        else
            call abort(0112186, "ERROR: discretization method in discretization is unknown")
        end if
    else
        !-----------------------------------------------------------------------
        ! 2D case 2D case 2D case 2D case
        !-----------------------------------------------------------------------
        if (discretization == "FD_2nd_central" ) then
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    fld_dx = (fld(ix+1,iy,1)-fld(ix-1,iy,1))*dx_inv*0.5_rk
                    fld_dy = (fld(ix,iy+1,1)-fld(ix,iy-1,1))*dy_inv*0.5_rk
                    grad_fld(ix,iy,1,1) = fld_dx
                    grad_fld(ix,iy,1,2) = fld_dy
                end do
            end do
        else if (discretization == "FD_4th_central_optimized") then
            do ix = g+1, Bs(1)+g
                do iy = g+1, Bs(2)+g
                    fld_dx = (a(-3)*fld(ix-3,iy,1) + a(-2)*fld(ix-2,iy,1) + &
                    a(-1)*fld(ix-1,iy,1) + a(0)*fld(ix,iy,1)&
                    +  a(+1)*fld(ix+1,iy,1) + a(+2)*fld(ix+2,iy,1) + a(+3)*fld(ix+3,iy,1))*dx_inv

                    fld_dy = (a(-3)*fld(ix,iy-3,1) + a(-2)*fld(ix,iy-2,1) + &
                    a(-1)*fld(ix,iy-1,1) + a(0)*fld(ix,iy,1)&
                    +  a(+1)*fld(ix,iy+1,1) + a(+2)*fld(ix,iy+2,1) + a(+3)*fld(ix,iy+3,1))*dy_inv

                    grad_fld(ix,iy,1,1) = fld_dx
                    grad_fld(ix,iy,1,2) = fld_dy
                end do
            end do
        else
            call abort(0112187, "ERROR: discretization method in discretization is unknown")
        end if
    end if

end subroutine gradient
