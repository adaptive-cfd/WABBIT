!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name sphere.f90
!> \version 0.5
!> \author sm
!
!> \brief  creates sphere geometry for the volume penalization method \n
!
!> \details
!! input:    - params, mask, center, spacing of block \n
!! output:   - mask term for every grid point of this block
!!
!!
!! = log ======================================================================================
!! \n
!! 21/11/17 - create
!
! ********************************************************************************************
subroutine sphere(params, mask, x0, dx, Bs, g )

    use module_params
    use module_precision

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                                    :: params
    !> mask term for every grid point of this block
    real(kind=rk), dimension(2*g+Bs, 2*g+Bs, 2*g+Bs), intent(inout)   :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(3), intent(in)                           :: x0, dx
    ! grid
    integer(kind=ik), intent(in)                                      :: Bs, g

    ! auxiliary variables
    real(kind=rk)                                                     :: x, y, z, R_sph, cx, cy, cz, r, h
    ! loop variables
    integer(kind=ik)                                                  :: ix, iy, iz

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset mask array
    mask = 0.0_rk

!---------------------------------------------------------------------------------------------
! main body

    ! place cylinder in the center of the domain
    cx = 0.5_rk * params%Lx
    cy = 0.5_rk * params%Ly
    cz = 0.5_rk * params%Lz
    ! radius of the sphere
    R_sph = params%inicond_width
    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))
    do ix=1, Bs+2*g
        x = dble(ix-(g+1)) * dx(1) + x0(1) - cx
        do iy=1, Bs+2*g
            y = dble(iy-(g+1)) * dx(2) + x0(2) - cy
            do iz = 1, Bs+2*g
                z = dble(iz-(g+1)) * dx(3) + x0(3) - cz
                ! distance from center of cylinder
                r = dsqrt(x**2 + y**2 + z**2)
                if (params%smooth_mask) then
                    call smoothstep(mask(ix,iy,iz), r, R_sph, h)
                else
                    ! if point is inside the cylinder, set mask to 1
                    if (r <= R_sph) then
                        mask(ix,iy,iz) = 1.0_rk
                    else
                        mask(ix,iy,iz) = 0.0_rk
                    end if
                end if
            end do
        end do
     end do

end subroutine sphere
