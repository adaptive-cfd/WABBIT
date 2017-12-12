!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name cylinder.f90
!> \version 0.5
!> \author sm
!
!> \brief  geometry of one cylinder for volume penalization method \n
!> \note 
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

subroutine cylinder(params, mask, x0, dx, Bs, g )

    use module_params
    use module_precision

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                            :: params
    !> mask term for every grid point of this block
    real(kind=rk), dimension(2*g+Bs, 2*g+Bs), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx
    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g

    ! auxiliary variables
    real(kind=rk)                                             :: x, y, R_cyl, cx, cy, r, h
    ! loop variables
    integer(kind=ik)                                          :: ix, iy

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset mask array
    mask = 0.0_rk

!---------------------------------------------------------------------------------------------
! main body

    ! place cylinder in the center of the domain
    cx = 0.5_rk * params%Lx
!    cx = 4.7_rk
    cy = 0.5_rk * params%Ly

    ! radius of the cylinder
    R_cyl = params%inicond_width
    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do ix=1, Bs+2*g
       x = dble(ix-(g+1)) * dx(1) + x0(1) - cx
       do iy=1, Bs+2*g
           y = dble(iy-(g+1)) * dx(2) + x0(2) - cy
           ! distance from center of cylinder
           r = dsqrt(x*x + y*y)
           if (params%smooth_mask) then
               call smoothstep(mask(ix,iy), r, R_cyl, h)
           else
               ! if point is inside the cylinder, set mask to 1
               if (r <= R_cyl) then
                   mask(ix,iy) = 1.0_rk
               else
                   mask(ix,iy) = 0.0_rk
               end if
           end if
       end do
    end do

end subroutine cylinder
