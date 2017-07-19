!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name create_mask.f90
!> \version 0.5
!> \author engels, sm
!
!> \brief  \n
!> \note 
!
!> \details
!! input:    - params, mask, center, spacing of block \n
!! output:   - mask term for every grid point of this block
!!
!!
!! = log ======================================================================================
!! \n
!
! ********************************************************************************************

subroutine create_mask( params, mask, x0, dx, Bs, g )

    use module_params
    use module_precision

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                            :: params
    !> mask term for every grid point of this block
    real(kind=rk), dimension(2*g+Bs, 2*g+Bs), intent(inout)   :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(3), intent(in)                   :: x0, dx
    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g

    ! auxiliary variable for gauss pulse
    real(kind=rk)                                             :: x, y, R_cyl, cx, cy, r, h
    ! loop variables
    integer(kind=ik)                                          :: ix, iy

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! reset mask array
    mask = 0.0_rk

    ! place cylinder in the center of the domain
    cx = 0.5_rk * params%Lx
    cy = 0.5_rk * params%Ly

    ! radius of the cylinder
    R_cyl = params%inicond_width * params%Lx * params%Ly
    ! parameter for smoothing function (width)
    h = 1e-1_rk*R_cyl
    
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


end subroutine create_mask

subroutine smoothstep(f,x,t,h)
!-------------------------------------------------------------------------------
!> This subroutine returns the value f of a smooth step function \n
!> The sharp step function would be 1 if x<=t and 0 if x>t \n
!> h is the semi-size of the smoothing area, so \n
!> f is 1 if x<=t-h \n
!> f is 0 if x>t+h \n
!> f is variable (smooth) in between
!-------------------------------------------------------------------------------
    use module_precision
    
    implicit none
    real(kind=rk), intent(out) :: f
    real(kind=rk), intent(in)  :: x,t,h

        !-------------------------------------------------
        ! cos shaped smoothing (compact in phys.space)
        !-------------------------------------------------
        if (x<=t-h) then
          f = 1.0_rk
        elseif (((t-h)<x).and.(x<(t+h))) then
          f = 0.5_rk * (1.0_rk + dcos((x-t+h) * pi / (2.0_rk*h)) )
        else
          f = 0.0_rk
        endif
  
end subroutine smoothstep
