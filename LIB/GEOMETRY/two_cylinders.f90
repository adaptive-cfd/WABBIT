!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name two_cylinders.f90
!> \version 0.5
!> \author sm
!
!> \brief  creates geometry of two cylinders for volume penalization method \n
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
! ********************************************************************************************
subroutine two_cylinders(params, mask, x0, dx, Bs, g)

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
    real(kind=rk)                                             :: x1, x2, y1, y2, R, cx1, cx2, cy1,&
                                                                cy2, r_1, r_2, h, mask1, mask2
    ! loop variables
    integer(kind=ik)                                          :: ix, iy

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset mask array
    mask = 0.0_rk
    mask1 = 0.0_rk
    mask2 = 0.0_rk

!---------------------------------------------------------------------------------------------
! main body

    ! center of the first cylinder
    cx1 = 0.5884_rk*params%Lx
    cy1 = 0.4116_rk*params%Ly

    ! center of the second cylinder
    cx2 = 0.4116_rk*params%Lx
    cy2 = 0.5884_rk*params%Ly

    ! radius of the cylinders
    R = params%inicond_width
    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do ix=1, Bs+2*g
       x1 = dble(ix-(g+1)) * dx(1) + x0(1) - cx1
       x2 = dble(ix-(g+1)) * dx(1) + x0(1) - cx2
       do iy=1, Bs+2*g
           y1 = dble(iy-(g+1)) * dx(2) + x0(2) - cy1
           y2 = dble(iy-(g+1)) * dx(2) + x0(2) - cy2
           ! distance from center of cylinder 1
           r_1 = dsqrt(x1*x1 + y1*y1)
           ! distance from center of cylinder 2
           r_2 = dsqrt(x2*x2 + y2*y2) 
           if (params%smooth_mask) then
               call smoothstep(mask1, r_1, R, h)
               call smoothstep(mask2, r_2, R, h)
               mask(ix,iy) = mask1 + mask2
           else
               ! if point is inside one of the cylinders, set mask to 1
               if (r_1 <= R) then
                   mask(ix,iy) = 1.0_rk
               elseif ( r_2 <= R) then
                   mask(ix,iy) = 1.0_rk
               else
                   mask(ix,iy) = 0.0_rk
               end if
           end if
       end do
    end do


end subroutine two_cylinders
