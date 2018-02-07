!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name sponge.f90
!> \version 0.5
!> \author sm
!
!> \brief
!
!>
!! input:    - params, origin and spacing of the block, grid parameters \n
!! output:   - sponge term \n
!!
!!
!! = log ======================================================================================
!! \n
!! 12/18 - create
!*********************************************************************************************
subroutine sponge_2D_NEW(sponge, x0, dx, Bs, g)

    implicit none

    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> sponge term for every grid point of this block
    real(kind=rk), dimension(2*g+Bs, 2*g+Bs), intent(out)     :: sponge
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx

    ! auxiliary variables
    real(kind=rk)                                             :: x, ddx
    ! loop variables
    integer(kind=ik)                                          :: ix, iy

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset sponge array
    sponge = 0.0_rk
!---------------------------------------------------------------------------------------------
! main body
    ddx = 0.1_rk*params_acm%Lx

    do ix=1, Bs+2*g
       x = dble(ix-(g+1)) * dx(1) + x0(1)
       do iy=1, Bs+2*g
           if ((params_acm%Lx-x) <= ddx) then
               sponge(ix,iy) = (x-(params_acm%Lx-ddx))**2
           elseif (x <= ddx) then
               sponge(ix,iy) = (x-ddx)**2
           else
               sponge(ix,iy) = 0.0_rk
           end if
       end do
    end do

end subroutine sponge_2D_NEW
