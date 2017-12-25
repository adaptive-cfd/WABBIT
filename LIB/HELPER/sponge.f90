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
!! input:    - velocity, volume integral, grid parameters \n
!! output:   - forcing term \n
!!
!!
!! = log ======================================================================================
!! \n
!! 21/07/17 - create \n
!! 19/09/17 - add 3D
!*********************************************************************************************
subroutine sponge_2D(params, mask, x0, dx, Bs, g)
    use module_params
    use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> user defined parameter structure
    type (type_params), intent(in)                            :: params
    !> mask term for every grid point of this block
    real(kind=rk), dimension(2*g+Bs, 2*g+Bs), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx

    ! auxiliary variables
    real(kind=rk)                                             :: x, ddx
    ! loop variables
    integer(kind=ik)                                          :: ix, iy

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset mask array
    mask = 0.0_rk

!---------------------------------------------------------------------------------------------
! main body

    ddx = 0.1_rk*params%Lx


    do ix=1, Bs+2*g
       x = dble(ix-(g+1)) * dx(1) + x0(1)
       do iy=1, Bs+2*g
           if ((params%Lx-x) <= ddx) then
               mask(ix,iy) = (x-(params%Lx-ddx))**2
           elseif (x <= ddx) then
               mask(ix,iy) = (x-ddx)**2
           else
               mask(ix,iy) = 0.0_rk
           end if
       end do
    end do

end subroutine sponge_2D
