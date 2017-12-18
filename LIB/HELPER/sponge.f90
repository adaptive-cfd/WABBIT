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
!! input:    - velocity,  grid parameters \n
!! output:   - sponge term \n
!!
!!
!! = log ======================================================================================
!! \n
!! 27/11/17 - create \n

!*********************************************************************************************
subroutine sponge_2D(params, sponge, x0, dx, Bs, g)
    use module_params
    use module_precision

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                            :: params
    !> sponge term for every grid point of this block
    real(kind=rk), dimension(2*g+Bs, 2*g+Bs), intent(out)     :: sponge
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx
    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g

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

    ddx = 0.1_rk*params%Lx
    
    
    do ix=1, Bs+2*g
       x = dble(ix-(g+1)) * dx(1) + x0(1)
       do iy=1, Bs+2*g    
           if ((params%Lx-x) <= ddx) then
               sponge(ix,iy) = (x-(params%Lx-ddx))**2
           elseif (x <= ddx) then
               sponge(ix,iy) = (x-ddx)**2
           else
               sponge(ix,iy) = 0.0_rk
           end if
       end do
    end do

end subroutine sponge_2D


    
