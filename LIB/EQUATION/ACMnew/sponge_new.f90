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
    real(kind=rk)                                             :: x, y, tmp
    ! loop variables
    integer(kind=ik)                                          :: ix, iy


    ! reset sponge array
    sponge = 0.0_rk

    do iy=1, Bs+2*g
      y = dble(iy-(g+1)) * dx(2) + x0(2)
       do ix=1, Bs+2*g
           x = dble(ix-(g+1)) * dx(1) + x0(1)
           ! distance to borders of domain
           tmp = minval( (/x,y,-(x-params_acm%Lx),-(y-params_acm%Ly)/) )

           call smoothstep(sponge(ix,iy), tmp, 0.5_rk*params_acm%L_sponge, 0.5_rk*params_acm%L_sponge)
       end do
    end do

end subroutine sponge_2D_NEW
