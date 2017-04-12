! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: inicond_sinus_2D.f90
! version: 0.5
! author: msr
!
! initialize sinus for 2D case
!
! input:    - params
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 21/03/17 - create
!
! ********************************************************************************************

subroutine inicond_sinus_2D( params, u, x0, dx )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(inout)       :: params

    ! actual block data (note this routine acts only on one block)
    real(kind=rk), intent(inout)            :: u(:,:,:,:)

    ! spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3),dx(1:3)

    ! auxiliary variable for gauss pulse
    real(kind=rk)                           :: x ,y
    ! loop variables
    integer(kind=ik)                        :: ix, iy
    ! grid
    integer(kind=ik)                        :: Bs, g

!---------------------------------------------------------------------------------------------
! variables initialization
    ! set parameters for readability
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

!---------------------------------------------------------------------------------------------
! main body

    ! create sinus2d
    do ix = g+1, Bs+g
      do iy = g+1, Bs+g
        ! compute x,y coordinates from spacing and origin
        x = dble(ix-(g+1)) * dx(1) + x0(1)
        y = dble(iy-(g+1)) * dx(2) + x0(2)
        ! set actual inicond sinus2d
        ! FIXME: df=2 ...
        u(ix,iy,1,2) = sin(1.0_rk*2.0_rk*pi*x/params%Lx)*sin(2.0_rk*pi*y/params%Ly)
      end do
    end do

end subroutine inicond_sinus_2D
