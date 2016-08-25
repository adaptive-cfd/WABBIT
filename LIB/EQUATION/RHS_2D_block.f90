! ********************************
! 2D AMR prototype
! --------------------------------
! 
! RHS, 2D, convection-diffusion equation
!
! name: RHS_2D_block.f90
! date: 05.08.2016
! author: msr
! version: 0.1
! 
! ********************************

subroutine RHS_2D_block(phi, dx, dy, g, N)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)					        :: g, N
    real(kind=rk), intent(in)                               :: dx, dy
    real(kind=rk), dimension(N+2*g, N+2*g), intent(inout)	:: phi

    real(kind=rk), dimension(N+2*g, N+2*g)			        :: grad_phi, laplace_phi

    grad_phi 		= 0.0_rk
    laplace_phi		= 0.0_rk

    grad_phi 		= params%u0(1) * matmul(blocks_params%D1, phi) / dx &
                      + params%u0(2) * matmul(phi, transpose(blocks_params%D1)) / dy

    laplace_phi 	= matmul(blocks_params%D2, phi) / (dx*dx) &
                    + matmul(phi, transpose(blocks_params%D2)) / (dy*dy)

    phi 			= grad_phi + params%nu * laplace_phi

end subroutine RHS_2D_block
