subroutine initial_condition_dense_field()
      use module_params
      use module_blocks

      implicit none

      integer(kind=ik) :: i, j, nx
      real(kind=rk) :: mux, muy, x ,y, lx
      ! gauss parameter
      real(kind=rk) :: sigma

      nx = blocks_params%size_domain
      lx = real(nx, kind=rk)
      mux     = 0.5_rk * ( real(nx,kind=rk) - 1.0_rk )
      muy     = 0.5_rk * ( real(nx,kind=rk) - 1.0_rk )
      sigma   = 1.0_rk
      allocate( blocks_params%phi(1:nx, 1:nx) )
      blocks_params%phi = 0.0_rk

      write(*,*) mux,muy,sigma

      do i = 1, nx
          do j = 1, nx
            x = real(i+1,kind=rk)
            y = real(j+1,kind=rk)
            blocks_params%phi(i,j) = exp( -((x-mux)**2 + (y-muy)**2) / sigma )
            blocks_params%phi(i,j) = dsin(6.283185307179586e+00_rk*x/lx) * dsin(6.283185307179586e+00_rk*y/lx)
          end do
      end do

      ! it sometimes causes bizarre effects not to delete extremely small numbers:
      ! so we do that now.
      ! where ( blocks_params%phi<1.0e-13_rk )
      !   blocks_params%phi = 0.0_rk
      ! end where
end subroutine
