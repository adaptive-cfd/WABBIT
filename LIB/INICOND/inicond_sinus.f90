subroutine inicond_sinus()
      use module_params
      use module_blocks

      implicit none

      integer(kind=ik) :: i, j, nx
      real(kind=rk) :: mux, muy, x ,y, lx

      nx = blocks_params%size_domain
      allocate( blocks_params%phi(1:nx, 1:nx) )
      blocks_params%phi = 0.0_rk

      lx = real(nx, kind=rk)
      mux     = 0.5_rk
      muy     = 0.5_rk

      do i = 1, nx
          do j = 1, nx
            ! x and y are normalized to the domain length
            x = real(i-1,kind=rk) / (lx-1.0_rk)
            y = real(j-1,kind=rk) / (lx-1.0_rk)
            blocks_params%phi(i,j) = dsin(2.0_rk*6.283185307179586e+00_rk*x) * dsin(2.0_rk*6.283185307179586e+00_rk*y)
          end do
      end do

end subroutine
