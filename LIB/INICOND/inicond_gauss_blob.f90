subroutine inicond_gauss_blob()
      use module_params
      use module_blocks

      implicit none

      integer(kind=ik) :: i, j, nx
      real(kind=rk) :: mux, muy, x ,y, lx
      ! gauss parameter
      real(kind=rk) :: sigma

      nx = blocks_params%size_domain
      lx = real(nx, kind=rk)
      mux     = 0.5_rk !* ( real(nx,kind=rk) - 1.0_rk )
      muy     = 0.5_rk !* ( real(nx,kind=rk) - 1.0_rk )
      sigma   = 2.0e-3_rk
      allocate( blocks_params%phi(1:nx, 1:nx) )
      blocks_params%phi = 0.0_rk

      write(*,'(a,f4.2,a,f4.2)') "initialize gauss pulse at x=", mux, " y=", muy
      write(*,'(a,e10.4)') "gauss sigma=", sigma
      write(*,'(80("-"))')

      do i = 1, nx
          do j = 1, nx
            ! x and y are normalized to the domain length
            x = real(i-1,kind=rk) / (lx-1.0_rk)
            y = real(j-1,kind=rk) / (lx-1.0_rk)
            blocks_params%phi(i,j) = exp( -((x-mux)**2 + (y-muy)**2) / sigma )
          end do
      end do

      ! it sometimes causes bizarre effects not to delete extremely small numbers:
      ! so we do that now.
      where ( blocks_params%phi<1.0e-13_rk )
        blocks_params%phi = 0.0_rk
      end where

end subroutine
