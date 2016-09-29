! ********************************
! 2D AMR prototype
! --------------------------------
!
! initialize gauss pulse
!
! name: inicond_gauss_blob.f90
! date: 31.08.2016
! author: engels, msr
! version: 0.1
!
! ********************************

subroutine inicond_gauss_blob()

      use module_params
      use module_blocks

      implicit none

      integer(kind=ik)      :: i, j, Ds
      real(kind=rk)         :: mux, muy, x ,y, Lx, Ly
      ! gauss parameter
      real(kind=rk)         :: sigma

      Ds        = blocks_params%size_domain
      Lx        = params%Lx
      Ly        = params%Ly

      ! place pulse in the center of the domain
      mux       = 0.5_rk * real(Ds+1, kind=rk)
      muy       = 0.5_rk * real(Ds+1, kind=rk)
      !muy       = 0.9_rk * real(Ds+1, kind=rk)

      !sigma     = 1.0e2_rk
      sigma     = 0.1e2_rk

      mux = 0.5_rk * Lx;
      muy = 0.5_rk * Ly;

      allocate( blocks_params%phi(1:Ds, 1:Ds) )
      blocks_params%phi = 0.0_rk

      write(*,'(a,f6.2,a,f6.2)') "initialize gauss pulse at x=", mux, " y=", muy
      write(*,'(a,e12.4)') "gauss sigma=", sigma
      write(*,'(80("-"))')

      do i = 1, Ds
          do j = 1, Ds
            x = real(i-1, kind=rk)
            y = real(j-1, kind=rk)

            x = Lx / real(Ds-1, kind=rk) * x
            y = Ly / real(Ds-1, kind=rk) * y

            blocks_params%phi(i,j) = dexp( -( (x-mux)**2 + (y-muy)**2 ) / sigma )
          end do
      end do

      ! it sometimes causes bizarre effects not to delete extremely small numbers:
      ! so we do that now.
      where ( blocks_params%phi<1.0e-13_rk )
        blocks_params%phi = 0.0_rk
      end where

end subroutine inicond_gauss_blob
