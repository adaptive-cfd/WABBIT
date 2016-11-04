! ********************************
! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: inicond_gauss_blob.f90
! version: 0.4
! author: engels, msr
!
! initialize gauss pulse
!
! input:    - grid parameter
! output:   - initial data field phi
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4
! ********************************************************************************************

subroutine inicond_gauss_blob(phi, Ds, Lx, Ly)

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! initial data field
    real(kind=rk), allocatable, intent(out)     :: phi(:, :)
    ! grid parameter, blocksize (Bs)
    integer(kind=ik), intent(in)                :: Ds
    ! domain length
    real(kind=rk), intent(in)                   :: Lx, Ly

    ! MPI error variable
    integer(kind=ik)                            :: ierr
    ! process rank
    integer(kind=ik)                            :: rank

    ! allocation error variable
    integer(kind=ik)                            :: allocate_error

    ! auxiliary variable for gauss pulse
    real(kind=rk)                               :: mux, muy, x ,y, sigma
    ! loop variables
    integer(kind=ik)                            :: i, j

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! allocate memory
    allocate( phi( Ds, Ds), stat=allocate_error )

    ! place pulse in the center of the domain
    mux = 0.5_rk * Lx;
    muy = 0.5_rk * Ly;

    ! pulse width
    sigma     = 0.1e2_rk

    ! create gauss pulse
    do i = 1, Ds
      do j = 1, Ds
        x = real(i-1, kind=rk)
        y = real(j-1, kind=rk)

        x = Lx / real(Ds-1, kind=rk) * x
        y = Ly / real(Ds-1, kind=rk) * y

        phi(i,j) = dexp( -( (x-mux)**2 + (y-muy)**2 ) / sigma )
      end do
    end do

    ! it sometimes causes bizarre effects not to delete extremely small numbers:
    ! so we do that now.
    where ( phi<1.0e-13_rk )
        phi = 0.0_rk
    end where

    ! output
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("INIT: initialize gauss pulse at x= ",f6.2," y= ", f6.2)') mux, muy
        write(*,'("INIT: with sigma= ",f6.2)') sigma
    end if

end subroutine inicond_gauss_blob
