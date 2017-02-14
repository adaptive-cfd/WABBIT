! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: inicond_gauss_blob.f90
! version: 0.5
! author: engels, msr
!
! initialize gauss pulse for 2D case
! note: field phi is 3D, but third dimension is not used
!
! input:    - params
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4
! 26/01/17 - use process rank from params struct
!          - use v0.5 hvy data array
!
! ********************************************************************************************

subroutine inicond_gauss_blob( params, lgt_block, hvy_block )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(inout)       :: params

    ! light data array
    integer(kind=ik), intent(inout)         :: lgt_block(:, :)

    ! heavy data array - block data
    real(kind=rk), intent(inout)            :: hvy_block(:, :, :, :, :)

    ! initial data field
    real(kind=rk), allocatable              :: phi(:, :, :)

    ! grid parameter, domainsize (Ds)
    integer(kind=ik)                        :: Ds
    ! number of datafields
    integer(kind=ik)                        :: dF
    ! domain length
    real(kind=rk)                           :: Lx, Ly

    ! process rank
    integer(kind=ik)                        :: rank

    ! allocation error variable
    integer(kind=ik)                        :: allocate_error

    ! auxiliary variable for gauss pulse
    real(kind=rk)                           :: mux, muy, x ,y, sigma
    ! loop variables
    integer(kind=ik)                        :: i, j, k

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank         = params%rank

    ! set parameters for readability
    Ds           = params%number_domain_nodes
    Lx           = params%Lx
    Ly           = params%Ly
    dF           = params%number_data_fields

!---------------------------------------------------------------------------------------------
! main body

    ! first: create field phi
    !-----------------------------------------------------------------------------------------

    ! allocate memory
    allocate( phi( Ds, Ds, 1), stat=allocate_error )
    call check_allocation(allocate_error)

    ! place pulse in the center of the domain
    mux = 0.5_rk * Lx;
    muy = 0.5_rk * Ly;

    ! pulse width
    sigma     = 0.00001e2_rk
    !sigma     = 0.1e2_rk

    ! create gauss pulse
    do i = 1, Ds
      do j = 1, Ds
        x = real(i-1, kind=rk)
        y = real(j-1, kind=rk)

        x = Lx / real(Ds-1, kind=rk) * x
        y = Ly / real(Ds-1, kind=rk) * y

        phi(i,j,1) = dexp( -( (x-mux)**2 + (y-muy)**2 ) / sigma )
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

    ! second: decompose init field phi to block data
    !-----------------------------------------------------------------------------------------
    ! first: init light and heavy data for datafield 1, create starting block distribution
    ! note: subroutine use allways a simple equal distribution
    ! for other distributions: use balance_load subroutine after initial distribution
    call initial_block_distribution_2D( params, lgt_block, hvy_block, phi )

    ! second: write heavy data for other datafields
    do k = 3, dF+1
        hvy_block( :, :, :, k, : ) = hvy_block( :, :, :, 2, : )
    end do

    ! navier stokes physics:
    ! set gauss blob + 1[bar] in datafield 5 (assume pressure)
    ! set velocity to zero and density to 1
    if ( params%physics_type == '2D_navier_stokes' ) then
        hvy_block( :, :, :, 2, : ) = 1.0_rk
        hvy_block( :, :, :, 3, : ) = 0.0_rk
        hvy_block( :, :, :, 4, : ) = 0.0_rk
        hvy_block( :, :, :, 5, : ) = hvy_block( :, :, :, 5, : ) + 1e5_rk
    end if

    ! clean up
    deallocate( phi, stat=allocate_error )

end subroutine inicond_gauss_blob
