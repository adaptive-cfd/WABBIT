!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_sinus_3D.f90
!> \version 0.5
!> \author msr
!
!> \brief initialize sinus for 3D case
!
!>
!! input:    - params \n
!! output:   - light and heavy data arrays \n
!!
!!
!! = log ======================================================================================
!! \n
!! 21/03/17 - create
!
! ********************************************************************************************

subroutine inicond_sinus_3D( params, lgt_block, hvy_block )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)       :: params

    !> light data array
    integer(kind=ik), intent(inout)         :: lgt_block(:, :)

    !> heavy data array - block data
    real(kind=rk), intent(inout)            :: hvy_block(:, :, :, :, :)

    ! initial data field
    real(kind=rk), allocatable              :: phi(:, :, :)

    ! grid parameter, domainsize (Ds)
    integer(kind=ik)                        :: Ds
    ! number of datafields
    integer(kind=ik)                        :: dF

    ! process rank
    integer(kind=ik)                        :: rank

    ! allocation error variable

    ! auxiliary variable for gauss pulse
    real(kind=rk)                           :: x ,y, z
    ! loop variables
    integer(kind=ik)                        :: i, j, k

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank         = params%rank

    ! set parameters for readability
    Ds           = params%number_domain_nodes
    dF           = params%number_data_fields

!---------------------------------------------------------------------------------------------
! main body

    ! first: create field phi
    !-----------------------------------------------------------------------------------------

    ! allocate memory
    allocate( phi( Ds, Ds, Ds)  )

    ! create sin
    do i = 1, Ds
      do j = 1, Ds
        do k = 1, Ds
          x = real(i-1, kind=rk)
          y = real(j-1, kind=rk)
          z = real(k-1, kind=rk)

          x = 10.0_rk*pi / real(Ds-1, kind=rk) * x
          y = 10.0_rk*pi / real(Ds-1, kind=rk) * y
          z = 10.0_rk*pi / real(Ds-1, kind=rk) * z

          phi(i,j,k) = sin(x)*sin(y)*sin(z)
        end do
      end do
    end do

    ! output
    if (rank==0) then
        if ( params%unit_test .eqv. .false. ) then
            write(*,'(80("_"))')
            write(*,'("INIT: initialize sin(x)*sin(y) ")')
        end if
    end if

    ! second: decompose init field phi to block data
    !-----------------------------------------------------------------------------------------
    ! init light and heavy data for datafield 1, create starting block distribution
    ! note: subroutine use allways a simple equal distribution
    ! for other distributions: use balance_load subroutine after initial distribution
    call initial_block_distribution_3D( params, lgt_block, hvy_block, phi )

    ! clean up
    deallocate( phi  )

end subroutine inicond_sinus_3D
