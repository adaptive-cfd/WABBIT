!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_richtmyer_meshkov.f90
!> \version 0.5
!> \author msr
!
!> \brief initialize richtmyer meshkov setup
!
!>
!! input:    - params \n
!! output:   - light and heavy data arrays \n
!!
!!
!! = log ======================================================================================
!! \n
!! 16/02/17 - create
!
! ********************************************************************************************

subroutine inicond_richtmyer_meshkov( params, lgt_block, hvy_block )

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
    ! domain length
    real(kind=rk)                           :: Lx, Ly, Lz

    ! process rank
    integer(kind=ik)                        :: rank

    ! allocation error variable
    integer(kind=ik)                        :: allocate_error

    ! loop variables
    integer(kind=ik)                        :: k

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank         = params%rank

    ! set parameters for readability
    Ds           = params%number_domain_nodes

    Lx           = params%Lx
    Ly           = params%Ly
    Lz           = params%Lz

    dF           = params%number_data_fields

!---------------------------------------------------------------------------------------------
! main body

    ! first: create field phi
    !-----------------------------------------------------------------------------------------

    ! allocate memory
    allocate( phi( Ds, Ds, Ds), stat=allocate_error )
    call check_allocation(allocate_error)

    phi = 0.0_rk

    ! it sometimes causes bizarre effects not to delete extremely small numbers:
    ! so we do that now.
    where ( phi<1.0e-13_rk )
       phi = 0.0_rk
    end where

    ! output
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("INIT: initialize richtmyer meshkov ")')
    end if

    ! second: decompose init field phi to block data
    !-----------------------------------------------------------------------------------------
    ! note: subroutine use allways a simple equal distribution
    ! for other distributions: use balance_load subroutine after initial distribution
    if ( params%threeD_case ) then
        ! 3D:
        call initial_block_distribution_3D( params, lgt_block, hvy_block, phi )
    else
        ! 2D:
        write(*,'(80("_"))')
        write(*,*) "ERROR: no richtmyer meshkov for 2D case"
        stop
    end if

    ! second: write heavy data for other datafields
    do k = 2, dF
        hvy_block( :, :, :, k, : ) = hvy_block( :, :, :, 1, : )
    end do

    ! navier stokes physics:
    ! set gauss blob + 1[bar] in datafield 6 (assume pressure)
    ! set velocity to zero and density to 1
    if ( params%physics_type == '3D_navier_stokes' ) then
        hvy_block( :, :, :, 1, : ) = 1.0_rk
        hvy_block( :, :, :, 2, : ) = 0.0_rk
        hvy_block( :, :, :, 3, : ) = 0.0_rk
        hvy_block( :, :, :, 4, : ) = 0.0_rk
        hvy_block( :, :, :, 5, : ) = 1.0_rk
    end if

    ! clean up
    deallocate( phi, stat=allocate_error )

end subroutine inicond_richtmyer_meshkov
