!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_vorticity_filaments.f90
!> \version 0.4
!> \author msr
!
!> \brief initialize the vorticity_filaments testcase \n
!! load velocity fields from file
!
!>
!! input:    - params, light and heavy data \n
!! output:   - params, light and heavy data \n
!!
!!
!! = log ======================================================================================
!! \n
!! 04/11/16 - switch to v0.4
! ********************************************************************************************

subroutine inicond_vorticity_filaments(params, lgt_block, hvy_block)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)       :: params

    !> light data array
    integer(kind=ik), intent(inout)         :: lgt_block(:, :)

    !> heavy data array - block data
    real(kind=rk), intent(inout)            :: hvy_block(:, :, :, :)

    ! field from file
    real(kind=rk), allocatable              :: phi(:, :), Ux(:, :, :)

    ! allocation error variable
    integer(kind=ik)                        :: allocate_error

    ! file IO error variable
    integer(kind=ik)                        :: io_error

!---------------------------------------------------------------------------------------------
! variables initialization

    ! allocate memory
    allocate( phi( params%number_domain_nodes, params%number_domain_nodes), stat=allocate_error )
    allocate( Ux( size(hvy_block,1), size(hvy_block,2), size(hvy_block,4)), stat=allocate_error )

!---------------------------------------------------------------------------------------------
! main body

    ! reset domain size
    params%Lx = 2.0_rk*pi
    params%Ly = 2.0_rk*pi

    ! ----------------------------------------------------------------------------------------
    ! read Ux velocity
    open(unit=99, file="Ux.start", status='old', action='read', iostat=io_error)
    read(unit=99, fmt=*) phi

    ! it sometimes causes bizarre effects not to delete extremely small numbers:
    ! so we do that now.
    where ( phi<1.0e-13_rk )
        phi = 0.0_rk
    end where

    ! decompose init field phi to block data
    ! first: init light and heavy data for datafield 1, create starting block distribution
    !call initial_block_distribution( params, lgt_block, hvy_block, phi )

    ! phi is on first datafield, save heavy data
    Ux = hvy_block( :, :, 2, : )

    ! ----------------------------------------------------------------------------------------
    ! reset data
    lgt_block = -1
    hvy_block = 0.0_rk

    ! read Uy velocity
    open(unit=99, file="Uy.start", status='old', action='read', iostat=io_error)
    read(unit=99, fmt=*) phi

    ! it sometimes causes bizarre effects not to delete extremely small numbers:
    ! so we do that now.
    where ( phi<1.0e-13_rk )
        phi = 0.0_rk
    end where

    ! decompose init field phi to block data
    ! first: init light and heavy data for datafield 1, create starting block distribution
    !call initial_block_distribution( params, lgt_block, hvy_block, phi )

    ! phi is now on first field, write to right position
    hvy_block( :, :, 4, : ) = hvy_block( :, :, 2, : )

    ! write Ux velocity
    hvy_block( :, :, 3, : ) = Ux

    ! clean up
    deallocate( phi, stat=allocate_error )
    deallocate( Ux, stat=allocate_error )

end subroutine inicond_vorticity_filaments
