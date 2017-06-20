!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_zeros.f90
!> \version 0.5
!> \author msr
!
!> \brief initialize zero for all fields \n
!! works for 2D and 3D data
!
!>
!! input:    - params \n
!! output:   - light and heavy data arrays \n
!!
!!
!! = log ======================================================================================
!! \n
!! 26/01/17 - create
!
! ********************************************************************************************

subroutine inicond_zeros( params, lgt_block, hvy_block )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)       :: params

    !> light data array
    integer(kind=ik), intent(inout)         :: lgt_block(:, :)

    !> heavy data array - block data
    real(kind=rk), intent(inout)            :: hvy_block(:, :, :, :, :)

    ! process rank
    integer(kind=ik)                        :: rank

    ! initial data field
    real(kind=rk), allocatable              :: phi(:, :, :)

    ! grid parameter, domainsize (Ds)
    integer(kind=ik)                        :: Ds

    ! number of datafields
    integer(kind=ik)                        :: dF


    ! loop variables
    integer(kind=ik)                        :: k

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank         = params%rank

    ! set parameters for readability
    Ds           = params%number_domain_nodes
    dF           = params%number_data_fields

!---------------------------------------------------------------------------------------------
! main body

    ! first: create start field
    !-----------------------------------------------------------------------------------------
    ! allocate memory
    if ( params%threeD_case ) then
        ! 3D:
        allocate( phi( Ds, Ds, Ds)  )
    else
        ! 2D:
        allocate( phi( Ds, Ds, 1)  )
    end if


    ! set phi
    phi = 0.0_rk

    ! output
    if (rank==0) then
        if ( params%unit_test .eqv. .false. ) then
            write(*,'(80("_"))')
            write(*,'("INIT: initialize zero for all fields ")')
        end if
    end if

    ! second: init light and heavy data for datafield 1, create starting block distribution
    !-----------------------------------------------------------------------------------------
    ! note: subroutine use allways a simple equal distribution
    ! for other distributions: use balance_load subroutine after initial distribution
    if ( params%threeD_case ) then
        ! 3D:
        call initial_block_distribution_3D( params, lgt_block, hvy_block, phi )
    else
        ! 2D:
        call initial_block_distribution_2D( params, lgt_block, hvy_block, phi )
    end if

    ! write heavy data for other datafields, copy first datafield
    do k = 2, dF
        hvy_block( :, :, :, k, : ) = hvy_block( :, :, :, 1, : )
    end do

    ! clean up
    deallocate( phi  )


end subroutine inicond_zeros
