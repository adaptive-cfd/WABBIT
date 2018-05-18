!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name get_inicond_from_file.f90
!> \version 0.5
!> \author sm
!
!> \brief call subroutines that read mesh and fields as initial condition from files
!
!>
!! input:
!!           - parameter array
!!
!! output:
!!           - light data array
!!           - heavy data array
!!           - number of active blocks (light and heavy)
!!           - time and iteration
!!
!!
!! = log ======================================================================================
!! \n
!! 28/09/17 - create
!
! ********************************************************************************************

subroutine get_inicond_from_file(params, lgt_block, hvy_block, hvy_n, lgt_n, time, iteration)
!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)        :: params
    !> light data array
    integer(kind=ik), intent(inout)       :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)          :: hvy_block(:, :, :, :, :)
    !> number of heavy and light active blocks
    integer(kind=ik), intent(inout)       :: hvy_n, lgt_n
    !> time loop variables
    real(kind=rk), intent(inout)          :: time
    integer(kind=ik), intent(inout)       :: iteration

    ! number of files to read from
    integer(kind=ik)                      :: N_files
    ! loop variable
    integer(kind=ik)                      :: k, dF
!---------------------------------------------------------------------------------------------
! variables initialization

    ! number of files to read from
    N_files = params%number_data_fields

!---------------------------------------------------------------------------------------------
! main body

    ! read treecode, time, iteration and total number of blocks from first input file
    call read_mesh_and_attributes(params%input_files(1), params, lgt_n, hvy_n, lgt_block, time, iteration)

    ! read datafields from files into hvy_block array
    do dF = 1, N_files
        call read_field(params%input_files(dF), dF, params, hvy_block, hvy_n )
    end do

end subroutine get_inicond_from_file
