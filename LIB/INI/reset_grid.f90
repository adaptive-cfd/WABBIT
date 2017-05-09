!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name init_data.f90
!> \version 0.5
!> \author msr
!
!> \brief reset grid, set all blocks to empty
!
!>
!! input:    
!!           - parameter array
!!           - light data array
!!           - heavy data array
!!           - neighbor data array
!!           - light and heavy active block list
!!
!! output:   
!!           - filled user defined data structure for global params
!!           - initialized light and heavy data arrays
!!
!! = log ======================================================================================
!! \n
!! 04/11/16 - switch to v0.4, now run complete initialization within these subroutine and return
!            initialized block data to main program \n
!! 07/12/16 - now uses heavy work data array \n
!! 25/01/17 - switch to 3D, v0.5
!
! ********************************************************************************************
subroutine reset_grid(params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, verbosity )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)       :: params
    !> light data array
    integer(kind=ik),  intent(inout)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)
    !> heavy work array  )
    real(kind=rk),  intent(inout)           :: hvy_work(:, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik),  intent(inout)        :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)         :: lgt_n
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)         :: hvy_n
    !> write output
    logical, intent(in)                     :: verbosity

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    if ( (params%rank == 0) .and. verbosity ) then
      write(*,'(80("_"))')
      write(*,'(A)') "RESET: resetting grid to empty (deactivate all blocks)."
    endif

    ! reset data:
    ! all blocks are inactive, reset treecode
    lgt_block(:, 1:params%max_treelevel) = -1
    ! all blocks are inactive, reset treelevel
    lgt_block(:, params%max_treelevel+1) = -1
    ! set refinement to 0
    lgt_block(:, params%max_treelevel+2) = 0


    ! reset data
    hvy_block = 9.99e99_rk
    hvy_work = 9.99e99_rk
    hvy_neighbor = -1
    ! as the grid has changed (we deleted it here), we now update the heavy and light
    ! active lists
    call create_lgt_active_list( lgt_block, lgt_active, lgt_n )
    call create_hvy_active_list( lgt_block, hvy_active, hvy_n )

end subroutine reset_grid
