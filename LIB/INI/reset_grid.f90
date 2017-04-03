! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: init_data.f90
! version: 0.5
! author: msr
!
! This routine initializes the block data, i.e. it evaluates the initial condition on the grid
!
! input:    - parameter array
!           - light data array
!           - heavy data array
!           - neighbor data array
!           - light and heavy active block list
! output:   - filled user defined data structure for global params
!           - initialized light and heavy data arrays
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4, now run complete initialization within these subroutine and return
!            initialized block data to main program
! 07/12/16 - now uses heavy work data array
! 25/01/17 - switch to 3D, v0.5
!
! ********************************************************************************************
subroutine reset_grid(params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(inout)     :: params
    ! light data array
    integer(kind=ik),  intent(inout)      :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk),  intent(inout)         :: hvy_block(:, :, :, :, :)
    ! heavy work array  )
    real(kind=rk),  intent(inout)         :: hvy_work(:, :, :, :, :)
    ! neighbor array (heavy data)
    integer(kind=ik),  intent(inout)      :: hvy_neighbor(:,:)
    ! list of active blocks (light data)
    integer(kind=ik),  intent(inout)      :: lgt_active(:)
    ! list of active blocks (light data)
    integer(kind=ik),  intent(inout)      :: hvy_active(:)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    if (params%rank == 0) then
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
    lgt_active = -1
    hvy_active = -1

end subroutine reset_grid
