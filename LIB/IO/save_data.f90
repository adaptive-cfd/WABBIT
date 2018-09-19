!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name save_data.f90
!> \version 0.5
!> \author msr
!
!> \brief save data main function, call write data routine
!
!>
!! input:
!!           - time loop parameter
!!           - parameter array
!!           - light data array
!!           - heavy data array
!!
!! output:
!!           -
!!
!!
!! = log ======================================================================================
!! \n
!! 07/11/16 - switch to v0.4 \n
!! 26/01/17 - switch to 3D, v0.5
!
! ********************************************************************************************

subroutine save_data(iteration, time, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n, hvy_work, hvy_active )

    !---------------------------------------------------------------------------------------------
    ! modules


    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> time loop parameters
    real(kind=rk), intent(in)                       :: time
    integer(kind=ik), intent(in)                    :: iteration

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> light data array
    integer(kind=ik), intent(inout)                 :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)                 :: lgt_active(:)
    !> number of active blocks (light/heavy data)
    integer(kind=ik), intent(inout)                 :: lgt_n, hvy_n

    ! loop variable
    integer(kind=ik)                                :: k, lgt_id
    ! file name
    character(len=80)                               :: fname, tmp
    ! cpu time variables for running time calculation
    real(kind=rk)                                   :: t0, x0(1:3), dx(1:3)
  
    !---------------------------------------------------------------------------------------------
    ! variables initialization
    ! start time
    t0 = MPI_Wtime()
    if (params%rank == 0) then
        write(*,'("IO: Saving data triggered, time=",g15.8)')  time
    endif

    !---------------------------------------------------------------------------------------------
    ! main body

    ! preparatory step. The physics modules have to copy everything they want to
    ! save to disk to the work array. missing qty's shall be computed.
    do k = 1, hvy_n
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        ! call preparatory routines. this routine saves the variables to be stored
        ! to disk in the work array. This way, we can also store derived variables
        ! such as the vorticity. Note in most cases, this copies just the state vector
        ! to work.
        call PREPARE_SAVE_DATA(params%physics_type, time, hvy_block(:,:,:,:,hvy_active(k)), &
        params%n_ghosts, x0, dx, hvy_work(:,:,:,:,hvy_active(k)))

    enddo

    ! actual saving step. one file per component.
    ! loop over components/qty's:
    do k = 1, params%N_fields_saved

        ! physics modules shall provide an interface for wabbit to know how to label
        ! the components to be stored to hard disk (in the work array)
        call FIELD_NAMES(params%physics_type, k, tmp)
        ! create filename
        write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(tmp)), nint(time * 1.0e6_rk)
        ! actual writing
        call write_field( fname, time, iteration, k, params, lgt_block, hvy_WORK, lgt_active, lgt_n, hvy_n, hvy_active)

    enddo

    call toc( params, "save_data", MPI_wtime()-t0 )
end subroutine save_data
