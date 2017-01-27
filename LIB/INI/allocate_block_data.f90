! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: allocate_block_data.f90
! version: 0.5
! author: msr, engels
!
! allocate and reset heavy block data
!
!                   -> dim 1: x coord   ( 1:number_block_nodes+2*number_ghost_nodes )
!                   -> dim 2: y coord   ( 1:number_block_nodes+2*number_ghost_nodes )
!                   -> dim 3: z coord   ( 1:number_block_nodes+2*number_ghost_nodes )
!                   -> dim 4: data type ( field_1, 2:number_data_fields+1)
! heavy data array  -> dim 5: block id  ( 1:number_blocks )
!           field_1 (to save mixed data):   line 1: x coordinates
!                                           line 2: y coordinates
!
! input:    - maximal number of blocks per process
!           - grid parameter
!           - number of data fields
! output:   - empty heavy data array
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4
! 26/01/17 - use process rank from params struct
!          - 3D hvy data array, for 2D cases: use dim_size=1 for third dimension
!
! ********************************************************************************************

subroutine  allocate_block_data( params, hvy_block )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)              :: params

    ! heavy data array
    real(kind=rk), allocatable, intent(out)     :: hvy_block(:, :, :, :, :)

    ! number of heavy data
    integer(kind=ik)                            :: number_blocks
    ! grid parameter, blocksize (Bs), ghostnodes (g), number of fields (F)
    integer(kind=ik)                            :: Bs, g, dF

    ! allocation error variable
    integer(kind=ik)                            :: allocate_error

    ! process rank
    integer(kind=ik)                            :: rank

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank         = params%rank

    ! set parameters for readability
    number_blocks   = params%number_blocks
    Bs              = params%number_block_nodes
    g               = params%number_ghost_nodes
    dF              = params%number_data_fields

!---------------------------------------------------------------------------------------------
! main body

    ! allocate memory
    if ( params%threeD_case ) then
        ! 3D:
        ! datafields + 1 -> first field for coordinates, ...
        allocate( hvy_block( Bs+2*g, Bs+2*g, Bs+2*g, dF+1, number_blocks ), stat=allocate_error )
        call check_allocation(allocate_error)

    else
        ! 2D:
        ! datafields + 1 -> first field for coordinates, ...
        allocate( hvy_block( Bs+2*g, Bs+2*g, 1, dF+1, number_blocks ), stat=allocate_error )
        call check_allocation(allocate_error)

    end if

    ! reset data
    hvy_block = 0.0_rk

    ! output
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("INIT: System is allocating heavy data for ",i7," blocks and ", i3, " fields" )') number_blocks, dF
    end if

end subroutine allocate_block_data
