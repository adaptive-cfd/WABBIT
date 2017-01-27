! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: allocate_work_data.f90
! version: 0.5
! author: msr
!
! allocate and reset heavy work data
!
!                   -> dim 1: x coord   ( 1:number_block_nodes+2*number_ghost_nodes )
!                   -> dim 2: y coord   ( 1:number_block_nodes+2*number_ghost_nodes )
!                   -> dim 3: z coord   ( 1:number_block_nodes+2*number_ghost_nodes )
!                   -> dim 4: data type ( old data, k1, k2, k3, k4 )
! heavy work array  -> dim 5: block id  ( 1:number_blocks )
!
! input:    - params
! output:   - empty heavy data array
!
! = log ======================================================================================
!
! 07/12/16 - create
! 26/01/17 - use process rank from params struct
!          - 3D hvy data array, for 2D cases: use dim_size=1 for third dimension
!
! ********************************************************************************************

subroutine  allocate_work_data( params, hvy_work )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)              :: params

    ! heavy work data array
    real(kind=rk), allocatable, intent(out)     :: hvy_work(:, :, :, :, :)

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
        allocate( hvy_work( Bs+2*g, Bs+2*g, Bs+2*g, dF*5, number_blocks ), stat=allocate_error )
        call check_allocation(allocate_error)

    else
        ! 2D:
        allocate( hvy_work( Bs+2*g, Bs+2*g, 1, dF*5, number_blocks ), stat=allocate_error )
        call check_allocation(allocate_error)

    end if

    ! reset data
    hvy_work = 0.0_rk

    ! output
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("INIT: System is allocating heavy work data for ",i7," blocks " )') number_blocks
    end if

end subroutine allocate_work_data
