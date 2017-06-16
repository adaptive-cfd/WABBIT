!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name allocate_com_arrays.f90
!> \version 0.5
!> \author msr
!
!> \brief Allocate communications arrays, to store com_matrix and com_lists between
!! synchronization steps with fixed grid
!
!>
!! input:
!!           - parameter array
!!           - com arrays
!!
!! output:
!!           - allocated com arrays
!!
!! = log ======================================================================================
!! \n
!! 13/06/17 - create
!
! ********************************************************************************************
subroutine allocate_com_arrays(params, com_lists, com_matrix)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)                   :: params

    !> communication lists:
    integer(kind=ik), allocatable, intent(out)          :: com_lists(:, :, :, :)

    !> communications matrix:
    integer(kind=ik), allocatable, intent(out)          :: com_matrix(:,:,:)

    ! local shortcuts:
    integer(kind=ik)                                    :: max_neighbors, number_procs, N

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    number_procs    = params%number_procs
    N               = params%number_blocks

    if ( params%threeD_case ) then
        ! 3D
        max_neighbors = 74
    else
        ! 2D
        max_neighbors = 12
    end if

!---------------------------------------------------------------------------------------------
! main body

    ! deallocate old arrays
    if (allocated(com_lists)) deallocate(com_lists)
    if (allocated(com_matrix)) deallocate(com_matrix)

    ! com lists
    allocate( com_lists( N*max_neighbors, 6, number_procs, 4) )
    ! allocate com matrix
    allocate( com_matrix(number_procs, number_procs, 4) )

end subroutine allocate_com_arrays
