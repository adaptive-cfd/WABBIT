!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name update_neighbors.f90
!> \version 0.4
!> \author engels
!
!> \brief functional wrapper for the 2D and 3D version of update neighbors.
!
!>
!! input:
!!            - light data array
!!            - params struct
!!
!! output:
!!            - neighbor list array
!!
!! = log ======================================================================================
!! \n
!! 31/03/17 - create (tired of if-clause everywhere...)
!
! ********************************************************************************************

subroutine update_neighbors(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n)

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(out)       :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(in)     :: lgt_sortednumlist(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)     :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    logical :: error = .false., error2 = .false.
    integer(kind=ik ) :: mpierror

    if ( params%dim == 3 ) then
        ! 3D:
        call update_neighbors_3D(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, error)
    else
        ! 2D:
        call update_neighbors_2D(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, error)
    end if


    call MPI_Allreduce(error, error2, 1, MPI_LOGICAL, MPI_LOR, WABBIT_COMM, mpierror )

    if (error2) then
        write(*,*) "DUMPING DEBUG DATA TO *_ERROR.dat"
        call write_neighbors(params, hvy_active, hvy_n, hvy_neighbor, 'neighbors_ERROR.dat')
        call write_lgt_data(params, lgt_block, 'lgt_block_ERROR.dat')

        call abort(31375162, "Grid error. This is very nasty: some neighbor-updates failed. the specific error message above &
        & is probably not very useful. We dump the entire light data to *.dat in the hope this helps you find the problem.")
    endif

    if (error) then
        call abort(31375162, "Grid error. This is very nasty: some neighbor-updates failed. the specific error message above &
        & is probably not very useful. We dump the entire light data to *.dat in the hope this helps you find the problem.")
    endif

end subroutine update_neighbors
