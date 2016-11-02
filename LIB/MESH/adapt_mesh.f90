! ********************************
! WABBIT
! --------------------------------
!
! mesh adapting main function
!
! name: adapt_mesh.f90
! date: 28.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine adapt_mesh()

    use mpi

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik)    :: i

    integer :: rank, ierr
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! maximal number of loops to coarsen the mesh == one block go down from max_treelevel to min_treelevel
    do i = 1, (params%max_treelevel - params%min_treelevel)

        ! check where to coarsen (refinement done with safety zone)
        call threshold_block()

        ! unmark blocks that cannot be coarsened due to gradedness
        call ensure_gradedness()

        ! ensure completeness
        call ensure_completeness()

        ! adapt the mesh
        call interpolate_mesh()

        ! update the neighbor relations
        call update_neighbors()

    end do

end subroutine adapt_mesh
