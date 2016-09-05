! ********************************
! 2D AMR prototype
! --------------------------------
!
! mesh adapting main function
!
! name: adapt_mesh.f90
! date: 02.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine adapt_mesh()

    use module_params
    use module_blocks

    implicit none

    ! check where to coarsen /refine
    call local_refinement_status()

    ! unmark blocks that cannot be coarsened due to gradedness
    call ensure_gradedness()

    ! ensure completeness
    call ensure_completeness()

    ! adapt the mesh
    call interpolate_mesh()
    call active_blocks_list()

    call block_check()

end subroutine adapt_mesh
