! While the state vector and many work variables (such as the mask function for penalization)
! are explicitly time dependent, some other quantities are not. They are rather grid-dependent
! but need not to be updated in every RK or krylov substep. Hence, those quantities are updated
! after the mesh is changed (i.e. after refine_mesh) and then kept constant during the evolution
! time step.
! An example for such a quantity would be geometry factors on non-cartesian grids, but also the
! body of an insect in tethered (=fixed) flight. In the latter example, only the wings need to be
! generated at every time t. This example generalizes to any combination of stationary and moving
! obstacle, i.e. insect behind fractal tree.
! Updating those grid-depend quantities is a task for the physics modules: they should provide interfaces, if they require such qantities. In many cases, the grid_qtys are probably not used.
! Please noe that in the current implementation, hvy_tmp also plays the role of a work array

subroutine update_grid_qtys_ACM(u, g, x0, dx )
    implicit none

    ! the grid-dependent qtys that are computed in this routine:
    real(kind=rk), intent(in) :: u(1:,1:,1:,1:)

    ! set the qty only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)




end subroutine update_grid_qtys_ACM
