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
! Please note that in the current implementation, hvy_tmp also plays the role of a work array
subroutine update_grid_qtys_ACM( time, field, g, x0, dx, stage )
    implicit none

    !> even though it is a bit odd, since those qtys shall be TIME INDEPENDENT, we pass time for debugging
    real(kind=rk), intent(in) :: time

    ! the grid-dependent qtys that are computed in this routine:
    real(kind=rk), intent(inout) :: field(1:,1:,1:,1:)

    ! set the qty only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)
    character(len=*), intent(in) :: stage

    integer(kind=2), allocatable, save :: mask_color(:,:,:)

    integer :: Bs
    ! compute the size of blocks
    Bs = size(field,1) - 2*g

    select case(stage)
    case("init_stage")
        ! even though this is strictly speaking not required, we update the insect here
        ! because if we don't we get an warning message that "time" and InsecT%time do not
        ! match
        if ( params_acm%geometry == "Insect" .and. Insect%body_moves == "no" .and. params_acm%dim==3 ) then
            call Update_Insect(time, Insect)
        endif

    case("main_stage")
        ! are we using insects? (3D only)
        if ( params_acm%geometry == "Insect" .and. Insect%body_moves == "no" .and. params_acm%dim==3 ) then
            ! in the tethered case, construct the body in the first 4 registers of
            ! grid_qty: mask, usx, usy, usz, color

            if (size(field,4) < 4) call abort(12121802,"[update_grid_qtys_ACM.f90]::not enough work arrays")

            if (.not. allocated(mask_color)) allocate(mask_color(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g))

            ! note the shift in origin: we pass the coordinates of point (1,1,1) since the insect module cannot
            ! know that the first g points are in fact ghost nodes...
            call Draw_Insect( time, Insect, x0, dx, field(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g,1), &
            mask_color, field(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g,2:4), with_body = .true., &
            with_wings = .false., delete_before_drawing = .true. )

            ! copy mask color array as well
            ! NOTE: I am not yet sure if I need this.
            ! NOTE: even us field is just zero, so maybe we shall not bother.
            field(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g,5) = dble( mask_color(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g) )
        endif

        ! are we using the sponge ?
        if ( params_acm%use_sponge ) then
            if ( params_acm%dim==3 ) then
                call sponge_3D( field(:, :, :, 6), x0, dx, Bs, g )
            else
                call sponge_2D( field(:, :, 1, 6), x0, dx, Bs, g )
            endif
        endif
    case default
        call abort(17121801,"[update_grid_qtys_ACM.f90]::unknown stage.")
    end select

end subroutine update_grid_qtys_ACM
