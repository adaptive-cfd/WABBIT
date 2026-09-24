! This file is included into module_mesh.



! Check if the query-point xq is in the extends of the block
! given by x0 and dx.
logical function pointInBlock_block(params, xq, x0, dx)
    implicit none
    type(type_params), intent(in) :: params
    ! query point
    real(kind=rk), intent(in) :: xq(3)
    ! origin and spacing of the block
    real(kind=rk), intent(in) :: x0(3)
    real(kind=rk), intent(in) :: dx(3)

    real(kind=rk) :: xmin, xmax
    integer(kind=ik) :: d

    pointInBlock_block = .true.

    ! loop over all directions
    do d = 1, params%dim

        ! extends of the blocks interior
        xmin = x0(d)
        ! NOTE: with the uniqueGrid definition, the block look likes this:
        !
        ! g g i i i i i i g g
        !     |           |
        !    x0          x0+Bs*dx
        !
        ! This is correct - the block can interpolate between the last interior
        ! point and the ghost nodes layer (hence the Bs*dx and not (Bs-1)*dx).
        xmax = x0(d) + dx(d) * real(params%Bs(d), rk)

        if (xq(d) < xmin - 1.0e-14_rk .or. xq(d) > xmax + 1.0e-14_rk) then
            pointInBlock_block = .false.
            return
        endif
    enddo
end function pointInBlock_block