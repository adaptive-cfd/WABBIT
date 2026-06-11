!-------------------------------------------------------------------------------
! 2D geometry primitives (loop-free, pointwise)
!
! This file is included into `module_geometry` after the `contains` statement.
! Therefore it must only contain procedure definitions.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Exact euclidean signed distance to a (filled) 2D circle.
!!
!! Inputs:
!! - `p` query point
!! - `center` circle center
!! - `radius` circle radius
!!
!! Return value:
!! - signed dist with dist<0 inside the boundary
!!
!! Notes:
!! - Exact Euclidean distance is computed as `sqrt((p-center)^2) - radius`.
function signed_distance_circle_2D(p, center, radius) result(dist)
    use module_globals, only: rk
    implicit none
    real(kind=rk), intent(in) :: p(1:2), center(1:2), radius
    real(kind=rk) :: dist
    real(kind=rk) :: d_center

    d_center = norm2_2d(p - center)
    dist = d_center - radius
end function signed_distance_circle_2D


!-------------------------------------------------------------------------------
!> Exact euclidean signed distance to a (filled) 2D rectangle.
!!
!! Inputs:
!! - `p` query point
!! - `center` rectangle center
!! - `half_size` half size of the rectangle (extent in x and y direction)
!! - `angle` rotation angle of the rectangle in radians (counter-clockwise)
!!
!! Return value:
!! - signed dist with dist<0 inside the boundary
!!
!! Notes:
!! - Projects the points into the rectangle's local coordinate system, then computes the distance to the axis-aligned rectangle.
function signed_distance_rectangle_2D(p, center, half_size, angle, cos_angle, sin_angle) result(dist)
    use module_globals, only: rk
    implicit none
    real(kind=rk), intent(in) :: p(1:2), center(1:2), half_size(1:2), angle
    real(kind=rk), intent(in), optional :: cos_angle, sin_angle
    real(kind=rk) :: dist, cos_a, sin_a, local_p(1:2), d(1:2)

    if (abs(angle) <= 1.0e-12_rk .or. abs(angle - 2*pi) <= 1.0e-12_rk) then
        ! No rotation: simple case with less computation.
        local_p = p - center
        ! there are more easy cases (90, 180, 270 degree), but we ignore them for now because who would even use them?
    else
        ! cos and sin are the most expensive, they could be precomputed and just passed each time
        if (present(cos_angle) .and. present(sin_angle)) then
            cos_a = cos_angle
            sin_a = sin_angle
        else
            cos_a = cos(angle)
            sin_a = sin(angle)
        end if
        local_p(1) = cos_a * (p(1)-center(1)) + sin_a * (p(2)-center(2))
        local_p(2) = -sin_a * (p(1)-center(1)) + cos_a * (p(2)-center(2))
    end if

    if (abs(local_p(1)) > half_size(1) .and. abs(local_p(2)) > half_size(2)) then
        ! Corners: distance to the closest corner.
        dist = norm2_2d(abs(local_p) - half_size)
        return
    else if (abs(local_p(1)) > half_size(1)) then
        ! Left or right edge: distance to the vertical edge.
        dist = abs(local_p(1)) - half_size(1)
        return
    else if (abs(local_p(2)) > half_size(2)) then
        ! Top or bottom edge: distance to the horizontal edge.
        dist = abs(local_p(2)) - half_size(2)
        return
    else
        ! Inside the rectangle: distance to the closest edge (negative).
        dist = -min(half_size(1)-abs(local_p(1)), half_size(2)-abs(local_p(2)))
        return
    end if
end function signed_distance_rectangle_2D


!-------------------------------------------------------------------------------
!> Exact euclidean signed distance to a (filled) 2D triangle.
!!
!! Inputs:
!! - `p` query point
!! - `a`,`b`,`c` triangle vertices
!!
!! Return value:
!! - signed dist with dist<0 inside the boundary
!!
!! Notes:
!! - Exact Euclidean distance is computed as the minimum distance to the three
!!   line segments; sign is determined by an orientation-consistent half-plane
!!   test (works for CW and CCW vertex order).
!! - Degenerate triangles (area ~ 0) are not catched - assuming everything is alright
!-------------------------------------------------------------------------------
function signed_distance_triangle_2D(p, a, b, c) result(dist)
    use module_globals, only: rk
    implicit none
    real(kind=rk), intent(in) :: p(1:2), a(1:2), b(1:2), c(1:2)
    real(kind=rk) :: dist
    real(kind=rk) :: d2_ab, d2_bc, d2_ca, d2_min
    real(kind=rk) :: area2
    logical :: inside

    d2_ab = distance2_point_segment_2D(p, a, b)
    d2_bc = distance2_point_segment_2D(p, b, c)
    d2_ca = distance2_point_segment_2D(p, c, a)
    d2_min = min(d2_ab, min(d2_bc, d2_ca))
    dist = sqrt(d2_min)

    ! ! Twice the signed area (orientation); used to detect degeneracy.
    ! area2 = cross2_2D(b-a, c-a)

    ! if (abs(area2) <= 1.0e-14_rk) then
    !     ! Degenerate triangle: treat it as three segments; signedness is undefined.
    !     return
    ! end if

    inside = point_in_triangle_halfplane_2D(p, a, b, c)
    if (inside) dist = -dist
end function signed_distance_triangle_2D


!===============================================================================
! Private helpers
!===============================================================================

!-------------------------------------------------------------------------------
!> Squared Euclidean distance from point `p` to the 2D segment `[a,b]`.
!! Returns squared distance (no `sqrt`), which is cheaper in tight loops.
!-------------------------------------------------------------------------------
pure function distance2_point_segment_2D(p, a, b) result(d2)
    use module_globals, only: rk
    implicit none
    real(kind=rk), intent(in) :: p(1:2), a(1:2), b(1:2)
    real(kind=rk) :: d2
    real(kind=rk) :: ab(1:2), ap(1:2)
    real(kind=rk) :: t, denom
    real(kind=rk) :: closest(1:2), diff(1:2)

    ab = b - a
    ap = p - a
    denom = dot_product(ab, ab)

    if (denom <= 1.0e-12_rk) then
        ! Degenerate segment: treat as point distance.
        diff = p - a
        d2 = dot_product(diff, diff)
        return
    end if

    t = dot_product(ap, ab) / denom
    t = max(0.0_rk, min(1.0_rk, t))
    closest = a + t * ab
    diff = p - closest
    d2 = dot_product(diff, diff)
end function distance2_point_segment_2D


!-------------------------------------------------------------------------------
!> 2D cross product (scalar): `u_x*v_y - u_y*v_x`.
!! Used for orientation tests and signed area (`2*area`).
!-------------------------------------------------------------------------------
pure function cross2_2D(u, v) result(c)
    use module_globals, only: rk
    implicit none
    real(kind=rk), intent(in) :: u(1:2), v(1:2)
    real(kind=rk) :: c

    c = u(1) * v(2) - u(2) * v(1)
end function cross2_2D


!-------------------------------------------------------------------------------
!> Half-plane point-in-triangle test (includes edges).
!! Works for both CW and CCW vertex order.
!-------------------------------------------------------------------------------
pure function point_in_triangle_halfplane_2D(p, a, b, c) result(inside)
    use module_globals, only: rk
    implicit none
    real(kind=rk), intent(in) :: p(1:2), a(1:2), b(1:2), c(1:2)
    logical :: inside
    real(kind=rk) :: s1, s2, s3

    s1 = cross2_2D(b-a, p-a)
    s2 = cross2_2D(c-b, p-b)
    s3 = cross2_2D(a-c, p-c)

    inside = ((s1 >= 0.0_rk .and. s2 >= 0.0_rk .and. s3 >= 0.0_rk) .or. &
              (s1 <= 0.0_rk .and. s2 <= 0.0_rk .and. s3 <= 0.0_rk))
end function point_in_triangle_halfplane_2D

