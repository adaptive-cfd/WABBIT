!-------------------------------------------------------------------------------
! 3D geometry primitives (loop-free, pointwise)
!
! This file is included into `module_geometry` after the `contains` statement.
! Therefore it must only contain procedure definitions.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Exact euclidean signed distance to a (filled) 3D sphere.
!!
!! Inputs:
!! - `p` query point
!! - `center` sphere center
!! - `radius` sphere radius
!!
!! Return value:
!! - signed dist with dist<0 inside the boundary
!!
!! Notes:
!! - Exact Euclidean distance is computed as `sqrt((p-center)^2) - radius`.
function signed_distance_sphere_3D(p, center, radius) result(dist)
    use module_globals, only: rk
    implicit none
    real(kind=rk), intent(in) :: p(1:3), center(1:3), radius
    real(kind=rk) :: dist
    real(kind=rk) :: d_center

    d_center = norm2_3d(p - center)
    dist = d_center - radius
end function signed_distance_sphere_3D

!-----------------------------------------------------------------------------
!> Exact signed distance to a finite cylinder segment (flat end caps).
!!
!! Geometry:
!! - Cylinder axis segment from `x1` to `x2`
!! - Radius `radius`
!!
!! Return value:
!! - signed dist with dist<0 inside the boundary
!!
!! Implementation:
!! - Analytic SDF for a capped cylinder segment (https://iquilezles.org/articles/distfunctions/ , https://www.shadertoy.com/view/wdXGDr)
!! - Uses a form that postpones dividing by |x2-x1|^2 (`axis_len2`) until the very end, which avoids repeated normalization.
!!
!! Optional precomputed inputs for performance:
!! - `ba_in`   : axis vector (x2 - x1)
!! - `baba_in` : squared axis length |x2-x1|^2
!! Passing these avoids re-computing `x2-x1` and the dot product in tight loops.
!!
!! Assumptions (for performance):
!! - `x1` and `x2` are distinct (so `axis_len2 > 0`). No degenerate check.
!-----------------------------------------------------------------------------
function signed_distance_cylinder_3D(point, x1, x2, radius, ba_in, baba_in) result(dist)
    implicit none
    real(kind=rk), dimension(3), intent(in) :: point, x1, x2
    real(kind=rk), intent(in) :: radius
    real(kind=rk), dimension(3), intent(in), optional :: ba_in
    real(kind=rk), intent(in), optional :: baba_in
    real(kind=rk) :: dist
    real(kind=rk), dimension(3) :: axis, point_rel, tmp_vec
    real(kind=rk) :: axis_len2, proj, radial_term, axial_term
    real(kind=rk) :: radial_term2, axial_term2

    if (present(ba_in)) then
        axis = ba_in
    else
        axis = x2 - x1
    end if
    if (present(baba_in)) then
        axis_len2 = baba_in
    else
        axis_len2 = dot_product(axis, axis)
    end if
    point_rel = point - x1
    proj = dot_product(point_rel, axis)
    
    tmp_vec = point_rel * axis_len2 - axis * proj
    radial_term = norm2_3d(tmp_vec) - radius * axis_len2
    axial_term = abs(proj - axis_len2 * 0.5_rk) - axis_len2 * 0.5_rk

    radial_term2 = radial_term * radial_term
    axial_term2 = axial_term * axial_term * axis_len2

    if (max(radial_term, axial_term) < 0.0_rk) then
        ! dist = -sqrt(min(radial_term2, axial_term2))
        dist = -min(radial_term2, axial_term2)
    else
        ! dist = sqrt((merge(radial_term2, 0.0_rk, radial_term > 0.0_rk) + merge(axial_term2, 0.0_rk, axial_term > 0.0_rk)))
        dist = 0.0_rk
        if (radial_term > 0.0_rk) then
            dist = dist + radial_term2
        end if
        if (axial_term > 0.0_rk) then
            dist = dist + axial_term2
        end if
    end if

    dist = sign(1.0_rk, dist) * sqrt(abs(dist)) / axis_len2
    ! dist = sign(1.0_rk, dist) * dist / axis_len2
end function signed_distance_cylinder_3D

!-------------------------------------------------------------------------------
!> Exact euclidean signed distance to a (filled) cylinder with spherical end caps
!!
!! Geometry:
!! - Cylinder axis segment from `x1` to `x2`
!! - Radius `radius`
!!
!! Return value:
!! - signed dist with dist<0 inside the boundary
!!
!! Implementation:
!! - Find closest point on line axis segment, then compute distance to that point minus radius.
!!
!! Optional precomputed inputs for performance:
!! - `ba_in`   : axis vector (x2 - x1)
!! - `baba_in` : squared axis length |x2-x1|^2
!! Passing these avoids re-computing `x2-x1` and the dot product in tight loops.
!!
!! Assumptions (for performance):
!! - `x1` and `x2` are distinct (so `axis_len2 > 0`). No degenerate check.
!-------------------------------------------------------------------------------
function signed_distance_cylinder_rounded_3D(point, x1, x2, radius, ba_in, baba_in) result(dist)
    implicit none
    real(kind=rk), dimension(3), intent(in) :: point, x1, x2
    real(kind=rk), intent(in) :: radius
    real(kind=rk), dimension(3), intent(in), optional :: ba_in
    real(kind=rk), intent(in), optional :: baba_in
    real(kind=rk) :: dist
    real(kind=rk), dimension(3) :: axis, point_rel
    real(kind=rk) :: axis_len2, proj

    if (present(ba_in)) then
        axis = ba_in
    else
        axis = x2 - x1
    end if
    if (present(baba_in)) then
        axis_len2 = baba_in
    else
        axis_len2 = dot_product(axis, axis)
    end if
    point_rel = point - x1
    proj = dot_product(point_rel, axis) / axis_len2
    proj = max(0.0_rk, min(1.0_rk, proj))  ! clamp to [0,1] to stay within the segment
    dist = norm2_3d(point_rel - proj * axis) - radius
end function signed_distance_cylinder_rounded_3D