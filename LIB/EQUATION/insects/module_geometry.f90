module module_geometry
    use module_globals
    use module_helpers
    use module_ini_files_parser_mpi

    implicit none

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    ! functions
    PUBLIC :: signed_distance_triangle_2D, signed_distance_circle_2D, signed_distance_rectangle_2D, signed_distance_sphere_3D, signed_distance_cylinder_3D, signed_distance_cylinder_rounded_3D, draw_sphere, draw_cylinder, draw_cylinder_rounded, draw_circle, draw_rectangle, draw_triangle, init_primitives_collection, draw_primitives_collection

    ! for 2d wing section optimization
    type :: primitives_collection
        logical :: initialized = .false.
        real(kind=rk), allocatable :: geometric_data(:,:), geometric_boundingbox(:,:)
        character(len=clong), allocatable :: geometric_type(:)
        integer(kind=ik) :: color_set, n_objects
    end type

    type(primitives_collection) :: geom_array(20)  ! we can have up to 20 primitives-collections, which should be enough for most cases. We can make this allocatable if we want to be fancy, but I don't see the need for now.

contains

!---------------------------------------
! note these include files also have to be specified as dependencies in the
! Makefile for make to check if one of them changed
#include "geometry_primitive_2D.f90"
#include "geometry_primitive_3D.f90"
#include "geometry_primitives_collection.f90"
!---------------------------------------

!---------------------------------------
!> Draw a sphere using the exact SDF
!!
!! This is a pointwise implementation that loops over all grid points in the bounding box of the sphere and applies the mask function based on the signed distance to the sphere.
!!
!! Inputs:
!! - `x0`, `dx`, `g` origin and spacing, ghost point size of the block
!! - `center`, `radius` geometry parameters
!! - `color_set` which color to set in the mask
!! - `smoothing_type_int`, `smoothing_width`, `smoothing_safety` parameters for the step function
!! - `bounding_box` optional pre-computed bounding box for the geometry
!!
!! Outputs:
!! - `mask`, `color` updated with the geometry
subroutine draw_sphere(mask, color, x0, dx, g, center, radius, color_set, smoothing_type_int, smoothing_width, smoothing_safety, bounding_box, x0_indices)
    use module_globals

    implicit none

    !> mask and color term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(out)     :: mask, color
    !> spacing and origin of block, being located at position g+1 in the mask array
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)
    !> grid information
    integer(kind=ik), intent(in) :: g
    real(kind=rk), dimension(3), intent(in) :: center  !< sphere center
    real(kind=rk), intent(in) :: radius                !< sphere radius
    integer(kind=ik), intent(in) :: color_set               !< which color to set the mask
    integer(kind=ik), intent(in) :: smoothing_type_int     !< which mask do we use?
    real(kind=rk), intent(in) :: smoothing_width       !< width of smoothing region
    real(kind=rk), optional, intent(in) :: smoothing_safety      !< safety margin for smoothing
    real(kind=rk),optional,intent(in) :: bounding_box(1:6)
    integer, optional, intent(in) :: x0_indices(1:3)  !< For WABBIT, x0 is defined at the first interior point g+1, for FLUSI, it is actually at index 0, so we can provide here the indices of x0, defaults to WABBIT formulation

    ! auxiliary variables
    real(kind=rk)  :: xyz(1:3), dist, safety, tmp
    integer, dimension(1:3) :: lbounds, ubounds
    integer(kind=ik) :: iz, ix, iy, Nsafety, bound_min(1:3), bound_max(1:3), x0_offset(1:3)

    ! set default value for x0_indices
    x0_offset = g+1
    if (present(x0_indices)) x0_offset = x0_indices

    ! Mask array is not resetted, as we may have other objects in the vicinity
    ! ! reset mask array
    ! mask = 0.0_rk

    ! safety for smoothing function
    safety = 3.0_rk * smoothing_width
    if (present(smoothing_safety)) safety = smoothing_safety
    Nsafety = ceiling(safety / minval(dx))

    ! bounds of the current patch of data
    lbounds = g+1
    ubounds = (/size(mask,1), size(mask,2), size(mask,3)/) -g+1

    ! bounding box check
    if (present(bounding_box)) then
        ! use pre-computed bounding box for sphere (can be faster if many spheres are to be drawn)
        bound_min = int((bounding_box(1:3)-x0)/dx) - Nsafety + g+1
        bound_max = int(ceiling((bounding_box(4:6)-x0)/dx)) + Nsafety + g+1
    else
        bound_min = int((center - radius - x0)/dx) - Nsafety + g+1
        bound_max = int(ceiling((center + radius - x0)/dx)) + Nsafety + g+1
    endif

    do iz = max(bound_min(3),lbounds(3)), min(bound_max(3),ubounds(3))
        xyz(3) = dble(iz-x0_offset(3)) * dx(3) + x0(3)
        do iy = max(bound_min(2),lbounds(2)), min(bound_max(2),ubounds(2))
            xyz(2) = dble(iy-x0_offset(2)) * dx(2) + x0(2)
            do ix = max(bound_min(1),lbounds(1)), min(bound_max(1),ubounds(1))
                xyz(1) = dble(ix-x0_offset(1)) * dx(1) + x0(1)
                ! distance from center of sphere
                dist = signed_distance_sphere_3D( xyz, center, radius )
                ! apply mask function
                tmp = step(dist, 0.0_rk, smoothing_width, safety, smoothing_type_int)
                if (tmp >= mask(ix,iy,iz)) then
                    mask(ix,iy,iz) = tmp    ! mask function
                    color(ix,iy,iz) = real(color_set, kind=rk) ! color
                endif
    end do; end do; end do

end subroutine draw_sphere

!---------------------------------------
!> Draw a cylinder with flat end points using the exact SDF
!!
!! This is a pointwise implementation that loops over all grid points in the bounding box of the cylinder and applies the mask function based on the signed distance to the cylinder.
!!
!! Inputs:
!! - `x0`, `dx`, `g` origin and spacing, ghost point size of the block
!! - `center`, `radius` geometry parameters
!! - `color_set` which color to set in the mask
!! - `smoothing_type_int`, `smoothing_width`, `smoothing_safety` parameters for the step function
!! - `bounding_box` optional pre-computed bounding box for the geometry
!!
!! Outputs:
!! - `mask`, `color` updated with the geometry
subroutine draw_cylinder(mask, color, x0, dx, g, endpoint_1, endpoint_2, radius, color_set, smoothing_type_int, smoothing_width, smoothing_safety, bounding_box, x0_indices)
    use module_globals

    implicit none

    !> mask and color term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(out)     :: mask, color
    !> spacing and origin of block, being located at position g+1 in the mask array
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)
    !> grid information
    integer(kind=ik), intent(in) :: g
    real(kind=rk), dimension(3), intent(in) :: endpoint_1, endpoint_2  !< cylinder endpoints
    real(kind=rk), intent(in) :: radius                !< cylinder radius
    integer(kind=ik), intent(in) :: color_set               !< which color to set the mask
    integer(kind=ik), intent(in) :: smoothing_type_int     !< which mask do we use?
    real(kind=rk), intent(in) :: smoothing_width       !< width of smoothing region
    real(kind=rk), optional, intent(in) :: smoothing_safety      !< safety margin for smoothing
    real(kind=rk),optional,intent(in) :: bounding_box(1:6)
    integer, optional, intent(in) :: x0_indices(1:3)  !< For WABBIT, x0 is defined at the first interior point g+1, for FLUSI, it is actually at index 0, so we can provide here the indices of x0, defaults to WABBIT formulation

    ! auxiliary variables
    real(kind=rk)  :: xyz(1:3), dist, safety, tmp, ba(1:3), baba, er(1:3), et(1:3)
    integer, dimension(1:3) :: lbounds, ubounds
    integer(kind=ik) :: iz, ix, iy, Nsafety, bound_min(1:3), bound_max(1:3), x0_offset(1:3)

    ! set default value for x0_indices
    x0_offset = g+1
    if (present(x0_indices)) x0_offset = x0_indices

    ! Mask array is not resetted, as we may have other objects in the vicinity
    ! ! reset mask array
    ! mask = 0.0_rk

    ! safety for smoothing function
    safety = 3.0_rk * smoothing_width
    if (present(smoothing_safety)) safety = smoothing_safety
    Nsafety = ceiling(safety / minval(dx))

    ! bounds of the current patch of data
    lbounds = g+1
    ubounds = (/size(mask,1), size(mask,2), size(mask,3)/) -g+1

    ! soem auxiliary variables for the cylinder, speeds up point-wise computations
    ba = endpoint_2 - endpoint_1
    baba = dot_product(ba, ba)
    er = ba / sqrt(baba)
    et = sqrt(1.0_rk - er**2)

    ! bounding box check for cyulinder
    ! for cylinder, we can go at most e_t * radius away from the endpoints. So we compute e_r and then use the relation e_t = sqrt(1-e_r^2) and do this for every component
    ! for angle=0 (parallel to x-axis) for example, e_r = (1,0,0) and e_t = (0,1,1)
    if (present(bounding_box)) then
        ! use pre-computed bounding box for sphere (can be faster if many spheres are to be drawn)
        bound_min = int((bounding_box(1:3)-x0)/dx) - Nsafety + g+1
        bound_max = int(ceiling((bounding_box(4:6)-x0)/dx)) + Nsafety + g+1
    else
        bound_min = int((min(endpoint_1, endpoint_2)-radius*et - x0)/dx) - Nsafety + g+1
        bound_max = int(ceiling((max(endpoint_1, endpoint_2)+radius*et - x0)/dx)) + Nsafety + g+1
    endif

    do iz = max(bound_min(3),lbounds(3)), min(bound_max(3),ubounds(3))
        xyz(3) = dble(iz-x0_offset(3)) * dx(3) + x0(3)
        do iy = max(bound_min(2),lbounds(2)), min(bound_max(2),ubounds(2))
            xyz(2) = dble(iy-x0_offset(2)) * dx(2) + x0(2)
            do ix = max(bound_min(1),lbounds(1)), min(bound_max(1),ubounds(1))
                xyz(1) = dble(ix-x0_offset(1)) * dx(1) + x0(1)
                ! distance from center of cylinder
                dist = signed_distance_cylinder_3D( xyz, endpoint_1, endpoint_2, radius, ba_in=ba, baba_in=baba )
                ! apply mask function
                tmp = step(dist, 0.0_rk, smoothing_width, safety, smoothing_type_int)
                if (tmp >= mask(ix,iy,iz)) then
                    mask(ix,iy,iz) = tmp    ! mask function
                    color(ix,iy,iz) = real(color_set, kind=rk) ! color
                endif
    end do; end do; end do

end subroutine draw_cylinder

!---------------------------------------
!> Draw a cylinder-rounded / cylinder with spherical end points using the exact SDF
!!
!! This is a pointwise implementation that loops over all grid points in the bounding box of the cylinder and applies the mask function based on the signed distance to the cylinder.
!!
!! Inputs:
!! - `x0`, `dx`, `g` origin and spacing, ghost point size of the block
!! - `center`, `radius` geometry parameters
!! - `color_set` which color to set in the mask
!! - `smoothing_type_int`, `smoothing_width`, `smoothing_safety` parameters for the step function
!! - `bounding_box` optional pre-computed bounding box for the geometry
!!
!! Outputs:
!! - `mask`, `color` updated with the geometry
subroutine draw_cylinder_rounded(mask, color, x0, dx, g, endpoint_1, endpoint_2, radius, color_set, smoothing_type_int, smoothing_width, smoothing_safety, bounding_box, x0_indices )
    use module_globals

    implicit none

    !> mask and color term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(out)     :: mask, color
    !> spacing and origin of block, being located at position g+1 in the mask array
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)
    !> grid information
    integer(kind=ik), intent(in) :: g
    real(kind=rk), dimension(3), intent(in) :: endpoint_1, endpoint_2  !< cylinder-rounded endpoints
    real(kind=rk), intent(in) :: radius                !< cylinder-rounded radius
    integer(kind=ik), intent(in) :: color_set               !< which color to set the mask
    integer(kind=ik), intent(in) :: smoothing_type_int     !< which mask do we use?
    real(kind=rk), intent(in) :: smoothing_width       !< width of smoothing region
    real(kind=rk), optional, intent(in) :: smoothing_safety      !< safety margin for smoothing
    real(kind=rk),optional,intent(in) :: bounding_box(1:6)
    integer, optional, intent(in) :: x0_indices(1:3)  !< For WABBIT, x0 is defined at the first interior point g+1, for FLUSI, it is actually at index 0, so we can provide here the indices of x0, defaults to WABBIT formulation

    ! auxiliary variables
    real(kind=rk)  :: xyz(1:3), dist, safety, tmp, ba(1:3), baba
    integer, dimension(1:3) :: lbounds, ubounds
    integer(kind=ik) :: iz, ix, iy, Nsafety, bound_min(1:3), bound_max(1:3), x0_offset(1:3)

    ! set default value for x0_indices
    x0_offset = g+1
    if (present(x0_indices)) x0_offset = x0_indices

    ! Mask array is not resetted, as we may have other objects in the vicinity
    ! ! reset mask array
    ! mask = 0.0_rk

    ! safety for smoothing function
    safety = 3.0_rk * smoothing_width
    if (present(smoothing_safety)) safety = smoothing_safety
    Nsafety = ceiling(safety / minval(dx))

    ! bounds of the current patch of data
    lbounds = g+1
    ubounds = (/size(mask,1), size(mask,2), size(mask,3)/) -g+1

    ! soem auxiliary variables for the cylinder-rounded, speeds up point-wise computations
    ba = endpoint_2 - endpoint_1
    baba = dot_product(ba, ba)

    ! bounding box check for cylinder-rounded
    ! for cylinder-rounded, we can go at most radius away from the endpoints.
    if (present(bounding_box)) then
        ! use pre-computed bounding box for sphere (can be faster if many spheres are to be drawn)
        bound_min = int((bounding_box(1:3)-x0)/dx) - Nsafety + g+1
        bound_max = int(ceiling((bounding_box(4:6)-x0)/dx)) + Nsafety + g+1
    else
        bound_min = int((min(endpoint_1, endpoint_2)-radius - x0)/dx) - Nsafety + g+1
        bound_max = int(ceiling((max(endpoint_1, endpoint_2)+radius - x0)/dx)) + Nsafety + g+1
    endif

    do iz = max(bound_min(3),lbounds(3)), min(bound_max(3),ubounds(3))
        xyz(3) = dble(iz-x0_offset(3)) * dx(3) + x0(3)
        do iy = max(bound_min(2),lbounds(2)), min(bound_max(2),ubounds(2))
            xyz(2) = dble(iy-x0_offset(2)) * dx(2) + x0(2)
            do ix = max(bound_min(1),lbounds(1)), min(bound_max(1),ubounds(1))
                xyz(1) = dble(ix-x0_offset(1)) * dx(1) + x0(1)
                ! distance from center of cylinder-rounded
                dist = signed_distance_cylinder_rounded_3D( xyz, endpoint_1, endpoint_2, radius, ba_in=ba, baba_in=baba )
                ! apply mask function
                tmp = step(dist, 0.0_rk, smoothing_width, safety, smoothing_type_int)
                if (tmp >= mask(ix,iy,iz)) then
                    mask(ix,iy,iz) = tmp    ! mask function
                    color(ix,iy,iz) = real(color_set, kind=rk) ! color
                endif
    end do; end do; end do

end subroutine draw_cylinder_rounded


!---------------------------------------
!> Draw a circle using the exact SDF
!!
!! This is a pointwise implementation that loops over all grid points in the bounding box of the circle and applies the mask function based on the signed distance to the circle.
!!
!! Inputs:
!! - `x0`, `dx`, `g` origin and spacing, ghost point size of the block
!! - `center`, `radius` geometry parameters
!! - `color_set` which color to set in the mask
!! - `smoothing_type_int`, `smoothing_width`, `smoothing_safety` parameters for the step function
!! - `bounding_box` optional pre-computed bounding box for the geometry
!!
!! Outputs:
!! - `mask`, `color` updated with the geometry
subroutine draw_circle(mask, color, x0, dx, g, center, radius, color_set, smoothing_type_int, smoothing_width, smoothing_safety, bounding_box, x0_indices )
    use module_globals

    implicit none

    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)     :: mask, color
    !> spacing and origin of block, being located at position g+1 in the mask array
    real(kind=rk), intent(in) :: x0(1:2), dx(1:2)
    !> grid information
    integer(kind=ik), intent(in) :: g
    real(kind=rk), dimension(2), intent(in) :: center  !< circle center
    real(kind=rk), intent(in) :: radius                !< circle radius
    integer(kind=ik), intent(in) :: color_set               !< which color to set the mask
    integer(kind=ik), intent(in) :: smoothing_type_int     !< which mask do we use?
    real(kind=rk), intent(in) :: smoothing_width       !< width of smoothing region
    real(kind=rk), optional, intent(in) :: smoothing_safety      !< safety margin for smoothing
    real(kind=rk), optional,intent(in) :: bounding_box(1:4)
    integer, optional, intent(in) :: x0_indices(1:2)  !< For WABBIT, x0 is defined at the first interior point g+1, for FLUSI, it is actually at index 0, so we can provide here the indices of x0, defaults to WABBIT formulation

    ! auxiliary variables
    real(kind=rk)  :: xyz(1:2), dist, safety, tmp
    integer, dimension(1:2) :: lbounds, ubounds
    integer(kind=ik) :: ix, iy, Nsafety, bound_min(1:2), bound_max(1:2), x0_offset(1:2)

    ! set default value for x0_indices
    x0_offset = g+1
    if (present(x0_indices)) x0_offset = x0_indices

    ! Mask array is not resetted, as we may have other objects in the vicinity
    ! ! reset mask array
    ! mask = 0.0_rk

    ! safety for smoothing function
    safety = 3.0_rk * smoothing_width
    if (present(smoothing_safety)) safety = smoothing_safety
    Nsafety = ceiling(safety / minval(dx))

    ! bounds of the current patch of data
    lbounds = g+1
    ubounds = (/size(mask,1), size(mask,2)/) -g+1

    ! bounding box check
    if (present(bounding_box)) then
        ! use pre-computed bounding box for sphere (can be faster if many spheres are to be drawn)
        bound_min = int((bounding_box(1:2)-x0)/dx) - Nsafety + g+1
        bound_max = int(ceiling((bounding_box(3:4)-x0)/dx)) + Nsafety + g+1
    else
        bound_min = int((center - radius - x0)/dx) - Nsafety + g+1
        bound_max = int(ceiling((center + radius - x0)/dx)) + Nsafety + g+1
    endif

    do iy = max(bound_min(2),lbounds(2)), min(bound_max(2),ubounds(2))
        xyz(2) = dble(iy-(x0_offset(2))) * dx(2) + x0(2)
        do ix = max(bound_min(1),lbounds(1)), min(bound_max(1),ubounds(1))
            xyz(1) = dble(ix-(x0_offset(1))) * dx(1) + x0(1)
            ! distance from center of sphere
            dist = signed_distance_circle_2D( xyz, center, radius )
            ! apply mask function
            tmp = step(dist, 0.0_rk, smoothing_width, safety, smoothing_type_int)
            if (tmp >= mask(ix,iy)) then
                mask(ix,iy) = tmp    ! mask function
                color(ix,iy) = real(color_set, kind=rk) ! color
            endif
    end do; end do

end subroutine draw_circle


!---------------------------------------
!> Draw a rectangle using the exact SDF
!!
!! This is a pointwise implementation that loops over all grid points in the bounding box of the rectangle and applies the mask function based on the signed distance to the rectangle.
!!
!! Inputs:
!! - `x0`, `dx`, `Bs`, `g` origin and spacing, Block size, ghost point size of the block
!! - `center`, `half_size`, `angle` geometry parameters
!! - `color_set` which color to set in the mask
!! - `smoothing_type_int`, `smoothing_width`, `smoothing_safety` parameters for the step function
!! - `bounding_box` optional pre-computed bounding box for the geometry
!!
!! Outputs:
!! - `mask`, `color` updated with the geometry
subroutine draw_rectangle(mask, color, x0, dx, g, center, half_size, angle, color_set, smoothing_type_int, smoothing_width, smoothing_safety, bounding_box, x0_indices )
    use module_globals

    implicit none

    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)     :: mask, color
    !> spacing and origin of block, being located at position g+1 in the mask array
    real(kind=rk), intent(in) :: x0(1:2), dx(1:2)
    !> grid information
    integer(kind=ik), intent(in) :: g
    real(kind=rk), dimension(2), intent(in) :: center  !< rectangle center
    real(kind=rk), intent(in) :: half_size(1:2)                !< rectangle half size
    real(kind=rk), intent(in) :: angle                   !< rectangle angle in radians
    integer(kind=ik), intent(in) :: color_set               !< which color to set the mask
    integer(kind=ik), intent(in) :: smoothing_type_int     !< which mask do we use?
    real(kind=rk), intent(in) :: smoothing_width       !< width of smoothing region
    real(kind=rk), optional, intent(in) :: smoothing_safety      !< safety margin for smoothing
    real(kind=rk),optional,intent(in) :: bounding_box(1:4)
    integer, optional, intent(in) :: x0_indices(1:2)  !< For WABBIT, x0 is defined at the first interior point g+1, for FLUSI, it is actually at index 0, so we can provide here the indices of x0, defaults to WABBIT formulation

    ! auxiliary variables
    real(kind=rk)  :: xyz(1:2), dist, safety, tmp, sin_a, cos_a
    integer, dimension(1:2) :: lbounds, ubounds
    integer(kind=ik) :: ix, iy, Nsafety, bound_min(1:2), bound_max(1:2), x0_offset(1:2)

    ! set default value for x0_indices
    x0_offset = g+1
    if (present(x0_indices)) x0_offset = x0_indices

    ! Mask array is not resetted, as we may have other objects in the vicinity
    ! ! reset mask array
    ! mask = 0.0_rk

    ! safety for smoothing function
    safety = 3.0_rk * smoothing_width
    if (present(smoothing_safety)) safety = smoothing_safety
    Nsafety = ceiling(safety / minval(dx))

    ! bounds of the current patch of data
    lbounds = g+1
    ubounds = (/size(mask,1), size(mask,2)/) -g+1

    ! precompute sin and cos of angle for rotation, speeds up point-wise computations
    sin_a = sin(angle)
    cos_a = cos(angle)

    if (present(bounding_box)) then
        ! use pre-computed bounding box for rectangle (can be faster if many rectangles are to be drawn)
        bound_min = int((bounding_box(1:2)-x0)/dx) - Nsafety + g+1
        bound_max = int(ceiling((bounding_box(3:4)-x0)/dx)) + Nsafety + g+1
    else
        bound_min = int((center - abs(half_size*cos_a) - abs(half_size(2:1:-1) * sin_a) - x0)/dx) - Nsafety + g+1
        bound_max = int(ceiling((center + abs(half_size*cos_a) + abs(half_size(2:1:-1) * sin_a) - x0)/dx)) + Nsafety + g+1
    endif

    do iy = max(bound_min(2),lbounds(2)), min(bound_max(2),ubounds(2))
        xyz(2) = dble(iy-(x0_offset(2))) * dx(2) + x0(2)
        do ix = max(bound_min(1),lbounds(1)), min(bound_max(1),ubounds(1))
            xyz(1) = dble(ix-(x0_offset(1))) * dx(1) + x0(1)
            ! distance from center of rectangle
            dist = signed_distance_rectangle_2D( xyz, center, half_size, angle, sin_angle=sin_a, cos_angle=cos_a )
            ! apply mask function
            tmp = step(dist, 0.0_rk, smoothing_width, safety, smoothing_type_int)
            if (tmp >= mask(ix,iy)) then
                mask(ix,iy) = tmp    ! mask function
                color(ix,iy) = real(color_set, kind=rk) ! color
            endif
    end do; end do

end subroutine draw_rectangle


!---------------------------------------
!> Draw a rectangle using the exact SDF
!!
!! This is a pointwise implementation that loops over all grid points in the bounding box of the rectangle and applies the mask function based on the signed distance to the rectangle.
!!
!! Inputs:
!! - `x0`, `dx`, `g` origin and spacing, ghost point size of the block
!! - `center`, `half_size`, `angle` geometry parameters
!! - `color_set` which color to set in the mask
!! - `smoothing_type_int`, `smoothing_width`, `smoothing_safety` parameters for the step function
!! - `bounding_box` optional pre-computed bounding box for the geometry
!!
!! Outputs:
!! - `mask`, `color` updated with the geometry
subroutine draw_triangle(mask, color, x0, dx, g, vertex1, vertex2, vertex3, color_set, smoothing_type_int, smoothing_width, smoothing_safety, bounding_box, x0_indices )
    use module_globals

    implicit none

    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)     :: mask, color
    !> spacing and origin of block, being located at position g+1 in the mask array
    real(kind=rk), intent(in) :: x0(1:2), dx(1:2)
    !> grid information
    integer(kind=ik), intent(in) :: g
    real(kind=rk), dimension(2), intent(in) :: vertex1, vertex2, vertex3
    integer(kind=ik), intent(in) :: color_set               !< which color to set the mask
    integer(kind=ik), intent(in) :: smoothing_type_int     !< which mask do we use?
    real(kind=rk), intent(in) :: smoothing_width       !< width of smoothing region
    real(kind=rk), optional, intent(in) :: smoothing_safety      !< safety margin for smoothing
    real(kind=rk),optional,intent(in) :: bounding_box(1:4)
    integer, optional, intent(in) :: x0_indices(1:2)  !< For WABBIT, x0 is defined at the first interior point g+1, for FLUSI, it is actually at index 0, so we can provide here the indices of x0, defaults to WABBIT formulation

    ! auxiliary variables
    real(kind=rk)  :: xyz(1:2), dist, safety, tmp
    integer, dimension(1:2) :: lbounds, ubounds
    integer(kind=ik) :: ix, iy, Nsafety, bound_min(1:2), bound_max(1:2), x0_offset(1:2)

    ! set default value for x0_indices
    x0_offset = g+1
    if (present(x0_indices)) x0_offset = x0_indices

    ! Mask array is not resetted, as we may have other objects in the vicinity
    ! ! reset mask array
    ! mask = 0.0_rk

    ! safety for smoothing function
    safety = 3.0_rk * smoothing_width
    if (present(smoothing_safety)) safety = smoothing_safety
    Nsafety = ceiling(safety / minval(dx))

    ! bounds of the current patch of data
    lbounds = g+1
    ubounds = (/size(mask,1), size(mask,2)/) -g+1

    if (present(bounding_box)) then
        ! use pre-computed bounding box for rectangle (can be faster if many rectangles are to be drawn)
        bound_min = int((bounding_box(1:2)-x0)/dx) - Nsafety + g+1
        bound_max = int(ceiling((bounding_box(3:4)-x0)/dx)) + Nsafety + g+1
    else
        bound_min = int((min(min(vertex1, vertex2), vertex3) - x0)/dx) - Nsafety + g+1
        bound_max = int(ceiling((max(max(vertex1, vertex2), vertex3) - x0)/dx)) + Nsafety + g+1
    endif

    do iy = max(bound_min(2),lbounds(2)), min(bound_max(2),ubounds(2))
        xyz(2) = dble(iy-(x0_offset(2))) * dx(2) + x0(2)
        do ix = max(bound_min(1),lbounds(1)), min(bound_max(1),ubounds(1))
            xyz(1) = dble(ix-(x0_offset(1))) * dx(1) + x0(1)
            ! distance from center of rectangle
            dist = signed_distance_triangle_2D( xyz, vertex1, vertex2, vertex3 )
            ! apply mask function
            tmp = step(dist, 0.0_rk, smoothing_width, safety, smoothing_type_int)
            if (tmp >= mask(ix,iy)) then
                mask(ix,iy) = tmp    ! mask function
                color(ix,iy) = real(color_set, kind=rk) ! color
            endif
    end do; end do

end subroutine draw_triangle

end module module_geometry