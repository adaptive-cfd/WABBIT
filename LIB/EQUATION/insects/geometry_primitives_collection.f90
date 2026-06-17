!---------------------------------------
!> Init a primitives collection of several primitive shapes
!!
!! The function first reads in the geometric data in the init stage call and computes all bounding boxes.
!!
!! Inputs:
!! - `center` and `scale` for the geometry - this is used to shift and scale the geometry, so that we can work with normalized geometry data
!! - `color_set` which color to set in the mask
!! - `string`, `file`, `stage` geometry parameters
!! - `i_collection` which collection we are setting up, this is used to store the geometric data in the collection array for later use in the draw stage
!!
!! Outputs:
!! - Sets up private type primitives_collection of module_geometry for the given index
subroutine init_primitives_collection(center, scale, color_set, string, file, i_collection)
    use module_globals

    implicit none

    !> center and scale for the geometry - this is used to shift and scale the geometry, so that we can work with normalized geometry data
    real(kind=rk), intent(in) :: center(1:), scale
    integer(kind=ik), intent(in) :: color_set               !< which color to set the mask
    character(len=*), intent(in), optional :: string, file  !< geometry parameters (e.g. file name for geometry data, or string with geometry parameters)
    integer(kind=ik), intent(in) :: i_collection  !< which collection we are setting up

    ! stuff for geometric primitve collections
    character(len=clong), allocatable :: collection_string(:,:), collection_string_line(:)
    character(len=clong) :: line_parameters

    integer :: i_pos, i_geom, mpirank, mpicode
    real(kind=rk) :: e12(1:3), et(1:3)
    logical :: root

    ! for debug printing
    call MPI_COMM_RANK (MPI_COMM_WORLD, mpirank, mpicode)
    root = (mpirank == 0)

    ! this function is called only once, not for each block.
    ! performs initializations in the RHS module, such as resetting integrals
    if (geom_array(i_collection)%initialized) return  ! already initialized, nothing to do

    ! write basic parameters
    geom_array(i_collection)%color_set = color_set

    ! we cannot pass string and file or none
    if ((.not. present(string) .and. .not. present(file)) .or. (present(string) .and. present(file))) then
        call abort(260604, "init_primitives_collection: exactly one of string or file must be passed in the init stage!")
    end if

    if (root) then
        write(*,'(A, i0, A)') "INIT: Initializing geometric collection ", i_collection, " for mask generation..."
    end if

    ! at first, process the file and read everything in
    if (present(file)) then
        call check_file_exists(file)

        !*****************************************************************************
        ! phase one: read the number of lines (which is the number of geometries)
        !*****************************************************************************
        call count_lines_in_ascii_file_mpi(file, geom_array(i_collection)%n_objects, 0)

        !*****************************************************************************
        ! phase two: read all lines, here we cheat and read them in as a char array with one column per line (so the whole line is one string)
        !*****************************************************************************
        allocate(collection_string(geom_array(i_collection)%n_objects,1))
        call read_chararray_from_ascii_file_mpi(file, collection_string, n_header=0)
    elseif (present(string)) then
        !*****************************************************************************
        ! phase one: read the number of geometries
        ! count how many delimiters we have in the string to determine how many objects we have
        ! Thomas already used ; for comments, so here we have to be creative. I am using : as delimiter                !*****************************************************************************
        call count_entries(string, geom_array(i_collection)%n_objects, ":")
        allocate(collection_string_line(geom_array(i_collection)%n_objects))
        allocate(collection_string(geom_array(i_collection)%n_objects,1))
        !*****************************************************************************
        ! phase two: split the string into lines
        !*****************************************************************************
        call split_string(string, collection_string_line, separator_optional=":")
        collection_string(:,1) = collection_string_line(:)
        ! for inlines, we use comma as separators (since semicolons are comments), so we need to translate them to semicolons
        do i_geom = 1, geom_array(i_collection)%n_objects
            collection_string(i_geom,1) = str_replace_text(collection_string(i_geom,1), ",", ";")
        end do

        deallocate(collection_string_line)
    end if

    if (root) then
        write(*,'(A,I0)') "INIT: Number of objects in collection: ", geom_array(i_collection)%n_objects
    end if

    !*****************************************************************************
    ! phase 3: parse what object we have and extract the geometric parameters, then shift and scale them directly
    !*****************************************************************************
    allocate(geom_array(i_collection)%geometric_data(geom_array(i_collection)%n_objects, 20))  ! we allocate a fixed number of columns for the geometric parameters, which should be enough for all supported geometries. Unused columns are set to zero.
    allocate(geom_array(i_collection)%geometric_boundingbox(geom_array(i_collection)%n_objects, 6))
    allocate(geom_array(i_collection)%geometric_type(geom_array(i_collection)%n_objects))
    do i_geom = 1, geom_array(i_collection)%n_objects
        ! split at the first ';' - this is the object type, the rest are parameters
        i_pos = index(collection_string(i_geom,1), ";")
        if (i_pos == 0) call abort(260605, "init_primitives_collection: invalid collection string, no ';' found in object definition: "//trim(collection_string(i_geom,1)))
        geom_array(i_collection)%geometric_type(i_geom) = trim(adjustl(collection_string(i_geom,1)(1:i_pos-1)))
        line_parameters = str_replace_text(str_replace_text(collection_string(i_geom,1)(i_pos+1:),",","."),";",",")
        select case (geom_array(i_collection)%geometric_type(i_geom))
        case ("sphere")
            ! sphere: center(1:3), radius
            read(line_parameters,*) geom_array(i_collection)%geometric_data(i_geom,1:4)
            geom_array(i_collection)%geometric_data(i_geom,:) = geom_array(i_collection)%geometric_data(i_geom,:) * scale  ! scale everything
            geom_array(i_collection)%geometric_data(i_geom,1:3) = geom_array(i_collection)%geometric_data(i_geom,1:3) + center(1:3)  ! shift center
        case ("cylinder")
            ! cylinder: endpoint1(1:3), endpoint2(1:3), radius
            read(line_parameters,*) geom_array(i_collection)%geometric_data(i_geom,1:7)
            geom_array(i_collection)%geometric_data(i_geom,:) = geom_array(i_collection)%geometric_data(i_geom,:) * scale  ! scale everything
            geom_array(i_collection)%geometric_data(i_geom,1:3) = geom_array(i_collection)%geometric_data(i_geom,1:3) + center(1:3)  ! shift endpoint 1
            geom_array(i_collection)%geometric_data(i_geom,4:6) = geom_array(i_collection)%geometric_data(i_geom,4:6) + center(1:3)  ! shift endpoint 2
        case ("cylinder-rounded", "capsule")
            ! cylinder-rounded: endpoint1(1:3), endpoint2(1:3), radius
            read(line_parameters,*) geom_array(i_collection)%geometric_data(i_geom,1:7)
            geom_array(i_collection)%geometric_data(i_geom,:) = geom_array(i_collection)%geometric_data(i_geom,:) * scale  ! scale everything
            geom_array(i_collection)%geometric_data(i_geom,1:3) = geom_array(i_collection)%geometric_data(i_geom,1:3) + center(1:3)  ! shift endpoint 1
            geom_array(i_collection)%geometric_data(i_geom,4:6) = geom_array(i_collection)%geometric_data(i_geom,4:6) + center(1:3)  ! shift endpoint 2
        case ("circle")
            ! circle: center(1:2), radius
            read(line_parameters,*) geom_array(i_collection)%geometric_data(i_geom,1:3)
            geom_array(i_collection)%geometric_data(i_geom,:) = geom_array(i_collection)%geometric_data(i_geom,:) * scale  ! scale everything
            geom_array(i_collection)%geometric_data(i_geom,1:2) = geom_array(i_collection)%geometric_data(i_geom,1:2) + center(1:2)  ! shift center
        case ("rectangle")
            ! rectangle: center(1:2), half_size(1:2), angle
            read(line_parameters,*) geom_array(i_collection)%geometric_data(i_geom,1:5)
            geom_array(i_collection)%geometric_data(i_geom,1:4) = geom_array(i_collection)%geometric_data(i_geom,1:4) * scale  ! scale everything but not the angle
            geom_array(i_collection)%geometric_data(i_geom,1:2) = geom_array(i_collection)%geometric_data(i_geom,1:2) + center(1:2)  ! shift center
        case ("triangle")
            ! triangle: vertex1(1:2), vertex2(1:2), vertex3(1:2)
            read(line_parameters,*) geom_array(i_collection)%geometric_data(i_geom,1:6)
            geom_array(i_collection)%geometric_data(i_geom,:) = geom_array(i_collection)%geometric_data(i_geom,:) * scale  ! scale everything
            geom_array(i_collection)%geometric_data(i_geom,1:2) = geom_array(i_collection)%geometric_data(i_geom,1:2) + center(1:2)  ! shift vertex 1
            geom_array(i_collection)%geometric_data(i_geom,3:4) = geom_array(i_collection)%geometric_data(i_geom,3:4) + center(1:2)  ! shift vertex 2
            geom_array(i_collection)%geometric_data(i_geom,5:6) = geom_array(i_collection)%geometric_data(i_geom,5:6) + center(1:2)  ! shift vertex 3
        case default
            call abort(260606, "draw_collection: invalid collection string, unknown geometry type: "//trim(adjustl(geom_array(i_collection)%geometric_type(i_geom))))
        end select
    end do

    ! ! debug print to check everything
    ! if (root) then
    !     do i_geom = 1, geom_array(i_collection)%n_objects
    !         write(*,'(A,I0,3(A),20(f0.2, 1x))') "Object ", i_geom, ": ", trim(adjustl(geom_array(i_collection)%geometric_type(i_geom))), " , params: ", geom_array(i_collection)%geometric_data(i_geom,1:20)
    !     enddo
    ! endif

    !*****************************************************************************
    ! phase 4: compute the bounding box for each geometry
    !*****************************************************************************
    do i_geom = 1, geom_array(i_collection)%n_objects
        select case (trim(adjustl(geom_array(i_collection)%geometric_type(i_geom))))
        case ("sphere")
            geom_array(i_collection)%geometric_boundingbox(i_geom,1:3) = geom_array(i_collection)%geometric_data(i_geom,1:3) - geom_array(i_collection)%geometric_data(i_geom,4)  ! center - radius
            geom_array(i_collection)%geometric_boundingbox(i_geom,4:6) = geom_array(i_collection)%geometric_data(i_geom,1:3) + geom_array(i_collection)%geometric_data(i_geom,4)  ! center + radius
        case ("cylinder")
            e12 = geom_array(i_collection)%geometric_data(i_geom,1:3) - geom_array(i_collection)%geometric_data(i_geom,4:6)  ! endpoint1 - endpoint2
            et = sqrt(1.0_rk - (e12 / sqrt(dot_product(e12, e12)))**2)  ! normed tangential direction
            ! bounding box is the box around the endpoints plus the radius in the tangential direction
            geom_array(i_collection)%geometric_boundingbox(i_geom,1:3) = min(geom_array(i_collection)%geometric_data(i_geom,1:3), geom_array(i_collection)%geometric_data(i_geom,4:6))-geom_array(i_collection)%geometric_data(i_geom,7)*et
            geom_array(i_collection)%geometric_boundingbox(i_geom,4:6) = max(geom_array(i_collection)%geometric_data(i_geom,1:3), geom_array(i_collection)%geometric_data(i_geom,4:6))+geom_array(i_collection)%geometric_data(i_geom,7)*et
        case ("cylinder-rounded", "capsule")
            ! bounding box is the box around the endpoints plus/minus the radius in all directions
            geom_array(i_collection)%geometric_boundingbox(i_geom,1:3) = min(geom_array(i_collection)%geometric_data(i_geom,1:3), geom_array(i_collection)%geometric_data(i_geom,4:6))-geom_array(i_collection)%geometric_data(i_geom,7)
            geom_array(i_collection)%geometric_boundingbox(i_geom,4:6) = max(geom_array(i_collection)%geometric_data(i_geom,1:3), geom_array(i_collection)%geometric_data(i_geom,4:6))+geom_array(i_collection)%geometric_data(i_geom,7)
        case ("circle")
            geom_array(i_collection)%geometric_boundingbox(i_geom,1:2) = geom_array(i_collection)%geometric_data(i_geom,1:2) - geom_array(i_collection)%geometric_data(i_geom,3)  ! center - radius
            geom_array(i_collection)%geometric_boundingbox(i_geom,3:4) = geom_array(i_collection)%geometric_data(i_geom,1:2) + geom_array(i_collection)%geometric_data(i_geom,3)  ! center + radius
        case ("rectangle")
            ! center +- half_size in the direction of the angle
            geom_array(i_collection)%geometric_boundingbox(i_geom,1:2) = geom_array(i_collection)%geometric_data(i_geom,1:2) - abs(geom_array(i_collection)%geometric_data(i_geom,3:4))*cos(geom_array(i_collection)%geometric_data(i_geom,5)) - abs(geom_array(i_collection)%geometric_data(i_geom,4:3:-1))*sin(geom_array(i_collection)%geometric_data(i_geom,5))
            geom_array(i_collection)%geometric_boundingbox(i_geom,3:4) = geom_array(i_collection)%geometric_data(i_geom,1:2) + abs(geom_array(i_collection)%geometric_data(i_geom,3:4))*cos(geom_array(i_collection)%geometric_data(i_geom,5)) + abs(geom_array(i_collection)%geometric_data(i_geom,4:3:-1))*sin(geom_array(i_collection)%geometric_data(i_geom,5))
        case ("triangle")
            geom_array(i_collection)%geometric_boundingbox(i_geom,1) = minval(geom_array(i_collection)%geometric_data(i_geom,1:6:2))  ! min of x-coordinates of vertices
            geom_array(i_collection)%geometric_boundingbox(i_geom,2) = maxval(geom_array(i_collection)%geometric_data(i_geom,1:6:2))  ! max of x-coordinates of vertices
            geom_array(i_collection)%geometric_boundingbox(i_geom,3) = minval(geom_array(i_collection)%geometric_data(i_geom,2:6:2))  ! min of y-coordinates of vertices
            geom_array(i_collection)%geometric_boundingbox(i_geom,4) = maxval(geom_array(i_collection)%geometric_data(i_geom,2:6:2))  ! max of y-coordinates of vertices
        end select
    end do

    ! ! debug print
    ! if (root) then
    !     do i_geom = 1, geom_array(i_collection)%n_objects
    !         write(*,'(A,I0,3(A),6(f0.2, 1x))') "Object ", i_geom, ": ", trim(adjustl(geom_array(i_collection)%geometric_type(i_geom))), " , bounding box: ", geom_array(i_collection)%geometric_boundingbox(i_geom,1:6)
    !     enddo
    ! endif

    ! at the end of this stage, we return
    geom_array(i_collection)%initialized = .true.
    deallocate(collection_string)

end subroutine init_primitives_collection


!! Every subsequent call then draws every object individually with a pointwise implementation that loops over all grid points in the bounding box of the geometry and applies the mask function based on the signed distance to the sphere.
!! Inputs:
!! - `mask` mask array of a block to be updated
!! - `x0`, `dx`, `Bs`, `g` origin and spacing, Block size, ghost point size of the block
!! - `smoothing_type`, `smoothing_width`, `smoothing_safety` parameters for the step function
!! - `i_collection` which collection we are drawing
!!
!! Outputs:
!! - `mask` updated with the geometry, ID-1 is mask array, ID-5 is color
subroutine draw_primitives_collection(mask, x0, dx, Bs, g, smoothing_type_int, smoothing_width, smoothing_safety, i_collection)
    use module_globals
    implicit none

    !> mask array of a block to be updated
    real(kind=rk), intent(inout) :: mask(:,:,:,:)
    real(kind=rk), intent(in) :: x0(:), dx(:)    !< origin and spacing of the block
    integer(kind=ik), intent(in) :: Bs(1:3), g  !< Block size, ghost point size of the block
    integer(kind=ik), intent(in) :: smoothing_type_int
    real(kind=rk), intent(in) :: smoothing_width
    real(kind=rk), intent(in), optional :: smoothing_safety
    integer(kind=ik), intent(in) :: i_collection  !< which collection we are drawing

    integer :: i_geom

    !-------------------------------------------------------------------------
    ! let's draw every object one after another
    do i_geom = 1, geom_array(i_collection)%n_objects
        ! names have to be standardized here, meaning all lower-case and no "_"
        select case (trim(standardize_string(geom_array(i_collection)%geometric_type(i_geom))))
        case ("sphere")
            call draw_sphere(mask(:,:,:,1),mask(:,:,:,5), x0, dx, g, geom_array(i_collection)%geometric_data(i_geom,1:3), geom_array(i_collection)%geometric_data(i_geom,4), color_set=geom_array(i_collection)%color_set, smoothing_type_int=smoothing_type_int, smoothing_width=smoothing_width, smoothing_safety=smoothing_safety, bounding_box=geom_array(i_collection)%geometric_boundingbox(i_geom,:))
        case ("cylinder")
            call draw_cylinder(mask(:,:,:,1),mask(:,:,:,5), x0, dx, g, geom_array(i_collection)%geometric_data(i_geom,1:3), geom_array(i_collection)%geometric_data(i_geom,4:6), geom_array(i_collection)%geometric_data(i_geom,7), color_set=geom_array(i_collection)%color_set, smoothing_type_int=smoothing_type_int, smoothing_width=smoothing_width, smoothing_safety=smoothing_safety, bounding_box=geom_array(i_collection)%geometric_boundingbox(i_geom,:))
        case ("capsule", "cylinder-rounded")
            call draw_cylinder_rounded(mask(:,:,:,1),mask(:,:,:,5), x0, dx, g, geom_array(i_collection)%geometric_data(i_geom,1:3), geom_array(i_collection)%geometric_data(i_geom,4:6), geom_array(i_collection)%geometric_data(i_geom,7), color_set=geom_array(i_collection)%color_set, smoothing_type_int=smoothing_type_int, smoothing_width=smoothing_width, smoothing_safety=smoothing_safety, bounding_box=geom_array(i_collection)%geometric_boundingbox(i_geom,:))
        case ("circle")
            call draw_circle(mask(:,:,1,1),mask(:,:,1,5), x0(1:2), dx(1:2), g, geom_array(i_collection)%geometric_data(i_geom,1:2), geom_array(i_collection)%geometric_data(i_geom,3), color_set=geom_array(i_collection)%color_set, smoothing_type_int=smoothing_type_int, smoothing_width=smoothing_width, smoothing_safety=smoothing_safety, bounding_box=geom_array(i_collection)%geometric_boundingbox(i_geom,:))
        case ("rectangle")
            call draw_rectangle(mask(:,:,1,1),mask(:,:,1,5), x0(1:2), dx(1:2), g, geom_array(i_collection)%geometric_data(i_geom,1:2), geom_array(i_collection)%geometric_data(i_geom,3:4), geom_array(i_collection)%geometric_data(i_geom,5), color_set=geom_array(i_collection)%color_set, smoothing_type_int=smoothing_type_int, smoothing_width=smoothing_width, smoothing_safety=smoothing_safety, bounding_box=geom_array(i_collection)%geometric_boundingbox(i_geom,:))
        case ("triangle")
            call draw_triangle(mask(:,:,1,1),mask(:,:,1,5), x0(1:2), dx(1:2), g, geom_array(i_collection)%geometric_data(i_geom,1:2), geom_array(i_collection)%geometric_data(i_geom,3:4), geom_array(i_collection)%geometric_data(i_geom,5:6), color_set=geom_array(i_collection)%color_set, smoothing_type_int=smoothing_type_int, smoothing_width=smoothing_width, smoothing_safety=smoothing_safety, bounding_box=geom_array(i_collection)%geometric_boundingbox(i_geom,:))
        end select
    end do

end subroutine draw_primitives_collection


!> This function is called for every geometry in the collection and should check if the geometry is present in the block and set geometry_in_block accordingly
subroutine primitives_collection_geometry_indicator(time, i_collection, BS, g, x0, dx, dim, geometry_in_block)
    use module_globals
    implicit none

    real(kind=rk), intent(in) :: time
    integer(kind=ik), intent(in) :: i_collection
    integer, intent(in) :: BS(1:3), g
    real(kind=rk), intent(in) :: x0(:), dx(:)
    integer(kind=ik), intent(in) :: dim
    logical, intent(out) :: geometry_in_block

    integer(kind=ik) :: i_geom
    real(kind=rk) :: xend(1:dim), block_extent(1:dim), point_check(1:dim), vec_normed(1:dim)

    ! compute end of blocks
    block_extent(1:dim) = dx(1:dim) * real(BS(1:dim), kind=rk)
    xend(1:dim) = x0(1:dim) + block_extent(1:dim)

    ! It would be tempting just to point to the center locations, but we actually need to point to a point on the BOUNDARY of the geometry
    do i_geom = 1, geom_array(i_collection)%n_objects
        select case (trim(standardize_string(geom_array(i_collection)%geometric_type(i_geom))))
        case ("sphere", "circle")
            ! boundary + r in one direction
            point_check(1:dim) = geom_array(i_collection)%geometric_data(i_geom,1:dim)
            point_check(1) = point_check(1) + geom_array(i_collection)%geometric_data(i_geom,4)  ! center + radius in x-direction
        case ("cylinder", "triangle")
            ! endpoint_1 or vertex_1 is part of the boundary
            point_check(1:dim) = geom_array(i_collection)%geometric_data(i_geom,1:dim)
        case ("cylinder-rounded", "capsule")
            ! endpoint_1 is not included, we need to extend it by the radius in the negative direction of the 
            ! first construct direction vector
            vec_normed = geom_array(i_collection)%geometric_data(i_geom,4:6) - geom_array(i_collection)%geometric_data(i_geom,1:3)
            ! now norm it
            vec_normed = vec_normed / sqrt(dot_product(vec_normed, vec_normed))
            point_check(1:dim) = geom_array(i_collection)%geometric_data(i_geom,1:dim) - geom_array(i_collection)%geometric_data(i_geom,7)*vec_normed(1:dim)  ! endpoint_1 minus radius in the direction of the cylinder
        case ("rectangle")
            ! center + first half_size in the direction of the angle
            point_check(1:dim) = geom_array(i_collection)%geometric_data(i_geom,1:dim)
            point_check(1) = point_check(1) + abs(geom_array(i_collection)%geometric_data(i_geom,3))*cos(geom_array(i_collection)%geometric_data(i_geom,5))
            point_check(2) = point_check(2) + abs(geom_array(i_collection)%geometric_data(i_geom,3))*sin(geom_array(i_collection)%geometric_data(i_geom,5))
        end select
        if (all(point_check(1:dim) >= x0(1:dim)) .and. all(point_check(1:dim) <= xend(1:dim))) then
            geometry_in_block = .true.
            return
        end if
    end do

end subroutine primitives_collection_geometry_indicator