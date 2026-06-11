!---------------------------------------
!> Init a composition of several primitive shapes
!!
!! The function first reads in the geometric data in the init stage call and computes all bounding boxes.
!!
!! Inputs:
!! - `center` and `scale` for the geometry - this is used to shift and scale the geometry, so that we can work with normalized geometry data
!! - `color_set` which color to set in the mask
!! - `string`, `file`, `stage` geometry parameters
!! - `i_composition` which composition we are setting up, this is used to store the geometric data in the composition array for later use in the draw stage
!!
!! Outputs:
!! - Sets up private type composition of module_geometry for the given index
subroutine init_composition(center, scale, color_set, string, file, i_composition)
    use module_globals

    implicit none

    !> center and scale for the geometry - this is used to shift and scale the geometry, so that we can work with normalized geometry data
    real(kind=rk), intent(in) :: center(1:), scale
    integer(kind=ik), intent(in) :: color_set               !< which color to set the mask
    character(len=*), intent(in), optional :: string, file  !< geometry parameters (e.g. file name for geometry data, or string with geometry parameters)
    integer(kind=ik), intent(in) :: i_composition  !< which composition we are setting up

    ! stuff for geometric compositions
    character(len=clong), allocatable :: composition_string(:,:), composition_string_line(:)

    integer :: i_pos, i_geom, mpirank, mpicode
    real(kind=rk) :: e12(1:3), et(1:3)
    logical :: root

    ! for debug printing
    call MPI_COMM_RANK (MPI_COMM_WORLD, mpirank, mpicode)
    root = (mpirank == 0)

    ! this function is called only once, not for each block.
    ! performs initializations in the RHS module, such as resetting integrals
    if (geom_array(i_composition)%initialized) return  ! already initialized, nothing to do

    ! write basic parameters
    geom_array(i_composition)%color_set = color_set

    ! we cannot pass string and file or none
    if ((.not. present(string) .and. .not. present(file)) .or. (present(string) .and. present(file))) then
        call abort(260604, "draw_composition: exactly one of string or file must be passed in the init stage!")
    end if

    if (root) then
        write(*,'(A, i0, A)') "INIT: Initializing geometric composition ", i_composition, " for mask generation..."
    end if

    ! at first, process the file and read everything in
    if (present(file)) then
        call check_file_exists(file)

        !*****************************************************************************
        ! phase one: read the number of lines (which is the number of geometries)
        !*****************************************************************************
        call count_lines_in_ascii_file_mpi(file, geom_array(i_composition)%n_objects, 0)

        !*****************************************************************************
        ! phase two: read all lines, here we cheat and read them in as a char array with one column per line (so the whole line is one string)
        !*****************************************************************************
        allocate(composition_string(geom_array(i_composition)%n_objects,1))
        call read_chararray_from_ascii_file_mpi(file, composition_string, n_header=0)
    elseif (present(string)) then
        !*****************************************************************************
        ! phase one: read the number of geometries
        ! count how many delimiters we have in the string to determine how many objects we have
        ! Thomas already used ; for comments, so here we have to be creative. I am using : as delimiter                !*****************************************************************************
        call count_entries(string, geom_array(i_composition)%n_objects, ":")
        allocate(composition_string_line(geom_array(i_composition)%n_objects))
        allocate(composition_string(geom_array(i_composition)%n_objects,1))
        !*****************************************************************************
        ! phase two: split the string into lines
        !*****************************************************************************
        call split_string(string, composition_string_line, separator_optional=":")
        composition_string(:,1) = composition_string_line(:)
        deallocate(composition_string_line)
    end if

    if (root) then
        write(*,'(A,I0)') "INIT: Number of objects in composition: ", geom_array(i_composition)%n_objects
    end if

    !*****************************************************************************
    ! phase 3: parse what object we have and extract the geometric parameters, then shift and scale them directly
    !*****************************************************************************
    allocate(geom_array(i_composition)%geometric_data(geom_array(i_composition)%n_objects, 20))  ! we allocate a fixed number of columns for the geometric parameters, which should be enough for all supported geometries. Unused columns are set to zero.
    allocate(geom_array(i_composition)%geometric_boundingbox(geom_array(i_composition)%n_objects, 6))
    allocate(geom_array(i_composition)%geometric_type(geom_array(i_composition)%n_objects))
    do i_geom = 1, geom_array(i_composition)%n_objects
        ! split at the first ',' - this is the object type, the rest are parameters
        i_pos = index(composition_string(i_geom,1), ",")
        if (i_pos == 0) call abort(260605, "draw_composition: invalid composition string, no ',' found in object definition: "//trim(composition_string(i_geom,1)))
        geom_array(i_composition)%geometric_type(i_geom) = trim(adjustl(composition_string(i_geom,1)(1:i_pos-1)))
        select case (geom_array(i_composition)%geometric_type(i_geom))
        case ("sphere")
            ! sphere: center(1:3), radius
            read(composition_string(i_geom,1)(i_pos+1:),*) geom_array(i_composition)%geometric_data(i_geom,1:4)
            geom_array(i_composition)%geometric_data(i_geom,:) = geom_array(i_composition)%geometric_data(i_geom,:) * scale  ! scale everything
            geom_array(i_composition)%geometric_data(i_geom,1:3) = geom_array(i_composition)%geometric_data(i_geom,1:3) + center(1:3)  ! shift center
        case ("cylinder")
            ! cylinder: endpoint1(1:3), endpoint2(1:3), radius
            read(composition_string(i_geom,1)(i_pos+1:),*) geom_array(i_composition)%geometric_data(i_geom,1:7)
            geom_array(i_composition)%geometric_data(i_geom,:) = geom_array(i_composition)%geometric_data(i_geom,:) * scale  ! scale everything
            geom_array(i_composition)%geometric_data(i_geom,1:3) = geom_array(i_composition)%geometric_data(i_geom,1:3) + center(1:3)  ! shift endpoint 1
            geom_array(i_composition)%geometric_data(i_geom,4:6) = geom_array(i_composition)%geometric_data(i_geom,4:6) + center(1:3)  ! shift endpoint 2
        case ("capsule")
            ! capsule: endpoint1(1:3), endpoint2(1:3), radius
            read(composition_string(i_geom,1)(i_pos+1:),*) geom_array(i_composition)%geometric_data(i_geom,1:7)
            geom_array(i_composition)%geometric_data(i_geom,:) = geom_array(i_composition)%geometric_data(i_geom,:) * scale  ! scale everything
            geom_array(i_composition)%geometric_data(i_geom,1:3) = geom_array(i_composition)%geometric_data(i_geom,1:3) + center(1:3)  ! shift endpoint 1
            geom_array(i_composition)%geometric_data(i_geom,4:6) = geom_array(i_composition)%geometric_data(i_geom,4:6) + center(1:3)  ! shift endpoint 2
        case ("circle")
            ! circle: center(1:2), radius
            read(composition_string(i_geom,1)(i_pos+1:),*) geom_array(i_composition)%geometric_data(i_geom,1:3)
            geom_array(i_composition)%geometric_data(i_geom,:) = geom_array(i_composition)%geometric_data(i_geom,:) * scale  ! scale everything
            geom_array(i_composition)%geometric_data(i_geom,1:2) = geom_array(i_composition)%geometric_data(i_geom,1:2) + center(1:2)  ! shift center
        case ("rectangle")
            ! rectangle: center(1:2), half_size(1:2), angle
            read(composition_string(i_geom,1)(i_pos+1:),*) geom_array(i_composition)%geometric_data(i_geom,1:5)
            geom_array(i_composition)%geometric_data(i_geom,1:4) = geom_array(i_composition)%geometric_data(i_geom,1:4) * scale  ! scale everything but not the angle
            geom_array(i_composition)%geometric_data(i_geom,1:2) = geom_array(i_composition)%geometric_data(i_geom,1:2) + center(1:2)  ! shift center
        case ("triangle")
            ! triangle: vertex1(1:2), vertex2(1:2), vertex3(1:2)
            read(composition_string(i_geom,1)(i_pos+1:),*) geom_array(i_composition)%geometric_data(i_geom,1:6)
            geom_array(i_composition)%geometric_data(i_geom,:) = geom_array(i_composition)%geometric_data(i_geom,:) * scale  ! scale everything
            geom_array(i_composition)%geometric_data(i_geom,1:2) = geom_array(i_composition)%geometric_data(i_geom,1:2) + center(1:2)  ! shift vertex 1
            geom_array(i_composition)%geometric_data(i_geom,3:4) = geom_array(i_composition)%geometric_data(i_geom,3:4) + center(1:2)  ! shift vertex 2
            geom_array(i_composition)%geometric_data(i_geom,5:6) = geom_array(i_composition)%geometric_data(i_geom,5:6) + center(1:2)  ! shift vertex 3
        case default
            call abort(260606, "draw_composition: invalid composition string, unknown geometry type: "//geom_array(i_composition)%geometric_type(i_geom))
        end select
    end do

    ! ! debug print to check everything
    ! if (root) then
    !     do i_geom = 1, geom_array(i_composition)%n_objects
    !         write(*,'(A,I0,3(A),20(f0.2, 1x))') "Object ", i_geom, ": ", trim(adjustl(geom_array(i_composition)%geometric_type(i_geom))), " , params: ", geom_array(i_composition)%geometric_data(i_geom,1:20)
    !     enddo
    ! endif

    !*****************************************************************************
    ! phase 4: compute the bounding box for each geometry
    !*****************************************************************************
    do i_geom = 1, geom_array(i_composition)%n_objects
        select case (geom_array(i_composition)%geometric_type(i_geom))
        case ("sphere")
            geom_array(i_composition)%geometric_boundingbox(i_geom,1:3) = geom_array(i_composition)%geometric_data(i_geom,1:3) - geom_array(i_composition)%geometric_data(i_geom,4)  ! center - radius
            geom_array(i_composition)%geometric_boundingbox(i_geom,4:6) = geom_array(i_composition)%geometric_data(i_geom,1:3) + geom_array(i_composition)%geometric_data(i_geom,4)  ! center + radius
        case ("cylinder")
            e12 = geom_array(i_composition)%geometric_data(i_geom,1:3) - geom_array(i_composition)%geometric_data(i_geom,4:6)  ! endpoint1 - endpoint2
            et = sqrt(1.0_rk - (e12 / sqrt(dot_product(e12, e12)))**2)  ! normed tangential direction
            ! bounding box is the box around the endpoints plus the radius in the tangential direction
            geom_array(i_composition)%geometric_boundingbox(i_geom,1:3) = min(geom_array(i_composition)%geometric_data(i_geom,1:3), geom_array(i_composition)%geometric_data(i_geom,4:6))-geom_array(i_composition)%geometric_data(i_geom,7)*et
            geom_array(i_composition)%geometric_boundingbox(i_geom,4:6) = max(geom_array(i_composition)%geometric_data(i_geom,1:3), geom_array(i_composition)%geometric_data(i_geom,4:6))+geom_array(i_composition)%geometric_data(i_geom,7)*et
        case ("capsule")
            ! bounding box is the box around the endpoints plus/minus the radius in all directions
            geom_array(i_composition)%geometric_boundingbox(i_geom,1:3) = min(geom_array(i_composition)%geometric_data(i_geom,1:3), geom_array(i_composition)%geometric_data(i_geom,4:6))-geom_array(i_composition)%geometric_data(i_geom,7)
            geom_array(i_composition)%geometric_boundingbox(i_geom,4:6) = max(geom_array(i_composition)%geometric_data(i_geom,1:3), geom_array(i_composition)%geometric_data(i_geom,4:6))+geom_array(i_composition)%geometric_data(i_geom,7)
        case ("circle")
            geom_array(i_composition)%geometric_boundingbox(i_geom,1:2) = geom_array(i_composition)%geometric_data(i_geom,1:2) - geom_array(i_composition)%geometric_data(i_geom,3)  ! center - radius
            geom_array(i_composition)%geometric_boundingbox(i_geom,3:4) = geom_array(i_composition)%geometric_data(i_geom,1:2) + geom_array(i_composition)%geometric_data(i_geom,3)  ! center + radius
        case ("rectangle")
            ! center +- half_size in the direction of the angle
            geom_array(i_composition)%geometric_boundingbox(i_geom,1:2) = geom_array(i_composition)%geometric_data(i_geom,1:2) - abs(geom_array(i_composition)%geometric_data(i_geom,3:4))*cos(geom_array(i_composition)%geometric_data(i_geom,5)) - abs(geom_array(i_composition)%geometric_data(i_geom,4:3:-1))*sin(geom_array(i_composition)%geometric_data(i_geom,5))
            geom_array(i_composition)%geometric_boundingbox(i_geom,3:4) = geom_array(i_composition)%geometric_data(i_geom,1:2) + abs(geom_array(i_composition)%geometric_data(i_geom,3:4))*cos(geom_array(i_composition)%geometric_data(i_geom,5)) + abs(geom_array(i_composition)%geometric_data(i_geom,4:3:-1))*sin(geom_array(i_composition)%geometric_data(i_geom,5))
        case ("triangle")
            geom_array(i_composition)%geometric_boundingbox(i_geom,1) = minval(geom_array(i_composition)%geometric_data(i_geom,1:6:2))  ! min of x-coordinates of vertices
            geom_array(i_composition)%geometric_boundingbox(i_geom,2) = maxval(geom_array(i_composition)%geometric_data(i_geom,1:6:2))  ! max of x-coordinates of vertices
            geom_array(i_composition)%geometric_boundingbox(i_geom,3) = minval(geom_array(i_composition)%geometric_data(i_geom,2:6:2))  ! min of y-coordinates of vertices
            geom_array(i_composition)%geometric_boundingbox(i_geom,4) = maxval(geom_array(i_composition)%geometric_data(i_geom,2:6:2))  ! max of y-coordinates of vertices
        end select
    end do

    ! ! debug print
    ! if (root) then
    !     do i_geom = 1, geom_array(i_composition)%n_objects
    !         write(*,'(A,I0,3(A),6(f0.2, 1x))') "Object ", i_geom, ": ", trim(adjustl(geom_array(i_composition)%geometric_type(i_geom))), " , bounding box: ", geom_array(i_composition)%geometric_boundingbox(i_geom,1:6)
    !     enddo
    ! endif

    ! at the end of this stage, we return
    geom_array(i_composition)%initialized = .true.
    deallocate(composition_string)

end subroutine init_composition


!! Every subsequent call then draws every object individually with a pointwise implementation that loops over all grid points in the bounding box of the geometry and applies the mask function based on the signed distance to the sphere.
!! Inputs:
!! - `mask` mask array of a block to be updated
!! - `x0`, `dx`, `Bs`, `g` origin and spacing, Block size, ghost point size of the block
!! - `smoothing_type`, `smoothing_width`, `smoothing_safety` parameters for the step function
!! - `i_composition` which composition we are drawing
!!
!! Outputs:
!! - `mask` updated with the geometry, ID-1 is mask array, ID-5 is color
subroutine draw_composition(mask, x0, dx, Bs, g, smoothing_type, smoothing_width, smoothing_safety, i_composition)
    use module_globals
    implicit none

    !> mask array of a block to be updated
    real(kind=rk), intent(inout) :: mask(:,:,:,:)
    real(kind=rk), intent(in) :: x0(:), dx(:)    !< origin and spacing of the block
    integer(kind=ik), intent(in) :: Bs(1:3), g  !< Block size, ghost point size of the block
    character(len=*), intent(in) :: smoothing_type
    real(kind=rk), intent(in) :: smoothing_width
    real(kind=rk), intent(in), optional :: smoothing_safety
    integer(kind=ik), intent(in) :: i_composition  !< which composition we are drawing

    integer :: i_geom

    !-------------------------------------------------------------------------
    ! let's draw every object one after another
    do i_geom = 1, geom_array(i_composition)%n_objects
        select case (geom_array(i_composition)%geometric_type(i_geom))
        case ("sphere")
            call draw_sphere(mask, x0, dx, Bs, g, geom_array(i_composition)%geometric_data(i_geom,1:3), geom_array(i_composition)%geometric_data(i_geom,4), color_set=geom_array(i_composition)%color_set, smoothing_type=smoothing_type, smoothing_width=smoothing_width, smoothing_safety=smoothing_safety, bounding_box=geom_array(i_composition)%geometric_boundingbox(i_geom,:))
        case ("cylinder")
            call draw_cylinder(mask, x0, dx, Bs, g, geom_array(i_composition)%geometric_data(i_geom,1:3), geom_array(i_composition)%geometric_data(i_geom,4:6), geom_array(i_composition)%geometric_data(i_geom,7), color_set=geom_array(i_composition)%color_set, smoothing_type=smoothing_type, smoothing_width=smoothing_width, smoothing_safety=smoothing_safety, bounding_box=geom_array(i_composition)%geometric_boundingbox(i_geom,:))
        case ("capsule")
            call draw_capsule(mask, x0, dx, Bs, g, geom_array(i_composition)%geometric_data(i_geom,1:3), geom_array(i_composition)%geometric_data(i_geom,4:6), geom_array(i_composition)%geometric_data(i_geom,7), color_set=geom_array(i_composition)%color_set, smoothing_type=smoothing_type, smoothing_width=smoothing_width, smoothing_safety=smoothing_safety, bounding_box=geom_array(i_composition)%geometric_boundingbox(i_geom,:))
        case ("circle")
            call draw_circle(mask(:,:,1,:), x0(1:2), dx(1:2), Bs, g, geom_array(i_composition)%geometric_data(i_geom,1:2), geom_array(i_composition)%geometric_data(i_geom,3), color_set=geom_array(i_composition)%color_set, smoothing_type=smoothing_type, smoothing_width=smoothing_width, smoothing_safety=smoothing_safety, bounding_box=geom_array(i_composition)%geometric_boundingbox(i_geom,:))
        case ("rectangle")
            call draw_rectangle(mask(:,:,1,:), x0(1:2), dx(1:2), Bs, g, geom_array(i_composition)%geometric_data(i_geom,1:2), geom_array(i_composition)%geometric_data(i_geom,3:4), geom_array(i_composition)%geometric_data(i_geom,5), color_set=geom_array(i_composition)%color_set, smoothing_type=smoothing_type, smoothing_width=smoothing_width, smoothing_safety=smoothing_safety, bounding_box=geom_array(i_composition)%geometric_boundingbox(i_geom,:))
        case ("triangle")
            call draw_triangle(mask(:,:,1,:), x0(1:2), dx(1:2), Bs, g, geom_array(i_composition)%geometric_data(i_geom,1:2), geom_array(i_composition)%geometric_data(i_geom,3:4), geom_array(i_composition)%geometric_data(i_geom,5:6), color_set=geom_array(i_composition)%color_set, smoothing_type=smoothing_type, smoothing_width=smoothing_width, smoothing_safety=smoothing_safety, bounding_box=geom_array(i_composition)%geometric_boundingbox(i_geom,:))
        end select
    end do

end subroutine draw_composition