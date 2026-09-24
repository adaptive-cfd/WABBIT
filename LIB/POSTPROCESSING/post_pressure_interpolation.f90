subroutine post_pressure_interpolation(params)

    use mpi
    use module_helpers
    use module_MPI
    use module_params               ! global parameters
    use module_timing
    use module_mesh                 ! mesh manipulation subroutines
    use module_time_step
    use module_unit_test
    use module_bridge_interface     ! bridge implementation of wabbit
    use module_forestMetaData
    use module_insects

    implicit none

    type (type_params), intent(inout)   :: params           ! user defined parameter structure

    integer(kind=ik)                    :: number_procs     ! number of processes
    ! real(kind=rk)                       :: t0, t1, t2       ! cpu time variables for running time calculation

    real(kind=rk), allocatable          :: hvy_mask(:, :, :, :, :), hvy_tmp(:, :, :, :, :), hvy_block(:, :, :, :, :)
    real(kind=rk)                       :: time, xx, yy, zz, delx, dely, delz, tmp           ! time loop variables
    character(len=cshort)               :: pressure_filename, ini_filename, wing_fname, fname_out
    integer(kind=ik)                    :: k, lgt_id, Bs(1:3), g, hvy_id, iter, iteration, tree_ID=1, N_support, nlines, ncols, n_blocks
    real(kind=rk)                       :: x(1:3), x0(1:3), dx(1:3), x_wing_w(1:3), x_wing_b(1:3), x_wing_g(1:3), x_wing_normal(1:3), coeff(1:3)
    real(kind=rk)                       :: block_x_min(1:3), block_x_max(1:3)
    integer(kind=ik)                    :: ipoint, ix, iy, iz, ix0, iy0, iz0
    logical                             :: help1, help2
    integer(kind=ik)                    :: nz, dim, tc_length
    real(kind=rk), dimension(3)         :: domain
    type(inifile) :: FILE
    character(len=cshort)               :: wing_type
    integer(kind=ik)                    :: surface_type    ! 1=bottom, 2=middle 3=top surface
    real(kind=rk), allocatable          :: wing_points_w(:,:), pressure_data(:,:), wing_points_g(:,:,:), xq(:,:), pq(:,:)
    integer(kind=ik)                    :: mpierr, n_header, isurface, npoints, nsurfaces

    allocate( hvy_n(1) )

    !---------------------------------------------------------------------------
    ! If called with '--help' or '-h', print a help message and exit.
    !---------------------------------------------------------------------------
    call get_cmd_arg( "--help", help1, default=.false. )
    call get_cmd_arg( "-h", help2, default=.false. )

    if ((help1 .or. help2)) then
        if (params%rank==0) then
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Wabbit postprocessing: extract (interpolate) wing pressure"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Given a txt file with coordinates on the wing surface (the midline, i.e. zw==0)"
            write(*, '(A)') " this routine interpolates a given wabbit field (usually this will be the pressure)"
            write(*, '(A)') " at the top, bottom and middle surface of the. The result is stored to CSV file."
            write(*, '(A)') " The input file is a semicolon-separated file with two columns, xw and yw, with one header line."
            write(*, '(A)') " "
            write(*, '(A)') " "
            write(*, '(A)') " Call:"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') "./wabbit-post --wing-pressure-interpolation INPUT_FIELD.h5 PARAMS.ini wing_query_points.txt output.csv --wing=[right/left] --memory=64.0GB"
            write(*, '(A)') "-----------------------------------------------------------"
        end if
        return
    endif

    !---------------------------------------------------------------------------
    ! Initialize parameters,bridge and grid
    !---------------------------------------------------------------------------
    ! read in the parameter file to setup the case
    ! get the second command line argument: this should be the ini-file name
    call get_command_argument( 2, pressure_filename )
    call get_command_argument( 3, ini_filename )
    call get_command_argument( 4, wing_fname )
    call check_file_exists( pressure_filename )
    call check_file_exists( ini_filename )
    call check_file_exists( wing_fname )
    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, ini_filename )
    params%n_eqn = 1

    call get_cmd_arg( "--wing", wing_type, default="right" )

    ! setup the wavelet etc
    call setup_wavelet(params)

    Bs = params%Bs
    g  = params%g
    ! Jmax = params%Jmax
    ! Jmin = params%Jmin
    !//FIXME tree_n = 1
    tree_n = params%forest_size ! used only for resetting at this point

    ! initializes the communicator for Wabbit and creates a bridge if needed
    call initialize_communicator(params)
    ! have the pysics module read their own parameters
    call init_physics_modules( params, ini_filename, params%N_mask_components )

    ! allocate memory for heavy, light, work and neighbor data
    call allocate_forest(params, hvy_block)

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    call read_attributes(pressure_filename, n_blocks, time, iteration, domain, Bs, tc_length, dim, &
         periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)
    ! reset the grid: all blocks are inactive and empty
    ! call reset_tree(params, .true., tree_ID)

    ! check if blocksize matches (otherwise reading will fail)
    if ( any( Bs-params%Bs  > 0) ) then
        call abort(180925, "Block size in file and parameter file do not match!")
    endif

    ! read input data
    ! read in pressure data
    call readHDF5vct_tree( (/pressure_filename/), params, hvy_block, tree_ID)

    ! let's update all insects. If there is none, then this is just an empty loop, so no problemo
    call Update_All_Insects(time)

    ! BEFORE WE CAN INTERPOLATE THE GHOTS NODES NEED TO BE FILLED
    call sync_ghosts_tree( params, hvy_block, tree_ID=tree_ID )

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! prepare pointcloud for interpolation (one for each surface)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! read in the wing points
    ! in the wing system
    n_header = 1
    call count_lines_in_ascii_file_mpi(wing_fname, nlines, n_header)
    call count_cols_in_ascii_file_mpi(wing_fname, ncols, n_header)

    if (ncols .ne. 2) then
        call abort(202512081, "Input file should just contain list of points on the wing midplance (xw,yw), zw==0")
    endif

    npoints = nlines

    allocate(wing_points_w( 1:nlines, 1:3) )
    allocate(wing_points_g( 1:nlines, 1:3, 1:3) ) ! for the 3 surfaces
    allocate(pressure_data( 1:nlines, 1:3) ) ! result for all 3 surfaces

    ! wing_points_w: xw, yw
    call read_array_from_ascii_file_mpi(wing_fname, wing_points_w(:,1:2), n_header)

    nsurfaces = 3

    ! always interpolate all surfaces
    do surface_type = 1, nsurfaces
        ! transform point data to global system 
        if (wing_type == "right") then
            do ipoint = 1, npoints
                ! input is in wing system, midplane of the wing
                x_wing_w(1:3) = (/ wing_points_w(ipoint, 1:2), 0.0_rk /)
                x_wing_normal(1:3) = (/ 0.0_rk, 0.0_rk, 1.0_rk /)
                
                ! //NOTE the normal_vector direction is to outside the surface, should be reversed
                if (surface_type == 2) then                    
                    x_wing_w(1:3) = x_wing_w(1:3)
                else if (surface_type == 3) then
                    ! Wings: 1=left 2=right 3=left hind 4=right hind
                    x_wing_w(1:3) = x_wing_w(1:3) - 0.5_rk*Insects(1)%Wings(2)%WingThickness * x_wing_normal(1:3)
                else if (surface_type == 1) then
                    x_wing_w(1:3) = x_wing_w(1:3) + 0.5_rk*Insects(1)%Wings(2)%WingThickness * x_wing_normal(1:3)
                else
                    call abort(372936, "surface type must be 1=bottom, 2=middle 3=top surface")
                end if

                ! then bring it to body system
                x_wing_b = matmul( transpose(Insects(1)%Wings(2)%M_b2w), x_wing_w ) + Insects(1)%Wings(2)%x_pivot_b
                ! and finnaly to global system
                x_wing_g = matmul( transpose(Insects(1)%M_g2b), x_wing_b ) + Insects(1)%xc_body_g
                ! save point on the wing now in global system
                wing_points_g(ipoint, 1:3, surface_type) = x_wing_g(1:3)
            enddo
        else if (wing_type == "left") then
            do ipoint = 1, npoints
                ! input is in wing system, midplane of the wing
                x_wing_w(1:3) = (/ wing_points_w(ipoint, 1:2), 0.0_rk /)

                x_wing_normal(1:3) = (/ 0.0_rk, 0.0_rk, 1.0_rk /)

                ! //NOTE the normal_vector direction is to outside the surface, should be reversed
                if (surface_type == 2) then                    
                    x_wing_w(1:3) = x_wing_w(1:3)
                else if (surface_type == 3) then
                    ! Wings: 1=left 2=right 3=left hind 4=right hind
                    x_wing_w(1:3) = x_wing_w(1:3) + 0.5_rk*Insects(1)%Wings(1)%WingThickness * x_wing_normal(1:3)
                else if (surface_type == 1) then
                    x_wing_w(1:3) = x_wing_w(1:3) - 0.5_rk*Insects(1)%Wings(1)%WingThickness * x_wing_normal(1:3)
                else
                    call abort(372936, "surface type must be 1=bottom, 2=middle 3=top surface")
                end if

                ! then bring it to body system
                x_wing_b = matmul( transpose(Insects(1)%Wings(1)%M_b2w), x_wing_w ) + Insects(1)%Wings(1)%x_pivot_b
                ! and finnaly to global system
                x_wing_g = matmul( transpose(Insects(1)%M_g2b), x_wing_b ) + Insects(1)%xc_body_g
                ! save point on the wing now in global system
                wing_points_g(ipoint, 1:3, surface_type) = x_wing_g(1:3)
            enddo
        else
            call abort(08122521, "Either right or left!")
        endif
    end do

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! interpolate pressure
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    allocate(xq( 1:3, 1:npoints*nsurfaces))
    allocate(pq( 1:1, 1:npoints*nsurfaces)) ! interpolate just one component (pressure)

    ! reshape in the format expected by interpolation (simple list of points)
    do isurface = 1, nsurfaces
        do ipoint = 1, npoints
            k = (isurface-1)*npoints + ipoint
            xq(:, k) = wing_points_g(ipoint, :, isurface)
        end do
    end do

    ! actual interpolation (in global system of course)
    call interpolatePointCloud_tree( params, hvy_block, tree_ID, xq, pq, "linear", sync=.true. )

    ! reshape result back in the format used here.
    pressure_data = reshape(pq, [npoints, nsurfaces])

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! output the data
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (params%rank == 0) then
        ! write original cell_id, face_id, qpoint and new interpolated pressure data to disk
        call get_command_argument(5,fname_out)
        open(14,file=fname_out, status='replace')
        write(14,*) "xw;yw;interpolated_value (bottom);interpolated_value (middle);interpolated_value (top)"

        do ipoint = 1, size(wing_points_w, 1)
            ! write qpoints and corresponding pressure
            write(14, '(ES13.6,";",ES13.6,";",ES13.6,";",ES13.6,";",ES13.6)') wing_points_w(ipoint,1:2), pressure_data(ipoint,1), pressure_data(ipoint,2), pressure_data(ipoint,3)
        end do

        close(14) 
    endif

    deallocate(wing_points_w, wing_points_g, pressure_data, xq, pq)
    call deallocate_forest(params, hvy_block)

end subroutine post_pressure_interpolation
