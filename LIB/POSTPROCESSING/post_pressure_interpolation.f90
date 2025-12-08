subroutine post_pressure_interpolation

    use mpi
    use module_helpers
    use module_MPI
    use module_params               ! global parameters
    use module_timing
    use module_mesh                 ! mesh manipulation subroutines
    use module_time_step
    use module_unit_test
    use module_bridge_interface     ! bridge implementation of wabbit
    ! HACK.We should load only the metamodule, but we require WRITE_INSECT_DATA(time)
    ! to dump kinematics data.
    use module_ACM
    use module_forestMetaData
    use module_insects

    implicit none

    integer(kind=ik)                    :: number_procs     ! number of processes
    ! real(kind=rk)                       :: t0, t1, t2       ! cpu time variables for running time calculation
    type (type_params)                  :: params           ! user defined parameter structure

    real(kind=rk), allocatable          :: hvy_mask(:, :, :, :, :), hvy_tmp(:, :, :, :, :), hvy_block(:, :, :, :, :)
    real(kind=rk)                       :: time, xx, yy, zz, delx, dely, delz, tmp           ! time loop variables
    character(len=cshort)               :: pressure_filename, ini_filename, wing_fname, fname_out
    integer(kind=ik)                    :: k, lgt_id, Bs(1:3), g, hvy_id, iter, iteration, tree_ID=1, N_support, nlines, ncols
    real(kind=rk)                       :: x(1:3), x0(1:3), dx(1:3), x_wing_w(1:3), x_wing_b(1:3), x_wing_g(1:3), x_wing_normal(1:3), coeff(1:3)
    real(kind=rk)                       :: block_x_min(1:3), block_x_max(1:3)
    integer(kind=ik)                    :: ipoint, ix, iy, iz, ix0, iy0, iz0
    logical                             :: help1, help2
    integer(kind=ik)                    :: nz, dim, tc_length
    real(kind=rk), dimension(3)         :: domain
    type(inifile) :: FILE
    character(len=cshort)               :: wing_type
    integer(kind=ik)                    :: surface_type    ! 1=bottom, 2=middle 3=top surface
    real(kind=rk), allocatable          :: wing_points_w(:,:), pressure_data(:,:), wing_points_g(:,:,:)
    integer(kind=ik)                    :: mpierr, n_header, isurface
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
            write(*, '(A)') " The input file is a SPACE-separated file with two columns, xw and yw, with one header line."
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

    ! ! modifications to parameters (because we use hvy_block instead of hvy_mask, NEQN set
    ! ! in ini file is not correct)
    ! deallocate( params%butcher_tableau )
    ! allocate( params%butcher_tableau(1,1) )
    ! ! mask, usx,usy,usz, color, sponge = 6 components
    ! params%n_eqn = 6
    ! deallocate(params%threshold_state_vector_component)
    ! allocate(params%threshold_state_vector_component(1:params%n_eqn))
    ! params%threshold_state_vector_component = 0
    ! params%threshold_state_vector_component(1) = 1

    ! deallocate(params%symmetry_vector_component)
    ! allocate(params%symmetry_vector_component(1:params%n_eqn))
    ! params%symmetry_vector_component = "0"


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

    call read_attributes(pressure_filename, iteration, time, iteration, domain, Bs, tc_length, dim, &
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

    ! update kinematics at the right time
    call Update_Insect(time, Insect)

    ! BEFORE WE CAN INTERPOLATE THE GHOTS NODES NEED TO BE FILLED
    call sync_ghosts_tree( params, hvy_block, tree_ID=tree_ID )


    ! read in the wing points
    ! in the wing system
    n_header = 1
    call count_lines_in_ascii_file_mpi(wing_fname, nlines, n_header)
    call count_cols_in_ascii_file_mpi(wing_fname, ncols, n_header)

    if (ncols .ne. 2) then
        call abort(202512081, "Input file should just contain list of points on the wing midplance (xw,yw), zw==0")
    endif

    allocate(wing_points_w( 1:nlines, 1:3) )
    allocate(wing_points_g( 1:nlines, 1:3, 1:3) )
    allocate(pressure_data( 1:nlines, 1:3) ) ! result for all 3 surfaces

    ! wing_points_w: xw, yw
    call read_array_from_ascii_file_mpi(wing_fname, wing_points_w(:,1:2), n_header)

    ! always interpolate all points
    do surface_type = 1, 3
        ! transform point data to global system 
        if (wing_type == "right") then
            do ipoint = 1, size(wing_points_w, 1)
                ! input is in wing system, midplane of the wing
                x_wing_w(1:3) = (/ wing_points_w(ipoint, 1:2), 0.0_rk /)
                x_wing_normal(1:3) = (/ 0.0_rk, 0.0_rk, 1.0_rk /)
                
                ! //NOTE the normal_vector direction is to outside the surface, should be reversed
                if (surface_type == 2) then                    
                    x_wing_w(1:3) = x_wing_w(1:3)
                else if (surface_type == 3) then
                    x_wing_w(1:3) = x_wing_w(1:3) - 0.5_rk*Insect%WingThickness * x_wing_normal(1:3)
                else if (surface_type == 1) then
                    x_wing_w(1:3) = x_wing_w(1:3) + 0.5_rk*Insect%WingThickness * x_wing_normal(1:3)
                else
                    call abort(372936, "surface type must be 1=bottom, 2=middle 3=top surface")
                end if

                ! then bring it to body system
                x_wing_b = matmul( transpose(Insect%M_b2w_r), x_wing_w ) + Insect%x_pivot_r_b
                ! and finnaly to global system
                x_wing_g = matmul( transpose(Insect%M_g2b), x_wing_b ) + Insect%xc_body_g
                ! save point on the wing now in global system
                wing_points_g(ipoint, 1:3, surface_type) = x_wing_g(1:3)
            enddo
        else if (wing_type == "left") then
            do ipoint = 1, size(wing_points_w, 1)
                ! input is in wing system, midplane of the wing
                x_wing_w(1:3) = (/ wing_points_w(ipoint, 1:2), 0.0_rk /)

                x_wing_normal(1:3) = (/ 0.0_rk, 0.0_rk, 1.0_rk /)

                ! //NOTE the normal_vector direction is to outside the surface, should be reversed
                if (surface_type == 2) then                    
                    x_wing_w(1:3) = x_wing_w(1:3)
                else if (surface_type == 3) then
                    x_wing_w(1:3) = x_wing_w(1:3) + 0.5_rk*Insect%WingThickness * x_wing_normal(1:3)
                else if (surface_type == 1) then
                    x_wing_w(1:3) = x_wing_w(1:3) - 0.5_rk*Insect%WingThickness * x_wing_normal(1:3)
                else
                    call abort(372936, "surface type must be 1=bottom, 2=middle 3=top surface")
                end if

                ! then bring it to body system
                x_wing_b = matmul( transpose(Insect%M_b2w_l), x_wing_w ) + Insect%x_pivot_l_b
                ! and finnaly to global system
                x_wing_g = matmul( transpose(Insect%M_g2b), x_wing_b ) + Insect%xc_body_g
                ! save point on the wing now in global system
                wing_points_g(ipoint, 1:3, surface_type) = x_wing_g(1:3)
            enddo
        else
            call abort(08122521, "Either right or left!")
        endif
    end do



    N_support = 3

    if ( params%g < N_support ) then
        call abort(1809251,"Error: not enough ghostpoints for delta interpolation. Increase number_ghost_nodes in PARAMS file!")
    endif

    ! set the array to a very large number. Why? If a point is NOT interpolated by a CPU, then the CPU
    ! will not touch the data. So it remains that large negative number. IN other words, the array looks like this:
    ! CPU1 = (/12.2, 13.4, -9e9, -9e9/)
    ! CPU2 = (/-9e9, -9e9, 7.2,  -39.9/)
    ! now I can just take the MAXIMUM of all values across all CPU
    ! and the final result is:
    ! (/12.2, 13.4, 7.2, -39.9 )
    pressure_data = -9.9e9_rk

    ! now we can interpolate each point
    do ipoint = 1, size(wing_points_w, 1)
        do isurface = 1, 3
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k,tree_ID)
                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, x0, dx )

                block_x_min(1:3) = x0(1:3)
                block_x_max(1:3) = x0(1:3) + real(Bs(1:3), rk) * dx(1:3)

                x_wing_g(1:3) =  wing_points_g(ipoint, 1:3, isurface)

                
                ! !//COMMENT for test, trilinear interpolation
                ! tmp = trilinear_interpolation(x0, dx, hvy_block( params%g+1: bs(1)+params%g+1,  params%g+1: bs(2)+params%g+1,  params%g+1: bs(3)+params%g+1, 1, hvy_id ), x_wing_g, .false.)  
                ! if (tmp /= -9.9e9_rk) then
                !     pressure_data(ipoint, jpoint) = tmp
                ! endif

                ! //COMMENT delta interpolation
                ! not all blocks are relevant: only one single block contains the interpolation 
                ! point we are looking at. Find the block! 
                if ( x_wing_g(1) >= block_x_min(1) .and. x_wing_g(2) >= block_x_min(2) .and. x_wing_g(3) >= block_x_min(3) .and. x_wing_g(1) < block_x_max(1) .and. x_wing_g(2) < block_x_max(2) .and. x_wing_g(3) < block_x_max(3)) then               

                    ! convert interpolation point to integer, nearest integer
                    x = x_wing_g(1:3) - x0(1:3)
                    ix0 = floor( x(1) / dx(1)) + (params%g + 1)
                    iy0 = floor( x(2) / dx(2)) + (params%g + 1)
                    iz0 = floor( x(3) / dx(3)) + (params%g + 1)

                    pressure_data(ipoint, isurface) = 0.0_rk
                    do iz = iz0-N_support,iz0+N_support ! the box size around the point
                        zz = real(iz - (params%g + 1), rk) * dx(3)
                        delz = delta_interpolation(abs(zz - x(3)),dx(3))

                        do iy = iy0-N_support,iy0+N_support
                            yy = real(iy - (params%g + 1), rk) * dx(2)
                            dely = delta_interpolation(abs(yy - x(2)),dx(2))

                            do ix = ix0-N_support,ix0+N_support
                                xx = real(ix - (params%g + 1), rk) * dx(1)
                                delx = delta_interpolation(abs(xx - x(1)),dx(1))

                                ! tmp = hvy_block( ix, iy, iz, 1, hvy_id )
                                pressure_data(ipoint, isurface) = pressure_data(ipoint, isurface) + delx * dely * delz * hvy_block( ix, iy, iz, 1, hvy_id )
                            enddo
                        enddo
                    enddo   
                endif
            end do
        end do
    end do


    call MPI_allreduce( MPI_IN_PLACE, pressure_data, size(pressure_data), MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

    ! output the data

    ! call summarize_profiling(WABBIT_COMM)
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

    deallocate(wing_points_w, wing_points_g, pressure_data)
    call deallocate_forest(params, hvy_block)

end subroutine post_pressure_interpolation
