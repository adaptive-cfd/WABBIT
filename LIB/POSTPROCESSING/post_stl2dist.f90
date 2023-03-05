subroutine post_stl2dist(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_physics_metamodule
    use module_time_step
    use module_stl_file_reader
    use module_helpers
    use module_ini_files_parser_mpi
    use module_forestMetaData

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=cshort) :: fname_ini, fname_stl, fname_out, dummy
    integer :: i, Bs(1:3), g, ntri, k, iter, skips, a
    integer :: ix, iy, iz, ivertex, xmin, xmax, ymin, ymax, zmin, zmax, safety, mpicode
    real(kind=rk), dimension(1:3) :: vertex1, vertex2, vertex3, vertex1_normal, vertex2_normal, &
    vertex3_normal, face_normal, edge1_normal, edge2_normal, edge3_normal
    real(kind=rk) :: scale, origin(1:3), dist
    ! origin and spacing of blocks
    real(kind=rk) :: x0(1:3), dx(1:3), x,y,z,tmp

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: xyz_nxnynz(:, :)
    integer :: hvy_id, lgt_id, tree_ID=1
    integer :: c_plus, c_minus, res
    integer :: ix1,iy1,iz1
    logical :: done, array_compare_real, pruning

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(2, fname_ini)

    ! does the user need help?
    if (fname_ini=='--help' .or. fname_ini=='--h' .or. fname_ini=='-h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "./wabbit-post --stl2dist --x0 2.0,3.0,4.0 --superstl vertexnormals.sstl --params PARAMS.ini --no-pruning --memory=10.0gb"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) " Output file name:"
            write(*,*) " -o mask_00000.h5 "
            write(*,*) ""
            write(*,*) ""
            write(*,*) ""
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    ! defaults
    scale = 1.0_rk
    origin = 0.0_rk
    pruning = .true.


    ! fetch parameters from command line call
    do i = 1, COMMAND_ARGUMENT_COUNT()

        call get_command_argument(i,dummy)

        select case (dummy)
        case ("--params")
            call get_command_argument(i+1, fname_ini)
            call check_file_exists( fname_ini )

        case ("--superstl")
            call get_command_argument(i+1, fname_stl)
            call check_file_exists( fname_stl )

        case ("-o")
            call get_command_argument(i+1, fname_out)

        case ("--no-pruning")
            pruning = .false.

        case ("--scale")
            call get_command_argument(i+1, dummy)
            read(dummy,*) scale

        case ("--x0")
            call get_command_argument(i+1, dummy)
            read(dummy,*) origin

        end select
    enddo

    if (params%rank==0) then
        write(*,'("super-STL file is ",A)') trim(adjustl(fname_stl))
        write(*,'("INI file is ",A)') trim(adjustl(fname_ini))
        write(*,'("STL scaling factor is scale=",g12.3)') scale
        write(*,'("origin shift=",3(g12.3,1x))') origin
    endif

    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, fname_ini )
    ! have the pysics module read their own parameters
    call init_physics_modules( params, fname_ini, params%N_mask_components  )

    ! one field for the result, one field to tag error points where we have trouble determining the sign.
    params%n_eqn = 1
    N_MAX_COMPONENTS = params%n_eqn
    deallocate(params%threshold_state_vector_component)
    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component = .false.
    params%threshold_state_vector_component(1) = .true.
    ! params%g =

    ! in usual parameter files, RK4 (or some other RK) is used an requires a lot of memory
    ! here we do not need that, and hence pretent to use a basic scheme (EE1 maybe)
    deallocate(params%butcher_tableau)
    allocate(params%butcher_tableau(1,1))

    Bs = params%Bs
    g = params%g

    ! allocate data
    call allocate_forest(params, hvy_block)

    call reset_tree( params, .true., tree_ID)

    ! start with an equidistant grid on coarsest level
    call createEquidistantGrid_tree( params, params%Jmin, .true., tree_ID)

    ! reset grid to zeros
    do k = 1, hvy_n(1)
        hvy_block(:,:,:,:,hvy_active(k,1)) = 0.0_rk
    enddo


    call count_lines_in_ascii_file_mpi(fname_stl, ntri, 0)
    allocate( xyz_nxnynz(1:ntri,1:30) )

    call read_array_from_ascii_file_mpi(fname_stl, xyz_nxnynz, 0)
    xyz_nxnynz(:,1:9) = xyz_nxnynz(:,1:9) / scale
    xyz_nxnynz(:,1) = xyz_nxnynz(:,1) +  origin(1)
    xyz_nxnynz(:,4) = xyz_nxnynz(:,4) +  origin(1)
    xyz_nxnynz(:,7) = xyz_nxnynz(:,7) +  origin(1)

    xyz_nxnynz(:,2) = xyz_nxnynz(:,2) +  origin(2)
    xyz_nxnynz(:,5) = xyz_nxnynz(:,5) +  origin(2)
    xyz_nxnynz(:,8) = xyz_nxnynz(:,8) +  origin(2)

    xyz_nxnynz(:,3) = xyz_nxnynz(:,3) +  origin(3)
    xyz_nxnynz(:,6) = xyz_nxnynz(:,6) +  origin(3)
    xyz_nxnynz(:,9) = xyz_nxnynz(:,9) +  origin(3)


    safety = 6 !Bs(1)
    do iter = params%Jmin, params%Jmax

        ! refine the mesh where the mask function is interesting
        call abort(99999, "need to adapt refine_tree call to include hvy_tmp")
        ! call refine_tree( params, hvy_block, "mask-threshold", tree_ID )

        skips = 0

        do k = 1, hvy_n(1)
            ! hvy_id of the block we're looking at
            hvy_id = hvy_active(k,1)


            ! The trick here is that the STL does not move, and that we set it on the ghost
            ! nodes as well. Then, a refined block should always end up with nonzero values
            ! only if its mother had nonzero values.
            ! That implies: if the block is zero, we can skip it
            if ( maxval(hvy_block(:,:,:,1,hvy_id)) <= 1.0e-6 .and. iter>params%Jmin ) then
                hvy_block(:,:,:,1,hvy_id) = 9e8_rk ! distance as far away
                skips = skips + 1
                cycle
            endif

            hvy_block(:,:,:,1,hvy_id) = 9e8_rk ! distance as far away


            ! compute block spacing and origin from treecode
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            call get_block_spacing_origin( params, lgt_id, x0, dx )

            ! shift origin to take ghost nodes into account
            x0 = x0 - dble(g)*dx

            do i = 1, ntri
                vertex1        = xyz_nxnynz(i, 1:3)
                vertex2        = xyz_nxnynz(i, 4:6)
                vertex3        = xyz_nxnynz(i, 7:9)
                face_normal    = xyz_nxnynz(i, 10:12)
                vertex1_normal = xyz_nxnynz(i, 13:15)
                vertex2_normal = xyz_nxnynz(i, 16:18)
                vertex3_normal = xyz_nxnynz(i, 19:21)
                edge1_normal   = xyz_nxnynz(i, 22:24)
                edge2_normal   = xyz_nxnynz(i, 25:27)
                edge3_normal   = xyz_nxnynz(i, 28:30)

                xmin = floor( ( minval((/xyz_nxnynz(i,1),xyz_nxnynz(i,4),xyz_nxnynz(i,7)/)) - x0(1) ) / dx(1)) - safety
                ymin = floor( ( minval((/xyz_nxnynz(i,2),xyz_nxnynz(i,5),xyz_nxnynz(i,8)/)) - x0(2) ) / dx(2)) - safety
                zmin = floor( ( minval((/xyz_nxnynz(i,3),xyz_nxnynz(i,6),xyz_nxnynz(i,9)/)) - x0(3) ) / dx(3)) - safety

                xmax = ceiling( ( maxval((/xyz_nxnynz(i,1),xyz_nxnynz(i,4),xyz_nxnynz(i,7)/)) - x0(1) ) / dx(1)) + safety
                ymax = ceiling( ( maxval((/xyz_nxnynz(i,2),xyz_nxnynz(i,5),xyz_nxnynz(i,8)/)) - x0(2) ) / dx(2)) + safety
                zmax = ceiling( ( maxval((/xyz_nxnynz(i,3),xyz_nxnynz(i,6),xyz_nxnynz(i,9)/)) - x0(3) ) / dx(3)) + safety

                xmin = max(xmin, 1)
                ymin = max(ymin, 1)
                zmin = max(zmin, 1)

                xmax = min(xmax, Bs(1)+2*g)
                ymax = min(ymax, Bs(2)+2*g)
                zmax = min(zmax, Bs(3)+2*g)

                do iz = zmin, zmax
                    z = dx(3)*dble(iz) + x0(3)
                    do iy = ymin, ymax
                        y = dx(2)*dble(iy) + x0(2)
                        do ix = xmin, xmax
                            x = dx(1)*dble(ix) + x0(1)

                            ! the distance to the current triangle:
                            tmp = pointTriangleDistance( &
                            vertex1, &
                            vertex2, &
                            vertex3, &
                            (/x,y,z/), & ! query point
                            face_normal, &
                            vertex1_normal, &
                            vertex2_normal, &
                            vertex3_normal, &
                            edge1_normal, &
                            edge2_normal, &
                            edge3_normal)


                            ! if closer (in abs value!) then use this now
                            if ( abs(tmp) < abs(hvy_block(ix,iy,iz,1,hvy_id)) .and. abs(tmp)<=dble(safety)*dx(1) ) then
                                hvy_block(ix,iy,iz,1,hvy_id) = (tmp)
                            endif
                        enddo
                    enddo
                enddo

            enddo ! loop over triangles
        enddo ! loop blocks

        !=======================================================================
        ! convert block to mask function
        !=======================================================================
        dx(1) = 2.0_rk**(-iter) * (params%domain_size(1) / real(Bs(1)-1, kind=rk))

        do k = 1, hvy_n(1)
            hvy_id = hvy_active(k,1)
            do iz = 1, Bs(3)+2*g
                do iy = 1, Bs(2)+2*g
                    do ix = 1, Bs(1)+2*g
                        hvy_block(ix,iy,iz,1,hvy_id) = smoothstep( hvy_block(ix,iy,iz,1,hvy_id), 0.0_rk, 1.5_rk*dx(1) )
                    enddo
                enddo
            enddo
        enddo ! loop over blocks

        write(*,'("rank=",i4," skipped ",i6," of its ",i6," blocks")') params%rank, skips, hvy_n(1)

        call MPI_barrier( WABBIT_COMM, mpicode)


        if (params%rank==0) then
            write(*, '("Nb=",i6," Jmin=",i2," Jmax=",i2)') &
            lgt_n(1), minActiveLevel_tree(1), &
            maxActiveLevel_tree(1)
        endif

    enddo ! loop over level

    !=======================================================================
    ! coarsening of blocks with constant values
    !=======================================================================
    call adapt_tree( 0.0_rk, params, hvy_block, 1, params%coarsening_indicator, hvy_block )


    call createActiveSortedLists_forest(params)

    if (pruning) then
        if (params%rank==0) write(*,*) "now pruning!"

        call prune_tree( params, hvy_block, tree_ID=1)
    endif

    call saveHDF5_tree(fname_out, 0.0_rk, -1_ik, 1, params, hvy_block, 1)
end subroutine
