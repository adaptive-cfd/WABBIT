subroutine post_stl2dist(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_mpi
    use module_operators
    use module_physics_metamodule
    use module_time_step
    use module_stl_file_reader
    use module_helpers
    use module_ini_files_parser_mpi

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=80) :: fname_ini, fname_stl, fname_out, dummy
    integer :: i, Bs(1:3), g, ntri, k, iter, skips, a
    integer :: ix, iy, iz, ivertex, xmin, xmax, ymin, ymax, zmin, zmax, safety, mpicode
    real(kind=4), dimension(1:3) :: vertex1, vertex2, vertex3, vertex1_normal, vertex2_normal, &
    vertex3_normal, face_normal, edge1_normal, edge2_normal, edge3_normal
    real(kind=rk) :: scale, origin(1:3), dist
    ! origin and spacing of blocks
    real(kind=rk) :: x0(1:3), dx(1:3), x,y,z,tmp

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: xyz_nxnynz(:, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:), hvy_n(:), lgt_n(:)
    integer(kind=ik), allocatable      :: lgt_active(:,:), hvy_active(:,:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:,:)
    integer :: hvy_id, lgt_id
    integer :: c_plus, c_minus,res, tree_n
    integer :: ix1,iy1,iz1
    logical :: done, array_compare_real

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(2, fname_ini)

    ! does the user need help?
    if (fname_ini=='--help' .or. fname_ini=='--h' .or. fname_ini=='-h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "./wabbit-post --stl2dist --x0 2.0,3.0,4.0 --superstl vertexnormals.sstl --params PARAMS.ini"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) ""
            write(*,*) ""
            write(*,*) ""
            write(*,*) ""
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    ! defaults
    scale = 1.0_rk
    N_MAX_COMPONENTS = 1

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
    ! params%n_ghosts =

    ! in usual parameter files, RK4 (or some other RK) is used an requires a lot of memory
    ! here we do not need that, and hence pretent to use a basic scheme (EE1 maybe)
    deallocate(params%butcher_tableau)
    allocate(params%butcher_tableau(1,1))

    Bs = params%Bs
    g = params%n_ghosts

    ! allocate data
    call allocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    hvy_active, lgt_sortednumlist, hvy_n=hvy_n, lgt_n=lgt_n)

    call reset_tree( params, lgt_block, lgt_active(:,1), &
    lgt_n(1), hvy_active(:,1), hvy_n(1), lgt_sortednumlist(:,:,1), .true., tree_ID=1)

    ! start with an equidistant grid on coarsest level
    call create_equidistant_grid( params, lgt_block, hvy_neighbor, lgt_active(:,1), lgt_n(1), &
    lgt_sortednumlist(:,:,1), hvy_active(:,1), hvy_n(1), params%min_treelevel, .true., tree_ID=1 )

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


    safety = 6
    do iter = params%min_treelevel, params%max_treelevel

        ! refine the mesh where the mask function is interesting
        call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,1), lgt_n(1), &
        lgt_sortednumlist(:,:,1), hvy_active(:,1), hvy_n(1), "mask-threshold", tree_ID=1 )

        skips = 0

        do k = 1, hvy_n(1)
            ! hvy_id of the block we're looking at
            hvy_id = hvy_active(k,1)


            ! The trick here is that the STL does not move, and that we set it on the ghost
            ! nodes as well. Then, a refined block should always end up with nonzero values
            ! only if its mother had nonzero values.
            ! That implies: if the block is zero, we can skip it
            if ( maxval(hvy_block(:,:,:,1,hvy_id)) <= 1.0e-6 .and. iter>params%min_treelevel ) then
                hvy_block(:,:,:,1,hvy_id) = 9e8_rk ! distance as far away
                skips = skips + 1
                cycle
            endif

            hvy_block(:,:,:,1,hvy_id) = 9e8_rk ! distance as far away


            ! compute block spacing and origin from treecode
            call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
            call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

            ! shift origin to take ghost nodes into account
            x0 = x0 - dble(g)*dx

            do i = 1, ntri
                vertex1        = real( xyz_nxnynz(i, 1:3)  , kind=4)
                vertex2        = real( xyz_nxnynz(i, 4:6)  , kind=4)
                vertex3        = real( xyz_nxnynz(i, 7:9)  , kind=4)
                face_normal    = real( xyz_nxnynz(i, 10:12), kind=4)
                vertex1_normal = real( xyz_nxnynz(i, 13:15), kind=4)
                vertex2_normal = real( xyz_nxnynz(i, 16:18), kind=4)
                vertex3_normal = real( xyz_nxnynz(i, 19:21), kind=4)
                edge1_normal   = real( xyz_nxnynz(i, 22:24), kind=4)
                edge2_normal   = real( xyz_nxnynz(i, 25:27), kind=4)
                edge3_normal   = real( xyz_nxnynz(i, 28:30), kind=4)

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
                            tmp = pointTriangleDistance3( &
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
                            if ( abs(tmp) < abs(hvy_block(ix,iy,iz,1,hvy_id)) ) then
                                hvy_block(ix,iy,iz,1,hvy_id) = tmp
                            endif
                        enddo
                    enddo
                enddo

            enddo ! loop over triangles
        enddo ! loop blocks

        !=======================================================================
        ! convert block to mask function
        !=======================================================================
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
            lgt_n(1), min_active_level( lgt_block, lgt_active(:,1), lgt_n(1) ), &
            max_active_level( lgt_block, lgt_active(:,1), lgt_n(1) )
        endif
    enddo ! loop over level

    !=======================================================================
    ! coarsening of blocks with constant values
    !=======================================================================
    call adapt_mesh( 0.0_rk, params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,1), &
    lgt_n(1), lgt_sortednumlist(:,:,1), hvy_active(:,1), &
    hvy_n(1), 1, params%coarsening_indicator, hvy_block )


    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)

    if (params%rank==0) write(*,*) "now pruning!"
    call prune_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id=1)

    call write_tree_field(fname_out, params, lgt_block, lgt_active, hvy_block, &
    lgt_n, hvy_n, hvy_active, dF=1, tree_id=1, time=0.0_rk, iteration=-1 )
end subroutine

! subroutine post_stl2dist(params)
!     use module_precision
!     use module_mesh
!     use module_params
!     use module_IO
!     use module_mpi
!     use module_operators
!     use module_physics_metamodule
!     use module_time_step
!     use module_stl_file_reader
!     use module_helpers
!     use module_ini_files_parser_mpi
!
!     implicit none
!
!     type (type_params), intent(inout)  :: params
!     character(len=80) :: fname_ini, fname_stl, fname_out, dummy, fname_xyz
!     integer :: i, Bs(1:3), g, ntri, k, iter, skips, a
!     integer :: ix, iy, iz, ivertex, xmin, xmax, ymin, ymax, zmin, zmax, safety, mpicode
!     real(kind=4), allocatable, dimension(:,:) :: triangles, normals
!     real(kind=rk) :: scale, origin(1:3), dist
!     ! origin and spacing of blocks
!     real(kind=rk) :: x0(1:3), dx(1:3), x,y,z,tmp
!     real(kind=4), allocatable :: n1(:,:), n2(:,:), n3(:,:)
!
!     integer(kind=ik), allocatable      :: lgt_block(:, :)
!     real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
!     real(kind=rk), allocatable         :: xyz_nxnynz(:, :)
!     integer(kind=ik), allocatable      :: hvy_neighbor(:,:), hvy_n(:), lgt_n(:)
!     integer(kind=ik), allocatable      :: lgt_active(:,:), hvy_active(:,:)
!     integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:,:)
!     integer :: hvy_id, lgt_id
!     integer :: c_plus, c_minus,res, tree_n
!     integer :: ix1,iy1,iz1
!     logical :: done, array_compare_real
!
!     !-----------------------------------------------------------------------------------------------------
!     ! get values from command line (filename and level for interpolation)
!     call get_command_argument(2, fname_ini)
!
!     ! does the user need help?
!     if (fname_ini=='--help' .or. fname_ini=='--h' .or. fname_ini=='-h') then
!         if (params%rank==0) then
!             write(*,*) "------------------------------------------------------------------"
!             write(*,*) "./wabbit-post --stl2dist --x0 2.0,3.0,4.0 --stl some.stl --xyz vertexnormals.xyz --params PARAMS.ini"
!             write(*,*) "./wabbit-post --stl2dist --x0 2.0,3.0,4.0 --superstl vertexnormals.sstl --params PARAMS.ini"
!             write(*,*) "------------------------------------------------------------------"
!             write(*,*) ""
!             write(*,*) ""
!             write(*,*) ""
!             write(*,*) ""
!             write(*,*) "------------------------------------------------------------------"
!         end if
!         return
!     endif
!
!     ! defaults
!     scale = 1.0_rk
!     N_MAX_COMPONENTS = 1
!
!     ! fetch parameters from command line call
!     do i = 1, COMMAND_ARGUMENT_COUNT()
!
!         call get_command_argument(i,dummy)
!
!         select case (dummy)
!         case ("--params")
!             call get_command_argument(i+1, fname_ini)
!             call check_file_exists( fname_ini )
!
!         case ("--stl")
!             call get_command_argument(i+1, fname_stl)
!             call check_file_exists( fname_stl )
!
!         case ("--xyz")
!             call get_command_argument(i+1, fname_xyz)
!             call check_file_exists( fname_xyz )
!
!         case ("-o")
!             call get_command_argument(i+1, fname_out)
!
!         case ("--scale")
!             call get_command_argument(i+1, dummy)
!             read(dummy,*) scale
!
!         case ("--x0")
!             call get_command_argument(i+1, dummy)
!             read(dummy,*) origin
!
!         end select
!     enddo
!
!     if (params%rank==0) then
!         write(*,'("STL file is ",A)') trim(adjustl(fname_stl))
!         write(*,'("INI file is ",A)') trim(adjustl(fname_ini))
!         write(*,'("STL scaling factor is scale=",g12.3)') scale
!         write(*,'("origin shift=",3(g12.3,1x))') origin
!     endif
!
!     ! read ini-file and save parameters in struct
!     call ini_file_to_params( params, fname_ini )
!     ! have the pysics module read their own parameters
!     call init_physics_modules( params, fname_ini, params%N_mask_components  )
!
!     ! one field for the result, one field to tag error points where we have trouble determining the sign.
!     params%n_eqn = 1
!     ! params%n_ghosts =
!
!     ! in usual parameter files, RK4 (or some other RK) is used an requires a lot of memory
!     ! here we do not need that, and hence pretent to use a basic scheme (EE1 maybe)
!     deallocate(params%butcher_tableau)
!     allocate(params%butcher_tableau(1,1))
!
!     Bs = params%Bs
!     g = params%n_ghosts
!
!     ! allocate data
!     call allocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
!     hvy_active, lgt_sortednumlist, hvy_n=hvy_n, lgt_n=lgt_n)
!
!     call reset_tree( params, lgt_block, lgt_active(:,1), &
!     lgt_n(1), hvy_active(:,1), hvy_n(1), lgt_sortednumlist(:,:,1), .true., tree_ID=1)
!
!     ! start with an equidistant grid on coarsest level
!     call create_equidistant_grid( params, lgt_block, hvy_neighbor, lgt_active(:,1), lgt_n(1), &
!     lgt_sortednumlist(:,:,1), hvy_active(:,1), hvy_n(1), params%min_treelevel, .true., tree_ID=1 )
!
!     ! reset grid to zeros
!     do k = 1, hvy_n(1)
!         hvy_block(:,:,:,:,hvy_active(k,1)) = 0.0_rk
!     enddo
!
!
!     call read_stl_file(fname_stl, ntri, triangles, normals)
!     call normalize_stl_file( ntri, triangles, "--scale", scale, origin )
!
!     call count_lines_in_ascii_file_mpi(fname_xyz, i, 0)
!     allocate( xyz_nxnynz(1:i,1:6) )
!
!     call read_array_from_ascii_file_mpi(fname_xyz, xyz_nxnynz, 0)
!     xyz_nxnynz(:,1:3) = xyz_nxnynz(:,1:3)/scale
!     xyz_nxnynz(:,1) = xyz_nxnynz(:,1)+origin(1)
!     xyz_nxnynz(:,2) = xyz_nxnynz(:,2)+origin(2)
!     xyz_nxnynz(:,3) = xyz_nxnynz(:,3)+origin(3)
!
!     allocate( n1(1:ntri,1:3), n2(1:ntri,1:3), n3(1:ntri,1:3) )
!
!     do i = 1, ntri
!         ivertex = 3*i - 2
!
!         dist = 9e9_rk
!         res = -1
!         do a = 1, size(xyz_nxnynz, 1)
!             if ( norm2(triangles(1:3,ivertex)-xyz_nxnynz(a,1:3)) < dist ) then
!                 res = a
!                 dist = norm2(triangles(1:3,ivertex)-xyz_nxnynz(a,1:3))
!             endif
!         enddo
!         n1(i,1:3) = real(xyz_nxnynz(res,4:6), kind=4)
!         if (dist>1.0e-5_rk) call abort(25010192, "STL/XYZ: vertex not found in both files. Do *.stl and *.xyz describe the same data?")
!
!
!         dist = 9e9_rk
!         res = -1
!         do a = 1, size(xyz_nxnynz, 1)
!             if ( norm2(triangles(1:3,ivertex+1)-xyz_nxnynz(a,1:3)) < dist ) then
!                 res = a
!                 dist = norm2(triangles(1:3,ivertex+1)-xyz_nxnynz(a,1:3))
!             endif
!         enddo
!         n2(i,1:3) = real(xyz_nxnynz(res,4:6), kind=4)
!         if (dist>1.0e-5_rk) call abort(25010192, "STL/XYZ: vertex not found in both files. Do *.stl and *.xyz describe the same data?")
!
!
!         dist = 9e9_rk
!         res = -1
!         do a = 1, size(xyz_nxnynz, 1)
!             if ( norm2(triangles(1:3,ivertex+2)-xyz_nxnynz(a,1:3)) < dist ) then
!                 res = a
!                 dist = norm2(triangles(1:3,ivertex+2)-xyz_nxnynz(a,1:3))
!             endif
!         enddo
!         n3(i,1:3) = real(xyz_nxnynz(res,4:6), kind=4)
!         if (dist>1.0e-5_rk) call abort(25010192, "STL/XYZ: vertex not found in both files. Do *.stl and *.xyz describe the same data?")
!     enddo
!
!     deallocate( xyz_nxnynz )
!
!     safety = 6
!     do iter = params%min_treelevel, params%max_treelevel
!
!         ! refine the mesh where the mask function is interesting
!         call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,1), lgt_n(1), &
!         lgt_sortednumlist(:,:,1), hvy_active(:,1), hvy_n(1), "mask-threshold", tree_ID=1 )
!
!         skips = 0
!
!         do k = 1, hvy_n(1)
!             ! hvy_id of the block we're looking at
!             hvy_id = hvy_active(k,1)
!
!
!             ! The trick here is that the STL does not move, and that we set it on the ghost
!             ! nodes as well. Then, a refined block should always end up with nonzero values
!             ! only if its mother had nonzero values.
!             ! That implies: if the block is zero, we can skip it
!             if ( maxval(hvy_block(:,:,:,1,hvy_id)) <= 1.0e-6 .and. iter>params%min_treelevel ) then
!                 hvy_block(:,:,:,1,hvy_id) = 9e8_rk ! distance as far away
!                 skips = skips + 1
!                 cycle
!             endif
!
!             hvy_block(:,:,:,1,hvy_id) = 9e8_rk ! distance as far away
!
!             ! compute block spacing and origin from treecode
!             call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
!             call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
!
!             ! shift origin to take ghost nodes into account
!             x0 = x0 - dble(g)*dx
!
!             do i = 1, ntri
!                 ivertex = 3*i - 2
!
!                 xmin = floor( minval( (triangles(1,ivertex:ivertex+2)-x0(1)) ) / dx(1)) - safety
!                 ymin = floor( minval( (triangles(2,ivertex:ivertex+2)-x0(2)) ) / dx(2)) - safety
!                 zmin = floor( minval( (triangles(3,ivertex:ivertex+2)-x0(3)) ) / dx(3)) - safety
!
!                 xmax = ceiling( maxval( (triangles(1,ivertex:ivertex+2)-x0(1)) ) / dx(1)) + safety
!                 ymax = ceiling( maxval( (triangles(2,ivertex:ivertex+2)-x0(2)) ) / dx(2)) + safety
!                 zmax = ceiling( maxval( (triangles(3,ivertex:ivertex+2)-x0(3)) ) / dx(3)) + safety
!
!                 xmin = max(xmin, 1)
!                 ymin = max(ymin, 1)
!                 zmin = max(zmin, 1)
!
!                 xmax = min(xmax, Bs(1)+2*g)
!                 ymax = min(ymax, Bs(2)+2*g)
!                 zmax = min(zmax, Bs(3)+2*g)
!
!                 do iz = zmin, zmax
!                     z = dx(3)*dble(iz) + x0(3)
!                     do iy = ymin, ymax
!                         y = dx(2)*dble(iy) + x0(2)
!                         do ix = xmin, xmax
!                             x = dx(1)*dble(ix) + x0(1)
!
!                             ! the distance to the current triangle:
!                             tmp = pointTriangleDistance2( triangles(:,ivertex), triangles(:,ivertex+1), &
!                             triangles(:,ivertex+2), (/x,y,z/), normals(:,i), n1(i,1:3), n2(i,1:3), n3(i,1:3) )
!
!                             ! if closer (in abs value!) then use this now
!                             if ( abs(tmp) < abs(hvy_block(ix,iy,iz,1,hvy_id)) ) then
!                                 hvy_block(ix,iy,iz,1,hvy_id) = tmp
!                             endif
!                         enddo
!                     enddo
!                 enddo
!
!             enddo ! loop over triangles
!         enddo ! loop blocks
!
!         !=======================================================================
!         ! convert block to mask function
!         !=======================================================================
!         do k = 1, hvy_n(1)
!             hvy_id = hvy_active(k,1)
!             do iz = 1, Bs(3)+2*g
!                 do iy = 1, Bs(2)+2*g
!                     do ix = 1, Bs(1)+2*g
!                         hvy_block(ix,iy,iz,1,hvy_id) = smoothstep( hvy_block(ix,iy,iz,1,hvy_id), 0.0_rk, 1.5_rk*dx(1) )
!                     enddo
!                 enddo
!             enddo
!         enddo ! loop over blocks
!
!         write(*,'("rank=",i4," skipped ",i6," of its ",i6," blocks")') params%rank, skips, hvy_n(1)
!
!         call MPI_barrier( WABBIT_COMM, mpicode)
!
!
!         if (params%rank==0) then
!             write(*, '("Nb=",i6," Jmin=",i2," Jmax=",i2)') &
!             lgt_n(1), min_active_level( lgt_block, lgt_active(:,1), lgt_n(1) ), &
!             max_active_level( lgt_block, lgt_active(:,1), lgt_n(1) )
!         endif
!     enddo ! loop over level
!
!     !=======================================================================
!     ! coarsening of blocks with constant values
!     !=======================================================================
!     call adapt_mesh( 0.0_rk, params, lgt_block, hvy_block, hvy_neighbor, lgt_active(:,1), &
!     lgt_n(1), lgt_sortednumlist(:,:,1), hvy_active(:,1), &
!     hvy_n(1), 1, params%coarsening_indicator, hvy_block )
!
!
!     call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
!     lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)
!
!     if (params%rank==0) write(*,*) "now pruning!"
!     call prune_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
!     hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id=1)
!
!     call write_tree_field(fname_out, params, lgt_block, lgt_active, hvy_block, &
!     lgt_n, hvy_n, hvy_active, dF=1, tree_id=1, time=0.0_rk, iteration=-1 )
! end subroutine


! subroutine post_stl2dist(params)
!     use module_precision
!     use module_mesh
!     use module_params
!     use module_IO
!     use module_mpi
!     use module_operators
!     use module_physics_metamodule
!     use module_time_step
!     use module_stl_file_reader
!     use module_helpers
!
!     implicit none
!
!     type (type_params), intent(inout)  :: params
!     character(len=80) :: fname_ini, fname_stl, fname_out, dummy
!     integer :: i, Bs(1:3), g, ntri, k, iter, skips
!     integer :: ix, iy, iz, ivertex, xmin, xmax, ymin, ymax, zmin, zmax, safety, mpicode
!     real(kind=4), allocatable, dimension(:,:) :: triangles, normals
!     real(kind=rk) :: scale, origin(1:3)
!     ! origin and spacing of blocks
!     real(kind=rk) :: x0(1:3), dx(1:3), x,y,z,tmp
!
!     integer(kind=ik), allocatable      :: lgt_block(:, :)
!     real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
!     real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
!     real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
!     real(kind=rk), allocatable         :: normal_block(:, :, :, :)
!     integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
!     integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
!     integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
!     integer :: hvy_n, lgt_n, hvy_id, lgt_id
!     integer :: c_plus, c_minus
!     integer :: ix1,iy1,iz1
!     logical :: done
!
!     !-----------------------------------------------------------------------------------------------------
!     ! get values from command line (filename and level for interpolation)
!     call get_command_argument(2, fname_ini)
!
!     ! does the user need help?
!     if (fname_ini=='--help' .or. fname_ini=='--h' .or. fname_ini=='-h') then
!         if (params%rank==0) then
!             write(*,*) "------------------------------------------------------------------"
!             write(*,*) "./wabbit-post --stl2dist --x0 2.0,3.0,4.0 --stl some.stl --params PARAMS.ini"
!             write(*,*) "------------------------------------------------------------------"
!             write(*,*) ""
!             write(*,*) ""
!             write(*,*) ""
!             write(*,*) ""
!             write(*,*) "------------------------------------------------------------------"
!         end if
!         return
!     endif
!
!     ! defaults
!     scale = 1.0_rk
!
!     ! fetch parameters from command line call
!     do i = 1, COMMAND_ARGUMENT_COUNT()
!
!         call get_command_argument(i,dummy)
!
!         select case (dummy)
!         case ("--params")
!             call get_command_argument(i+1, fname_ini)
!             call check_file_exists( fname_ini )
!
!         case ("--stl")
!             call get_command_argument(i+1, fname_stl)
!             call check_file_exists( fname_stl )
!
!         case ("-o")
!             call get_command_argument(i+1, fname_out)
!
!         case ("--scale")
!             call get_command_argument(i+1, dummy)
!             read(dummy,*) scale
!
!         case ("--x0")
!             call get_command_argument(i+1, dummy)
!             read(dummy,*) origin
!
!         end select
!     enddo
!
!     if (params%rank==0) then
!         write(*,'("STL file is ",A)') trim(adjustl(fname_stl))
!         write(*,'("INI file is ",A)') trim(adjustl(fname_ini))
!         write(*,'("STL scaling factor is scale=",g12.3)') scale
!         write(*,'("origin shift=",3(g12.3,1x))') origin
!     endif
!
!     ! read ini-file and save parameters in struct
!     call ini_file_to_params( params, fname_ini )
!     ! have the pysics module read their own parameters
!     call init_physics_modules( params, fname_ini, params%N_mask_components  )
!
!     ! one field for the result, one field to tag error points where we have trouble determining the sign.
!     params%n_eqn = 2
!     ! params%n_ghosts =
!
!     ! in usual parameter files, RK4 (or some other RK) is used an requires a lot of memory
!     ! here we do not need that, and hence pretent to use a basic scheme (EE1 maybe)
!     deallocate(params%butcher_tableau)
!     allocate(params%butcher_tableau(1,1))
!
!     Bs = params%Bs
!     g = params%n_ghosts
!
!     ! allocate data
!     call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, &
!     lgt_active, hvy_active, lgt_sortednumlist)
!
!     call reset_tree( params, lgt_block, lgt_active, &
!     lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true., tree_ID=1)
!
!     ! start with an equidistant grid on coarsest level
!     call create_equidistant_grid( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
!     lgt_sortednumlist, hvy_active, hvy_n, params%min_treelevel, .true., tree_ID=1 )
!
!     ! reset grid to zeros
!     do k = 1, hvy_n
!         hvy_block(:,:,:,:,hvy_active(k)) = 0.0_rk
!     enddo
!
!
!     call read_stl_file(fname_stl, ntri, triangles, normals)
!     call normalize_stl_file( ntri, triangles, "--scale", scale, origin )
!
!     safety = 6!12
!
!     allocate( normal_block(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g, 1:3) )
!     normal_block = 0.0_rk
!
!     do iter = params%min_treelevel, params%max_treelevel
!
!         ! synchronization before refinement (because the interpolation takes place on the extended blocks
!         ! including the ghost nodes)
!         ! Note: at this point the grid is rather coarse (fewer blocks), and the sync step is rather cheap.
!         ! Snych'ing becomes much mor expensive one the grid is refined.
!         call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
!
!         ! refine the mesh. Note: afterwards, it can happen that two blocks on the same level differ
!         ! in their redundant nodes, but the ghost node sync'ing later on will correct these mistakes.
!         call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
!         lgt_sortednumlist, hvy_active, hvy_n, "mask-threshold", tree_ID=1 )
!
!         skips = 0
!
!         do k = 1, hvy_n
!             ! hvy_id of the block we're looking at
!             hvy_id = hvy_active(k)
!
!             ! light id of this block
!
!             hvy_block(:,:,:,2,hvy_id) = 0.0_rk ! no problem
!             normal_block = 0.0_rk
!
!             ! The trick here is that the STL does not move, and that we set it on the ghost
!             ! nodes as well. Then, a refined block should always end up with nonzero values
!             ! only if its mother had nonzero values.
!             ! That implies: if the block is zero, we can skip it
!             if ( maxval(hvy_block(:,:,:,1,hvy_id)) <= 1.0e-6 .and. iter>params%min_treelevel ) then
!                 hvy_block(:,:,:,1,hvy_id) = 9e8_rk ! distance as far away
!                 skips = skips + 1
!                 cycle
!             endif
!
!             hvy_block(:,:,:,1,hvy_id) = 9e8_rk ! distance as far away
!
!             ! compute block spacing and origin from treecode
!             call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
!             call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
!
!             ! shift origin to take ghost nodes into account
!             x0 = x0 - dble(g)*dx
!
!             do i = 1, ntri
!                 ivertex = 3*i - 2
!
!                 xmin = floor( minval( (triangles(1,ivertex:ivertex+2)-x0(1)) ) / dx(1)) - safety
!                 ymin = floor( minval( (triangles(2,ivertex:ivertex+2)-x0(2)) ) / dx(2)) - safety
!                 zmin = floor( minval( (triangles(3,ivertex:ivertex+2)-x0(3)) ) / dx(3)) - safety
!
!                 xmax = ceiling( maxval( (triangles(1,ivertex:ivertex+2)-x0(1)) ) / dx(1)) + safety
!                 ymax = ceiling( maxval( (triangles(2,ivertex:ivertex+2)-x0(2)) ) / dx(2)) + safety
!                 zmax = ceiling( maxval( (triangles(3,ivertex:ivertex+2)-x0(3)) ) / dx(3)) + safety
!
!                 xmin = max(xmin, 1)
!                 ymin = max(ymin, 1)
!                 zmin = max(zmin, 1)
!
!                 xmax = min(xmax, Bs(1)+2*g)
!                 ymax = min(ymax, Bs(2)+2*g)
!                 zmax = min(zmax, Bs(3)+2*g)
!
!                 do iz = zmin, zmax
!                     z = dx(3)*dble(iz) + x0(3)
!                     do iy = ymin, ymax
!                         y = dx(2)*dble(iy) + x0(2)
!                         do ix = xmin, xmax
!                             x = dx(1)*dble(ix) + x0(1)
!
!                             ! the distance to the current triangle:
!                             tmp = pointTriangleDistance( triangles(:,ivertex), triangles(:,ivertex+1), &
!                             triangles(:,ivertex+2), (/x,y,z/), normals(:,i) )
!
!                             ! same distance ?
!                             if (abs(abs(tmp) - abs(hvy_block(ix,iy,iz,1,hvy_id))) < 1.0e-11_rk) then
!                                 ! different sign
!                                 if (sign(1.0_rk, tmp) /= sign(1.0_rk, hvy_block(ix,iy,iz,1,hvy_id))) then
!                                     ! problem case
!                                     ! different normal
!                                     hvy_block(ix,iy,iz,1,hvy_id) = pointTriangleDistance( triangles(:,ivertex), triangles(:,ivertex+1), &
!                                     triangles(:,ivertex+2), (/x,y,z/), normals(:,i)+real(normal_block(ix,iy,iz,1:3),kind=4) )
!
!                                     normal_block(ix,iy,iz,1:3) = normals(:,i) + normal_block(ix,iy,iz,1:3)
!                                 endif
!                             else
!                                 ! if closer (in abs value!) then use this now
!                                 if ( abs(tmp) < abs(hvy_block(ix,iy,iz,1,hvy_id)) ) then
!                                     hvy_block(ix,iy,iz,1,hvy_id) = tmp
!                                     normal_block(ix,iy,iz,1:3) = normals(:,i)
!                                 endif
!                             endif
!                         enddo
!                     enddo
!                 enddo
!
!             enddo ! loop over triangles
!         enddo ! loop blocks
!
!
!
!         call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
!
!         ! convert block to mask function
!         do k = 1, hvy_n
!             hvy_id = hvy_active(k)
!             do iz = 1, Bs(3)+2*g
!                 do iy = 1, Bs(2)+2*g
!                     do ix = 1, Bs(1)+2*g
!                         hvy_block(ix,iy,iz,1,hvy_id) = smoothstep( hvy_block(ix,iy,iz,1,hvy_id), 0.0_rk, 1.5_rk*dx(1) )
!                     enddo
!                 enddo
!             enddo
!         enddo ! loop over blocks
!
!         write(*,'("rank=",i4," skipped ",i6," of its ",i6," blocks")') params%rank, skips, hvy_n
!
!         call MPI_barrier( WABBIT_COMM, mpicode)
!
!
!         if (params%rank==0) then
!             write(*, '("Nb=",i6," Jmin=",i2," Jmax=",i2)') &
!             lgt_n, min_active_level( lgt_block, lgt_active, lgt_n ), &
!             max_active_level( lgt_block, lgt_active, lgt_n )
!         endif
!     enddo ! loop over level
!
!     call write_field(fname_out, 0.0_rk, -99, 1, params, lgt_block, &
!     hvy_block(:,:,:,:,:), lgt_active, lgt_n, hvy_n, hvy_active )
!
! end subroutine




! subroutine post_stl2dist(params)
!     use module_precision
!     use module_mesh
!     use module_params
!     use module_IO
!     use module_mpi
!     use module_operators
!     use module_physics_metamodule
!     use module_time_step
!     use module_stl_file_reader
!     use module_helpers
!
!     implicit none
!
!     type (type_params), intent(inout)  :: params
!     character(len=80) :: fname_ini, fname_stl, fname_out, dummy
!     integer :: i, Bs(1:3), g, ntri, k, iter, skips
!     integer :: ix, iy, iz, ivertex, xmin, xmax, ymin, ymax, zmin, zmax, safety, mpicode
!     real(kind=4), allocatable, dimension(:,:) :: triangles, normals
!     real(kind=rk) :: scale, origin(1:3)
!     ! origin and spacing of blocks
!     real(kind=rk) :: x0(1:3), dx(1:3), x,y,z,tmp
!
!     integer(kind=ik), allocatable      :: lgt_block(:, :)
!     real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
!     real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
!     real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
!     integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
!     integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
!     integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
!     integer :: hvy_n, lgt_n, hvy_id, lgt_id
!     integer :: c_plus, c_minus, N_problems1, N_problems2
!     integer :: ix1,iy1,iz1
!     logical :: done
!
!     !-----------------------------------------------------------------------------------------------------
!     ! get values from command line (filename and level for interpolation)
!     call get_command_argument(2, fname_ini)
!
!     ! does the user need help?
!     if (fname_ini=='--help' .or. fname_ini=='--h' .or. fname_ini=='-h') then
!         if (params%rank==0) then
!             write(*,*) "------------------------------------------------------------------"
!             write(*,*) "./wabbit-post --stl2dist --x0 2.0,3.0,4.0 --stl some.stl --params PARAMS.ini"
!             write(*,*) "------------------------------------------------------------------"
!             write(*,*) ""
!             write(*,*) ""
!             write(*,*) ""
!             write(*,*) ""
!             write(*,*) "------------------------------------------------------------------"
!         end if
!         return
!     endif
!
!     ! defaults
!     scale = 1.0_rk
!
!     ! fetch parameters from command line call
!     do i = 1, COMMAND_ARGUMENT_COUNT()
!
!         call get_command_argument(i,dummy)
!
!         select case (dummy)
!         case ("--params")
!             call get_command_argument(i+1, fname_ini)
!             call check_file_exists( fname_ini )
!
!         case ("--stl")
!             call get_command_argument(i+1, fname_stl)
!             call check_file_exists( fname_stl )
!
!         case ("-o")
!             call get_command_argument(i+1, fname_out)
!
!         case ("--scale")
!             call get_command_argument(i+1, dummy)
!             read(dummy,*) scale
!
!         case ("--x0")
!             call get_command_argument(i+1, dummy)
!             read(dummy,*) origin
!
!         end select
!     enddo
!
!     if (params%rank==0) then
!         write(*,'("STL file is ",A)') trim(adjustl(fname_stl))
!         write(*,'("INI file is ",A)') trim(adjustl(fname_ini))
!         write(*,'("STL scaling factor is scale=",g12.3)') scale
!         write(*,'("origin shift=",3(g12.3,1x))') origin
!     endif
!
!     ! read ini-file and save parameters in struct
!     call ini_file_to_params( params, fname_ini )
!     ! have the pysics module read their own parameters
!     call init_physics_modules( params, fname_ini, params%N_mask_components  )
!
!     ! one field for the result, one field to tag error points where we have trouble determining the sign.
!     params%n_eqn = 2
!
!     ! in usual parameter files, RK4 (or some other RK) is used an requires a lot of memory
!     ! here we do not need that, and hence pretent to use a basic scheme (EE1 maybe)
!     deallocate(params%butcher_tableau)
!     allocate(params%butcher_tableau(1,1))
!
!     Bs = params%Bs
!     g = params%n_ghosts
!
!     ! allocate data
!     call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, &
!     lgt_active, hvy_active, lgt_sortednumlist)
!
!     call reset_tree( params, lgt_block, lgt_active, &
!     lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true., tree_ID=1)
!
!     ! start with an equidistant grid on coarsest level
!     call create_equidistant_grid( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
!     lgt_sortednumlist, hvy_active, hvy_n, params%min_treelevel, .true., tree_ID=1 )
!
!     ! reset grid to zeros
!     do k = 1, hvy_n
!         hvy_block(:,:,:,:,hvy_active(k)) = 0.0_rk
!     enddo
!
!     ! read the STL we want to tranform to a signed distance function.
!     call read_stl_file(fname_stl, ntri, triangles, normals)
!     call normalize_stl_file( ntri, triangles, "--scale", scale, origin )
!
!     safety = 6!12
!
!     !=====================================================================================================
!     ! start from coarsest level, work up to Jmax
!     !=====================================================================================================
!     call debug_header_barrier("Phase outer starting")
!     do iter = params%min_treelevel, params%max_treelevel
!         ! synchronization before refinement (because the interpolation takes place on the extended blocks
!         ! including the ghost nodes)
!         call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
!
!         ! refine the mesh. Note: afterwards, it can happen tha  t two blocks on the same level differ
!         ! in their redundant nodes, but the ghost node sync'ing later on will correct these mistakes.
!         call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
!         lgt_sortednumlist, hvy_active, hvy_n, "mask-threshold", tree_ID=1 )
!
!         skips = 0
!
!         !=====================================================================================================
!         ! Compute distance and sign where possiblem flag points where sign could not be determined
!         !=====================================================================================================
!         call debug_header_barrier("Phase 1 starting")
!         do k = 1, hvy_n
!             ! hvy_id of the block we're looking at
!             hvy_id = hvy_active(k)
!
!             hvy_block(:,:,:,2,hvy_id) = 0.0_rk ! no problem
!
!             ! The trick here is that the STL does not move, and that we set it on the ghost
!             ! nodes as well. Then, a refined block should always end up with nonzero values
!             ! only if its mother had nonzero values.
!             ! That implies: if the block is zero, we can skip it
!             if ( maxval(hvy_block(:,:,:,1,hvy_id)) <= 1.0e-6 .and. iter>params%min_treelevel ) then
!                 hvy_block(:,:,:,1,hvy_id) = 9e8_rk ! distance as far away
!                 skips = skips + 1
!                 cycle
!             endif
!
!             hvy_block(:,:,:,1,hvy_id) = 9e8_rk ! distance as far away
!
!             ! compute block spacing and origin from treecode
!             call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
!             call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
!
!             ! shift origin to take ghost nodes into account
!             x0 = x0 - dble(g)*dx
!
!             do i = 1, ntri
!                 ivertex = 3*i - 2
!
!                 xmin = floor( minval( (triangles(1,ivertex:ivertex+2)-x0(1)) ) / dx(1)) - safety
!                 ymin = floor( minval( (triangles(2,ivertex:ivertex+2)-x0(2)) ) / dx(2)) - safety
!                 zmin = floor( minval( (triangles(3,ivertex:ivertex+2)-x0(3)) ) / dx(3)) - safety
!
!                 xmax = ceiling( maxval( (triangles(1,ivertex:ivertex+2)-x0(1)) ) / dx(1)) + safety
!                 ymax = ceiling( maxval( (triangles(2,ivertex:ivertex+2)-x0(2)) ) / dx(2)) + safety
!                 zmax = ceiling( maxval( (triangles(3,ivertex:ivertex+2)-x0(3)) ) / dx(3)) + safety
!
!                 xmin = max(xmin, 1)
!                 ymin = max(ymin, 1)
!                 zmin = max(zmin, 1)
!
!                 xmax = min(xmax, Bs(1)+2*g)
!                 ymax = min(ymax, Bs(2)+2*g)
!                 zmax = min(zmax, Bs(3)+2*g)
!
!                 do iz = zmin, zmax
!                     z = dx(3)*dble(iz) + x0(3)
!                     do iy = ymin, ymax
!                         y = dx(2)*dble(iy) + x0(2)
!                         do ix = xmin, xmax
!                             x = dx(1)*dble(ix) + x0(1)
!
!                             ! the distance to the current triangle:
!                             tmp = pointTriangleDistance( triangles(:,ivertex), triangles(:,ivertex+1), &
!                             triangles(:,ivertex+2), (/x,y,z/), normals(:,i) )
!
!                             ! same distance ?
!                             if (abs(abs(tmp) - abs(hvy_block(ix,iy,iz,1,hvy_id))) < 1.0e-11_rk) then
!                                 ! different sign
!                                 if (sign(1.0_rk, tmp) /= sign(1.0_rk, hvy_block(ix,iy,iz,1,hvy_id))) then
!                                     ! problem case
!                                     hvy_block(ix,iy,iz,1,hvy_id) = abs(tmp) ! set only value of distance, not the sign
!                                     hvy_block(ix,iy,iz,2,hvy_id) = 1.0_rk ! problem
!                                 endif
!                             else
!                                 ! if closer (in abs value!) then use this now
!                                 if ( abs(tmp) < abs(hvy_block(ix,iy,iz,1,hvy_id)) ) then
!                                     hvy_block(ix,iy,iz,1,hvy_id) = tmp
!                                     hvy_block(ix,iy,iz,2,hvy_id) = 0.0_rk ! no problem
!                                 endif
!                             endif
!                         enddo
!                     enddo
!                 enddo
!
!             enddo ! loop over triangles
!         enddo ! loop blocks
!
!         call debug_header_barrier("Phase 1 completed")
!
!         !=====================================================================================================
!         ! fix problem signs by assigning the sign of nearest neighbors
!         !=====================================================================================================
!         call debug_header_barrier("Phase 2 starting")
!         done = .false.
!         do while (.not. done)
!             done = .true. ! tentatively
!
!             ! rarely required..but blocks that have "problems-only" require it
!             call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
!
!             do k = 1, hvy_n
!                 hvy_id = hvy_active(k)
!
!                 N_problems1 = sum(hvy_block(:,:,:,2,hvy_id))
!                 N_problems2 = 0
!
!                 if (N_problems1 == 0) cycle ! no problem here...
!
!                 do while (N_problems1 /= N_problems2)
!                     N_problems1 = sum(hvy_block(:,:,:,2,hvy_id))
!                     do iz = g, Bs(3)+g
!                         do iy = g, Bs(2)+g
!                             do ix = g, Bs(1)+g
!                                 if (hvy_block(ix,iy,iz,2,hvy_id) > 0.0_rk) then ! problem
!                                     xmin = max(ix-1, 1)
!                                     ymin = max(iy-1, 1)
!                                     zmin = max(iz-1, 1)
!
!                                     xmax = min(ix+1, Bs(1)+2*g)
!                                     ymax = min(iy+1, Bs(2)+2*g)
!                                     zmax = min(iz+1, Bs(3)+2*g)
!
!                                     c_plus = 0
!                                     c_minus = 0
!
!                                     do iz1 = zmin, zmax
!                                         do iy1 = ymin, ymax
!                                             do ix1 = xmin, xmax
!                                                 if (hvy_block(ix1,iy1,iz1,2,hvy_id) < 1.0_rk) then ! sign is okay, no problem at this neighboring point
!                                                     if ( sign(1.0_rk, hvy_block(ix1,iy1,iz1,1,hvy_id)) > 0.0_rk ) then
!                                                         c_plus = c_plus+1
!                                                     else
!                                                         c_minus = c_minus+1
!                                                     endif
!                                                 endif
!                                             enddo
!                                         enddo
!                                     enddo
!
!                                     if (c_plus >0 .or. c_minus>0) then
!                                         if (c_plus<c_minus) then
!                                             hvy_block(ix,iy,iz,1,hvy_id) = -hvy_block(ix,iy,iz,1,hvy_id) !sign inversion
!                                         endif
!                                         hvy_block(ix,iy,iz,2,hvy_id) = 0.0_rk ! problem solved
!                                     endif
!
!                                 endif
!                             enddo
!                         enddo
!                     enddo
!                     N_problems2 = sum(hvy_block(:,:,:,2,hvy_id))
!                 enddo
!
!                 N_problems1 = sum(hvy_block(:,:,:,2,hvy_id))
!                 if ( N_problems1 > 0.0_rk ) done = .false.
!             enddo ! loop over blocks
!
!
!             call MPI_Allreduce(MPI_IN_PLACE, done, 1, MPI_LOGICAL, MPI_LAND, WABBIT_COMM, k )
!         enddo ! while not done
!         call debug_header_barrier("Phase 2 completed")
!         call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
!
!         !=====================================================================================================
!         ! convert block to mask function (required to check if we can skip a block or not.)
!         !=====================================================================================================
!         do k = 1, hvy_n
!             hvy_id = hvy_active(k)
!             do iz = 1, Bs(3)+2*g
!                 do iy = 1, Bs(2)+2*g
!                     do ix = 1, Bs(1)+2*g
!                         hvy_block(ix,iy,iz,1,hvy_id) = smoothstep( hvy_block(ix,iy,iz,1,hvy_id), 0.0_rk, 1.5_rk*dx(1) )
!                     enddo
!                 enddo
!             enddo
!         enddo ! loop over blocks
!
!         write(*,'("rank=",i4," skipped ",i6," of its ",i6," blocks")') params%rank, skips, hvy_n
!         call MPI_barrier( WABBIT_COMM, mpicode)
!
!         if (params%rank==0) then
!             write(*, '("Nb=",i6," Jmin=",i2," Jmax=",i2)') &
!             lgt_n, min_active_level( lgt_block, lgt_active, lgt_n ), &
!             max_active_level( lgt_block, lgt_active, lgt_n )
!         endif
!     enddo ! loop over level
!
!     call write_field(fname_out, 0.0_rk, -99, 1, params, lgt_block, &
!     hvy_block(:,:,:,:,:), lgt_active, lgt_n, hvy_n, hvy_active )
!
! end subroutine
