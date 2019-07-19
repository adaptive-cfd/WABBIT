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

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=80) :: fname_ini, fname_stl, fname_out, dummy
    integer :: i, Bs(1:3), g, ntri, k, iter, skips
    integer :: ix, iy, iz, ivertex, xmin, xmax, ymin, ymax, zmin, zmax, safety, mpicode
    real(kind=4), allocatable, dimension(:,:) :: triangles, normals
    real(kind=rk) :: scale, origin(1:3)
    ! origin and spacing of blocks
    real(kind=rk) :: x0(1:3), dx(1:3), x,y,z,tmp

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
    integer :: hvy_n, lgt_n, hvy_id, lgt_id

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(2, fname_ini)

    ! does the user need help?
    if (fname_ini=='--help' .or. fname_ini=='--h' .or. fname_ini=='-h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "./wabbit-post --stl2dist --x0 2.0,3.0,4.0 --stl some.stl --params PARAMS.ini"
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

    ! fetch parameters from command line call
    do i = 1, COMMAND_ARGUMENT_COUNT()

        call get_command_argument(i,dummy)

        select case (dummy)
        case ("--params")
            call get_command_argument(i+1, fname_ini)
            call check_file_exists( fname_ini )

        case ("--stl")
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
        write(*,'("STL file is ",A)') trim(adjustl(fname_stl))
        write(*,'("INI file is ",A)') trim(adjustl(fname_ini))
        write(*,'("STL scaling factor is scale=",g12.3)') scale
        write(*,'("origin shift=",3(g12.3,1x))') origin
    endif

    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, fname_ini )
    ! have the pysics module read their own parameters
    call init_physics_modules( params, fname_ini, params%N_mask_components  )

    !
    params%n_eqn = 1
    ! params%n_ghosts =

    ! in usual parameter files, RK4 (or some other RK) is used an requires a lot of memory
    ! here we do not need that, and hence pretent to use a basic scheme (EE1 maybe)
    deallocate(params%butcher_tableau)
    allocate(params%butcher_tableau(1,1))

    Bs = params%Bs
    g = params%n_ghosts

    ! allocate data
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, &
    lgt_active, hvy_active, lgt_sortednumlist)

    call reset_tree( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true., tree_ID=1)

    ! start with an equidistant grid on coarsest level
    call create_equidistant_grid( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, params%min_treelevel, .true., tree_ID=1 )

    ! read the STL we want to tranform to a signed distance function.
    ! NOTE:
    call read_stl_file(fname_stl, ntri, triangles, normals)
    call normalize_stl_file( ntri, triangles, "--scale", scale, origin )

    safety = 6

    do iter = 1, params%max_treelevel - params%min_treelevel

        ! synchronization before refinement (because the interpolation takes place on the extended blocks
        ! including the ghost nodes)
        ! Note: at this point the grid is rather coarse (fewer blocks), and the sync step is rather cheap.
        ! Snych'ing becomes much mor expensive one the grid is refined.
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

        ! refine the mesh. Note: afterwards, it can happen that two blocks on the same level differ
        ! in their redundant nodes, but the ghost node sync'ing later on will correct these mistakes.
        call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
        lgt_sortednumlist, hvy_active, hvy_n, "mask-threshold", tree_ID=1 )

        skips = 0

        do k = 1, hvy_n
            ! hvy_id of the block we're looking at
            hvy_id = hvy_active(k)

            ! light id of this block
            call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )

            ! The trick here is that the STL does not move, and that we set it on the ghost
            ! nodes as well. Then, a refined block should always end up with nonzero values
            ! only if its mother had nonzero values.
            ! That implies: if the block is zero, we can skip it
            if ( maxval(hvy_block(:,:,:,:,hvy_id)) <= 1.0e-6 ) then
                skips = skips + 1
                cycle
            endif

            ! compute block spacing and origin from treecode
            call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

            hvy_block(:,:,:,:,hvy_id) = 9e8

            ! shift origin to take ghost nodes into account
            x0 = x0 - dble(g)*dx

            do i = 1, ntri
                ivertex = 3*i - 2

                xmin = floor( minval( (triangles(1,ivertex:ivertex+2)-x0(1)) ) / dx(1)) - safety
                ymin = floor( minval( (triangles(2,ivertex:ivertex+2)-x0(2)) ) / dx(2)) - safety
                zmin = floor( minval( (triangles(3,ivertex:ivertex+2)-x0(3)) ) / dx(3)) - safety

                xmax = ceiling( maxval( (triangles(1,ivertex:ivertex+2)-x0(1)) ) / dx(1)) + safety
                ymax = ceiling( maxval( (triangles(2,ivertex:ivertex+2)-x0(2)) ) / dx(2)) + safety
                zmax = ceiling( maxval( (triangles(3,ivertex:ivertex+2)-x0(3)) ) / dx(3)) + safety

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
                            tmp = pointTriangleDistance( triangles(:,ivertex), triangles(:,ivertex+1), &
                            triangles(:,ivertex+2), (/x,y,z/), normals(:,i) )


                            ! if closer (in abs value!) then use this now
                            if ( dabs(tmp) < dabs(hvy_block(ix,iy,iz,1,hvy_id))) then
                                hvy_block(ix,iy,iz,1,hvy_id) = tmp
                            endif
                        enddo
                    enddo
                enddo

            enddo ! loop over triangles

            ! convert block to mask function
            do iz = 1, Bs(3)+2*g
                do iy = 1, Bs(2)+2*g
                    do ix = 1, Bs(1)+2*g
                        hvy_block(ix,iy,iz,1,hvy_id) = smoothstep( hvy_block(ix,iy,iz,1,hvy_id), 0.0_rk, 1.5_rk*dx(1) )
                    enddo
                enddo
            enddo

        enddo ! loop over blocks

        write(*,'("rank=",i4," skipped ",i6," of its ",i6," blocks")') params%rank, skips, hvy_n

        call MPI_barrier( WABBIT_COMM, mpicode)


        if (params%rank==0) then
            write(*, '("Nb=",i6," Jmin=",i2," Jmax=",i2)') &
             lgt_n, min_active_level( lgt_block, lgt_active, lgt_n ), &
             max_active_level( lgt_block, lgt_active, lgt_n )
         endif
    enddo

    call write_field(fname_out, 0.0_rk, -99, 1, params, lgt_block, &
    hvy_block(:,:,:,:,:), lgt_active, lgt_n, hvy_n, hvy_active )

end subroutine
