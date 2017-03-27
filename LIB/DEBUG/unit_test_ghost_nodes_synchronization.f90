! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: unit_test_ghost_nodes_synchronization.f90
! version: 0.5
! author: msr
!
! unit test for ghost nodes synchronization
! note: input only params struct to this subroutine
!       create new light/heavy data arrays here and deallocate them after this function
!
! input:    - params
! output:   -
!
! = log ======================================================================================
!
! 21/01/17 - create
!
! ********************************************************************************************

subroutine unit_test_ghost_nodes_synchronization( params )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)          :: params

    ! local copy of params struct
    ! user defined parameter structure
    type (type_params)                      :: params_loc

    ! light data array
    integer(kind=ik), allocatable           :: lgt_block_loc(:, :)
    ! heavy data array - block data
    real(kind=rk), allocatable              :: hvy_block_loc(:, :, :, :, :), hvy_block_loc_exact(:, :, :, :, :)

    ! list of active blocks (heavy data)
    integer(kind=ik), allocatable           :: hvy_active(:)
    ! number of active blocks (heavy data)
    integer(kind=ik), allocatable           :: hvy_n

    ! list of active blocks (light data)
    integer(kind=ik), allocatable           :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik)                        :: lgt_n

    ! neighbor array
    integer(kind=ik), allocatable           :: hvy_neighbor_loc(:,:)

    ! allocation error variable
    integer(kind=ik)                        :: allocate_error

    ! loop variables
    integer(kind=ik)                        :: k, N, l, i, j

    ! process rank
    integer(kind=ik)                        :: rank

    ! number of blocks in one dimension
    integer(kind=ik)                        :: block_num

    ! coordinates vectors
    real(kind=rk), allocatable              :: coord_x(:), coord_y(:), coord_z(:)
    ! spacing
    real(kind=rk)                           :: dx, dy, dz

    ! grid parameter
    integer(kind=ik)                        :: Bs, g
    real(kind=rk)                           :: Lx, Ly, Lz

    ! error variable
    real(kind=rk)                           :: error, my_error

    ! MPI error variable
    integer(kind=ik)                        :: ierr

    ! chance for block refinement, random number
    real(kind=rk)                           :: ref_chance, r

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! local copy of params struct
    params_loc = params

    ! number of blocks per proc
    N = params_loc%number_blocks

    ! set MPI parameters
    rank = params_loc%rank

    ! grid parameter
    Lx = params_loc%Lx
    Ly = params_loc%Ly
    Lz = params_loc%Lz

!---------------------------------------------------------------------------------------------
! main body

    ! choose number of blocks here
    block_num = 4!4!2

    ! reset blocksize in params struct
    params_loc%number_block_nodes = (params_loc%number_domain_nodes + (block_num-1)) / block_num
    Bs = params_loc%number_block_nodes
    g = params_loc%number_ghost_nodes

    !---------------------------------------------------------------------------------------------
    ! first: initializing new block data
    ! allocate block_list
    call allocate_block_list( params_loc, lgt_block_loc )
    ! allocate heavy data
    call allocate_block_data( params_loc, hvy_block_loc )
    call allocate_block_data( params_loc, hvy_block_loc_exact )

    ! active heavy block list
    allocate( hvy_active( N ), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    ! reset active list
    hvy_active = -1
    ! list index
    hvy_n = 0

    ! init data array with zeros => after this: blocks have correct coordinate vectors
    call inicond_zeros( params_loc, lgt_block_loc, hvy_block_loc )

    ! list of active blocks, heavy data
    call create_hvy_active_list( lgt_block_loc, hvy_active, hvy_n )

    !---------------------------------------------------------------------------------------------
    ! second: refine some blocks (random)

    ! create active light block list
    allocate( lgt_active( size(lgt_block_loc, 1) ), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    ! reset active list
    lgt_active = -1
    ! list index
    lgt_n = 0

    ! list of active blocks, light data
    call create_lgt_active_list( lgt_block_loc, lgt_active, lgt_n )

    ! refinement chance
    ref_chance = 0.1_rk

    ! set random seed
    call init_random_seed()

    ! set refinement status
    do k = 1, lgt_n
        ! random number
        call random_number(r)
        ! set refinement status
        if ( r <= ref_chance ) then
            lgt_block_loc( lgt_active(k), params_loc%max_treelevel+2 ) = 1
        end if
    end do

    ! refine blocks
    ! check if block has reached maximal level
    call respect_min_max_treelevel( params_loc, lgt_block_loc, lgt_active, lgt_n )

    ! interpolate the new mesh
    if ( params%threeD_case ) then
        ! 3D:
        call refine_mesh_3D( params_loc, lgt_block_loc, hvy_block_loc, hvy_active, hvy_n )
    else
        ! 2D:
        call refine_mesh_2D( params_loc, lgt_block_loc, hvy_block_loc(:,:,1,:,:), hvy_active, hvy_n )
    end if

    ! update lists of active blocks (light and heavy data)
    call create_lgt_active_list( lgt_block_loc, lgt_active, lgt_n )
    call create_hvy_active_list( lgt_block_loc, hvy_active, hvy_n )

    !---------------------------------------------------------------------------------------------
    ! third: fill blocks with data

    ! allocate coord arrays
    allocate(coord_x( Bs + 2*g ), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if
    allocate(coord_y( Bs + 2*g ), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if
    allocate(coord_z( Bs + 2*g ), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    ! loop over all active blocks
    do k = 1, hvy_n

        ! write coord arrays
        coord_x(g+1:Bs+g) = hvy_block_loc(1, 1:Bs, 1, 1, hvy_active(k))
        coord_y(g+1:Bs+g) = hvy_block_loc(2, 1:Bs, 1, 1, hvy_active(k))
        coord_z(g+1:Bs+g) = hvy_block_loc(3, 1:Bs, 1, 1, hvy_active(k))

        ! spacing
        dx = abs( hvy_block_loc(1, 2, 1, 1, hvy_active(k)) - hvy_block_loc(1, 1, 1, 1, hvy_active(k)) )
        dy = abs( hvy_block_loc(2, 2, 1, 1, hvy_active(k)) - hvy_block_loc(2, 1, 1, 1, hvy_active(k)) )
        dz = abs( hvy_block_loc(3, 2, 1, 1, hvy_active(k)) - hvy_block_loc(3, 1, 1, 1, hvy_active(k)) )

        ! ghost nodes coordinates
        do l = 1, g

            ! minus direction
            coord_x(l) = hvy_block_loc(1, 1, 1, 1, hvy_active(k)) - (g+1-l)*dx
            coord_y(l) = hvy_block_loc(2, 1, 1, 1, hvy_active(k)) - (g+1-l)*dy
            coord_z(l) = hvy_block_loc(3, 1, 1, 1, hvy_active(k)) - (g+1-l)*dz

            ! plus direction
            coord_x(Bs+g+l) = hvy_block_loc(1, Bs, 1, 1, hvy_active(k)) + (l)*dx
            coord_y(Bs+g+l) = hvy_block_loc(2, Bs, 1, 1, hvy_active(k)) + (l)*dy
            coord_z(Bs+g+l) = hvy_block_loc(3, Bs, 1, 1, hvy_active(k)) + (l)*dz

        end do

        ! calculate f(x,y,z) for first datafield
        if ( params_loc%threeD_case ) then
            ! 3D:
            call f_xyz_3D( coord_x, coord_y, coord_z, hvy_block_loc(:, :, :, 2, hvy_active(k)), Bs, g, Lx, Ly, Lz )

        else
            ! 2D:
            call f_xy_2D( coord_x, coord_y, hvy_block_loc(:, :, 1, 2, hvy_active(k)), Bs, g, Lx, Ly )
        end if

    end do

    !---------------------------------------------------------------------------------------------
    ! fourth: synchronize ghost nodes

    ! save heavy data
    hvy_block_loc_exact = hvy_block_loc

    ! init neighbor data array
    ! 2D: maximal 16 neighbors per block
    ! 3D: maximal 74 neighbors per block
    if ( params%threeD_case ) then
        ! 3D:
        allocate( hvy_neighbor_loc( params_loc%number_blocks, 74 ), stat=allocate_error )
        !call check_allocation(allocate_error)
        if ( allocate_error /= 0 ) then
            write(*,'(80("_"))')
            write(*,*) "ERROR: memory allocation fails"
            stop
        end if

    else
        ! 2D:
        allocate( hvy_neighbor_loc( params_loc%number_blocks, 16 ), stat=allocate_error )
        !call check_allocation(allocate_error)
        if ( allocate_error /= 0 ) then
            write(*,'(80("_"))')
            write(*,*) "ERROR: memory allocation fails"
            stop
        end if

    end if

    ! update neighbors
    if ( params_loc%threeD_case ) then
        ! 3D:
        call update_neighbors_3D( params_loc, lgt_block_loc, hvy_neighbor_loc, lgt_active, lgt_n, hvy_active, hvy_n )
    else
        ! 2D:
        call update_neighbors_2D( params_loc, lgt_block_loc, hvy_neighbor_loc, lgt_active, lgt_n, hvy_active, hvy_n )
    end if

    ! synchronize ghost nodes
    call synchronize_ghosts( params_loc, lgt_block_loc, hvy_block_loc, hvy_neighbor_loc, hvy_active, hvy_n )

    !---------------------------------------------------------------------------------------------
    ! fifth: compare values

    ! reset error
    error = 0.0_rk
    my_error = 0.0_rk

    ! loop over all active blocks
    do k = 1, hvy_n

        if ( params_loc%threeD_case ) then
            ! 3D:
            do i = 1, Bs+2*g
                do j = 1, Bs+2*g
                    do l = 1, Bs+2*g

                        if ( sqrt( ( hvy_block_loc(i,j,l,2,hvy_active(k)) - hvy_block_loc_exact(i,j,l,2,hvy_active(k)) )**2  ) > 1e-6) then
                            my_error = my_error + sqrt( ( hvy_block_loc(i,j,l,2,hvy_active(k)) - hvy_block_loc_exact(i,j,l,2,hvy_active(k)) )**2  )
                        end if

                    end do
                end do
            end do
        else
            ! 2D:
            do i = 1, Bs+2*g
                do j = 1, Bs+2*g
                    my_error = my_error + sqrt( ( hvy_block_loc(i,j,1,2,hvy_active(k)) - hvy_block_loc_exact(i,j,1,2,hvy_active(k)) )**2  )
                end do
            end do

        end if

    end do

    ! synchronize errors
    call MPI_Allreduce(my_error, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! output
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("UNIT TEST: ghost nodes synchronization error = ",f16.8)')  error
    end if

    !---------------------------------------------------------------------------------------------
    ! last: clean up
    deallocate(lgt_block_loc, stat=allocate_error )
    deallocate(hvy_block_loc, stat=allocate_error )
    deallocate(hvy_block_loc_exact, stat=allocate_error )
    deallocate( hvy_active, stat=allocate_error )
    deallocate( lgt_active, stat=allocate_error )
    deallocate(coord_x, stat=allocate_error )
    deallocate(coord_y, stat=allocate_error )
    deallocate(coord_z, stat=allocate_error )
    deallocate( hvy_neighbor_loc, stat=allocate_error )

end subroutine unit_test_ghost_nodes_synchronization
