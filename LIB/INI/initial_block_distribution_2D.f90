! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: initial_block_distribution_2D.f90
! version: 0.5
! author: msr
!
! distribute blocks at start => create light data array
! note: field phi is 3D, but third dimension is not used
!
! input:    - parameters
!           - light block data array
!           - heavy block data array
!           - start field phi
! output:   - filled light and heavy data array
!
! = log ======================================================================================
!
! 07/11/16  - switch to v0.4
! 05/12/16  - add dummy space filling curve distribution
! 26/01/17  - use process rank from params struct
!           - remove sfc distribution: for sfc start dist, call balance_load after this routine
!           - use v0.5 hvy data array, but still work only for 2D data
!
! ********************************************************************************************

subroutine initial_block_distribution_2D( params, lgt_block, hvy_block, phi )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    ! initial data field
    real(kind=rk), intent(in)           :: phi(:, :, :)

    ! process rank
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! grid parameters (domain size, block size, number of ghost nodes, grid sizes)
    integer(kind=ik)                    :: Ds, Bs, g
    real(kind=rk)                       :: Lx, Ly

    ! block decomposition parameters
    integer(kind=ik)                    :: max_blocks, max_treelevel
    integer(kind=ik)                    :: num_blocks_x, num_blocks_y, num_blocks
    integer(kind=ik), allocatable       :: block_proc_list(:)

    ! loop variables
    integer(kind=ik)                    :: i, j, k

    ! domain coordinate vectors
    real(kind=rk), allocatable          :: coord_x(:), coord_y(:)

    ! heavy and light data id
    integer(kind=ik)                    :: heavy_id, light_id

    ! treecode variable and function to calculate size of treecode
    integer(kind=ik), allocatable       :: treecode(:)
    integer(kind=ik)                    :: treecode_size

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank            = params%rank
    number_procs    = params%number_procs

    ! set parameters for readability
    Ds              = params%number_domain_nodes
    Bs              = params%number_block_nodes
    g               = params%number_ghost_nodes

    Lx              = params%Lx
    Ly              = params%Ly

    max_blocks      = params%number_blocks
    max_treelevel   = params%max_treelevel

    ! allocate block to proc list
    allocate( block_proc_list( number_procs ), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    ! allocate domain coordinate vectors
    allocate( coord_x( Ds ), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    allocate( coord_y( Ds ), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    ! allocate treecode
    allocate( treecode( params%max_treelevel ), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

!---------------------------------------------------------------------------------------------
! main body

    ! simple uniformly distribution

    ! calculate starting block decomposition
    ! print block decomposition information
    ! every block has two more points in a single direction (from his neighbors)
    ! therefore the complete domain has also two additional points
    num_blocks_x        = (Ds-1) / (Bs-1)
    num_blocks_y        = (Ds-1) / (Bs-1)
    num_blocks          = num_blocks_x * num_blocks_y

    ! check given domain and block size
    if ( Ds /= (num_blocks_x-1)*(Bs-1) + Bs ) then
        print*, "ERROR: blocksize do not fit into domain size"
        stop
    end if

    ! output decomposition information
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("INIT: Field with res: ", i5, " x", i5, " gives: ", i5,  " x", i5, " (", i5, ") blocks of size: ", i5)') Ds, Ds, num_blocks_x, num_blocks_y, num_blocks, Bs
    end if

    ! decompose blocks to procs
    block_proc_list = (num_blocks - mod(num_blocks, number_procs))/number_procs
    ! distribute remaining blocks
    if (mod(num_blocks, number_procs) > 0) then
        block_proc_list(1:mod(num_blocks, number_procs)) = (num_blocks - mod(num_blocks, number_procs))/number_procs + 1
    end if

    ! calculate domain coordinate vectors
    do i = 1, Ds
        coord_x(i) = (i-1) * Lx / (Ds-1)
        coord_y(i) = Lx - (i-1) * Ly / (Ds-1)
    end do

    ! create block-tree
    k = 1
    do i = 1, num_blocks_x
        do j = 1, num_blocks_y
            ! ------------------------------------------------------------------------------------------------------
            ! write heavy data
            ! determine proc
            if (block_proc_list(k) == 0) then
                k = k + 1
                block_proc_list(k) = block_proc_list(k) - 1
            else
                block_proc_list(k) = block_proc_list(k) - 1
            end if

            ! find and set free heavy data id, note: look for free id in light data
            ! search routine only on corresponding light data -> so, returned id works directly on heavy data
            call get_free_light_id( heavy_id, lgt_block( (k-1)*max_blocks + 1 : ((k-1)+1)*max_blocks, 1 ), max_blocks )

            ! save data, write start field phi in first datafield
            ! note: third dimension of phi has size 1
            ! note: coordinates for z-dim are set to zero
            if (rank == (k-1)) then
                call new_block_heavy(params, &
                                     hvy_block, &
                                     heavy_id, &
                                     phi( (i-1)*(Bs-1) + 1 : i*(Bs-1) + 1 , (j-1)*(Bs-1) + 1 : j*(Bs-1) + 1, : ), &
                                     coord_x( (j-1)*(Bs-1) + 1 : j*(Bs-1) + 1 ), &
                                     coord_y( (i-1)*(Bs-1) + 1 : i*(Bs-1) + 1 ), &
                                     0.0_rk*coord_y( (i-1)*(Bs-1) + 1 : i*(Bs-1) + 1 ) )

            end if

            ! ------------------------------------------------------------------------------------------------------
            ! encoding treecode
            call encoding_2D(treecode, i, j, num_blocks_x, num_blocks_y, max_treelevel )

            ! ------------------------------------------------------------------------------------------------------
            ! write light data
            ! light data id is calculated from proc rank and heavy_id
            light_id = (k-1)*max_blocks + heavy_id
            ! write treecode
            lgt_block( light_id, 1 : max_treelevel ) = treecode
            ! treecode level (size)
            lgt_block( light_id, max_treelevel + 1 ) = treecode_size( treecode, max_treelevel )

        end do
    end do

    ! clean up
    deallocate( block_proc_list, stat=allocate_error )
    deallocate( coord_x, stat=allocate_error )
    deallocate( coord_y, stat=allocate_error )
    deallocate( treecode, stat=allocate_error )

end subroutine initial_block_distribution_2D
