! ********************************
! 2D AMR prototype
! --------------------------------
!
! initialize all data (params, fields, blocks, ...)
!
! name: init_data.f90
! date: 02.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine init_data()

    use module_params
    use module_blocks

    implicit none

    integer 								        :: allocate_error

    ! grid size
    integer(kind=ik)             					:: nx
    integer(kind=ik)							    :: i, j
    real(kind=rk) 							        :: mux, muy
    ! block size,
    integer(kind=ik)        						:: blocksize
    integer(kind=ik)                                :: max_blocks
    integer(kind=ik)        						:: ghosts
    ! derivatives
    real(kind=rk), dimension(:,:), allocatable      :: D1, D2
    ! gauss parameter
    real(kind=rk) 							        :: sigma
    ! 2D grid
    real(kind=rk), dimension(:,:), allocatable  	:: phi

    !------------------------------
    ! grid and block parameter
    nx          = 513 !512
    blocksize   = 17 !33 !65 !129 !257
    ! allocate twice as many blocks as required
    max_blocks = 2*((nx-1)/(blocksize-1))**2
    ghosts      = 4

    write(*,'(80("-"))')
    write(*,*) "INITIALIZATION PHASE"
    write(*,*) "we use a maximum number of blocks of",max_blocks
    write(*,*) "nghosts=", ghosts
    write(*,'(80("-"))')

    !------------------------------
    ! allocate memory for local variables
    allocate( phi(1:nx, 1:nx), stat=allocate_error )
    allocate( D1(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )
    allocate( D2(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )

    !------------------------------
    ! initial data field, memory allocation for module variables
    mux     = 0.5_rk * ( real(nx,kind=rk) - 1.0_rk )
    muy     = 0.5_rk * ( real(nx,kind=rk) - 1.0_rk )
    sigma   = 100.0_rk
    phi     = 0.0_rk

    do i = 1, nx
        do j = 1, nx
            phi(i,j) = dexp( -((real(i,kind=rk)-mux)**2 + (real(j,kind=rk)-muy)**2) / sigma)
        end do
    end do

    allocate( blocks_params%phi(nx, nx), stat=allocate_error )

    allocate( blocks(max_blocks), stat=allocate_error )
    ! dummy allocation
    allocate( blocks_params%active_list(1), stat=allocate_error )

    do i = 1, max_blocks
        allocate( blocks(i)%data1(blocksize, blocksize), stat=allocate_error )
        allocate( blocks(i)%data2(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )
        allocate( blocks(i)%data_old(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )
        allocate( blocks(i)%k1(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )
        allocate( blocks(i)%k2(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )
        allocate( blocks(i)%k3(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )
        allocate( blocks(i)%k4(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )
        allocate( blocks(i)%treecode(10), stat=allocate_error )
        allocate( blocks(i)%neighbor_treecode(8,10), stat=allocate_error )
        allocate( blocks(i)%neighbor2_treecode(4,10), stat=allocate_error )
    end do

    blocks_params%size_domain		    = nx
    blocks_params%size_block		    = blocksize
    blocks_params%number_max_blocks = max_blocks
    blocks_params%number_ghost_nodes	= ghosts

    !------------------------------
    ! start field, blocks
    blocks_params%phi 			            = phi

    do i = 1, max_blocks
        blocks(i)%data1(:,:)                = 0.0_rk
        blocks(i)%data2(:,:)                = 0.0_rk
        blocks(i)%data_old(:,:)             = 0.0_rk
        blocks(i)%k1(:,:)                   = 0.0_rk
        blocks(i)%k2(:,:)                   = 0.0_rk
        blocks(i)%k3(:,:)                   = 0.0_rk
        blocks(i)%k4(:,:)                   = 0.0_rk
        blocks(i)%active                    = .false.
        blocks(i)%treecode(:)               = -1
        blocks(i)%neighbor_treecode(:,:)    = -1
        blocks(i)%neighbor2_treecode(:,:)   = -1
        blocks(i)%neighbor_id(:)            = -1
        blocks(i)%neighbor2_id(:)           = -1
        blocks(i)%refinement                = 0
        blocks(i)%level                     = -1
        blocks(i)%neighbor_number(:)        = 1
    end do

    !------------------------------
    ! derivatives
    allocate( blocks_params%D1(blocksize+2*ghosts,blocksize+2*ghosts), stat=allocate_error )
    allocate( blocks_params%D2(blocksize+2*ghosts,blocksize+2*ghosts), stat=allocate_error )

    call D18j(D1, blocksize+2*ghosts, 1.0_rk)
    call D26p(D2, blocksize+2*ghosts, 1.0_rk)

    blocks_params%D1 			        = D1
    blocks_params%D2			        = D2

    !------------------------------
    ! time loop parameter
    params%time_max             = 100.0_rk
    params%CFL 		            = 0.5_rk

    !------------------------------
    ! convective velocity
    params%u0 		            = (/1.0_rk, 0.5_rk/)

    !------------------------------
    ! diffusion coeffcient
    params%nu 		            = 1e-2_rk

    !------------------------------
    ! domain size in [m]
    params%Lx 		            = 256.0_rk
    params%Ly                   = params%Lx

    !------------------------------
    ! workdir, case name, write frequency
    params%name_workdir 	    = "./data/"
    params%name_case 	        = "nghosts4"
    params%write_freq	        =  20

    !------------------------------
    ! spacing

    do i = 1, max_blocks
        blocks(i)%dx    	    = params%Lx / real(nx,8)
        blocks(i)%dy		    = params%Lx / real(nx,8)
    end do

    !------------------------------
    ! coordinates
    do i = 1, max_blocks
        allocate( blocks(i)%coord_x(blocksize) )
        allocate( blocks(i)%coord_y(blocksize) )
    end do

    do i = 1, max_blocks
        blocks(i)%coord_x(:)    = 0.0_rk
        blocks(i)%coord_y(:)    = 0.0_rk
    end do

    !------------------------------
    ! deallocate memory for local variables
    deallocate( phi, stat=allocate_error )
    deallocate( D1, stat=allocate_error )
    deallocate( D2, stat=allocate_error )

    !------------------------------
    ! eps for coarsen and refine the block
    params%eps_coarsen          = 1e-3_rk
    params%eps_refine           = 1e-2_rk

    !------------------------------
    ! set treelevel
    params%max_treelevel        = 5
    params%min_treelevel        = 1

end subroutine init_data
