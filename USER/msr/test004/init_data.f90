! ********************************
! 2D AMR prototype
! --------------------------------
!
! initialize all data (params, fields, blocks, ...)
!
! name: init_data.f90
! date: 02.08.2016
! author: msr, engels
! version: 0.1
!
! ********************************

subroutine init_data()

    use module_params
    use module_blocks

    implicit none

    integer :: allocate_error

    ! grid size
    integer(kind=ik) :: nx
    ! block size,
    integer(kind=ik) :: blocksize
    integer(kind=ik) :: max_blocks
    integer(kind=ik) :: ghosts

    !------------------------------
    ! grid and block parameter
    nx          = 513 !512
    blocksize   = 17 !33 !65 !129 !257
    ! allocate 4x as many blocks as required
    ! TODO: allocate max level full
    max_blocks = 4*((nx-1)/(blocksize-1))**2
    ghosts     = 4

    blocks_params%size_domain		    = nx
    blocks_params%size_block		    = blocksize
    blocks_params%number_max_blocks = max_blocks
    blocks_params%number_ghost_nodes	= ghosts

    write(*,'(80("-"))')
    write(*,*) "INITIALIZATION PHASE"
    write(*,*) "we use a maximum number of blocks of", max_blocks
    write(*,*) "nghosts=", ghosts
    write(*,'(80("-"))')

    !***************************************************************************
    ! GENERAL PARAMETERS (STRUCT "PARAMS")
    !***************************************************************************
    ! time loop parameter
    params%time_max = 160.0_rk
    params%CFL 		  = 0.5_rk
    ! convective velocity
    params%u0 		  = (/1.0_rk, 0.5_rk/) !(/0.0_rk, 0.0_rk/)
    ! diffusion coeffcient
    params%nu 		            = 0.0e-3_rk
    ! domain size
    params%Lx 		            = 256.0_rk
    params%Ly                   = params%Lx
    ! workdir, case name, write frequency
    params%name_workdir 	    = "./data/"
    params%name_case 	        = "eps1e-3_level6"
    params%write_freq	        =  25
    ! eps for coarsen and refine the block
    params%eps_coarsen          = 1e-3_rk
    params%eps_refine           = 5.0_rk * params%eps_coarsen
    ! set treelevel
    params%max_treelevel        = 6
    params%min_treelevel        = 1
    params%order_predictor      = "multiresolution_4th" !"multiresolution_2nd"  ! "multiresolution_4th"
    params%order_discretization = "FD_4th_central_optimized"!"FD_2nd_central" ! "FD_4th_central_optimized"

    allocate( blocks_params%active_list(1), stat=allocate_error )
    ! allocate the individual block's memory
    call allocate_block_memory()
    ! initial data field
    call initial_condition_dense_field()
end subroutine init_data
