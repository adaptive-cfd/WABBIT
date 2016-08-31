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

    !***************************************************************************
    ! BLOCK STRUCTURE PARAMETERS
    !***************************************************************************
    blocks_params%size_domain		    = 513 !512
    blocks_params%size_block		    = 17 !33 !65 !129 !257
    ! allocate 4x as many blocks as required
    blocks_params%number_max_blocks = 4*((blocks_params%size_domain-1)/(blocks_params%size_block-1))**2
    blocks_params%number_ghost_nodes	= 4

    !***************************************************************************
    ! GENERAL PARAMETERS (STRUCT "PARAMS")
    !***************************************************************************
    ! time loop parameter
    params%time_max = 200.0_rk
    params%CFL 		  = 0.5_rk
    ! convective velocity
    params%u0 		  = (/1.0_rk, 0.5_rk/) !(/0.0_rk, 0.0_rk/)
    ! diffusion coeffcient
    params%nu 		            = 0.0e-3_rk
    ! domain size
    params%Lx 		            = 256.0_rk
    params%Ly                 = params%Lx
    ! output write frequency
    params%write_freq	        =  25
    ! eps for coarsen and refine the block
    params%eps_coarsen          = 1e-3_rk
    params%eps_refine           = 5.0_rk * params%eps_coarsen
    ! set treelevel
    params%max_treelevel        = 6
    params%min_treelevel        = 1
    params%order_predictor      = "multiresolution_4th" !"multiresolution_2nd"  ! "multiresolution_4th"
    params%order_discretization = "FD_4th_central_optimized"!"FD_2nd_central" ! "FD_4th_central_optimized"
    params%inicond              = "gauss_blob" !"sinus" ! "gauss_blob"

    write(*,'(80("-"))')
    write(*,*) "INITIALIZATION PHASE"
    write(*,*) "we use a maximum number of blocks of:", blocks_params%number_max_blocks
    write(*,*) "nghosts:", blocks_params%number_ghost_nodes
    write(*,'(80("-"))')

    ! allocate the individual block's memory
    call allocate_block_memory()
    ! initial data field
    call initial_condition_dense_field()

end subroutine init_data
