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

    integer 								                        :: allocate_error

    ! grid size
    integer(kind=ik), parameter 						            :: nx = 513
    integer(kind=ik)							                    :: i, j
    real(kind=rk) 							                        :: mux, muy
    ! block size
    integer(kind=ik), parameter						                :: blocksize = 513
    integer(kind=ik), parameter						                :: ghosts = 4
    ! derivatives
    real(kind=rk), dimension(blocksize+2*ghosts,blocksize+2*ghosts)	:: D1, D2
    ! gauss parameter
    real(kind=rk) 							                        :: sigma = 300.0_rk
    ! 2D grid
    real(kind=rk), dimension(nx,nx) 					            :: phi

    !------------------------------
    ! initial data field
    mux = 0.5_rk * ( real(nx,8) - 1.0_rk )
    muy = 0.5_rk * ( real(nx,8) - 1.0_rk )

    do i = 1, nx
        do j = 1, nx
            phi(i,j) = exp( -( (real(i,8)-mux)*(real(i,8)-mux) + (real(j,8)-muy)*(real(j,8)-muy) )/sigma )
        end do
    end do

    allocate( blocks%test(nx,nx), stat=allocate_error )
    allocate( blocks%test2(nx+2*ghosts,nx+2*ghosts), stat=allocate_error )

    blocks%size_domain		    = nx
    blocks%size_block		    = nx
    blocks%number_ghost_nodes	= ghosts
    blocks%test 			    = phi

    !------------------------------
    ! derivatives
    allocate( blocks%D1(blocksize+2*ghosts,blocksize+2*ghosts), stat=allocate_error )
    allocate( blocks%D2(blocksize+2*ghosts,blocksize+2*ghosts), stat=allocate_error )

    call D18j(D1, blocksize+2*ghosts, 1.0_rk)
    call D26p(D2, blocksize+2*ghosts, 1.0_rk)

    blocks%D1 			        = D1
    blocks%D2			        = D2

    !------------------------------
    ! time loop parameter
    params%time_max             = 250.0_rk
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

    !------------------------------
    ! workdir, case name, write frequency
    params%name_workdir 	    = "/work/sroka/2D_AMR_Proto/"
    params%name_case 	        = "001Test"
    params%write_freq	        =  5

    !------------------------------
    ! spacing
    blocks%dx 			        = params%Lx / real(nx,8)
    blocks%dy			        = params%Lx / real(nx,8)

end
