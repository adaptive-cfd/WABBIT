!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================

!
!> = log ======================================================================================
!! \n
!! 04/04/17 - create
!
! ********************************************************************************************

subroutine unit_test_wavelet_compression( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active )

!---------------------------------------------------------------------------------------------
! modules
! IO module
use module_IO

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> user defined parameter structure
    type (type_params), intent(inout)     :: params
    !> light data array
    integer(kind=ik),  intent(inout)      :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk),  intent(inout)         :: hvy_block(:, :, :, :, :)
    !> heavy work array
    real(kind=rk),  intent(inout)         :: hvy_work (:, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik),  intent(inout)      :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)      :: lgt_active(:)
    !> list of active blocks (light data)
    integer(kind=ik),  intent(inout)      :: hvy_active(:)

    ! number of active blocks (heavy data)
    integer(kind=ik)                        :: hvy_n
    ! number of active blocks (light data)
    integer(kind=ik)                        :: lgt_n
    ! loop variables
    integer(kind=ik)                        :: k, l, lgt_id, hvy_id

    ! process rank
    integer(kind=ik)                        :: rank, ierr

    ! spacing
    real(kind=rk)                           :: eeps(1:6)

    real(kind=rk)                           :: Lx, Ly, Lz
    ! data dimensionality
    integer(kind=ik)                        :: d

    ! error variable
    real(kind=rk)                           :: error_l2, my_error_l2, norm_l2, my_norm_l2
    real(kind=rk)                           :: error_linf, my_error_linf, norm_linf, my_norm_linf

    ! MPI error variable
    integer(kind=ik)                        :: Jmax
    ! origin and spacing of blocks
    real(kind=rk)                           :: x0(1:3), dx(1:3)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank = params%rank

    ! grid parameter
    Lx = params%Lx
    Ly = params%Ly
    Lz = params%Lz

    ! set data dimension
    if ( params%threeD_case ) then
        d = 3
    else
        d = 2
    endif


!---------------------------------------------------------------------------------------------
! main body

    ! largest equidistant grid which we can allocate in the memory is
    Jmax = floor( log(dble(params%number_procs*params%number_blocks)**(1.0_rk/dble(d)))/log(2.0_rk) )
    Jmax = min( Jmax, params%max_treelevel)

    if (rank == 0) then
        write(*,'(80("_"))')
        write(*,'("UNIT TEST: Beginning test")')
        write(*,*) Jmax
    end if

    ! start on the finest level Jmax' which we fit in the memory  and initialize the gauss
    ! blob thereon.
    ! loop over eps
    !     coarsen until grid has reached steady state
    !     refine back to the original level Jmax'
    !     evaluate error on this level
eeps= (/1.0e-1_rk, 1.0e-2_rk, 1.0e-3_rk, 1.0e-4_rk, 1.0e-5_rk, 0.0_rk/)
    do l = 1,5

    call reset_grid( params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, .true. )
    ! setup the coarsest grid level with some data (we don't care what data, we'll erase it)
    call create_equidistant_base_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, Jmax, .true. )
    ! call create_equidistant_base_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, l )

    ! FIXME: always set gauss blob
    call set_inicond_all_blocks( params, lgt_block, hvy_block, hvy_active, hvy_n, "sin2d")

    params%eps = eeps(l)



    ! now, evaluate the refinement criterion on each block, and coarsen the grid where possible.
    ! adapt-mesh also performs neighbor and active lists updates. Adapt mesh also contains a
    ! loop to coarsen until no more coarsening is possible
    call adapt_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n , "threshold" )
    if (params%rank == 0) then
      write(*,'("Coarsening done. Nblocks=",i6, " Jmax=",i2)') lgt_n, maxval(lgt_block(:,params%max_treelevel+1))
    endif

! call save_data( 1, 1.0_rk, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n )

    ! refine until grid is full again
    do while (lgt_n < (2**Jmax)**d)
        ! push up the entire grid one level
        call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, "everywhere"  )

        if (params%rank == 0) then
          write(*,'("did one refinement for test Nblocks=",i6, " Jmax=",i2)') lgt_n, maxval(lgt_block(:,params%max_treelevel+1))
        endif
    end do

    ! error eval
    my_error_l2 = 0.0_rk
    my_error_linf = 0.0_rk
    my_norm_l2 = 0.0_rk
    my_norm_linf = 0.0_rk

! call save_data( 10, 10.0_rk, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n )

    ! loop over my active heavy data
    do k = 1, hvy_n
      ! hvy_id of the block we're looking at
      hvy_id = hvy_active(k)
      ! light id of this block
      call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
      ! compute block spacing and origin from treecode
      call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
      ! set the initial condition on this block
      ! copy result of upsampling
      hvy_work(:,:,:,1,hvy_id) = hvy_block(:,:,:,2,hvy_id)
      call initial_condition_on_block_wrapper( params, hvy_block(:,:,:,:,hvy_id), x0, dx, "sin2d" )

      my_error_l2 = my_error_l2 + sqrt(sum( (hvy_work(:,:,:,1,hvy_id) - hvy_block(:,:,:,2,hvy_id))**2  ))
      my_error_linf = max( my_error_linf,  maxval( (hvy_work(:,:,:,1,hvy_id) - hvy_block(:,:,:,2,hvy_id)) ) )


      my_norm_l2 = my_norm_l2 + sqrt( sum( (hvy_block(:,:,:,2,hvy_id))**2 ) )
      my_norm_linf = max( my_norm_linf, maxval( abs((hvy_block(:,:,:,2,hvy_id))) ) )

    enddo

! call save_data( 100, 100.0_rk, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n )

    ! synchronize errors
    call MPI_Allreduce(my_error_l2, error_l2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(my_error_linf, error_linf, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(my_norm_l2, norm_l2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(my_norm_linf, norm_linf, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (params%rank==0) then
      ! write(*,*) params%eps, error_l2 / norm_l2, error_linf / norm_linf
      write(*,*) params%eps, error_l2 , error_linf, norm_l2, norm_linf
    endif

  enddo

end subroutine unit_test_wavelet_compression
