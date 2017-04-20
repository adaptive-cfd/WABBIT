! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: init_data.f90
! version: 0.5
! author: msr
!
! This routine initializes the block data, i.e. it evaluates the initial condition on the grid
!
! input:    - parameter array
!           - light data array
!           - heavy data array
!           - neighbor data array
!           - light and heavy active block list
! output:   - filled user defined data structure for global params
!           - initialized light and heavy data arrays
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4, now run complete initialization within these subroutine and return
!            initialized block data to main program
! 07/12/16 - now uses heavy work data array
! 25/01/17 - switch to 3D, v0.5
!
! ********************************************************************************************
subroutine set_blocks_initial_condition(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, lgt_n, hvy_n, adapt)

  !---------------------------------------------------------------------------------------------
  ! variables

  implicit none

  ! user defined parameter structure
  type (type_params), intent(inout)    :: params
  ! light data array
  integer(kind=ik), intent(inout)      :: lgt_block(:, :)
  ! heavy data array - block data
  real(kind=rk), intent(inout)         :: hvy_block(:, :, :, :, :)
  ! neighbor array (heavy data)
  integer(kind=ik), intent(inout)      :: hvy_neighbor(:,:)
  ! list of active blocks light data)
  integer(kind=ik), intent(inout)      :: lgt_active(:)
  ! list of active blocks (light data)
  integer(kind=ik), intent(inout)      :: hvy_active(:)
  ! number of heavy and light active blocks
  integer(kind=ik), intent(inout)      :: hvy_n, lgt_n
  ! if .false. the code initializes on the coarsest grid, if .true. iterations
  ! are performed and the mesh is refined to gurantee the error eps
  logical, intent(in) :: adapt
  integer(kind=ik)                     :: lgt_n_old

  !---------------------------------------------------------------------------------------------
  ! interfaces

  !---------------------------------------------------------------------------------------------
  ! variables initialization
    lgt_n_old = 9999999

  !---------------------------------------------------------------------------------------------
  ! main body

    !---------------------------------------------------------------------------
    ! Create the first mesh on the coarsest treelevel
    !---------------------------------------------------------------------------
    call create_equidistant_base_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, params%min_treelevel, .true. )

    !---------------------------------------------------------------------------
    ! on the grid, evaluate the initial condition
    !---------------------------------------------------------------------------
    call set_inicond_all_blocks(params, lgt_block, hvy_block, hvy_active, hvy_n, params%initial_cond)

    !---------------------------------------------------------------------------
    ! grid adaptation
    !---------------------------------------------------------------------------
    if (adapt) then
      ! we have to repeat the adapation process until the grid has reached a final
      ! state. Since we start on the coarsest level, in each iteration we cannot loose
      ! blocks, but only gain or no change. Therefore, iterate until lgt_n is constant.
      do while ( lgt_n /= lgt_n_old)
        lgt_n_old = lgt_n
        ! push up the entire grid one level. TODO: It would be better to selectively
        ! go up one level where a refinement indicator tells us to do so, but in teh current code
        ! versions it is easier to use everywhere
        call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, "everywhere"  )

        ! It may seem surprising, but we now have to re-set the inicond on the blocks. if
        ! not, the detail coefficients for all blocks are zero. In the time stepper, this
        ! corresponds to advancing the solution in time, it's just that here we know the exact
        ! solution (the inicond)
        call set_inicond_all_blocks(params, lgt_block, hvy_block, hvy_active, hvy_n, params%initial_cond)

        ! now, evaluate the refinement criterion on each block, and coarsen the grid where possible.
        ! adapt-mesh also performs neighbor and active lists updates
        call adapt_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, "threshold" )

        if (params%rank == 0) then
          write(*,'(" did one mesh adaptation for the initial condition. Nblocks=",i6, " Jmax=",i2)') lgt_n, maxval(lgt_block(:,params%max_treelevel+1))
        endif
      enddo
    endif


  !     ! start time
  !     sub_t0 = MPI_Wtime()
  !
  !     ! initial data field
  !     select case( params%initial_cond )
  !         case ("gauss-blob","gauss_blob")
  !
  !             ! 2D: gaussblob, 3D: zero fields
  !             if ( params%threeD_case ) then
  !                 ! 3D:
  !                 call inicond_sphere( params, lgt_block, hvy_block )
  !             else
  !                 ! 2D:
  !                 call inicond_gauss_blob( params, lgt_block, hvy_block )
  !             end if
  !
  !         case ("vorticity_filaments")
  ! !            call inicond_vorticity_filaments(params, lgt_block, hvy_block)
  ! !            ! set density and pressure
  ! !            hvy_block( :, :, 2, : ) = 1.0_rk
  ! !            hvy_block( :, :, 5, : ) = 1e5_rk
  !
  !         case ("richtmyer_meshkov")
  ! !            call inicond_richtmyer_meshkov(params, lgt_block, hvy_block)
  !
  !         case ("shear_layer")
  !             call inicond_shear_layer(params, lgt_block, hvy_block)
  !
  !         case default
  !             write(*,'(80("_"))')
  !             write(*,*) "ERROR: initial condition is unknown"
  !             write(*,*) params%initial_cond
  !             stop
  !
  !     end select
  !
  !
  !     ! end time
  !     sub_t1 = MPI_Wtime()
  !     ! write time
  !     if ( params%debug ) then
  !         ! find free or corresponding line
  !         k = 1
  !         do while ( debug%name_comp_time(k) /= "---" )
  !             ! entry for current subroutine exists
  !             if ( debug%name_comp_time(k) == "init_data" ) exit
  !             k = k + 1
  !         end do
  !         ! write time
  !         debug%name_comp_time(k) = "init_data"
  !         debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
  !         debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
  !
  !     end if

end subroutine set_blocks_initial_condition
