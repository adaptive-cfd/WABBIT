!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name coarsening_indicator.f90
!> \version 0.5
!> \author engels
!> \brief Set coarsening status for active blocks, different methods possible
!
!> \details This routine sets the coarsening flag for all blocks. We allow for different
!! mathematical methods (everywhere / random) currently not very compley, but expected to grow
!! in the future.
!! \n
!! ------------------ \n
!! Refinement status: \n
!! ------------------ \n
!! +1 refine \n
!! 0 do nothing \n
!! -1 block wants to refine (ignoring other constraints, such as gradedness) \n
!! -2 block will refine and be merged with her sisters \n
!! ------------------ \n
!! \n
!! = log ======================================================================================
!! \n
!! 23/05/2017 create
! ********************************************************************************************

subroutine coarsening_indicator( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n, indicator, iteration )
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n
    !> how to choose blocks for refinement
    character(len=*), intent(in)        :: indicator
    !> coarsening iteration index
    integer(kind=ik), intent(in)        :: iteration

    ! local variables
    integer(kind=ik) :: k, Jmax, d, ierr, j
    ! chance for block refinement, random number
    real(kind=rk) :: ref_chance, r

!---------------------------------------------------------------------------------------------
! variables initialization

        Jmax = params%max_treelevel
        ! set data dimension
        if ( params%threeD_case ) then
            d = 3
        else
            d = 2
        endif

!---------------------------------------------------------------------------------------------
! main body

    ! reset refinement status to "stay"
    do k = 1, lgt_n
      lgt_block( lgt_active(k), Jmax+2 ) = 0
    enddo


    !> (a) loop over the blocks and set their refinement status.
    !! NOTE: refinement is an absolute statement, that means once set, the block will be refined
    !! (which is not the case in block coarsening), it may even entrail other blocks in
    !! its vicinity to be refined as well.
    select case (indicator)
      case ("threshold")
        ! use wavelet indicator to check where to coarsen. threshold_block performs
        ! the required ghost node sync and loops over all active blocks.
        call threshold_block( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )

      case ("random")
        ! randomly coarse some blocks. used for testing. note we tag for coarsening
        ! only once in the first iteration.
        if (iteration == 0) then
          call init_random_seed()
          ref_chance = 0.25_rk
          ! unset all refinement flags
          lgt_block( :,Jmax+2 ) = 0
          ! only root rank sets the flag, then we sync. It is messy if all procs set a
          ! random value which is not sync'ed
          if (params%rank == 0) then
            do j = 1, lgt_n
              ! random number
              call random_number(r)
              ! set refinement status to coarsen
              if ( r <= ref_chance ) then
                  lgt_block( lgt_active(j), Jmax+2 ) = -1
              end if
            end do
          endif
          ! sync light data, as only root sets random coarsening
          call MPI_BCAST( lgt_block(:,params%max_treelevel+2), size(lgt_block,1), MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr )
        endif

      case default
          call error_msg("ERROR: unknown coarsening operator")

    end select

end subroutine coarsening_indicator
