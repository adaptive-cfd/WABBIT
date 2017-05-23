!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name refinement_indicator.f90
!> \version 0.5
!> \author engels
!> \brief Set refinement status for active blocks, different methods possible
!
!> \details This routine sets the refinement flag for all blocks. We allow for different
!! mathematical methods (everywhere / random) currently not very compley, but expected to grow
!! in the future.
!! \n
!! ------------------
!! Refinement status:
!! ------------------
!! +1 refine
!! 0 do nothing
!! -1 block wants to refine (ignoring other constraints, such as gradedness)
!! -2 block will refine and be merged with her sisters
!! ------------------
!! \n
!! = log ======================================================================================
!! \n
!! 23/05/2017 create
! ********************************************************************************************

subroutine refinement_indicator( params, lgt_block, hvy_block, lgt_active, lgt_n, indicator )
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> how to choose blocks for refinement
    character(len=*), intent(in)        :: indicator

    ! local variables
    integer(kind=ik) :: k, Jmax, max_blocks, d, ierr
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
        case ("everywhere")
          ! set status "refine" for all active blocks, which is just setting the
          ! last index in the light data block list to +1. This indicator is used
          ! to refine the entire mesh at the beginning of a time step, if error
          ! control is desired.
          do k = 1, lgt_n
              lgt_block( lgt_active(k), params%max_treelevel+2 ) = +1
          end do

      case ("random")
          ! randomized refinement. This can be used to generate debug meshes for
          ! testing purposes. For example the unit tests use that
          ref_chance = 0.25_rk
          ! random refinement can set at most this many blocks to refine (avoid errors
          ! sue to insufficient memory) (since we already have lgt_n blocks we can set the status
          ! at most for Nmax-lgt_n blocks, whcih genrate ech 2**d new blocks)
          max_blocks = (size(lgt_block,1)-lgt_n-10) / 2**d
          ! set random seed
          call init_random_seed()
          ! unset all refinement flags
          lgt_block( :,Jmax+2 ) = 0
          ! only root rank sets the flag, then we sync. It is messy if all procs set a
          ! random value which is not sync'ed
          if (params%rank == 0) then
            do k = 1, lgt_n
              ! random number
              call random_number(r)
              ! set refinement status to refine
              if ( r <= ref_chance .and. sum(lgt_block( :,Jmax+2 )) <= max_blocks) then
                  lgt_block( lgt_active(k), Jmax+2 ) = 1
              else
                  lgt_block( lgt_active(k), Jmax+2 ) = 0
              end if
            end do
          endif
          ! sync light data, as only root sets random refinement
          call MPI_BCAST( lgt_block(:,params%max_treelevel+2), size(lgt_block,1), MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr )

      case default
          call error_msg("ERROR: refine_mesh: the refinement indicator is unkown")

    end select

end subroutine refinement_indicator
