! ********************************
! 2D AMR prototype
! --------------------------------
!
! check the gradedness after new
! refinement status
!
! name: ensure_gradedness.f90
! date: 18.08.2016
! author: msr, engels
! version: 0.2
!
! ********************************

subroutine ensure_gradedness()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik) :: myid, l, neighbor_id, mylevel, neighbor_level, neighbor2_id
    logical :: grid_changed
    integer :: counter

    ! we repeat the ensure_gradedness procedure until this flag is .false. since as long
    ! as the grid changes due to gradedness requirements, we have to check it again
    grid_changed = .true. ! set true to tigger the loop
    counter = 0

    do while ( grid_changed )
      ! we hope not to set the flag to .true. again in this iteration
      grid_changed = .false.

      !-------------------------------------------------------------------------
      ! loop over all active blocks
      !-------------------------------------------------------------------------
      do myid = 1, blocks_params%number_max_blocks
          ! is this block in use?
          if (blocks(myid)%active) then

          !-----------------------------------------------------------------------
          ! block wants to refine
          !-----------------------------------------------------------------------
          if (blocks(myid)%refinement == 1) then

              ! loop over all directions
              do l = 1, 8

                  ! check number of neighbors, first case: exactly 1 neighbor
                  if (blocks(myid)%neighbor_number(l) == 1) then

                      mylevel     = blocks(myid)%level
                      ! check neighbor treelevel
                      neighbor_id = blocks(myid)%neighbor_id(l)

                      if (neighbor_id == -1) then
                          print*, 'error: gradedness error, can not refine - neighbor does not exists'
                          stop
                      end if

                      neighbor_level = blocks(neighbor_id)%level

                      if (mylevel == neighbor_level) then
                          ! neighbor on same level
                          ! neighbor can not coarsen, because I want to refine
                          if (blocks(neighbor_id)%refinement == -1) then
                            blocks(neighbor_id)%refinement = 0
                            ! if the neighbor wants to stay or refine, that is fine with us
                            grid_changed = .true. ! we changed something in the grid
                          end if

                      elseif (mylevel - neighbor_level == 1) then
                          ! neighbor one level lower (=coarser)
                          ! neighbor has to refine, because I want to refine
                          if (blocks(neighbor_id)%refinement /= 1) then
                            ! note this branch takes only place if the neighbor does not want
                            ! to refine anyways already
                            blocks(neighbor_id)%refinement = 1

                            ! this inherited (involuntary) refinement has important consequences: we
                            ! have to be sure that this block also respects gradedness
                            ! after its refinement. therefore, we have to repeat the entire process
                            grid_changed = .true.
                          endif

                      elseif (mylevel - neighbor_level == -1) then
                          ! neighbor one level higher (=finer), and I want to refine
                          ! myself. My neighbor can
                          ! a) refine, b) stay, c) coarsen but in any of these cases
                          ! the level difference will be at most 1, so we have nothing
                          ! to do (this case is active for edges)
                      else
                          ! error
                          print*, 'error: neighbor more than one level up/down in tree'
                          stop

                      end if

                  ! second case: 2 neighbors
                  elseif (blocks(myid)%neighbor_number(l) == 2) then
                      ! two neighbors imply that they are one level finer than I am
                      ! neighbor one level higher (=finer), and I want to refine
                      ! myself. My neighbor can
                      ! a) refine, b) stay, c) coarsen but in any of these cases
                      ! the level difference will be at most 1, so we have nothing
                      ! to do

                  ! error case
                  else
                      print*, "error: block has wrong number of neighbors!"
                      stop

                  end if
              end do
          !-----------------------------------------------------------------------
          ! block wants to coarsen
          !-----------------------------------------------------------------------
          elseif (blocks(myid)%refinement == -1) then

              ! loop over all directions
              do l = 1, 8

                  ! check number of neighbors, first case: exactly 1 neighbor
                  if (blocks(myid)%neighbor_number(l) == 1) then

                      ! check neighbor treelevel
                      neighbor_id = blocks(myid)%neighbor_id(l)
                      mylevel     = blocks(myid)%level

                      if (neighbor_id == -1) then
                        ! error case
                        write(*,*) myid, "myid"
                        blocks(myid)%data2 = 9e9_rk
                        blocks(myid)%data1 = 9e9_rk
                        call save_data(9999, 9.999_rk)
                        print*, 'error: gradedness error, can not coarsen - neighbor does not exists'
                        stop
                      end if

                      neighbor_level  = blocks(neighbor_id)%level

                      if (mylevel == neighbor_level) then
                          ! neighbor on same level
                          ! block can not coarsen, because neighbor wants to refine
                          if (blocks(neighbor_id)%refinement == 1) then
                            blocks(myid)%refinement = max(0, blocks(myid)%refinement)
                            grid_changed = .true. ! we changed something in the grid
                          end if

                      elseif (mylevel - neighbor_level == 1) then
                          ! neighbor on lower level
                          ! nothing to do

                      elseif (neighbor_level - mylevel == 1) then
                          ! neighbor on higher level
                          ! neighbor want to refine, ...
                          if (blocks(neighbor_id)%refinement == 1) then
                              ! ... so I also have to refine
                              blocks(myid)%refinement = 1
                              grid_changed = .true. ! we changed something in the grid
                          else
                              ! neighbor stays on his level or want to coarsen,
                              ! so block can not coarsen
                              ! TODO: block can coarsen, if neighbor also want to coarsen
                              blocks(myid)%refinement = max(0, blocks(myid)%refinement)
                              grid_changed = .true. ! we changed something in the grid
                          end if

                      else
                          ! error case
                          print*, "error: block has wrong number of neighbors!"
                          stop

                      end if

                  ! second case: 2 neighbors
                  elseif (blocks(myid)%neighbor_number(l) == 2) then
                      ! two neighbor, both are one level higher in tree
                      neighbor_id     = blocks(myid)%neighbor_id(l)

                      if (neighbor_id == -1) then
                          print*, 'error: gradedness error, can not coarsen - neighbor 1 does not exists'
                          stop
                      end if

                      ! look for second neighbor
                      neighbor2_id    = blocks(myid)%neighbor2_id(l)

                      if (neighbor_id == -1) then
                          print*, 'error: gradedness error, can not coarsen - neighbor 2 does not exists'
                          stop
                      end if

                      if ( (blocks(neighbor_id)%refinement == 1) .or. (blocks(neighbor2_id)%refinement == 1) ) then
                          ! one of the neighbors want to refine,
                          ! so block have also refine (and this is sure, so definetly set +1)
                          blocks(myid)%refinement = +1
                          grid_changed = .true. ! we changed something in the grid
                      else
                          ! neighbors stay on tree level or want to coarsen,
                          ! so block can not coarsen
                          ! TODO: block can coarsen, if both neighbors also want to coarsen
                          blocks(myid)%refinement = max( 0,blocks(myid)%refinement )
                          grid_changed = .true. ! we changed something in the grid
                      end if

                  ! error case
                  else
                      print*, "error: block has wrong number of neighbors!"
                      stop
                  end if

              end do

          end if
          end if
      end do
      counter = counter + 1
    ! end do of repeat procedure until grid_changed==.false.
    end do
end subroutine ensure_gradedness
