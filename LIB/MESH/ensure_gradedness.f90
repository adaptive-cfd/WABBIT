! ********************************
! WABBIT
! --------------------------------
!
! check the gradedness after new
! refinement status
!
! name: ensure_gradedness.f90
! date: 28.10.2016
! author: msr, engels
! version: 0.3
!
! ********************************

subroutine ensure_gradedness()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik) :: myid, l, neighbor_id, mylevel, neighbor_level
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
      ! loop over all blocks, use only active ones
      !-------------------------------------------------------------------------
      do myid = 1, blocks_params%number_max_blocks
          ! is this block in use?
          if (blocks(myid)%active) then

              !-----------------------------------------------------------------------
              ! block wants to coarsen, refinement only in safety zone step
              !-----------------------------------------------------------------------
              if (blocks(myid)%refinement == -1) then

                  ! loop over all directions
                  do l = 1, 16

                  if ( blocks(myid)%neighbor_id(l) /= -1 ) then

                      ! active neighbor
                      ! check neighbor treelevel
                      neighbor_id = blocks(myid)%neighbor_id(l)
                      mylevel     = blocks(myid)%level

                      if (neighbor_id == -1) then
                        ! error case
                        print*, 'error: gradedness error, can not coarsen - neighbor does not exists'
                        stop
                      end if

                      neighbor_level  = blocks(neighbor_id)%level

                      if (mylevel == neighbor_level) then
                          ! neighbor on same level
                          ! block can not coarsen, if neighbor wants to refine
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

                  end if

                  end do

              end if

          end if

      end do
      counter = counter + 1
      if (counter == 99) then
        print*, "error: unable to build a graded mesh"
        stop
      end if

    ! end do of repeat procedure until grid_changed==.false.
    end do

end subroutine ensure_gradedness
