
subroutine unit_test_treecode( params )

  use module_treelib

  implicit none
  type (type_params), intent(inout) :: params                   !> user defined parameter structure
  integer(kind=ik), ALLOCATABLE     :: treecode(:), n(:)
  integer(kind=ik)                  :: treeN, ix, iy, iz, k, kcheck, kk, j
  integer(kind=ik)                  :: ixy(3)
  character(len=3), dimension(4)    :: dir_2d
  character(len=3), dimension(6)    :: dir_3d
  integer(kind=tsize)               :: newtreecode, neighbor
  real(kind=rk)                     :: t, rand1, rand2
  logical                           :: array_compare

  Kcheck = 50
  KK = 1000000
  ! KK = 1

  treeN = params%Jmax
  allocate (treecode(1:treeN), n(1:treeN))

  ! define primary sides to be investigated for 2d
  dir_2d(1) = "__N"
  dir_2d(2) = "__S"
  dir_2d(3) = "__W"
  dir_2d(4) = "__E"

  ! define primary sides to be investigated for 3d
  dir_3d(1) = "__1"
  dir_3d(2) = "__2"
  dir_3d(3) = "__3"
  dir_3d(4) = "__4"
  dir_3d(5) = "__5"
  dir_3d(6) = "__6"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (params%rank == 0) then
    write(*,'(80("_"))')
    write(*,'("UNIT TEST: Beginning treecode test")')
    write(*,*) "UNIT TEST: treecode - checking all directions in 2D"

    ! test several times if indices increase correctly
    do k = 1, kcheck
      !> start coordinates random between 20 and 100
      call random_number(rand1)
      call random_number(rand2)
      ixy(1) = 20 + floor(rand1 * 80)
      ixy(2) = 20 + floor(rand2 * 80)
      call encoding_n( ixy, 2, newtreecode)
      call encoding(treecode, ixy, 2, 4**treeN, treeN)

      ! tests for PK version
      ! loop over all directions - in if condition we set the correct neighbouring index change
      do j= 1,4
        call adjacent4_NESW( newtreecode, dir_2d(j), neighbor)
        call decoding_n(neighbor, ix,iy,iz)
        if ( ix /= ixy(2) + (1 - (j-1)/2) * (-1 + 2*modulo((j-1), 2)) .or. iy /= ixy(1) + (j-1)/2 * (-1 + 2*modulo((j-1), 2))) then
          write(*,*) "UNIT TEST FAILED: treecode num ", dir_2d(j), " - orig ", newtreecode,", sx=", ixy(2), "sy=", ixy(1), " ", neighbor,", ix=", ix, "iy=", iy
          call abort(123980)
        endif
      end do

      ! tests for JB version
      ! loop over all directions - in if condition we set the correct neighbouring index change
      do j= 1,4
        call adjacent_JB_NESW( newtreecode, neighbor, dir_2d(j), params%Jmax, params%Jmax)
        call decoding_n(neighbor, ix,iy,iz)
        if ( ix /= ixy(2) + (1 - (j-1)/2) * (-1 + 2*modulo((j-1), 2)) .or. iy /= ixy(1) + (j-1)/2 * (-1 + 2*modulo((j-1), 2))) then
          write(*,*) "UNIT TEST FAILED: treecode num JB ", dir_2d(j), " - orig ", newtreecode,", sx=", ixy(2), "sy=", ixy(1), " ", neighbor,", ix=", ix, "iy=", iy
          call abort(123980)
        endif
      end do

      ! tests for old version
      ! loop over all directions - in if condition we set the correct neighbouring index change
      do j= 1,4
        call adjacent_block_2D( treecode, n, dir_2d(j), params%Jmax, params%Jmax)
        call decoding(n, ix,iy,iz, treeN)
        if ( ix /= ixy(2) + (1 - (j-1)/2) * (-1 + 2*modulo((j-1), 2)) .or. iy /= ixy(1) + (j-1)/2 * (-1 + 2*modulo((j-1), 2))) then
          write(*,*) "UNIT TEST FAILED: treecode num old ", dir_2d(j), " - orig sx=", ixy(2), "sy=", ixy(1), "neighbour ix=", ix, "iy=", iy
          call abort(123980)
        endif
      end do
    end do
    write(*,*) "UNIT TEST: treecode - Bravooo, 2D leaves don't have to be lonely and can find their neighbours"

    ! write(*,*) "UNIT TEST: treecode - checking all directions in 3D"

    ! ! test several times if indices increase correctly
    ! do k = 1, kcheck
    !   !> start coordinates random between 20 and 100
    !   call random_number(rand1)
    !   call random_number(rand2)
    !   ixy(1) = 20 + floor(rand1 * 80)
    !   ixy(2) = 20 + floor(rand2 * 80)
    !   call encoding_n( ixy, 2, newtreecode)
    !   call encoding(treecode, ixy, 2, 4**treeN, treeN)

    !   ! tests for JB version
    !   ! loop over all directions - in if condition we set the correct neighbouring index change
    !   do j= 1,4
    !     call adjacent_JB_NESW( newtreecode, neighbor, dir_2d(j), params%Jmax, params%Jmax)
    !     call decoding_n(neighbor, ix,iy,iz)
    !     if ( ix /= ixy(2) + (1 - (j-1)/2) * (-1 + 2*modulo((j-1), 2)) .or. iy /= ixy(1) + (j-1)/2 * (-1 + 2*modulo((j-1), 2))) then
    !       write(*,*) "UNIT TEST FAILED: treecode num JB ", dir_2d(j), " - orig ", newtreecode,", sx=", ixy(2), "sy=", ixy(1), " ", neighbor,", ix=", ix, "iy=", iy
    !       call abort(123980)
    !     endif
    !   end do

    !   ! tests for old version
    !   ! loop over all directions - in if condition we set the correct neighbouring index change
    !   do j= 1,4
    !     call adjacent_block_3D( treecode, n, dir_2d(j), params%Jmax, params%Jmax)
    !     call decoding(n, ix,iy,iz, treeN)
    !     if ( ix /= ixy(2) + (1 - (j-1)/2)*(2 - (j-1)/2)/2 * (-1 + 2*modulo((j-1), 2))
    !     .or. iy /= ixy(1) + (j-1)/2*(2 - (j-1)/2)/2 * (-1 + 2*modulo((j-1), 2))
    !     .or. iz /= ixy(3) + (j-1)/2*(1 - (j-1)/2) * (-1 + 2*modulo((j-1), 2))) then
    !       write(*,*) "UNIT TEST FAILED: treecode num old ", dir_2d(j), " - orig sx=", ixy(2), "sy=", ixy(1), "neighbour ix=", ix, "iy=", iy
    !       call abort(123980)
    !     endif
    !   end do
    ! end do
    ! write(*,*) "UNIT TEST: treecode - Bravooo, 3D leaves don't have to be lonely and can find their neighbours"


    write(*,*) "===new lib==="
    write(*,*) "UNIT TEST: treecode - measuring time"
  end if

  !> reset number to set equal conditions
  ixy(1) = 77
  ixy(2) = 42
  call encoding_n( ixy, 2, newtreecode)
  call encoding(treecode, ixy, 2, 4**treeN, treeN)
  
  t = MPI_wtime()
  do k = 1, kk
    call encoding_n(ixy, 2, newtreecode)
    ! call encoding4( ix, iy, newtreecode)
    call decoding_n(newtreecode, ix,iy,iz)
  enddo
  write(*,*) "treecode numerical encoding: rank=", params%rank, " elapsed=", MPI_wtime()-t

  t = MPI_wtime()
  do k = 1, kk
    ! call adjacent4_NESW( newtreecode, "__N", neighbor)
    ! call adjacent4_NESW( neighbor, "__S", newtreecode)
    ! call adjacent4_NESW( newtreecode, "__W", neighbor)
    ! call adjacent4_NESW( neighbor, "__E", newtreecode)

    call adjacent_JB_NESW(newtreecode, neighbor, "__N", params%Jmax, params%Jmax)
    call adjacent_JB_NESW(neighbor, newtreecode, "__S", params%Jmax, params%Jmax)
    call adjacent_JB_NESW(newtreecode, neighbor, "__W", params%Jmax, params%Jmax)
    call adjacent_JB_NESW(neighbor, newtreecode, "__E", params%Jmax, params%Jmax)

    ! loop of 1000 simulates searching loop for neighbour
    ! do iz = 1, 1000
    if ( newtreecode == neighbor ) then
      call abort(123980)
    endif
    ! enddo
  enddo
  write(*,*) "treecode numerical adjacent: rank=", params%rank, " elapsed=", MPI_wtime()-t


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (params%rank == 0) then
    write(*,'(80("-"))')
    write(*,*) "===old lib==="
    write(*,*) "UNIT TEST: treecode - measuring time"
  end if

  t = MPI_wtime()
  do k = 1, kk
    call encoding(treecode, ixy, 2, 4**treeN, treeN)
    call decoding(treecode, ix, iy, iz, treeN)
  enddo
  write(*,*) "treecode array encoding: rank=", params%rank, " elapsed=", MPI_wtime()-t

  t = MPI_wtime()
  do k = 1, kk
    call adjacent_block_2D(treecode, n, "__N", params%Jmax, params%Jmax)
    call adjacent_block_2D(n, treecode, "__S", params%Jmax, params%Jmax)
    call adjacent_block_2D(treecode, n, "__W", params%Jmax, params%Jmax)
    call adjacent_block_2D(n, treecode, "__E", params%Jmax, params%Jmax)

    ! loop of 1000 simulates searching loop for neighbour
    ! do iz = 1,1000
    if (array_compare(treecode, n, params%Jmax)) then
      call abort(123980)
    endif
    ! enddo
  enddo
  write(*,*) "treecode array adjacent: rank=", params%rank, " elapsed=", MPI_wtime()-t

end subroutine unit_test_treecode
