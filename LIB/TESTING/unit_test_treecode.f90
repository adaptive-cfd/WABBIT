
subroutine unit_test_treecode( params, hvy_block, hvy_work, hvy_tmp, tree_ID, abort_on_fail)

  implicit none
  !> user defined parameter structure
  type (type_params), intent(inout) :: params                   
  !> heavy data array - block data
  real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)   
  !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
  real(kind=rk), intent(out)              :: hvy_tmp(:, :, :, :, :)
  !> heavy work array: used for RHS evaluation in multistep methods (like RK4: u0, k1, k2 etc)
  real(kind=rk), intent(out)              :: hvy_work(:, :, :, :, :, :)
  !> Tree_id for test, should be flow_id
  integer(kind=ik), intent(in)            :: tree_ID
  logical, intent(in)                     :: abort_on_fail

  integer(kind=ik), ALLOCATABLE     :: treecode(:), n(:)
  integer(kind=ik)                  :: treeN, ix, iy, iz, k, kcheck, kk, j, level
  integer(kind=ik)                  :: ixy(3), oxy(3)
  character(len=3), dimension(4)    :: dir_2d
  character(len=7), dimension(6)    :: dir_3d
  integer(kind=ik), dimension(4)    :: dir_2d_int
  integer(kind=ik), dimension(6)    :: dir_3d_int
  integer(kind=tsize)               :: newtreecode, tc_b, tc_b_n, neighbor
  real(kind=rk)                     :: t, rand(3)
  logical                           :: array_compare, do_timing
  character(len=params%Jmax)        :: tc_str

  if (params%Jmax == 0) then
    if (params%rank == 0) then
      write(*,'(20("_/¯\"))')
      write(*,'("UNIT TEST: Skipping treecode test, as JMax=0. You should consider going back to Matlab for whatever you are trying to do.")')
    endif
    return
  endif

  Kcheck = 50  ! how many random times the direction should be checked
  KK = 1000000  ! How many iterations for measuring time
  do_timing = .false.

  treeN = params%Jmax
  allocate (treecode(1:treeN), n(1:treeN))

  ! define primary sides to be investigated for 2d
  dir_2d(1) = "__N"
  dir_2d(2) = "__S"
  dir_2d(3) = "__W"
  dir_2d(4) = "__E"
  dir_2d_int(1) = 900  ! N - x-1
  dir_2d_int(2) = 100  ! S - x+1
  dir_2d_int(3) = 090  ! W - y-1
  dir_2d_int(4) = 010  ! E - y+1

  ! define primary sides to be investigated for 3d
  dir_3d(6) = "__1/___"  ! top z+1
  dir_3d(3) = "__2/___"  ! front y-1
  dir_3d(2) = "__3/___"  ! right x+1
  dir_3d(4) = "__4/___"  ! back y+1
  dir_3d(1) = "__5/___"  ! left x-1
  dir_3d(5) = "__6/___"  ! bottom z-1
  dir_3d_int(1) = 900  ! left x-1
  dir_3d_int(2) = 100  ! right x+1
  dir_3d_int(3) = 090  ! front y-1
  dir_3d_int(4) = 010  ! back y+1
  dir_3d_int(5) = 009  ! bottom z-1
  dir_3d_int(6) = 001  ! top z+1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (params%rank == 0) then
    write(*,'(20("_/¯\"))')
    write(*,'("UNIT TEST: Beginning treecode test")')
    write(*,'("UNIT TEST: It tests for treecode functions to find the correct neighbors.")')

    ! ! space for initial testing
    ! !> start coordinates random between 20 and 100
    ! call random_number(rand)
    ! ! ixy(:) = 20 + floor(rand(:) * 80)
    ! ixy(:) = (/ 1536, 4069, 6738/)
    ! ixy(:) = 2 + floor(rand(:) * (2**params%Jmax-2))  ! ignore edge cases

    ! call encoding_b(ixy, tc_b, dim=3, level=-1, max_level=params%Jmax)
    ! call decoding_b(oxy, tc_b, dim=3, level=-1, max_level=params%Jmax)
    ! call tc_to_str(tc_b, tc_str, dim=3, level=-1, max_level=params%Jmax)
    ! write(*,'("UNIT TEST: treecode num bin - orig=", 3(i0, 1x), " TC dec=", a, 3(1x, i0), " TC=", b64.64)') ixy, tc_str, oxy, tc_b
    ! call encoding_n(ixy, newtreecode, dim=3, level=-1, max_level=params%Jmax)
    ! call decoding_n(oxy, newtreecode, dim=3, level=-1, max_level=params%Jmax)
    ! write(*,'("UNIT TEST: treecode num dec - orig=", 3(i0, 1x), " TC dec=", i0, 3(1x, i0))') ixy, newtreecode, oxy

    ! level = 4
    ! call adjacent_wrapper_b(tc_b, tc_b_n, "__N", dim=3, level=-1, max_level=params%Jmax)
    ! call decoding_b(oxy, tc_b_n, dim=3, level=-1, max_level=params%Jmax)
    ! call tc_to_str(tc_b_n, tc_str, dim=3, level=-1, max_level=params%Jmax)
    ! write(*,'("UNIT TEST: treecode num bin - orig=", 3(i0, 1x), " TC dec=", a, 3(1x, i0), " N TC=", b64.64)') ixy, tc_str, oxy, tc_b_n

    ! call adjacent_wrapper_n(newtreecode, neighbor, "__N", level=-1, max_level=params%Jmax)
    ! call decoding_n(oxy, neighbor, dim=3, level=-1, max_level=params%Jmax)
    ! write(*,'("UNIT TEST: treecode num dec - orig=", 3(i0, 1x), " TC dec=", i0, 3(1x, i0))') ixy, neighbor, oxy

    ! call abort(123980)


    write(*,'("UNIT TEST: Checking encoding and decoding")')
    ! test several times if indices increase correctly
    do k = 1, kcheck
      !> start coordinates random between 20 and 100
      call random_number(rand)
      ixy(:) = 1 + floor(rand(:) * 2**params%Jmax)

      call encoding_b( ixy, tc_b, 3)
      call decoding_b( oxy, tc_b, 3)
      call tc_to_str(tc_b, tc_str, dim=3, level=-1, max_level=params%Jmax)
      if (ixy(1) /= oxy(1) .or. ixy(2) /= oxy(2) .or. ixy(3) /= oxy(3)) then
        write(*,'("UNIT TEST FAILED: treecode num bin - orig=", 3(i0, 1x), "now=", 3(i0, 1x), "TC=", b64.64, " TC dec=", a)') ixy, oxy, tc_b, tc_str
        if (abort_on_fail) call abort(123980)
      end if

      ! call encoding_n( ixy, newtreecode, 3)
      ! call decoding_n( oxy, newtreecode, 3)
      ! if (ixy(1) /= oxy(1) .or. ixy(2) /= oxy(2) .or. ixy(3) /= oxy(3)) then
      !   write(*,'("UNIT TEST FAILED: treecode dec - orig=", 3(i0, 1x), "now=", 3(i0, 1x), "TC=", i0)') ixy, oxy, newtreecode
      !   if (abort_on_fail) call abort(123980)
      ! end if

      ! call encoding_revised(treecode, ixy, 3, params%Jmax)
      ! call decoding(treecode, ix, iy, iz, treeN)
      ! if (ixy(1) /= iy .or. ixy(2) /= ix .or. ixy(3) /= iz) then
      !   write(*,'("UNIT TEST FAILED: treecode old - orig=", 3(i0, 1x), "now=", 3(i0, 1x), "TC=", 13(i0))') ixy, iy, ix, iz, treecode
      !   if (abort_on_fail) call abort(123980)
      ! end if
    end do
    write(*, '("UNIT TEST: -> Bravooo, Treecodes can be en- and decoded from block positions")')
    write(*,'("UNIT TEST: Checking all cardinal directions in 2D")')

    ! test several times if indices increase correctly
    do k = 1, kcheck
      !> start coordinates random between 20 and 100
      call random_number(rand)
      ixy(:) = 2 + floor(rand(:) * (2**params%Jmax-2))  ! ignore edge cases

      call encoding_b( ixy, tc_b, 2)
      call encoding_d( ixy, newtreecode, 2)
      ! call encoding(treecode, ixy, 2, 4**treeN, treeN)

      ! tests for JB version - binary
      ! loop over all directions - in if condition we set the correct neighbouring index change
      do j= 1,4
        call adjacent_wrapper_b( tc_b, tc_b_n, dir_2d_int(j), -1, 2)
        call decoding_b(oxy, tc_b_n, 2)

        if ( (oxy(1) /= ixy(1) + (1 - (j-1)/2) * (-1 + 2*modulo((j-1), 2))) &
        .or. (oxy(2) /= ixy(2) + (j-1)/2       * (-1 + 2*modulo((j-1), 2)))) then
          write(*,'("UNIT TEST FAILED: treecode num bin ", i3.3, " - orig", 2(1x, i0), " neighbour ", 2(1x, i0))') dir_2d_int(j), ixy(1:2), oxy(1:2)
          if (abort_on_fail) call abort(123980)
        endif
      end do

      ! tests for JB version - decimal
      ! loop over all directions - in if condition we set the correct neighbouring index change
      do j= 1,4
        call adjacent_wrapper_d( newtreecode, neighbor, dir_2d_int(j), params%Jmax, params%Jmax)
        call decoding_d(oxy, neighbor, 2)
        if ( (oxy(1) /= ixy(1) + (1 - (j-1)/2) * (-1 + 2*modulo((j-1), 2))) &
        .or. (oxy(2) /= ixy(2) + (j-1)/2       * (-1 + 2*modulo((j-1), 2)))) then
          write(*,'("UNIT TEST FAILED: treecode num dec ", i3.3, " - orig", 2(1x, i0), " neighbour ", 2(1x, i0))') dir_2d_int(j), ixy(1:2), oxy(1:2)
          if (abort_on_fail) call abort(123980)
        endif
      end do

      ! ! tests for old version
      ! ! loop over all directions - in if condition we set the correct neighbouring index change
      ! do j= 1,4
      !   call adjacent_block_2D( treecode, n, dir_2d(j), params%Jmax, params%Jmax)
      !   call decoding(n, ix,iy,iz, treeN)
      !   if ( ix /= ixy(2) + (1 - (j-1)/2) * (-1 + 2*modulo((j-1), 2)) &
      !   .or. iy /= ixy(1) + (j-1)/2       * (-1 + 2*modulo((j-1), 2))) then
      !     write(*,'("UNIT TEST FAILED: treecode num old ", a3, " - orig", 2(1x, i0), " neighbour ", 2(1x, i0))') dir_2d(j), ixy(1:2), iy, ix
      !     if (abort_on_fail) call abort(123980)
      !   endif
      ! end do
    end do
    write(*,'("UNIT TEST: -> Bravooo, 2D leaves don''t have to be lonely and can find their neighbours")')

    write(*,'("UNIT TEST: Checking all cardinal directions in 3D")')
    ! test several times if indices increase correctly
    do k = 1, kcheck
      !> start coordinates random between 20 and 100
      call random_number(rand)
      ixy(:) = 2 + floor(rand(:) * (2**params%Jmax-2))  ! ignore edge cases

      call encoding_b( ixy, tc_b, 3)
      call encoding_d( ixy, newtreecode, 3)
      call encoding_revised(treecode, ixy, 3, params%Jmax)

      ! tests for numerical binary treecode
      ! loop over all directions - in if condition we set the correct neighbouring index change
      do j= 1,6
        call adjacent_wrapper_b( tc_b, tc_b_n, dir_3d_int(j), -1, 3)
        call decoding_b(oxy, tc_b_n, 3)

        if ( (oxy(1) /= ixy(1) + (1 - (j-1)/2)*(2 - (j-1)/2)/2 * (-1 + 2*modulo((j-1), 2))) &
        .or. (oxy(2) /= ixy(2) + (j-1)/2*(2 - (j-1)/2)         * (-1 + 2*modulo((j-1), 2))) &
        .or. (oxy(3) /= ixy(3) + (j-1)/4*(1 - (j-1)/2)*-1      * (-1 + 2*modulo((j-1), 2)))) then
          write(*,'("UNIT TEST FAILED: treecode num bin 3D ", i3.3, " - orig", 3(1x, i0), " neighbour ", 3(1x, i0))') dir_3d_int(j), ixy, oxy
          if (abort_on_fail) call abort(123980)
      endif
      end do

      ! tests for JB version - decimal
      ! loop over all directions - in if condition we set the correct neighbouring index change
      do j= 1,6
        call adjacent_wrapper_d( newtreecode, neighbor, dir_3d_int(j), params%Jmax, params%Jmax)
        call decoding_d(oxy, neighbor, 3)

        if ( (oxy(1) /= ixy(1) + (1 - (j-1)/2)*(2 - (j-1)/2)/2 * (-1 + 2*modulo((j-1), 2))) &
        .or. (oxy(2) /= ixy(2) + (j-1)/2*(2 - (j-1)/2)         * (-1 + 2*modulo((j-1), 2))) &
        .or. (oxy(3) /= ixy(3) + (j-1)/4*(1 - (j-1)/2)*-1      * (-1 + 2*modulo((j-1), 2)))) then
          write(*,'("UNIT TEST FAILED: treecode num dec 3D ", i3.3, " - orig", 3(1x, i0), " neighbour ", 3(1x, i0))') dir_3d_int(j), ixy, oxy
          if (abort_on_fail) call abort(123980)
        endif
      end do

      ! ! tests for old version
      ! ! loop over all directions - in if condition we set the correct neighbouring index change
      ! do j= 1,6
      !   call adjacent_block_3D( treecode, n, dir_3d(j), params%Jmax, params%Jmax)
      !   call decoding(n, ix,iy,iz, treeN)

      !   if ( (ix /= ixy(2) + (1 - (j-1)/2)*(2 - (j-1)/2)/2 * (-1 + 2*modulo((j-1), 2))) &
      !   .or. (iy /= ixy(1) + (j-1)/2*(2 - (j-1)/2)         * (-1 + 2*modulo((j-1), 2))) &
      !   .or. (iz /= ixy(3) + (j-1)/4*(1 - (j-1)/2)*-1      * (-1 + 2*modulo((j-1), 2)))) then
      !     write(*,'("UNIT TEST FAILED: treecode num old 3D ", a7, " - orig", 3(1x, i0), " neighbour ", 3(1x, i0))') dir_3d(j), ixy, iy, ix, iz
      !     if (abort_on_fail) call abort(123980)
      !   endif
      ! end do
    end do
    write(*,'("UNIT TEST: -> Bravooo, 3D leaves don''t have to be lonely and can find their neighbours")')
  end if

  if (do_timing) then
    if (params%rank == 0) then
      write(*,'(80("─"))')
      write(*,'("UNIT TEST: Measuring time for different treecode implementations with ", i0, " iterations")') kk
      write(*,'("===num bin=== measuring time")')
    end if
    ! reset number to set equal conditions
    if (params%Jmax >= 8) then
      ixy(1) = 79
      ixy(2) = 42
      ixy(3) = 27
    else
      call random_number(rand)
      ixy(:) = 1 + floor(rand(:) * 2**params%Jmax)
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    t = MPI_wtime()
    do k = 1, kk
      call encoding_b(ixy, tc_b, 2)
      call decoding_b(oxy, tc_b, 2)
    enddo
    if (params%rank == 0) then
      write(*,'("2D bin enc: rank=", i0, " elapsed= ", f6.4, " s")') params%rank, MPI_wtime()-t
    end if

    t = MPI_wtime()
    do k = 1, kk
      call adjacent_wrapper_b(tc_b, tc_b_n, dir_2d_int(1), -1, 2)
      call adjacent_wrapper_b(tc_b_n, tc_b, dir_2d_int(2), -1, 2)
      call adjacent_wrapper_b(tc_b, tc_b_n, dir_2d_int(3), -1, 2)
      call adjacent_wrapper_b(tc_b_n, tc_b, dir_2d_int(4), -1, 2)

      ! loop of 1000 simulates searching loop for neighbour
      ! do iz = 1, 1000
      if ( tc_b == tc_b_n ) then
        call abort(123980)
      endif
      ! enddo
    enddo
    if (params%rank == 0) then
      write(*,'("2D bin adj: rank=", i0, " elapsed= ", f6.4, " s")') params%rank, MPI_wtime()-t
    end if

    t = MPI_wtime()
    do k = 1, kk
      call encoding_b(ixy, tc_b, 3)
      call decoding_b(oxy, tc_b, 3)
    enddo
    if (params%rank == 0) then
      write(*,'("3D bin enc: rank=", i0, " elapsed= ", f6.4, " s")') params%rank, MPI_wtime()-t
    end if

    t = MPI_wtime()
    do k = 1, kk
      call adjacent_wrapper_b(tc_b, tc_b_n, dir_3d_int(1), -1, 3)
      call adjacent_wrapper_b(tc_b_n, tc_b, dir_3d_int(2), -1, 3)
      call adjacent_wrapper_b(tc_b, tc_b_n, dir_3d_int(3), -1, 3)
      call adjacent_wrapper_b(tc_b_n, tc_b, dir_3d_int(4), -1, 3)
      call adjacent_wrapper_b(tc_b, tc_b_n, dir_3d_int(5), -1, 3)
      call adjacent_wrapper_b(tc_b_n, tc_b, dir_3d_int(6), -1, 3)

      ! loop of 1000 simulates searching loop for neighbour
      ! do iz = 1, 1000
      if ( tc_b == tc_b_n ) then
        call abort(123980)
      endif
      ! enddo
    enddo
    if (params%rank == 0) then
      write(*,'("3D bin adj: rank=", i0, " elapsed= ", f6.4, " s")') params%rank, MPI_wtime()-t
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (params%rank == 0) then
      write(*,'("===num dec=== measuring time")')
    end if

    t = MPI_wtime()
    do k = 1, kk
      call encoding_d(ixy, newtreecode, 2)
      call decoding_d(oxy, newtreecode, 2)
    enddo
    if (params%rank == 0) then
      write(*,'("2D dec enc: rank=", i0, " elapsed= ", f6.4, " s")') params%rank, MPI_wtime()-t
    end if

    t = MPI_wtime()
    do k = 1, kk
      call adjacent_wrapper_d(newtreecode, neighbor, dir_2d_int(1), params%Jmax, params%Jmax)
      call adjacent_wrapper_d(neighbor, newtreecode, dir_2d_int(2), params%Jmax, params%Jmax)
      call adjacent_wrapper_d(newtreecode, neighbor, dir_2d_int(3), params%Jmax, params%Jmax)
      call adjacent_wrapper_d(neighbor, newtreecode, dir_2d_int(4), params%Jmax, params%Jmax)

      ! loop of 1000 simulates searching loop for neighbour
      ! do iz = 1, 1000
      if ( newtreecode == neighbor ) then
        call abort(123980)
      endif
      ! enddo
    enddo
    if (params%rank == 0) then
      write(*,'("2D dec adj: rank=", i0, " elapsed= ", f6.4, " s")') params%rank, MPI_wtime()-t
    end if

    t = MPI_wtime()
    do k = 1, kk
      call encoding_d(ixy, newtreecode, 3)
      call decoding_d(oxy, newtreecode, 3)
    enddo
    if (params%rank == 0) then
      write(*,'("3D dec enc: rank=", i0, " elapsed= ", f6.4, " s")') params%rank, MPI_wtime()-t
    end if

    t = MPI_wtime()
    do k = 1, kk
      call adjacent_wrapper_d(newtreecode, neighbor, dir_3d_int(1), params%Jmax, params%Jmax)
      call adjacent_wrapper_d(neighbor, newtreecode, dir_3d_int(2), params%Jmax, params%Jmax)
      call adjacent_wrapper_d(newtreecode, neighbor, dir_3d_int(3), params%Jmax, params%Jmax)
      call adjacent_wrapper_d(neighbor, newtreecode, dir_3d_int(4), params%Jmax, params%Jmax)
      call adjacent_wrapper_d(newtreecode, neighbor, dir_3d_int(5), params%Jmax, params%Jmax)
      call adjacent_wrapper_d(neighbor, newtreecode, dir_3d_int(6), params%Jmax, params%Jmax)

      ! loop of 1000 simulates searching loop for neighbour
      ! do iz = 1, 1000
      if ( newtreecode == neighbor ) then
        call abort(123980)
      endif
      ! enddo
    enddo
    if (params%rank == 0) then
      write(*,'("3D dec adj: rank=", i0, " elapsed= ", f6.4, " s")') params%rank, MPI_wtime()-t
    end if


    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! if (params%rank == 0) then
    !   write(*,'("=== array === measuring time")')
    ! end if

    ! t = MPI_wtime()
    ! do k = 1, kk
    !   call encoding(treecode, ixy, 2, 4**treeN, treeN)
    !   call decoding(treecode, ix, iy, iz, treeN)
    ! enddo
    ! if (params%rank == 0) then
    !   write(*,'("2D arr enc: rank=", i0, " elapsed= ", f6.4, " s")') params%rank, MPI_wtime()-t
    ! end if

    ! t = MPI_wtime()
    ! do k = 1, kk
    !   call adjacent_block_2D(treecode, n, dir_2d(1), params%Jmax, params%Jmax)
    !   call adjacent_block_2D(n, treecode, dir_2d(2), params%Jmax, params%Jmax)
    !   call adjacent_block_2D(treecode, n, dir_2d(3), params%Jmax, params%Jmax)
    !   call adjacent_block_2D(n, treecode, dir_2d(4), params%Jmax, params%Jmax)

    !   ! loop of 1000 simulates searching loop for neighbour
    !   ! do iz = 1,1000
    !   if (array_compare(treecode, n, params%Jmax)) then
    !     call abort(123980)
    !   endif
    !   ! enddo
    ! enddo
    ! if (params%rank == 0) then
    !   write(*,'("2D arr adj: rank=", i0, " elapsed= ", f6.4, " s")') params%rank, MPI_wtime()-t
    ! end if

    ! t = MPI_wtime()
    ! do k = 1, kk
    !   call encoding_revised(treecode, ixy, 3, params%Jmax)
    !   call decoding(treecode, ix, iy, iz, treeN)
    ! enddo
    ! if (params%rank == 0) then
    !   write(*,'("3D arr enc: rank=", i0, " elapsed= ", f6.4, " s")') params%rank, MPI_wtime()-t
    ! end if

    ! t = MPI_wtime()
    ! do k = 1, kk
    !   call adjacent_block_3D(treecode, n, dir_3d(1), params%Jmax, params%Jmax)
    !   call adjacent_block_3D(n, treecode, dir_3d(2), params%Jmax, params%Jmax)
    !   call adjacent_block_3D(treecode, n, dir_3d(3), params%Jmax, params%Jmax)
    !   call adjacent_block_3D(n, treecode, dir_3d(4), params%Jmax, params%Jmax)
    !   call adjacent_block_3D(treecode, n, dir_3d(5), params%Jmax, params%Jmax)
    !   call adjacent_block_3D(n, treecode, dir_3d(6), params%Jmax, params%Jmax)

    !   ! loop of 1000 simulates searching loop for neighbour
    !   ! do iz = 1,1000
    !   if (array_compare(treecode, n, params%Jmax)) then
    !     call abort(123980)
    !   endif
    !   ! enddo
    ! enddo
    ! if (params%rank == 0) then
    !   write(*,'("3D arr adj: rank=", i0, " elapsed= ", f6.4, " s")') params%rank, MPI_wtime()-t
    ! end if
  end if


  !> reset number to set equal conditions
  ixy(1) = 79
  ixy(2) = 42
  ixy(3) = 27
  level = -3

  ! call encoding_n( ixy, newtreecode, dim=3, level=level)
  ! call tcd2array(newtreecode, treecode, level=level, max_level=params%Jmax)
  ! write(*,'("treecode num dec", 4(1x, i0), 1x, 13(i0))') ixy, newtreecode, treecode
  ! call array2tcd(newtreecode, treecode, level=level, max_level=params%Jmax)
  ! write(*,'("treecode num dec", 4(1x, i0), 1x, 13(i0))') ixy, newtreecode, treecode

  ! call encoding_b( ixy, tc_b, dim=3, level=level)
  ! call tc_to_str(tc_b, tc_str, dim=3, level=params%Jmax)
  ! call tcb2array(tc_b, treecode, dim=3, level=level, max_level=params%Jmax)
  ! tc_b_n = treearray2bid(treecode, tree_ID=3, dim=3, level=level, max_level=params%Jmax)
  ! write(*,'("treecode num bin", 3(1x, i0), 2(1x, b64.64))') ixy, tc_b, tc_b_n
  ! write(*,'("treecode num bin", 3(1x, i0), 1x, a, 1x, 13(i0))') ixy, tc_str, treecode
  ! call array2tcb(tc_b, treecode, dim=3, level=level, max_level=params%Jmax)
  ! write(*,'("treecode num bin", 3(1x, i0), 1x, a, 1x, 13(i0))') ixy, tc_str, treecode

  ! call tc_to_str(tc_b, tc_str, dim=3, level=params%Jmax)
  ! call decoding_b( oxy, tc_b, dim=3, level=level)
  ! write(*,'("treecode num bin", 3(1x, i0), 1(1x, a), 3(1x, i0))') ixy, tc_str, oxy
  ! call encoding_n( ixy, newtreecode, dim=3, level=level)
  ! call decoding_n( oxy, newtreecode, dim=3, level=level)
  ! write(*,'("treecode num dec", 7(1x, i0))') ixy, newtreecode, oxy

  ! write(*,'("treecode num bin", 3(1x, i0), 2(1x, b64.64))') ixy, tc_b, tcb2id(tc_b, dim=3, tree_ID=23, level=level, max_level=params%Jmax)
  ! write(*,'("treecode num dec", 5(1x, i0))') ixy, newtreecode, tcd2id(newtreecode, tree_ID=23, level=level, max_level=params%Jmax)
  

  if (params%rank == 0) then
    WRITE(*,'("UNIT TEST: Finished all treecode tests")')
  endif

end subroutine unit_test_treecode
