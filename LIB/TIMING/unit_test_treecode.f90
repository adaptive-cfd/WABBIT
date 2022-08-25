subroutine unit_test_treecode( params )

use module_treelib

    implicit none
    type (type_params), intent(inout) :: params                   !> user defined parameter structure
    integer(kind=ik), ALLOCATABLE     ::  treecode(:), n(:)
    integer(kind=ik)                  :: treeN, ix, iy, iz, k, kk
    integer(kind=tsize)               :: newtreecode, neighbor
    real(kind=rk)                     :: t
    logical                           :: array_compare

    ix = 77
    iy = 69

    KK = 100000

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "===new lib==="
    t = MPI_wtime()
    do k = 1, kk
            call encoding4( ix, iy, newtreecode)
            call decoding4(newtreecode, ix,iy,iz)

            call adjacent4_NESW( newtreecode, "__N", neighbor)
            call adjacent4_NESW( neighbor, "__S", newtreecode)
            call adjacent4_NESW( newtreecode, "__W", neighbor)
            call adjacent4_NESW( neighbor, "__E", newtreecode)

            do iz = 1, 1000
            if ( newtreecode == neighbor ) then
              call abort(123980)
            endif
            enddo
    enddo
    write(*,*) "elapsed=", MPI_wtime()-t
    write(*,*) "new neighbor", neighbor, "org:", newtreecode


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "===old lib==="
    treeN = params%max_treelevel
    allocate (treecode(1:treeN), n(1:treeN))

    t = MPI_wtime()
    do k = 1, kk

            ! call encoding_2D(treecode, ix, iy, 2**treeN, 2**treeN, treeN)
            ! call decoding(treecode, ix, iy, iz, treeN)
            !
            ! call adjacent_block_2D(treecode, n, "__N", params%max_treelevel, params%max_treelevel)
            ! call adjacent_block_2D(n, treecode, "__S", params%max_treelevel, params%max_treelevel)
            ! call adjacent_block_2D(treecode, n, "__W", params%max_treelevel, params%max_treelevel)
            ! call adjacent_block_2D(n, treecode, "__E", params%max_treelevel, params%max_treelevel)

            do iz = 1,1000
            if (array_compare(treecode, n, params%max_treelevel)) then
              call abort(123980)
            endif
          enddo
    enddo
    write(*,*) "elapsed=", MPI_wtime()-t
    write(*,'("old nachbar=",20(i1))') n

end subroutine unit_test_treecode
