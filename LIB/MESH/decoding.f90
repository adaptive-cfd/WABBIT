!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name decoding.f90
!> \version 0.5
!
!> \brief from the treecode (quad/octtree), go back to cartesian coordinates i,j,k
!
!>
!! a word on normalization:
!!
!! the pair i,j is to be understood on the level the quaterny code k
!! brings us to. That means a code XXX is within [1,8],[1,8] while a
!! code XXXX is within [1,16],[1,16]
!
!>
!!
!! tested with [1] and other figures procuded by "Encoding"
!! [1] Garagantini. An effective Way to Represent Quadtrees, Graphics and
!! Image Processing (1982)
!
!> \note Subroutine works with both 2D and 3D (quad/octtrees)
!
! ********************************************************************************************

subroutine decoding(treecode, i, j, k, treeN)
    ! global parameters
    use module_params

    implicit none

    !> block position coordinates
    integer(kind=ik), intent(out)    :: i, j, k
    !> treecode size
    integer(kind=ik), intent(in)    :: treeN
    !> treecode
    integer(kind=ik), intent(in)   :: treecode(treeN)
    integer(kind=ik) :: nx, step, l

    ! this is the maximum index possible (the last one on the finest grid)
    nx = 2**treeN

    ! NOTE: one-based indexing
    i = 1
    j = 1
    k = 1

    ! stepping for j=1 level (not j=0, where only one block exists) is nx/2, since
    ! there are two blocks in each direction. For subsequent level, the displacement
    ! get smaller by factors of 2
    step = nx / 2

    ! first entry treecode(1) is indeed level one (and not the root, which would be 0 and is excluded)
    do l = 1, treeN
        select case (treecode(l))
          case (0)
            ! nothing ( origin of does not change, as finer block origin indeed coincides
            ! with her mothers one )
          case (1)
            j = j + step
          case (2)
            i = i + step
          case (3)
            j = j + step
            i = i + step
          case (4)
            k = k + step
          case (5)
            k = k + step
            j = j + step
          case (6)
            k = k + step
            i = i + step
          case (7)
            k = k + step
            j = j + step
            i = i + step
        end select
        step = step / 2
    end do
end subroutine decoding
