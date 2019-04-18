!> \file
module module_treelib

  use module_globals


          
contains

    include "get_neighbor_treecode.f90"

  !===============================================================================
  !> \brief from a tree id and treecode we make a tree-identifieer, which is 
  !! an integer made of the tree_id and the treecode
  !! Note: tree_id s start from 1
  function treecode2int(treearray, tree_id) 
    implicit none
    integer(kind=ik), intent(in) :: treearray(:)
    integer(kind=ik),optional, intent(in) :: tree_id
    integer(kind=ik) :: N,i, potency
    integer(kind=tsize) :: treecode2int
    N = size(treearray,1)

    if (present(tree_id) .and. tree_id>0) then
        treecode2int = tree_id
        potency = floor(log10(real(tree_id))) + 1
        ! the +1 is necessary because it seperates the treecode from the tree_id with a 0
    else
    ! For the rest of wabbit not using tree_ids we have to make sure that
    ! the results stay the same
        potency = -1
        treecode2int = 0 
    endif

    do i = 1, N
      if (treearray(i) >= 0) then
        ! note zero is a bit tedious for comparison, as 0002 and 2 are the same
        ! therefore, we shift all treecodes by 1, thus 0012 gives 1123 as integer
        treecode2int = treecode2int + (10**(i+potency)) * ( treearray(i) + 1 )
      endif
    enddo
  end function

  !===============================================================================
  !> \brief from an integer, return the first (rightmost) digit and remove it from
  !! the number
  elemental subroutine pop( number, element )
    implicit none
    integer(kind=tsize), intent(inout) :: number, element
    ! extract digit
    element = mod( number, 10_tsize)
    ! remove it from number
    number = number / 10
  end subroutine
!===============================================================================

!===============================================================================
!> \brief convert treecodenumber to treearray
subroutine treecodenumber2array( number, level, array )
  implicit none
  integer(kind=tsize),intent(in):: number
  integer(kind=ik),intent(in)      :: level
  integer(kind=ik),intent(inout)   :: array(:)
  integer(kind=ik)                 :: i
  integer(kind=tsize)              :: element,tmp


  tmp=number
  array=0
  i=0
  do while (tmp .ne. 0 )
    call pop(tmp,element)
    array(level-i)=int(element,kind=ik)
    i=i+1
  enddo

end subroutine
!===============================================================================




subroutine adjacent4( treecode, direction, treecode_neighbor)
    implicit none

    integer(kind=tsize), intent(in) :: treecode
    integer(kind=tsize), intent(out) :: treecode_neighbor
    character(len=3), intent(in) :: direction

    integer(kind=tsize) :: treecode_tmp

    select case(direction)
      case('__N','__S','__E','__W')
          call adjacent4_NESW( treecode, direction, treecode_neighbor)
      case('_NE')
          call adjacent4_NESW( treecode, '__N', treecode_tmp)
          call adjacent4_NESW( treecode_tmp, '__E', treecode_neighbor)
      case('_NW')
          call adjacent4_NESW( treecode, '__N', treecode_tmp)
          call adjacent4_NESW( treecode_tmp, '__E', treecode_neighbor)
      case('_SE')
          call adjacent4_NESW( treecode, '__N', treecode_tmp)
          call adjacent4_NESW( treecode_tmp, '__E', treecode_neighbor)
      case('_SW')
          call adjacent4_NESW( treecode, '__N', treecode_tmp)
          call adjacent4_NESW( treecode_tmp, '__E', treecode_neighbor)
      case default
          call abort(118118, "Lord vader, the treelib does not know the direction")
    end select

end subroutine adjacent4

!===============================================================================

  subroutine adjacent4_NESW( treecode, direction, treecode_neighbor)
    implicit none
    integer(kind=tsize), intent(in) :: treecode
    integer(kind=tsize) :: treecode_this, tmp
    integer(kind=tsize), intent(out) :: treecode_neighbor
    character(len=3), intent(in) :: direction
    integer(kind=ik) :: i, j
    logical :: go

    ! copy treecode, as we modify it, but not return this modified value
    treecode_this = treecode
    ! this is the neighbors treecode we're looking for
    treecode_neighbor = 0
    go = .true.

    select case(direction)
      case('__N')
          !*********************************************************************
          ! NORTH
          !*********************************************************************
          call pop(treecode_this, tmp)
          treecode_neighbor = modulo( tmp + 2_tsize, 4_tsize)
          i = 1 ! i=0 is done above
          do while (i < maxdigits .and. go)
            ! last element is used to distinguish
            if (tmp == 2 .or. tmp == 3) then
              ! neighbor on same quadrant: copy, be happy.
              do j = i, maxdigits
                call pop(treecode_this, tmp)
                treecode_neighbor = treecode_neighbor + tmp * 10**j
              enddo
              ! done, so escape loop
              go = .false.
            else
              ! transition to different quadrant.
              call pop(treecode_this, tmp)
              treecode_neighbor = treecode_neighbor + (modulo(tmp+2,4_tsize)) * 10**i
              i = i + 1
            endif
          enddo

        case('__S')
          !*********************************************************************
          ! SOUTH
          !*********************************************************************
          call pop(treecode_this, tmp)
          treecode_neighbor = modulo( tmp + 2_tsize, 4_tsize)
          i = 1 ! i=0 is done above
          do while (i < maxdigits .and. go)
            ! last element is used to distinguish
            if (tmp == 0 .or. tmp == 1) then
              ! neighbor on same quadrant: copy, be happy.
              do j = i, maxdigits
                call pop(treecode_this, tmp)
                treecode_neighbor = treecode_neighbor + tmp * 10**j
              enddo
              ! done, so escape loop
              go = .false.
            else
              ! transition to different quadrant.
              call pop(treecode_this, tmp)
              treecode_neighbor = treecode_neighbor + (modulo(tmp+2_tsize,4_tsize)) * 10**i
              i = i + 1
            endif
          enddo
        case('__E')
          !*********************************************************************
          ! EAST
          !*********************************************************************
          call pop(treecode_this, tmp)
          if (tmp == 0 .or. tmp == 2) then
            treecode_neighbor = modulo( tmp+1_tsize , 4_tsize)
          else
            treecode_neighbor = modulo( tmp-1_tsize , 4_tsize)
          endif

          i = 1 ! i=0 is done above
          do while (i < maxdigits .and. go)
            ! last element is used to distinguish
            if (tmp == 0 .or. tmp == 2) then
              ! neighbor on same quadrant: copy, be happy.
              do j = i, maxdigits
                call pop(treecode_this, tmp)
                treecode_neighbor = treecode_neighbor + tmp * 10**j
              enddo
              ! done, so escape loop
              go = .false.
            else
              ! transition to different quadrant.
              call pop(treecode_this, tmp)
              if ( tmp == 1 .or. tmp == 3) then
                treecode_neighbor = treecode_neighbor + (modulo(tmp-1_tsize,4_tsize)) * 10**i
              else
                treecode_neighbor = treecode_neighbor + (modulo(tmp+1_tsize,4_tsize)) * 10**i
              endif
              i = i + 1
            endif
          enddo
        case('__W')
          !*********************************************************************
          ! WEST
          !*********************************************************************
          call pop(treecode_this, tmp)
          if (tmp == 0 .or. tmp == 2) then
            treecode_neighbor = modulo( tmp+1_tsize , 4_tsize)
          else
            treecode_neighbor = modulo( tmp-1_tsize , 4_tsize)
          endif

          i = 1 ! i=0 is done above
          do while (i < maxdigits .and. go)
            ! last element is used to distinguish
            if (tmp == 1 .or. tmp == 3) then
              ! neighbor on same quadrant: copy, be happy.
              do j = i, maxdigits
                call pop(treecode_this, tmp)
                treecode_neighbor = treecode_neighbor + tmp * 10**j
              enddo
              ! done, so escape loop
              go = .false.
            else
              ! transition to different quadrant.
              call pop(treecode_this, tmp)
              if ( tmp == 1 .or. tmp == 3) then
                treecode_neighbor = treecode_neighbor + (modulo(tmp-1_tsize,4_tsize)) * 10**i
              else
                treecode_neighbor = treecode_neighbor + (modulo(tmp+1_tsize,4_tsize)) * 10**i
              endif
              i = i + 1
            endif
          enddo
      case default
          call abort(118118, "Lord vader, the treelib does not know the direction")
    end select

  end subroutine


!===============================================================================


  subroutine decoding4(treecode1, i, j, k)
      implicit none

      !> block position coordinates
      integer(kind=ik), intent(out)    :: i, j, k
      !> treecode
      integer(kind=tsize), intent(in) :: treecode1
      integer(kind=tsize) :: nx, step, l, treeL, treecode, ix, iy, iz

      ! copy+flip treecode (we modify it but to do give caller the modification back)
      ! gargantini gives the example (I,J)=(6,5) which gives K=321. To reproduce it, set (ix,iy)=7,6
      ! which gives you at this point
      ! treecode=1230000000000000 this means first index (rightmost, "0") is COARSEST
      ! now we reverse the direction (flipint) and end up with
      ! treecode=000000000000321, this means first index (rightmost, "1") is FINEST
      ! Then, we have K[0] (which is the rightmost entry of the code) = 1 as
      ! gargantini has. in the decoding prodecure, we flip the treecode again and start from the coarsest
      ! to finest level
      treecode = flipint( treecode1 )

      ! this is the maximum index possible (the last one on the finest grid)
      nx = 2**maxdigits

      ! NOTE: one-based indexing
      ix = 1
      iy = 1
      iz = 1

      ! stepping for j=1 level (not j=0, where only one block exists) is nx/2, since
      ! there are two blocks in each direction. For subsequent level, the displacement
      ! get smaller by factors of 2
      step = nx / 2

      ! first entry treecode(1) is indeed level one (and not the root, which would be 0 and is excluded)
      do l = 1, maxdigits
          treeL = mod(treecode,10_tsize)
          treecode = treecode / 10
          select case (treeL)
            case (0)
              ! nothing ( origin of does not change, as finer block origin indeed coincides
              ! with her mothers one )
            case (1)
              iy = iy + step
            case (2)
              ix = ix + step
            case (3)
              iy = iy + step
              ix = ix + step
            case (4)
              iz = iz + step
            case (5)
              iz = iz + step
              iy = iy + step
            case (6)
              iz = iz + step
              ix = ix + step
            case (7)
              iz = iz + step
              iy = iy + step
              ix = ix + step
          end select
          step = step / 2
      end do

      i = int(ix, kind=ik)
      j = int(iy, kind=ik)
      k = int(iz, kind=ik)

  end subroutine decoding4



  subroutine encoding4( ix, iy, treecode ) !, Jmax )
    implicit none
    integer(kind=ik), intent(in) :: ix, iy!, Jmax
    integer(kind=tsize), intent(out) :: treecode
    integer(kind=tsize) :: c, d, cl, dl, tl, j
    !integer(kind=ik) :: b(Jmax)

    ! following gargantini, we first require the binary representation of ix,iy
    ! note algorithm ENCODING in her paper requires us to loop down from the highest
    ! level (the last entry of the binary), so here we FLIP the number.
    ! NOTE: gargantini uses 0-based indexing, which is a source of errors. we use 1-based
    c = flipint( toBinary( ix-1 ) )
    d = flipint( toBinary( iy-1 ) )
    treecode = 0_tsize

    ! this is the encoding part from gargantinis paper. Note since we reversed the
    ! binary representations, we loop from 1:end and not from end:-1:1
    do j = 1, maxdigits
      ! pop current digit:
      cl = modulo(c, 10_tsize)
      dl = modulo(d, 10_tsize)
      ! remove digit:
      c = c / 10_tsize
      d = d / 10_tsize
      tl = 9_tsize

      if (cl==0_tsize .and. dl==0_tsize) then
        tl = 0_tsize
      end if

      if (cl==0_tsize .and. dl==1_tsize) then
        tl = 1_tsize
      end if

      if (cl==1_tsize .and. dl==0_tsize) then
        tl = 2_tsize
      end if

      if (cl==1_tsize .and. dl==1_tsize) then
        tl = 3_tsize
      end if

      treecode = treecode + tl * 10_tsize**(j-1)
    end do

    ! gargantini gives the example (I,J)=(6,5) which gives K=321. To reproduce it, set (ix,iy)=7,6
    ! which gives you at this point
    ! treecode=1230000000000000 this means first index (rightmost, "0") is COARSEST
    ! now we reverse the direction (flipint) and end up with
    ! treecode=000000000000321, this means first index (rightmost, "1") is FINEST
    ! Then, we have K[0] (which is the rightmost entry of the code) = 1 as
    ! gargantini has. in the decoding prodecure, we flip the treecode again and start from the coarsest
    ! to finest level
    treecode = flipint(treecode)
  end subroutine encoding4


  ! convert given integger decnum to binary representation
  ! textbook code
  INTEGER(kind=tsize) FUNCTION toBinary(decNum)
    IMPLICIT NONE

    INTEGER(kind=ik), INTENT(IN) :: decNum !!ARGUMENT COMING IN TO FUNCTION
    INTEGER(kind=tsize) :: bineq !!THE BINARY EQUIVALENT OF decNum
    INTEGER(kind=tsize) :: power !!THE POWER OF 10 WE ARE CURRENTLY DEALING WITH
    INTEGER(kind=tsize) :: temp !!temp VALUE TO HOLD A COPY OF THE ARGUMENT

    bineq = 0
    power = 0
    temp = 0

    temp = decNum

    DO WHILE(temp > 0)
      bineq = bineq + (MOD(temp, 2_tsize) * (MOD(temp, 2_tsize) * (10 ** power)))
      power = power + 1
      temp = temp / 2
    END DO

    toBinary = bineq
  END FUNCTION toBinary


  ! flip an integer number, eg given 523 return 325
  ! or 000021 gives 120000, so we have to use zero-padding!
  ! The length of returned (zero padded) reverse int is specified by maxdigits
  integer(kind=tsize) function flipint(i)
    implicit none
    integer(kind=tsize), intent(in) :: i
    integer(kind=tsize) ::  tmp, j
    flipint = 0
    tmp = i
    j = 1
    do while ( j <= maxdigits )
      flipint = MOD(tmp,10_tsize) + flipint*10
      tmp = tmp / 10
      j = j+1
    end do
    return
  end function



  !> \callgraph
  ! ********************************************************************************************
  ! WABBIT
  ! ============================================================================================
  !> \name    adjacent_block_2D.f90
  !> \version  0.4
  !> \author   msr
  !
  !> \brief give treecode for adjacent block
  !
  !>
  !! input:
  !!                    - treecode for block N
  !!                    - direction for neighbor search
  !!                    - max treelevel
  !!
  !! output:
  !!                    - neighbor treecode, for neighbor on same level
  !!
  !!
  !! = log ======================================================================================
  !! \n
  !! 07/11/16 - switch to v0.4
  ! ********************************************************************************************

  recursive subroutine adjacent_block_2D(me, neighbor, direction, level, max_treelevel)

  !---------------------------------------------------------------------------------------------
  ! modules

      ! global parameters
      use module_params

  !---------------------------------------------------------------------------------------------
  ! variables

      implicit none

      !> max treelevel
      integer(kind=ik), intent(in)        :: max_treelevel
      !> mesh level
      integer(kind=ik), intent(in)        :: level
      !> block treecode
      integer(kind=ik), intent(in)        :: me(max_treelevel)
      !> direction for neighbor search
      character(len=3), intent(in)        :: direction

      !> neighbor treecode
      integer(kind=ik), intent(out)       :: neighbor(max_treelevel)

      ! treecode variable
      integer(kind=ik)                    :: neighbor2(max_treelevel)
      ! loop variable
      integer(kind=ik)                    :: i

  !---------------------------------------------------------------------------------------------
  ! variables initialization

      neighbor = -1

  !---------------------------------------------------------------------------------------------
  ! main body

      select case(direction)
          case('__N')
          ! north
              neighbor(level) = modulo(me(level)+2, 4)
              i = level - 1

              do while (i /= 0)
                  if ( (me(i+1)==2) .or. (me(i+1)==3) ) then
                      neighbor(1:i) = me(1:i)
                      i = 0
                  else
                      neighbor(i) = modulo(me(i)+2, 4)
                      i = i - 1
                  end if
              end do

          case('__S')
          ! south
              neighbor(level) = modulo(me(level)+2, 4)
              i = level - 1

              do while (i /= 0)
                  if ( (me(i+1)==0) .or. (me(i+1)==1) ) then
                      neighbor(1:i) = me(1:i)
                      i = 0
                  else
                      neighbor(i) = modulo(me(i)+2, 4)
                      i = i - 1
                  end if
              end do

          case('__E')
          ! east
              if ( (me(level)==0) .or. (me(level)==2) ) then
                  neighbor(level) = modulo(me(level)+1, 4)
              else
                  neighbor(level) = modulo(me(level)-1, 4)
              end if

              i = level - 1

              do while (i /= 0)
                  if ( (me(i+1)==0) .or. (me(i+1)==2) ) then
                      neighbor(1:i) = me(1:i)
                      i = 0
                  else
                      if ( (me(i)==1) .or. (me(i)==3) ) then
                          neighbor(i) = modulo(me(i)-1, 4)
                      else
                          neighbor(i) = modulo(me(i)+1, 4)
                      end if
                      i = i - 1
                  end if
              end do

          case('__W')
          ! west
              if ( (me(level)==0) .or. (me(level)==2) ) then
                  neighbor(level) = modulo(me(level)+1, 4)
              else
                  neighbor(level) = modulo(me(level)-1, 4)
              end if

              i = level - 1

              do while (i /= 0)
                  if ( (me(i+1)==1) .or. (me(i+1)==3) ) then
                      neighbor(1:i) = me(1:i)
                      i = 0
                  else
                      if ( (me(i)==1) .or. (me(i)==3) ) then
                          neighbor(i) = modulo(me(i)-1, 4)
                      else
                          neighbor(i) = modulo(me(i)+1, 4)
                      end if
                      i = i - 1
                  end if
              end do

          case('_NE')
          ! northeast
              call adjacent_block_2D(me, neighbor2, '__N', level, max_treelevel)
              call adjacent_block_2D(neighbor2, neighbor, '__E', level, max_treelevel)

          case('_NW')
          ! northwest
              call adjacent_block_2D(me, neighbor2, '__N', level, max_treelevel)
              call adjacent_block_2D(neighbor2, neighbor, '__W', level, max_treelevel)

          case('_SE')
          ! southeast
              call adjacent_block_2D(me, neighbor2, '__S', level, max_treelevel)
              call adjacent_block_2D(neighbor2, neighbor, '__E', level, max_treelevel)

          case('_SW')
          ! southwest
              call adjacent_block_2D(me, neighbor2, '__S', level, max_treelevel)
              call adjacent_block_2D(neighbor2, neighbor, '__W', level, max_treelevel)

          case default
              call abort(118119, "Lord vader, the treelib does not know the direction")
      end select
  end subroutine adjacent_block_2D


  !> \file
  !> \callgraph
  ! ********************************************************************************************
  ! WABBIT
  ! ============================================================================================
  !> \name adjacent_block_3D.f90
  !> \version 0.5
  !> \author msr
  !
  !> \brief give treecode for adjacent block in 3D \n
  !
  !>
  !! input:
  !!                    - treecode for block N
  !!                    - direction for neighbor search
  !!                    - max treelevel
  !!
  !!output:
  !!                    - neighbor treecode, for neighbor on same level
  !!
  !! \n
  !!
  !  --------------------------------------------------------------------------------------------
  !> neighbor codes: \n
  !  ---------------
  !> for imagination:
  !!                   - 6-sided dice with '1'-side on top, '6'-side on bottom, '2'-side in front
  !!                   - edge: boundary between two sides - use sides numbers for coding
  !!                   - corner: between three sides - so use all three sides numbers
  !!                   - block on higher/lower level: block shares face/edge and one unique corner,
  !!                     so use this corner code in second part of neighbor code
  !!
  !! faces:  '__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___' \n
  !! edges:  '_12/___', '_13/___', '_14/___', '_15/___' \n
  !!         '_62/___', '_63/___', '_64/___', '_65/___' \n
  !!         '_23/___', '_25/___', '_43/___', '_45/___' \n
  !! corner: '123/___', '134/___', '145/___', '152/___' \n
  !!         '623/___', '634/___', '645/___', '652/___' \n
  !! \n
  !! complete neighbor code array, 74 possible neighbor relations \n
  !! neighbors = (/'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___', '_12/___', '_13/___', '_14/___', '_15/___',
  !!               '_62/___', '_63/___', '_64/___', '_65/___', '_23/___', '_25/___', '_43/___', '_45/___', '123/___', '134/___',
  !!               '145/___', '152/___', '623/___', '634/___', '645/___', '652/___', '__1/123', '__1/134', '__1/145', '__1/152',
  !!               '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634',
  !!               '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652',
  !!               '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652',
  !!               '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652',
  !!               '_43/134', '_43/634', '_45/145', '_45/645' /) \n
  ! --------------------------------------------------------------------------------------------
  !
  !> \details = log ======================================================================================
  !! \n
  !! 27/01/17 - start
  !
  ! ********************************************************************************************

  recursive subroutine adjacent_block_3D(me, neighbor, direction, level, max_treelevel)

  !---------------------------------------------------------------------------------------------
  ! modules

      ! global parameters
      use module_params

  !---------------------------------------------------------------------------------------------
  ! variables

      implicit none

      !> max treelevel
      integer(kind=ik), intent(in)        :: max_treelevel
      !> mesh level
      integer(kind=ik), intent(in)        :: level
      !> block treecode
      integer(kind=ik), intent(in)        :: me(max_treelevel)
      !> direction for neighbor search
      character(len=7), intent(in)        :: direction

      !> neighbor treecode
      integer(kind=ik), intent(out)       :: neighbor(max_treelevel)
      ! treecode variable
      integer(kind=ik)                    :: neighbor2(max_treelevel), neighbor3(max_treelevel)

      ! loop variable
      integer(kind=ik)                    :: i

  !---------------------------------------------------------------------------------------------
  ! variables initialization

      neighbor = me

  !---------------------------------------------------------------------------------------------
  ! main body

      select case(direction)

          case('__1/___')
          ! '1'-side
              neighbor(level) = modulo(me(level) + 4, 8)
              i = level - 1

              do while (i /= 0)
                  select case( me(i+1) )
                      case (0,1,2,3)
                          ! nothing to do, leave loop
                          i = 0

                      case (4,5,6,7)
                          neighbor(i) = modulo(me(i) + 4, 8)
                          i = i - 1

                  end select
              end do

          case('__6/___')
          ! '6'-side
              neighbor(level) = modulo(me(level) - 4, 8)
              i = level - 1

              do while (i /= 0)
                  select case( me(i+1) )
                      case (4,5,6,7)
                          ! nothing to do, leave loop
                          i = 0

                      case (0,1,2,3)
                          neighbor(i) = modulo(me(i) - 4, 8)
                          i = i - 1

                  end select
              end do

          case('__3/___')
          ! '3'-side
              select case(me(level))
                  case(0,1,4,5)
                      neighbor(level) = modulo(me(level) + 2, 8)
                  case(2,3,6,7)
                      neighbor(level) = modulo(me(level) + 6, 8)
              end select

              i = level - 1

              do while (i /= 0)
                  select case( me(i+1) )
                      case (0,1,4,5)
                          ! nothing to do, leave loop
                          i = 0

                      case (2,3,6,7)
                          select case(me(i))
                              case (0,1,4,5)
                                  neighbor(i) = modulo(me(i) + 2, 8)
                              case (2,3,6,7)
                                  neighbor(i) = modulo(me(i) + 6, 8)
                          end select
                          i = i - 1

                  end select
              end do

          case('__5/___')
          ! '5'-side
              select case(me(level))
                  case(2,3,6,7)
                      neighbor(level) = modulo(me(level) - 2, 8)
                  case(0,1,4,5)
                      neighbor(level) = modulo(me(level) - 6, 8)
              end select

              i = level - 1

              do while (i /= 0)
                  select case( me(i+1) )
                      case (2,3,6,7)
                          ! nothing to do, leave loop
                          i = 0

                      case (0,1,4,5)
                          select case(me(i))
                              case (2,3,6,7)
                                  neighbor(i) = modulo(me(i) - 2, 8)
                              case (0,1,4,5)
                                  neighbor(i) = modulo(me(i) - 6, 8)
                          end select
                          i = i - 1

                  end select
              end do

          case('__2/___')
          ! '2'-side
              select case(me(level))
                  case(1,3,5,7)
                      neighbor(level) = modulo(me(level) - 1, 8)
                  case(0,2,4,6)
                      neighbor(level) = modulo(me(level) - 7, 8)
              end select

              i = level - 1

              do while (i /= 0)
                  select case( me(i+1) )
                      case(1,3,5,7)
                          ! nothing to do, leave loop
                          i = 0

                      case(0,2,4,6)
                          select case(me(i))
                              case(1,3,5,7)
                                  neighbor(i) = modulo(me(i) - 1, 8)
                              case(0,2,4,6)
                                  neighbor(i) = modulo(me(i) - 7, 8)
                          end select
                          i = i - 1

                  end select
              end do

          case('__4/___')
          ! '4'-side
              select case(me(level))
                  case(0,2,4,6)
                      neighbor(level) = modulo(me(level) + 1, 8)
                  case(1,3,5,7)
                      neighbor(level) = modulo(me(level) + 7, 8)
              end select

              i = level - 1

              do while (i /= 0)
                  select case( me(i+1) )
                      case(0,2,4,6)
                          ! nothing to do, leave loop
                          i = 0

                      case(1,3,5,7)
                          select case(me(i))
                              case(0,2,4,6)
                                  neighbor(i) = modulo(me(i) + 1, 8)
                              case(1,3,5,7)
                                  neighbor(i) = modulo(me(i) + 7, 8)
                          end select
                          i = i - 1

                  end select
              end do

          case('_12/___')
          ! '12'-edge
              call adjacent_block_3D(me, neighbor2, '__1/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__2/___', level, max_treelevel)

          case('_13/___')
          ! '13'-edge
              call adjacent_block_3D(me, neighbor2, '__1/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__3/___', level, max_treelevel)

          case('_14/___')
          ! '14'-edge
              call adjacent_block_3D(me, neighbor2, '__1/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__4/___', level, max_treelevel)

          case('_15/___')
          ! '15'-edge
              call adjacent_block_3D(me, neighbor2, '__1/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__5/___', level, max_treelevel)

          case('_62/___')
          ! '62'-edge
              call adjacent_block_3D(me, neighbor2, '__6/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__2/___', level, max_treelevel)

          case('_63/___')
          ! '63'-edge
              call adjacent_block_3D(me, neighbor2, '__6/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__3/___', level, max_treelevel)

          case('_64/___')
          ! '64'-edge
              call adjacent_block_3D(me, neighbor2, '__6/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__4/___', level, max_treelevel)

          case('_65/___')
          ! '65'-edge
              call adjacent_block_3D(me, neighbor2, '__6/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__5/___', level, max_treelevel)

          case('_23/___')
          ! '23'-edge
              call adjacent_block_3D(me, neighbor2, '__2/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__3/___', level, max_treelevel)

          case('_25/___')
          ! '25'-edge
              call adjacent_block_3D(me, neighbor2, '__2/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__5/___', level, max_treelevel)

          case('_43/___')
          ! '43'-edge
              call adjacent_block_3D(me, neighbor2, '__4/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__3/___', level, max_treelevel)

          case('_45/___')
          ! '45'-edge
              call adjacent_block_3D(me, neighbor2, '__4/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__5/___', level, max_treelevel)

          case('123/___')
          ! '123'-corner
              call adjacent_block_3D(me, neighbor3, '__1/___', level, max_treelevel)
              call adjacent_block_3D(neighbor3, neighbor2, '__2/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__3/___', level, max_treelevel)

          case('134/___')
          ! '134'-corner
              call adjacent_block_3D(me, neighbor3, '__1/___', level, max_treelevel)
              call adjacent_block_3D(neighbor3, neighbor2, '__3/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__4/___', level, max_treelevel)

          case('145/___')
          ! '145'-corner
              call adjacent_block_3D(me, neighbor3, '__1/___', level, max_treelevel)
              call adjacent_block_3D(neighbor3, neighbor2, '__4/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__5/___', level, max_treelevel)

          case('152/___')
          ! '152'-corner
              call adjacent_block_3D(me, neighbor3, '__1/___', level, max_treelevel)
              call adjacent_block_3D(neighbor3, neighbor2, '__5/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__2/___', level, max_treelevel)

          case('623/___')
          ! '623'-corner
              call adjacent_block_3D(me, neighbor3, '__6/___', level, max_treelevel)
              call adjacent_block_3D(neighbor3, neighbor2, '__2/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__3/___', level, max_treelevel)

          case('634/___')
          ! '634'-corner
              call adjacent_block_3D(me, neighbor3, '__6/___', level, max_treelevel)
              call adjacent_block_3D(neighbor3, neighbor2, '__3/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__4/___', level, max_treelevel)

          case('645/___')
          ! '645'-corner
              call adjacent_block_3D(me, neighbor3, '__6/___', level, max_treelevel)
              call adjacent_block_3D(neighbor3, neighbor2, '__4/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__5/___', level, max_treelevel)

          case('652/___')
          ! '652'-corner
              call adjacent_block_3D(me, neighbor3, '__6/___', level, max_treelevel)
              call adjacent_block_3D(neighbor3, neighbor2, '__5/___', level, max_treelevel)
              call adjacent_block_3D(neighbor2, neighbor, '__2/___', level, max_treelevel)

          case default
              call abort(118112, "Lord vader, the treelib does not know the direction")
      end select

  end subroutine adjacent_block_3D

  !> \brief from the treecode (quad/octtree), go back to cartesian coordinates i,j,k
  !
  !> \details
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


  !> \author PKrah
  !> \brief encoding from block position (cartesian coordinates) to treecode
  !! \details
  !> \note Subroutine works with both 2D and 3D (quad/octtrees)
  !! \date 3.08.2018 - created from 2D/3D and simplified for d dimensions
! ix(1)->   0      1     2     3      4     5     6     7
! ix(2)  ---------------------------------------------------
!   0    || 000 | 001 | 010 | 011 || 100  | 101 | 110 | 111||
!        ---------------------------------------------------
!   1    || 002 | 003 | 012 | 013 || 102  | 103 | 112 | 113||
!        ---------------------------------------------------
!   2    || 020 | 021 | 030 | 031 || 120  | 121 | 130 | 131||
!        ---------------------------------------------------
!   3    || 022 | 023 | 032 | 033 || 122  | 123 | 132 | 133||
!        :::::::::::::::::::::::::::::::::::::::::::::::::::
!   4    || 200 | 201 | 210 | 211 || 300  | 301 | 310 | 311||
!        ---------------------------------------------------
!   5    || 202 | 203 | 212 | 213 || 302  | 303 | 312 | 313||
!       ----------------------------------------------------
!   6    || 220 | 221 | 230 | 231 || 320  | 321 | 330 | 331||
!       ----------------------------------------------------
!   7    || 222 | 223 | 232 | 233 || 322  | 323 | 332 | 333||
!       ----------------------------------------------------
! Coordinate Transformation:
! 1. transformation of the cartesian coordinates from dezimal(basis 10) to binary (basis 2)
!       ix      = 2^{N-1} c_{N-1} + ... +  2^2 c_2 + 2^1 c_1 + 2^0 c_0 =(c_{N-1},...,c_2,c_1,c_0)_2
!       example:
!       ix      =(5,6)_10                 =((101),(110))_2
! 2. coordinate transformation from cartesian to treecode
!       ix -> treecode
!       example for d=2
!       ((101),(011))_2 -> 2^0 (101)_2 + 2^1 (110)_2=(3 2 1)
!       compare to the quaternary codes in the quadrants above
!




  subroutine encoding(treearray, ix, dim , block_num, treeN)
  !---------------------------------------------------------------------------------------------
  ! modules

      ! global parameters
      use module_params

  !---------------------------------------------------------------------------------------------
  ! variables

      implicit none
      !> dimension (2 or 3)
      integer(kind=ik), intent(in)    :: dim
      !> block position coordinates
      integer(kind=ik), intent(in)    :: ix(dim)
      !> number of blocks
      integer(kind=ik), intent(in)    :: block_num
      !> treecode size
      integer(kind=ik), intent(in)    :: treeN
      !> treecode
      integer(kind=ik), intent(out)   :: treearray(treeN)

      ! variables for calculate real treecode length N
      real(kind=rk)                   :: Jn
      integer(kind=ik)                :: N, l,d


      ! auxiliary vectors
      integer(kind=ik), allocatable   :: ix_binary(:)
      ! real treecode length
      Jn = log(dble(block_num)) / log(2.0_rk**dim)
      N  = nint(Jn)
      ! set N to 1, for one block decomposition
      if (N==0) N=1
      ! reset output
      treearray = 0
      ! allocate auxiliary vectors
      allocate( ix_binary(N) )

      ! convert block coordinates into binary numbers
      do d=1, dim
          call int_to_binary(ix(d)-1, N, ix_binary)
          treearray  =     treearray+ix_binary(1:N)*2**(d-1)
      end do
      ! converts treecodearray to treecode
      ! if treearray is (00123) the treecode will be 00123+11111=11234
      !treecode=treecode2int(treearray)
      ! clean up
      deallocate( ix_binary )

  end subroutine encoding


  !> \brief convert a integer i to binary b \n
  !! binary return as vector with length N
  !! \details
  !!\date 07/11/16 - switch to v0.4
  subroutine int_to_binary(i, N, b)

      ! global parameters
      use module_params

      implicit none

      !> integer to convert into binary
      integer(kind=ik), intent(in)    :: i
      !> length of binary output vector
      integer(kind=ik), intent(in)    :: N
      !> output vector
      integer(kind=ik), intent(out)   :: b(N)
      ! loop variables
      integer(kind=ik)                :: j, k, tmp(N)


      j = 1
      b = 0
      tmp=0
      k = i

      do while (k > 0)
          tmp(j) = mod(k, 2)
          k = int(k/2)
          j = j + 1
      end do

      ! binary has to be flipped to the right order of the basis vectors (2^{N-1},...,2^1,2^0)
      do j = 1, N
       b(N+1-j) = tmp(j)
     end do
  end subroutine int_to_binary


end module
