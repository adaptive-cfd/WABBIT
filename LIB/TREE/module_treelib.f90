!===========================================================================
!> Module to deal with the tree-structure in quadtree or octree format
! *****************************************************************************
module module_treelib

  use module_globals

  contains

#include "neighborhood.f90"

  !> \author engels
  !> \brief Compute block spacing and origin from lgt_block
  !> \details For any block lgt_id this routine computes, from the treearray stored in
  !! lgt_block( lgt_id, : ), the block's origin and grid spacing. Note spacing
  !! and origin are 3D vectors, the third component being zero in a 2D case. \n
  subroutine get_block_spacing_origin_array( treecode, domain, Bs, dim, x0, dx )
      implicit none

      integer(kind=ik), intent(in)               :: dim
      real(kind=rk), dimension(1:3), intent(in)  :: domain
      integer(kind=ik), dimension(1:3)           :: Bs
      integer(kind=ik), intent(in)               :: treecode(:)
      !> output
      real(kind=rk), dimension(1:3), intent(out) :: x0, dx

      integer(kind=ik)                           :: ix, iy, iz, J

      ! fetch this blocks level:
      J = size(treecode)

      ! compute its coordinates in ijk space
      call decoding(treecode, ix, iy, iz, J)

      ! the spacing on a block is the basic spacing Lx/Bs of the coarsest block (if there
      ! is only one block, j=0) divided by 2 for each level, thus the 2^-j factor
      dx = 0.0_rk
      dx(1:dim) = 2.0_rk**(-J) * domain(1:dim) / real( Bs(1:dim), kind=rk )
      ! dx(1:dim) = 2.0_rk**(-J) * domain(1:dim) / real( Bs(1:dim)-1, kind=rk )

      ! note zero based indexing:
      ! x0 = real( ((/ix,iy,iz/) - 1)*(Bs-1), kind=rk) * dx
      x0 = real( ((/ix,iy,iz/) - 1)*(Bs), kind=rk) * dx

  end subroutine get_block_spacing_origin_array

  !===============================================================================
  !> \author JB
  !> \brief Compute block spacing and origin from lgt_block
  !> \details For any block lgt_id this routine computes, from the numerical
  !! treecode stored in lgt_block( lgt_id, : ), the block's origin and grid spacing.
  !! \note spacing and origin are 3D vectors, the third component being zero in a 2D case.
  subroutine get_block_spacing_origin_b( treecode, domain, Bs, x0, dx, dim, level, max_level, input_dec)
    implicit none

    !> Dimension, 2 or 3 and defaults to 3
    integer(kind=ik), optional, intent(in)     :: dim
    !> Maximum level to be encoded, should be params%Jmax, defaults to maxdigits
    integer(kind=ik), optional, intent(in)     :: max_level
    !> Level at which the block is situated, defaults to max_level.
    integer(kind=ik), optional, intent(in)     :: level
    !> Domain extent
    real(kind=rk), dimension(1:3), intent(in)  :: domain
    !> Block size / Points in block
    integer(kind=ik), dimension(1:3)           :: Bs
    !> Numerical binary
    integer(kind=tsize), intent(in)            :: treecode
    !> If this is a decimal treecode input we dont have to write a new function,
    !! defaults to false
    logical, optional, intent(in)              :: input_dec
    !> output
    real(kind=rk), dimension(1:3), intent(out) :: x0, dx

    integer(kind=ik)                           :: n_dim, n_level, max_tclevel
    logical                                    :: n_input_dec
    integer(kind=ik), dimension(1:3)           :: ixyz

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1
    ! Set default for input dec
    n_input_dec = .false.; if (present(input_dec)) n_input_dec = input_dec

    ! compute its coordinates in ijk space
    if (n_input_dec) then
      call decoding_d(ixyz, treecode, dim=n_dim, level=n_level, max_level=max_tclevel)
    else
      call decoding_b(ixyz, treecode, dim=n_dim, level=n_level, max_level=max_tclevel)
    end if

    ! the spacing on a block is the basic spacing Lx/Bs of the coarsest block (if there
    ! is only one block, j=0) divided by 2 for each level, thus the 2^-j factor
    dx = 0.0_rk
    dx(1:dim) = 2.0_rk**(-n_level) * domain(1:dim) / real( Bs(1:dim), kind=rk )
    ! dx(1:dim) = 2.0_rk**(-J) * domain(1:dim) / real( Bs(1:dim)-1, kind=rk )

    ! note zero based indexing:
    ! x0 = real( ((/ix,iy,iz/) - 1)*(Bs-1), kind=rk) * dx
    x0 = real( ((ixyz) - 1)*Bs * dx, kind=rk)

  end subroutine get_block_spacing_origin_b

  !===============================================================================
  !> \brief get digit at specific index and return digit
  function tc_get_digit_at_level_b( treecode, dim, level, max_level)
    implicit none
    !> Treecode in decimal numerical representation
    integer(kind=tsize),intent(in)   :: treecode
    !> Depth to read the digit
    integer(kind=ik), optional, intent(in)      :: level
    !> Dimension, 2 or 3, defaults to 3
    integer(kind=ik), optional, intent(in)       :: dim
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level

    integer(kind=ik)                 :: i_level, n_level, max_tclevel, n_dim

    integer(kind=ik)              :: tc_get_digit_at_level_b

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    ! now output specific level
    tc_get_digit_at_level_b = int(ibits(treecode, (max_tclevel - n_level)*n_dim, n_dim), kind=ik)
  end function
  !===============================================================================

  !===============================================================================
  !> \brief set number at specific index and return new altered treecode
  function tc_set_digit_at_level_b( treecode, digit, dim, level, max_level)
    implicit none
    !> Treecode in decimal numerical representation
    integer(kind=tsize),intent(in)   :: treecode
    !> Digit to be set at level, has to be 0-3 for 2D and 0-7 for 3D
    integer(kind=ik), intent(in)        :: digit
    !> Depth to set the digit
    integer(kind=ik), optional, intent(in)      :: level
    !> Dimension, 2 or 3, defaults to 3
    integer(kind=ik), optional, intent(in)       :: dim
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level

    integer(kind=ik)                 :: i_level, n_level, max_tclevel, n_dim

    integer(kind=tsize)              :: tc_set_digit_at_level_b

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    ! subtract specific level
    tc_set_digit_at_level_b = treecode - ishft(ibits(treecode, (max_tclevel - n_level)*n_dim, n_dim), (max_tclevel - n_level)*n_dim)
    ! now set specific level
    tc_set_digit_at_level_b = tc_set_digit_at_level_b + ishft(int(digit, kind=tsize), (max_tclevel - n_level)*n_dim)
  end function
  !===============================================================================

  !===============================================================================
  !> \brief Define a "<"-operator for two blocks based on their tree_id, level and tc
  !> \details In order to sort the array we want to uniquely arrange them to find
  !! blocks quickly. This function defines the "<"-operator based on following input:
  !!  - tree_ID
  !!  - level
  !!  - treecode in binary or numerical encoding
  function tc_id_lower( array1, array2, check_meta)
    implicit none
    !> Array1 of size 2-3 to be compared: (tc1, tc2, (level+100*tree_id))
    integer(kind=ik),intent(in)   :: array1(:)
    !> Array2 of size 2-3 to be compared: (tc1, tc2, (level+100*tree_id))
    integer(kind=ik),intent(in)   :: array2(:)
    !> Is tree_id and level included and should be checked?
    logical                       :: check_meta
    !> Output
    logical                       :: tc_id_lower

    ! check if tree_id or level are larger, they are combined in one number
    ! check if treecode is larger if tree_id and level are similar
    if (check_meta) then
      tc_id_lower = (array1(3) < array2(3) .or. &
      (array1(3) == array2(3) .and. get_tc(array1(1:2)) < get_tc(array2(1:2))))
    else
      tc_id_lower = (get_tc(array1(1:2)) < get_tc(array2(1:2)))
    endif
  end function
  !===============================================================================

  ! !===============================================================================
  ! !> \brief Define a "<"-operator for two blocks based on their tree_id, level and tc
  ! !> \details In order to sort the array we want to uniquely arrange them to find
  ! !! blocks quickly. This function defines the "<"-operator based on following input:
  ! !!  - tree_ID
  ! !!  - level
  ! !!  - treecode in binary or numerical encoding
  ! function tc_id_lower2( array1, array2)
  !   implicit none
  !   !> Array1 of size 4 to be compared: (tc1, tc2, level, tree_id)
  !   integer(kind=ik),intent(in)   :: array1(:)
  !   !> Array2 of size 4 to be compared: (tc1, tc2, level, tree_id)
  !   integer(kind=ik),intent(in)   :: array2(:)
  !   !> Output
  !   logical                       :: tc_id_lower2

  !   ! check if tree_id or level are larger, they are combined in one number
  !   ! check if treecode is larger if tree_id and level are similar
  !   if (array1(4) > array2(4) .or. &
  !     (array1(4) == array2(4) .and. array1(3) > array2(3)) .or. &
  !     (array1(4) == array2(4) .and. array1(3) == array2(3) &
  !     ! .and. array1(2) >= array2(2))) then
  !     .and. get_tc(array1(1:2)) >= get_tc(array2(1:2)))) then
  !       tc_id_lower2 = .false.
  !   else
  !     tc_id_lower2 = .true.
  !   endif
  ! end function
  ! !===============================================================================

  !===============================================================================
  !> \brief Extract from lgt_block from the specific ID the treecode
  !> \details As treecode needs two normal integers in size, this is used to decode
  function get_tc( tc_2)
    implicit none
    !> Array of two elements where we get / decoded ID
    integer(kind=ik),intent(in)   :: tc_2(:)
    !> Output, Treecode to be encoded
    integer(kind=tsize)           :: get_tc

    ! first number is shifted by amount of bits in second int
    ! We have to left-and-right shift the second int as when left-most bit is set to 1
    ! it is treated as negative number and left-extended with 1s, with the shift we clear that
    get_tc = ishft(int(tc_2(1), kind=tsize), bit_size(tc_2(2))) &
           + ishft(ishft(int(tc_2(2), kind=tsize), bit_size(tc_2(2))), -bit_size(tc_2(2)))
  end function
  !===============================================================================

  !===============================================================================
  !> \brief Encode in lgt_block with the specific ID the treecode
  !> \details As treecode needs two normal integers in size, this is used to encode
  subroutine set_tc( tc_2, treecode)
    implicit none
    !> Array of two elements where we set / encode the ID
    integer(kind=ik),intent(out)     :: tc_2(:)
    !> Treecode to be encoded
    integer(kind=tsize),intent(in)   :: treecode

    ! select first and second number and store them
    tc_2(1) = int(ibits(treecode, bit_size(tc_2(2)), tc_2(1)), kind=ik)
    tc_2(2) = int(ibits(treecode, 0, tc_2(2)), kind=ik)
  end subroutine
  !===============================================================================

  !===============================================================================
  !> \brief set all levels higher as level to 0 to ensure coarsening can be comparable
  function tc_clear_until_level_b( treecode, dim, level, max_level)
    implicit none
    !> Treecode in decimal numerical representation
    integer(kind=tsize),intent(in)   :: treecode
    !> Levels to keep
    integer(kind=ik), optional, intent(in)      :: level
    !> Dimension, 2 or 3, defaults to 3
    integer(kind=ik), optional, intent(in)       :: dim
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level

    integer(kind=ik)                 :: i_level, n_level, max_tclevel, n_dim

    integer(kind=tsize)              :: tc_clear_until_level_b

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    ! clearing by double shifting
    tc_clear_until_level_b = ishft(ishft(treecode, -(max_tclevel - n_level)*n_dim), (max_tclevel - n_level)*n_dim)
  end function
  !===============================================================================

  !===============================================================================
  !> \brief from an integer, return the first (rightmost) digit and remove it from
  !! number
  elemental subroutine pop( number, element )
    implicit none
    integer(kind=tsize), intent(inout) :: number, element
    ! extract digit
    element = mod( number, 10_tsize)
    ! remove it from number
    number = number / 10_tsize
  end subroutine
  !===============================================================================

  !===============================================================================
  !> \brief convert treearray to numerical decimal treecodenumber
  subroutine array2tcd( treecode, treearray, level, max_level)
    implicit none
    !> Treecode in decimal numerical representation
    integer(kind=tsize),intent(out)   :: treecode
    !> Depth to be iterated over
    integer(kind=ik), optional, intent(in)      :: level
    !> Treearray
    integer(kind=ik),intent(in)   :: treearray(:)
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level

    integer(kind=ik)                 :: i_level, n_level, max_tclevel
    integer(kind=tsize)              :: element,tmp

    ! Set defaults for dimension, level and max_level
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    treecode = 0_tsize
    do i_level = 1, n_level
      if (treearray(i_level) >= 0) then
        treecode = treecode + treearray(i_level) * 10_tsize**(max_tclevel-i_level)
      end if
    end do
  end subroutine
  !===============================================================================

  !===============================================================================
  !> \brief convert treearray to numerical binary treecodenumber
  subroutine array2tcb( treecode, treearray, dim, level, max_level)
    implicit none
    !> Treecode in binary numerical representation
    integer(kind=tsize),intent(out)   :: treecode
    !> Depth to be iterated over
    integer(kind=ik), optional, intent(in)      :: level
    !> Treearray
    integer(kind=ik),intent(in)   :: treearray(:)
    !> Dimension, 2 or 3, defaults to 3
    integer(kind=ik), optional, intent(in)       :: dim
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level

    integer(kind=ik)                 :: i_level, n_level, max_tclevel, n_dim

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    treecode = 0_tsize
    do i_level = 1, n_level
      if (treearray(i_level) >= 0) then
        treecode = treecode + ishft(int(treearray(i_level), kind=tsize), n_dim*(max_tclevel-i_level))
      end if
    end do
  end subroutine
  !===============================================================================

  !===============================================================================
  !> \brief convert numerical decimal treecodenumber to treearray
  subroutine tcd2array( treecode, array, level, max_level)
    implicit none
    !> Treecode in decimal numerical representation
    integer(kind=tsize),intent(in)   :: treecode
    !> Depth to be iterated over
    integer(kind=ik), optional, intent(in)      :: level
    !> Treearray
    integer(kind=ik),intent(out)   :: array(:)
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level

    integer(kind=ik)                 :: i_level, n_level, max_tclevel
    integer(kind=tsize)              :: element,tmp

    ! Set defaults for dimension, level and max_level
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    tmp=treecode
    ! init as -1 or unset
    array=-1
    ! skip some elements
    do i_level = n_level, max_tclevel-1
      call pop(tmp,element)
    end do
    ! repeatedly pop of numbers and insert them
    do i_level = 0, n_level-1
      call pop(tmp,element)
      array(n_level-i_level)=int(element,kind=ik)
    enddo
  end subroutine
  !===============================================================================

  !===============================================================================
  !> \brief convert numerical binary treecodenumber to treearray
  !> \author JB
  subroutine tcb2array( treecode, array, dim, level, max_level)
    implicit none
    !> Treecode in binary numerical representation
    integer(kind=tsize),intent(in)   :: treecode
    !> Depth to be iterated over
    integer(kind=ik), optional, intent(in)      :: level
    !> Treearray
    integer(kind=ik),intent(out)   :: array(:)
    !> Dimension, 2 or 3, defaults to 3
    integer(kind=ik), optional, intent(in)       :: dim
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level

    integer(kind=ik)                 :: i_level, n_dim, n_level, max_tclevel

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1
    
    ! init as -1 or unset
    array=-1
    do i_level = 0,n_level-1
        array(n_level-i_level) = int(ibits(treecode, (i_level + max_tclevel - n_level)*n_dim, n_dim), ik)
    end do
  end subroutine
  !===============================================================================

  !===============================================================================
  !> \brief Obtain neighbour for given direction with numerical binary treecode
  !> \details This function takes the directions and splits it into the desire cardinal directions and calls function adjacent
  !> \author JB
  subroutine adjacent_wrapper_b(treecode, treecode_neighbor, direction, level, dim, max_level)
    implicit none
    !> dimension (2 or 3), defaults to 3
    integer(kind=ik), optional    :: dim
    !> Level at which to encode, can be negative to set from max_level, defaults to max_level
    integer(kind=ik), optional, intent(in)    :: level
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level
    !> Numerical treecode in
    integer(kind=tsize), intent(in)     :: treecode
    !> Numerical treecoude out
    integer(kind=tsize), intent(out)    :: treecode_neighbor
    !> direction for neighbor search - number where each digit represents a cardinal direction
    !> XYZ, values are 1 for dir +, 9 for dir - and 0 for nothing
    integer(kind=ik), intent(in)        :: direction
    integer(kind=tsize)                 :: tc1
    integer(kind=ik)                    :: i_dim, n_dim, n_level, max_tclevel
    integer(kind=tsize)                 :: dir_now

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    treecode_neighbor = treecode
    ! loop over all letters in direction and call the cardinal directions
    do i_dim = 1,3
      tc1 = treecode_neighbor
      dir_now = mod(direction/(10**(3-i_dim)),10)

      if (dir_now == 0) then
        cycle  ! do nothing and cycle
      elseif (dir_now == 9) then
        dir_now = 0
      endif
      dir_now = dir_now + 2*(i_dim-1) + 1
      ! dir_now is now between 1 and 6, being index of [+x, -x, +y, -y, +z, -z]
      call adjacent_faces_b(tc1, treecode_neighbor, dir_now, level=n_level, dim=n_dim, max_level=max_tclevel)
    end do
  end subroutine
  !===============================================================================

  !===============================================================================
  !> \author JB
  !> \brief Obtain neighbour in 3D for given direction with numerical binary treecode
  !> \details Use math representation and loop over digits
  !  --------------------------------------------------------------------------------------------
  !> neighbor codes: \n
  !  ---------------
  !> Each digit in variable direction represents one of the dimensions. Digits can take following values:
  !!    - 9, we look in negative direction in this dimension
  !!    - 1, we look in positive direction in this dimension
  !!    - 0, this dimension does not change
  !! Dependend on the number of 0s, we either look at a face, edge or corner
  ! ********************************************************************************************
  subroutine adjacent_faces_b( treecode, treecode_neighbor, direction, level, dim, max_level)
    implicit none

    !> dimension (2 or 3), defaults to 3
    integer(kind=ik), optional    :: dim  
    !> Level at which to encode, can be negative to set from max_level, defaults to max_level
    integer(kind=ik), optional, intent(in)    :: level
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level
    !> Numerical treecode in
    integer(kind=tsize), intent(in)     :: treecode
    !> Numerical treecoude out
    integer(kind=tsize), intent(out)    :: treecode_neighbor
    !> direction for neighbor search - number with one digit representing a cardinal direction
    integer(kind=tsize), intent(in)        :: direction
    integer(kind=tsize)                 :: tc_reduce, digit_last, dir_sign, dir_fac
    integer(kind=ik)                    :: i, n_dim, max_tclevel, n_level

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    ! We need direction and level from each direction
    select case(direction)
      ! At first we assign the main cardinal directions
      case(1)  ! N - left side x-1
        dir_sign = -1
        dir_fac = 2
      case(2)  ! S - right side x+1
        dir_sign = 1
        dir_fac = 2
      case(3)  ! W - front side y-1
        dir_sign = -1
        dir_fac = 1
      case(4)  ! E - back side y+1
        dir_sign = 1
        dir_fac = 1
      case(5)  ! B - bottom side z-1
        dir_sign = -1
        dir_fac = 4
      case(6)  ! T - top side z+1
        dir_sign = 1
        dir_fac = 4
      case default
        call abort(118118, "Lord vader, the treelib does not know the direction")
    end select

    ! copy treecode, as we modify it, but not return this modified value
    tc_reduce = treecode

    ! this is the neighbors treecode we're looking for
    treecode_neighbor = 0_tsize

    ! scales finer as neighbour search - keep as 0 and pop of numbers
    do i = max_tclevel-n_level+1, max_tclevel
      ! extract binary-duplet YX or triplet ZYX
      digit_last = ibits(treecode, (i-1)*n_dim, n_dim)

      ! add number with change from neighbour or overflow
      treecode_neighbor = treecode_neighbor + ishft(modulo(digit_last + dir_sign*dir_fac, 2*dir_fac) + digit_last/(2*dir_fac)*2*dir_fac, n_dim*(i-1))

      ! compute overflow
      dir_sign = ibits(digit_last/dir_fac + (2 - dir_sign)/2, 0, 1)*dir_sign

      ! copy directly the rest and exit loop if numbers are not gonna change anymore
      ! shift back and forth in order to set lower values to zero
      if (dir_sign == 0) then
        treecode_neighbor = treecode_neighbor + ishft(ishft(treecode, -n_dim*i), n_dim*i)
        exit
      end if
    end do

  end subroutine
  !===============================================================================

  !===============================================================================
  !> \brief Obtain neighbour for given direction with numerical decimal treecode
  !> \details This function takes the directions and splits it into the desire cardinal directions and calls function adjacent
  !> \author JB
  subroutine adjacent_wrapper_d(treecode, treecode_neighbor, direction, level, max_level)
    implicit none
    !> Level at which to search the neighbour, this is from coarse to fine
    integer(kind=ik), optional, intent(in)        :: level
    !> Max treelevel, needed to loop correctly
    integer(kind=ik), optional, intent(in)        :: max_level
    !> Numerical treecode in
    integer(kind=tsize), intent(in)     :: treecode
    !> Numerical treecoude out
    integer(kind=tsize), intent(out)    :: treecode_neighbor
    !> direction for neighbor search - number where each digit represents a cardinal direction
    !> XYZ, values are 1 for dir +, 9 for dir - and 0 for nothing
    integer(kind=ik), intent(in)        :: direction
    integer(kind=tsize)                 :: tc1
    integer(kind=ik)                    :: i_dim, n_level, max_tclevel
    integer(kind=tsize)                 :: dir_now

    ! Set defaults for dimension, level and max_level
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    treecode_neighbor = treecode
    ! loop over all letters in direction and call the cardinal directions
    do i_dim = 1,3
      tc1 = treecode_neighbor
      dir_now = mod(direction/(10**(3-i_dim)),10)

      if (dir_now == 0) then
        cycle  ! do nothing and cycle
      elseif (dir_now == 9) then
        dir_now = 0
      endif
      dir_now = dir_now + 2*(i_dim-1) + 1
      ! dir_now is now between 1 and 6, being index of [+x, -x, +y, -y, +z, -z]
      call adjacent_faces_d(tc1, treecode_neighbor, dir_now, level=n_level, max_level=max_tclevel)
    end do
  end subroutine
  !===============================================================================

  !===============================================================================
  !> \author JB
  !> \brief Obtain neighbour in 3D for given direction with numerical decimal treecode
  !> \details Use math representation and loop over digits
  !  --------------------------------------------------------------------------------------------
  !> neighbor codes: \n
  !  ---------------
  !> Each digit in variable direction represents one of the dimensions. Digits can take following values:
  !!    - 9, we look in negative direction in this dimension
  !!    - 1, we look in positive direction in this dimension
  !!    - 0, this dimension does not change
  !! Dependend on the number of 0s, we either look at a face, edge or corner
  ! ********************************************************************************************
  subroutine adjacent_faces_d( treecode, treecode_neighbor, direction, level, max_level)
    implicit none
    !> Level at which to search the neighbour, this is from coarse to fine
    integer(kind=ik), optional, intent(in)        :: level
    !> Max treelevel, needed to loop correctly
    integer(kind=ik), optional, intent(in)        :: max_level
    !> Numerical treecode in
    integer(kind=tsize), intent(in)     :: treecode
    !> Numerical treecoude out
    integer(kind=tsize), intent(out)    :: treecode_neighbor
    !> direction for neighbor search - number with one digit representing a cardinal direction
    integer(kind=tsize), intent(in)        :: direction
    integer(kind=tsize)                 :: tc_reduce, digit_last, dir_sign, dir_fac
    integer(kind=ik)                    :: i, n_level, max_tclevel

    ! Set defaults for dimension, level and max_level
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    ! We need direction and level from each direction
    select case(direction)
      ! At first we assign the main cardinal directions
      case(1)  ! N - left side x-1
        dir_sign = -1
        dir_fac = 2
      case(2)  ! S - right side x+1
        dir_sign = 1
        dir_fac = 2
      case(3)  ! W - front side y-1
        dir_sign = -1
        dir_fac = 1
      case(4)  ! E - back side y+1
        dir_sign = 1
        dir_fac = 1
      case(5)  ! B - bottom side z-1
        dir_sign = -1
        dir_fac = 4
      case(6)  ! T - top side z+1
        dir_sign = 1
        dir_fac = 4
      case default
        call abort(118118, "Lord vader, the treelib does not know the direction")
    end select

    ! copy treecode, as we modify it, but not return this modified value
    tc_reduce = treecode

    ! this is the neighbors treecode we're looking for
    treecode_neighbor = 0_tsize

    ! scales finer as neighbour search - keep as 0 and pop of numbers
    do i = 1, max_tclevel-n_level
      call pop(tc_reduce, digit_last)
    end do
    do i = max_tclevel-n_level+1, maxdigits
      ! pop last digit
      call pop(tc_reduce, digit_last)
      ! add number with change from neighbour or overflow
      treecode_neighbor = treecode_neighbor + (modulo(digit_last + dir_sign*dir_fac, 2*dir_fac) + digit_last/(2*dir_fac)*2*dir_fac)* 10_tsize**(i-1)
      
      ! compute overflow as -1 or 1
      dir_sign = modulo(digit_last/dir_fac + (2 - dir_sign)/2, 2_tsize)*dir_sign

      ! copy directly the rest and exit loop if numbers are not gonna change anymore
      if (dir_sign == 0) then
        treecode_neighbor = treecode_neighbor + tc_reduce* 10_tsize**i
        exit
      end if
    end do

  end subroutine
  !===============================================================================

  !===============================================================================
  !> \brief Obtain block position coordinates from numerical decimal treecode
  !> \details Works for 2D and 3D. Considers each digit and adds their level-shift to each coordinate
  subroutine decoding_d(ix, treecode, dim, level, max_level)
    implicit none

    !> dimension (2 or 3), defaults to 3
    integer(kind=ik), optional, intent(in)    :: dim              
    !> block position coordinates
    ! set to fixed size of 3, so ensure that for dim=2 the third will not be accessed or we get oob
    integer(kind=ik), intent(out)    :: ix(3)
    !> Treecode
    integer(kind=tsize), intent(in) :: treecode
    !> Level at which to encode, can be negative to set from max_level, defaults to max_level
    integer(kind=ik), optional, intent(in)    :: level
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level

    integer(kind=tsize) :: tc_digit, tc_temp
    integer(kind=ik) :: i_dim, n_dim, max_tclevel, i_level, n_level

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    tc_temp = treecode

    ! NOTE: one-based indexing
    do i_dim = 1,n_dim
      ix(i_dim) = 1
    end do

    ! skip first digits
    do i_level = 0, max_tclevel - n_level - 1
      call pop(tc_temp, tc_digit)
    end do
    ! extract corresponding bit from bit-duplet YX or -triplet ZYX of level i_nx, then multiply by 10**(i_nx-1)
    do i_level = 0, n_level-1
      call pop(tc_temp, tc_digit)
      do i_dim = 1,n_dim
        ix(i_dim) = ix(i_dim) + int(ibits(tc_digit, i_dim-1, 1) * 2_tsize**i_level, kind=ik)
      end do
    end do

    ! IY is encoded in first entry and IX in second so we need to swap them
    i_level = ix(1); ix(1) = ix(2); ix(2) = i_level

  end subroutine decoding_d
  !===============================================================================

  !===============================================================================
  !> \brief Obtain numerical decimal treecode from block position coordinates for 2 or 3 dimensions
  !> \author JB
  subroutine encoding_d( ix, treecode, dim, level, max_level )
    implicit none
    !> dimension (2 or 3), defaults to 3
    integer(kind=ik), optional, intent(in)    :: dim
    !> Level at which to encode, can be negative to set from max_level, defaults to max_level
    integer(kind=ik), optional, intent(in)    :: level
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level
    !> block position coordinates
    ! set to fixed size of 3, so ensure that for dim=2 the third will not be accessed or we get oob
    integer(kind=ik), intent(in)    :: ix(3)
    !> treecode
    integer(kind=tsize), intent(out) :: treecode
    integer(kind=ik) :: i_dim, n_dim, i_level, n_level, max_tclevel, ix_p(3)
    !integer(kind=ik) :: b(Jmax)

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    ! IY is encoded in first entry and IX in second so we need to swap them
    ix_p(1:n_dim) = ix(1:n_dim)
    i_level = ix_p(1); ix_p(1) = ix_p(2); ix_p(2) = i_level

    ! NOTE: gargantini uses 0-based indexing, which is a source of errors. we use 1-based (-1)
    treecode = 0_tsize
    do i_dim = 1, n_dim
      ! loop over all bits in ix which are set
      do i_level = 0, bit_size(ix_p(i_dim))-leadz(ix_p(i_dim))
        treecode = treecode + ibits(ix_p(i_dim)-1, i_level, 1) * 2_tsize**(i_dim-1) * 10_tsize**(i_level + (max_tclevel - n_level))
      end do
    end do

  end subroutine encoding_d
  !===============================================================================

  !===============================================================================
  !> \brief Obtain block position coordinates from numerical binary treecode
  !> \details Works for 2D and 3D. Considers each digit and adds their level-shift to each coordinate
  subroutine decoding_b(ix, treecode, dim, level, max_level)
    implicit none

    !> dimension (2 or 3), defaults to 3
    integer(kind=ik), optional, intent(in)    :: dim              
    !> block position coordinates
    ! set to fixed size of 3, so ensure that for dim=2 the third will not be accessed or we get oob
    integer(kind=ik), intent(out)    :: ix(:)
    !> Level at which to encode, can be negative to set from max_level, defaults to max_level
    integer(kind=ik), optional, intent(in)    :: level
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level
    !> treecode
    integer(kind=tsize), intent(in) :: treecode
    integer(kind=ik) :: n_dim, i_dim, i_level, max_tclevel, n_level

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    ! NOTE: one-based indexing
    do i_dim = 1,n_dim
      ix(i_dim) = 1
    end do

    ! extract corresponding bit from bit-duplet YX or -triplet ZYX of level i_nx, then multiply by 2**(i_nx-1)
    ! convert from int8 to int4 is not necessary but gets rid of conversion warning
    do i_level = 0, n_level-1
      do i_dim = 1,n_dim
        ix(i_dim) = ix(i_dim) + int(ishft(ibits(treecode, (i_level + max_tclevel - n_level)*n_dim+(i_dim-1), 1), i_level), kind=ik)
      end do
    end do

    ! IY is encoded in first entry and IX in second so we need to swap them
    i_level = ix(1); ix(1) = ix(2); ix(2) = i_level

  end subroutine decoding_b
  !===============================================================================

  !===============================================================================
  !> \brief Obtain numerical binary treecode from block position coordinates for 2 or 3 dimensions
  !> \author JB
  subroutine encoding_b( ix, treecode, dim, level, max_level) !, Jmax )
    implicit none
    !> dimension (2 or 3), defaults to 3
    integer(kind=ik), optional    :: dim              
    !> block position coordinates
    ! set to fixed size of 3, so ensure that for dim=2 the third will not be accessed or we get oob
    integer(kind=ik), intent(in)    :: ix(:)
    !> Level at which to encode, can be negative to set from max_level, defaults to max_level
    integer(kind=ik), optional, intent(in)    :: level
    !> Max level possible, should be set after params%Jmax
    integer(kind=ik), optional, intent(in)    :: max_level
    !> treecode
    integer(kind=tsize), intent(out) :: treecode
    integer(kind=ik) :: n_dim, i_dim, i_level, max_tclevel, n_level, ix_p(3)

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    treecode = 0_tsize
    ! IY is encoded in first entry and IX in second so we need to swap them
    ix_p(1:n_dim) = ix(1:n_dim)
    i_level = ix_p(1); ix_p(1) = ix_p(2); ix_p(2) = i_level

    ! extract value as bit from each direction and add correspondend position in level-triplet or -duplet to treecode
    ! subtract 1 as it is one-based indexing
    ! loop over all bits set in index
    do i_dim = 1, n_dim
      do i_level = 0,  bit_size(ix_p(i_dim)) - leadz(ix_p(i_dim)) -1
        treecode = treecode + ishft(int(ibits(ix_p(i_dim)-1, i_level,1), kind=tsize), (i_level + (max_tclevel - n_level))*n_dim + (i_dim - 1))
      end do
    end do
  end subroutine encoding_b
  !===============================================================================

  !===============================================================================
  !> \brief Convert numerical binary treecode to str for readability where each digit is from 0-7
  !> \details Str representation is needed as max length of 31/21 exceeds maximum digits a decimal representation can have
  !> \author JB
  subroutine tc_to_str( treecode, tc_str, dim, level, max_level) !, Jmax )
    implicit none
    !> Treecode in binary numerical representation
    integer(kind=tsize), intent(in) :: treecode
    !> String where each digit is 0-7
    character(len=*), intent(out)   :: tc_str
    !> Dimension (2 or 3) - defaults to 3
    integer(kind=ik), optional, intent(in)    :: dim
    !> Level, cut trailing zeros - defaults to max_level
    integer(kind=ik), optional, intent(in)    :: level
    !> Maximum level, should be params%Jmax, defaults to maxdigits
    integer(kind=ik), optional, intent(in) :: max_level
    
    character(len=:), allocatable :: temp_str ! allocatable as variable length
    integer(kind=tsize) :: i_dim, i_nx, temp, n_level, n_dim, max_tclevel

    ! Set defaults for dimension, level and max_level
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = max_tclevel; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    ! Dynamically allocate temp_str based on n_level and init
    allocate(character(len=n_level) :: temp_str)
    temp_str = ""

    do i_nx = max_tclevel-n_level, max_tclevel-1
        ! extract bit-triplet ZYX on current level
        temp = ibits(treecode, 3*i_nx, 3)

        ! convert to str and add at correct position
        temp_str = achar(temp + ichar('0')) // temp_str
    end do
    ! append X for levels unset
    do i_nx = n_level, max_tclevel -1
      temp_str = temp_str // 'X'
    end do

    tc_str = temp_str
    ! allocated needs to be deallocated
    deallocate(temp_str)
  end subroutine tc_to_str
  !===============================================================================

  !===============================================================================
  !> \brief Convert given integer decnum to binary representation
  !> \details Each digit represents a power of 2 and can be either 1 or 0, output is as integer.
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
  !===============================================================================

  !===============================================================================
  !> \brief Flip an integer number
  !! \details Flip an integer number, eg given 523 return 325
  !! or 000021 gives 120000, so we have to use zero-padding! \n
  !! The length of returned (zero padded) reverse int is specified by maxdigits
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
  !===============================================================================



  !> \brief give treecode for adjacent block \n
  !! \details input:
  !!   - treecode for block N
  !!   - direction for neighbor search
  !!   - max treelevel
  !!
  !! output:
  !!   - neighbor treecode, for neighbor on same level
  ! ********************************************************************************************
  recursive subroutine adjacent_block_2D(tcBlock, tcNeighbor, direction, level, max_treelevel)

      implicit none
      integer(kind=ik), intent(in)        :: max_treelevel
      integer(kind=ik), intent(in)        :: level
      integer(kind=ik), intent(in)        :: tcBlock(max_treelevel)          !> block treecode
      character(len=3), intent(in)        :: direction                  !> direction for neighbor search
      integer(kind=ik), intent(out)       :: tcNeighbor(max_treelevel)    !> neighbor treecode
      integer(kind=ik)                    :: tcNeighbor2(max_treelevel)   ! treecode variable
      integer(kind=ik)                    :: i                          ! loop variable

      tcNeighbor = -1

      select case(direction)
      case('__N')
          ! north
          tcNeighbor(level) = modulo(tcBlock(level)+2, 4)
          i = level - 1

          do while (i /= 0)
              if ( (tcBlock(i+1)==2) .or. (tcBlock(i+1)==3) ) then
                  tcNeighbor(1:i) = tcBlock(1:i)
                  i = 0
              else
                  tcNeighbor(i) = modulo(tcBlock(i)+2, 4)
                  i = i - 1
              end if
          end do

      case('__S')
          ! south
          tcNeighbor(level) = modulo(tcBlock(level)+2, 4)
          i = level - 1

          do while (i /= 0)
              if ( (tcBlock(i+1)==0) .or. (tcBlock(i+1)==1) ) then
                  tcNeighbor(1:i) = tcBlock(1:i)
                  i = 0
              else
                  tcNeighbor(i) = modulo(tcBlock(i)+2, 4)
                  i = i - 1
              end if
          end do

      case('__E')
          ! east
          if ( (tcBlock(level)==0) .or. (tcBlock(level)==2) ) then
              tcNeighbor(level) = modulo(tcBlock(level)+1, 4)
          else
              tcNeighbor(level) = modulo(tcBlock(level)-1, 4)
          end if

          i = level - 1

          do while (i /= 0)
              if ( (tcBlock(i+1)==0) .or. (tcBlock(i+1)==2) ) then
                  tcNeighbor(1:i) = tcBlock(1:i)
                  i = 0
              else
                  if ( (tcBlock(i)==1) .or. (tcBlock(i)==3) ) then
                      tcNeighbor(i) = modulo(tcBlock(i)-1, 4)
                  else
                      tcNeighbor(i) = modulo(tcBlock(i)+1, 4)
                  end if
                  i = i - 1
              end if
          end do

      case('__W')
          ! west
          if ( (tcBlock(level)==0) .or. (tcBlock(level)==2) ) then
              tcNeighbor(level) = modulo(tcBlock(level)+1, 4)
          else
              tcNeighbor(level) = modulo(tcBlock(level)-1, 4)
          end if

          i = level - 1

          do while (i /= 0)
              if ( (tcBlock(i+1)==1) .or. (tcBlock(i+1)==3) ) then
                  tcNeighbor(1:i) = tcBlock(1:i)
                  i = 0
              else
                  if ( (tcBlock(i)==1) .or. (tcBlock(i)==3) ) then
                      tcNeighbor(i) = modulo(tcBlock(i)-1, 4)
                  else
                      tcNeighbor(i) = modulo(tcBlock(i)+1, 4)
                  end if
                  i = i - 1
              end if
          end do

      case('_NE')
          ! northeast
          call adjacent_block_2D(tcBlock, tcNeighbor2, '__N', level, max_treelevel)
          call adjacent_block_2D(tcNeighbor2, tcNeighbor, '__E', level, max_treelevel)

      case('_NW')
          ! northwest
          call adjacent_block_2D(tcBlock, tcNeighbor2, '__N', level, max_treelevel)
          call adjacent_block_2D(tcNeighbor2, tcNeighbor, '__W', level, max_treelevel)

      case('_SE')
          ! southeast
          call adjacent_block_2D(tcBlock, tcNeighbor2, '__S', level, max_treelevel)
          call adjacent_block_2D(tcNeighbor2, tcNeighbor, '__E', level, max_treelevel)

      case('_SW')
          ! southwest
          call adjacent_block_2D(tcBlock, tcNeighbor2, '__S', level, max_treelevel)
          call adjacent_block_2D(tcNeighbor2, tcNeighbor, '__W', level, max_treelevel)

      case default
          call abort(118119, "Lord vader, the treelib does not know the direction")
      end select
  end subroutine adjacent_block_2D


  !> \brief give treecode for adjacent block in 3D \n
  !! \details input:
  !!                    - treecode for block N
  !!                    - direction for neighbor search
  !!                    - max treelevel
  !! output:             
  !!
  !!                    - neighbor treecode, for neighbor on same level
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
  ! ********************************************************************************************
  recursive subroutine adjacent_block_3D(tcBlock, tcNeighbor, direction, level, max_treelevel)

      implicit none
      integer(kind=ik), intent(in)        :: max_treelevel
      integer(kind=ik), intent(in)        :: level
      integer(kind=ik), intent(in)        :: tcBlock(max_treelevel)      !> block treecode
      character(len=7), intent(in)        :: direction              !> direction for neighbor search
      integer(kind=ik), intent(out)       :: tcNeighbor(max_treelevel)!> neighbor treecode
      integer(kind=ik)                    :: tcNeighbor2(max_treelevel), tcNeighbor3(max_treelevel)   ! treecode variable
      integer(kind=ik)                    :: i                      ! loop variable

      tcNeighbor = tcBlock

      select case(direction)

      case('__1/___')
          ! '1'-side
          tcNeighbor(level) = modulo(tcBlock(level) + 4, 8)
          i = level - 1

          do while (i /= 0)
              select case( tcBlock(i+1) )
              case (0,1,2,3)
                  ! nothing to do, leave loop
                  i = 0

              case (4,5,6,7)
                  tcNeighbor(i) = modulo(tcBlock(i) + 4, 8)
                  i = i - 1

              end select
          end do

      case('__6/___')
          ! '6'-side
          tcNeighbor(level) = modulo(tcBlock(level) - 4, 8)
          i = level - 1

          do while (i /= 0)
              select case( tcBlock(i+1) )
              case (4,5,6,7)
                  ! nothing to do, leave loop
                  i = 0

              case (0,1,2,3)
                  tcNeighbor(i) = modulo(tcBlock(i) - 4, 8)
                  i = i - 1

              end select
          end do

      case('__3/___')
          ! '3'-side
          select case(tcBlock(level))
          case(0,1,4,5)
              tcNeighbor(level) = modulo(tcBlock(level) + 2, 8)
          case(2,3,6,7)
              tcNeighbor(level) = modulo(tcBlock(level) + 6, 8)
          end select

          i = level - 1

          do while (i /= 0)
              select case( tcBlock(i+1) )
              case (0,1,4,5)
                  ! nothing to do, leave loop
                  i = 0

              case (2,3,6,7)
                  select case(tcBlock(i))
                  case (0,1,4,5)
                      tcNeighbor(i) = modulo(tcBlock(i) + 2, 8)
                  case (2,3,6,7)
                      tcNeighbor(i) = modulo(tcBlock(i) + 6, 8)
                  end select
                  i = i - 1

              end select
          end do

      case('__5/___')
          ! '5'-side
          select case(tcBlock(level))
          case(2,3,6,7)
              tcNeighbor(level) = modulo(tcBlock(level) - 2, 8)
          case(0,1,4,5)
              tcNeighbor(level) = modulo(tcBlock(level) - 6, 8)
          end select

          i = level - 1

          do while (i /= 0)
              select case( tcBlock(i+1) )
              case (2,3,6,7)
                  ! nothing to do, leave loop
                  i = 0

              case (0,1,4,5)
                  select case(tcBlock(i))
                  case (2,3,6,7)
                      tcNeighbor(i) = modulo(tcBlock(i) - 2, 8)
                  case (0,1,4,5)
                      tcNeighbor(i) = modulo(tcBlock(i) - 6, 8)
                  end select
                  i = i - 1

              end select
          end do

      case('__2/___')
          ! '2'-side
          select case(tcBlock(level))
          case(1,3,5,7)
              tcNeighbor(level) = modulo(tcBlock(level) - 1, 8)
          case(0,2,4,6)
              tcNeighbor(level) = modulo(tcBlock(level) - 7, 8)
          end select

          i = level - 1

          do while (i /= 0)
              select case( tcBlock(i+1) )
              case(1,3,5,7)
                  ! nothing to do, leave loop
                  i = 0

              case(0,2,4,6)
                  select case(tcBlock(i))
                  case(1,3,5,7)
                      tcNeighbor(i) = modulo(tcBlock(i) - 1, 8)
                  case(0,2,4,6)
                      tcNeighbor(i) = modulo(tcBlock(i) - 7, 8)
                  end select
                  i = i - 1

              end select
          end do

      case('__4/___')
          ! '4'-side
          select case(tcBlock(level))
          case(0,2,4,6)
              tcNeighbor(level) = modulo(tcBlock(level) + 1, 8)
          case(1,3,5,7)
              tcNeighbor(level) = modulo(tcBlock(level) + 7, 8)
          end select

          i = level - 1

          do while (i /= 0)
              select case( tcBlock(i+1) )
              case(0,2,4,6)
                  ! nothing to do, leave loop
                  i = 0

              case(1,3,5,7)
                  select case(tcBlock(i))
                  case(0,2,4,6)
                      tcNeighbor(i) = modulo(tcBlock(i) + 1, 8)
                  case(1,3,5,7)
                      tcNeighbor(i) = modulo(tcBlock(i) + 7, 8)
                  end select
                  i = i - 1

              end select
          end do

      case('_12/___')
          ! '12'-edge
          call adjacent_block_3D(tcBlock, tcNeighbor2, '__1/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__2/___', level, max_treelevel)

      case('_13/___')
          ! '13'-edge
          call adjacent_block_3D(tcBlock, tcNeighbor2, '__1/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__3/___', level, max_treelevel)

      case('_14/___')
          ! '14'-edge
          call adjacent_block_3D(tcBlock, tcNeighbor2, '__1/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__4/___', level, max_treelevel)

      case('_15/___')
          ! '15'-edge
          call adjacent_block_3D(tcBlock, tcNeighbor2, '__1/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__5/___', level, max_treelevel)

      case('_62/___')
          ! '62'-edge
          call adjacent_block_3D(tcBlock, tcNeighbor2, '__6/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__2/___', level, max_treelevel)

      case('_63/___')
          ! '63'-edge
          call adjacent_block_3D(tcBlock, tcNeighbor2, '__6/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__3/___', level, max_treelevel)

      case('_64/___')
          ! '64'-edge
          call adjacent_block_3D(tcBlock, tcNeighbor2, '__6/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__4/___', level, max_treelevel)

      case('_65/___')
          ! '65'-edge
          call adjacent_block_3D(tcBlock, tcNeighbor2, '__6/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__5/___', level, max_treelevel)

      case('_23/___')
          ! '23'-edge
          call adjacent_block_3D(tcBlock, tcNeighbor2, '__2/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__3/___', level, max_treelevel)

      case('_25/___')
          ! '25'-edge
          call adjacent_block_3D(tcBlock, tcNeighbor2, '__2/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__5/___', level, max_treelevel)

      case('_43/___')
          ! '43'-edge
          call adjacent_block_3D(tcBlock, tcNeighbor2, '__4/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__3/___', level, max_treelevel)

      case('_45/___')
          ! '45'-edge
          call adjacent_block_3D(tcBlock, tcNeighbor2, '__4/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__5/___', level, max_treelevel)

      case('123/___')
          ! '123'-corner
          call adjacent_block_3D(tcBlock, tcNeighbor3, '__1/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor3, tcNeighbor2, '__2/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__3/___', level, max_treelevel)

      case('134/___')
          ! '134'-corner
          call adjacent_block_3D(tcBlock, tcNeighbor3, '__1/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor3, tcNeighbor2, '__3/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__4/___', level, max_treelevel)

      case('145/___')
          ! '145'-corner
          call adjacent_block_3D(tcBlock, tcNeighbor3, '__1/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor3, tcNeighbor2, '__4/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__5/___', level, max_treelevel)

      case('152/___')
          ! '152'-corner
          call adjacent_block_3D(tcBlock, tcNeighbor3, '__1/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor3, tcNeighbor2, '__5/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__2/___', level, max_treelevel)

      case('623/___')
          ! '623'-corner
          call adjacent_block_3D(tcBlock, tcNeighbor3, '__6/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor3, tcNeighbor2, '__2/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__3/___', level, max_treelevel)

      case('634/___')
          ! '634'-corner
          call adjacent_block_3D(tcBlock, tcNeighbor3, '__6/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor3, tcNeighbor2, '__3/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__4/___', level, max_treelevel)

      case('645/___')
          ! '645'-corner
          call adjacent_block_3D(tcBlock, tcNeighbor3, '__6/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor3, tcNeighbor2, '__4/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__5/___', level, max_treelevel)

      case('652/___')
          ! '652'-corner
          call adjacent_block_3D(tcBlock, tcNeighbor3, '__6/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor3, tcNeighbor2, '__5/___', level, max_treelevel)
          call adjacent_block_3D(tcNeighbor2, tcNeighbor, '__2/___', level, max_treelevel)

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

      implicit none

      integer(kind=ik), intent(out)  :: i, j, k          !> block position coordinates
      integer(kind=ik), intent(in)   :: treeN            !> treecode size
      integer(kind=ik), intent(in)   :: treecode(treeN)
      integer(kind=ik)               :: nx, step, l

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
      implicit none
      integer(kind=ik), intent(in)    :: dim              !> dimension (2 or 3)
      integer(kind=ik), intent(in)    :: ix(dim)          !> block position coordinates
      integer(kind=ik), intent(in)    :: block_num        !> number of blocks
      integer(kind=ik), intent(in)    :: treeN            !> treecode size
      integer(kind=ik), intent(out)   :: treearray(treeN) !> treecode

      ! variables for calculate real treecode length N
      real(kind=rk)                   :: Jn
      integer(kind=ik)                :: N, l,d
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

  !===============================================================================
  !> \brief Encoding from block position (cartesian coordinates) to treecode with different input
  !
  ! JB: Redundant and only used once in post_extract_slice
  subroutine encoding_revised(treecode, ix, dim, level)
      implicit none
      integer(kind=ik), intent(in)    :: dim            !> dimension (2 or 3)
      integer(kind=ik), intent(in)    :: ix(dim)        !> block position coordinates
      integer(kind=ik), intent(in)    :: level
      integer(kind=ik), intent(out)   :: treecode(1:)

      ! variables for calculate real treecode length N
      integer(kind=ik)                :: treeN
      real(kind=rk)                   :: Jn
      integer(kind=ik)                :: l,d
      integer(kind=ik), allocatable   :: ix_binary(:)

      treeN = size(treecode)

      ! reset output
      treecode = 0
      ! allocate auxiliary vectors
      allocate( ix_binary(level) )

      ! convert block coordinates into binary numbers
      do d = 1, dim
          call int_to_binary(ix(d)-1, level, ix_binary)
          treecode(1:level) = treecode(1:level) + ix_binary(1:level)*2**(d-1)
      end do

      deallocate( ix_binary )
  end subroutine encoding_revised
  !===============================================================================

  !> \brief convert a integer i to binary b \n
  !> \details binary return as vector with length N
  subroutine int_to_binary(i, N, b)
      implicit none
      integer(kind=ik), intent(in)    :: i            !> integer to convert into binary
      integer(kind=ik), intent(in)    :: N            !> length of binary output vector
      integer(kind=ik), intent(out)   :: b(N)         !> output vector
      integer(kind=ik)                :: j, k, tmp(N) ! loop variables

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
