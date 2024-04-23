! subroutine get_neighbor_treecode( tcBlock, neighbor_treecode, direction, &
!     leveldiff, J, dim, max_treelevel, level_down_neighbor_possible)
!     implicit none

!     integer(kind=ik), intent(in)  :: direction, leveldiff, J, dim, max_treelevel
!     integer(kind=ik), intent(in)  :: tcBlock(max_treelevel)
!     integer(kind=ik), intent(out) :: neighbor_treecode(max_treelevel)
!     logical, intent(out) :: level_down_neighbor_possible

!     integer(kind=ik) :: level
!     integer(kind=ik) :: virt_treecode(max_treelevel)
!     integer(kind=ik) :: virtual_treecode_finer_level(1:74) = -9
!     character(len=3) :: direction_names_2d(16) = "xxx"
!     character(len=7) :: direction_names_3d(74) = "xxxxxxx"
!     logical :: array_compare


!     level = J

!     !-------2D
!     direction_names_2d(1)    = '__N'
!     direction_names_2d(9:10) = '__N'
!     virtual_treecode_finer_level(10) = 0
!     virtual_treecode_finer_level(9) = 1

!     direction_names_2d(2)     = '__E'
!     direction_names_2d(13:14) = '__E'
!     virtual_treecode_finer_level(13) = 1
!     virtual_treecode_finer_level(14) = 3

!     direction_names_2d(3)     = '__S'
!     direction_names_2d(11:12) = '__S'
!     virtual_treecode_finer_level(12) = 2
!     virtual_treecode_finer_level(11) = 3

!     direction_names_2d(4)     = '__W'
!     direction_names_2d(15:16) = '__W'
!     virtual_treecode_finer_level(15) = 0
!     virtual_treecode_finer_level(16) = 2

!     direction_names_2d(5) = '_NE'
!     virtual_treecode_finer_level(5) = 1

!     direction_names_2d(6) = '_NW'
!     virtual_treecode_finer_level(6) = 0

!     direction_names_2d(7) = '_SE'
!     virtual_treecode_finer_level(7) = 3

!     direction_names_2d(8) = '_SW'
!     virtual_treecode_finer_level(8) = 2

!     !------- faces 3d
!     direction_names_3d(1) = '__1/___'
!     direction_names_3d(27:30) = '__1/___'
!     virtual_treecode_finer_level(30) = 4
!     virtual_treecode_finer_level(29) = 5
!     virtual_treecode_finer_level(27) = 6
!     virtual_treecode_finer_level(28) = 7

!     direction_names_3d(2) = '__2/___'
!     direction_names_3d(31:34) = '__2/___'
!     virtual_treecode_finer_level(34) = 0
!     virtual_treecode_finer_level(32) = 2
!     virtual_treecode_finer_level(33) = 4
!     virtual_treecode_finer_level(31) = 6

!     direction_names_3d(3) = '__3/___'
!     direction_names_3d(35:38) = '__3/___'
!     virtual_treecode_finer_level(36) = 2
!     virtual_treecode_finer_level(38) = 3
!     virtual_treecode_finer_level(35) = 6
!     virtual_treecode_finer_level(37) = 7

!     direction_names_3d(4) = '__4/___'
!     direction_names_3d(39:42) = '__4/___'
!     virtual_treecode_finer_level(42) = 1
!     virtual_treecode_finer_level(40) = 3
!     virtual_treecode_finer_level(41) = 5
!     virtual_treecode_finer_level(39) = 7

!     direction_names_3d(5) = '__5/___'
!     direction_names_3d(43:46) = '__5/___'
!     virtual_treecode_finer_level(46) = 0
!     virtual_treecode_finer_level(44) = 1
!     virtual_treecode_finer_level(45) = 4
!     virtual_treecode_finer_level(43) = 5

!     direction_names_3d(6) = '__6/___'
!     direction_names_3d(47:50) = '__6/___'
!     virtual_treecode_finer_level(50) = 0
!     virtual_treecode_finer_level(49) = 1
!     virtual_treecode_finer_level(47) = 2
!     virtual_treecode_finer_level(48) = 3

!     ! find_neighbor_edge_3D
!     direction_names_3d(7) = '_12/___'
!     virtual_treecode_finer_level(52) = 4
!     direction_names_3d(52) ='_12/___'
!     virtual_treecode_finer_level(51) = 6
!     direction_names_3d(51) = '_12/___'

!     direction_names_3d(8) = '_13/___'
!     virtual_treecode_finer_level(53) = 6
!     direction_names_3d(53) ='_13/___'
!     virtual_treecode_finer_level(54) = 7
!     direction_names_3d(54)= '_13/___'

!     direction_names_3d(9) = '_14/___'
!     virtual_treecode_finer_level(56) = 5
!     direction_names_3d(56) ='_14/___'
!     virtual_treecode_finer_level(55) = 7
!     direction_names_3d(55) = '_14/___'

!     direction_names_3d(10) = '_15/___'
!     virtual_treecode_finer_level(58) = 4
!     direction_names_3d(58) = '_15/___'
!     virtual_treecode_finer_level(57) = 5
!     direction_names_3d(57) = '_15/___'

!     direction_names_3d(11) = '_62/___'
!     virtual_treecode_finer_level(60) = 0
!     direction_names_3d(60) = '_62/___'
!     virtual_treecode_finer_level(59) = 2
!     direction_names_3d(59) = '_62/___'

!     direction_names_3d(12) = '_63/___'
!     virtual_treecode_finer_level(61) = 2
!     direction_names_3d(61) = '_63/___'
!     virtual_treecode_finer_level(62) = 3
!     direction_names_3d(62) = '_63/___'

!     direction_names_3d(13) = '_64/___'
!     virtual_treecode_finer_level(64) = 1
!     direction_names_3d(64) = '_64/___'
!     virtual_treecode_finer_level(63) = 3
!     direction_names_3d(63) = '_64/___'

!     direction_names_3d(14) = '_65/___'
!     virtual_treecode_finer_level(66) = 0
!     direction_names_3d(66) = '_65/___'
!     virtual_treecode_finer_level(65) = 1
!     direction_names_3d(65) = '_65/___'

!     direction_names_3d(15) = '_23/___'
!     virtual_treecode_finer_level(68) = 2
!     direction_names_3d(68) = '_23/___'
!     virtual_treecode_finer_level(67) = 6
!     direction_names_3d(67) = '_23/___'

!     direction_names_3d(16) = '_25/___'
!     virtual_treecode_finer_level(70) = 0
!     direction_names_3d(70) = '_25/___'
!     virtual_treecode_finer_level(69) = 4
!     direction_names_3d(69) = '_25/___'

!     direction_names_3d(17) = '_43/___'
!     virtual_treecode_finer_level(72) = 3
!     direction_names_3d(72) = '_43/___'
!     virtual_treecode_finer_level(71) = 7
!     direction_names_3d(71) = '_43/___'

!     direction_names_3d(18) = '_45/___'
!     virtual_treecode_finer_level(74) = 1
!     direction_names_3d(74) = '_45/___'
!     virtual_treecode_finer_level(73) = 5
!     direction_names_3d(73) = '_45/___'

!     ! corners 3d
!     direction_names_3d(19) = '123/___'
!     virtual_treecode_finer_level(19) = 6
!     direction_names_3d(20) = '134/___'
!     virtual_treecode_finer_level(20) = 7
!     direction_names_3d(21) = '145/___'
!     virtual_treecode_finer_level(21) = 5
!     direction_names_3d(22) = '152/___'
!     virtual_treecode_finer_level(22) = 4
!     direction_names_3d(23) = '623/___'
!     virtual_treecode_finer_level(23) = 2
!     direction_names_3d(24) = '634/___'
!     virtual_treecode_finer_level(24) = 3
!     direction_names_3d(25) = '645/___'
!     virtual_treecode_finer_level(25) = 1
!     direction_names_3d(26) = '652/___'
!     virtual_treecode_finer_level(26) = 0

!     level = J


!     level_down_neighbor_possible = .true.

!     if ( leveldiff == 0 .or. leveldiff == +1 ) then

!         ! calculate treecode for neighbor on same level
!         if (dim==3) then
!             call adjacent_block_3D( tcBlock, neighbor_treecode, direction_names_3d(direction), &
!             level, max_treelevel)
!         else
!             call adjacent_block_2D( tcBlock, neighbor_treecode, direction_names_2d(direction), &
!             level, max_treelevel)
!         endif

!         ! coarser level
!         if (leveldiff == +1) then
!             ! just remove the last digit: go one level up (coarser)
!             neighbor_treecode( level ) = -1
!             ! is it possible to have a neighbor on that level at all?
!             ! in the coarser neighbor case, finding a valid treecode can be tricky. not all treecodes
!             ! have a VALID coarser neighbor in the specified direction.
!             if (array_compare( tcBlock(1:J-leveldiff), neighbor_treecode(1:J-leveldiff), J-leveldiff) ) then
!                 level_down_neighbor_possible = .false.
!             else
!                 level_down_neighbor_possible = .true.
!             endif
!         endif

!     elseif (leveldiff==-1) then
!         ! first neighbor virtual treecode, one level up
!         virt_treecode = tcBlock
!         ! append last digit
!         virt_treecode( level+1 ) = virtual_treecode_finer_level(direction)

!         if (virt_treecode( level+1 ) < 0 .or. virt_treecode( level+1 ) > 7) then
!             call abort(587361912, "Senator Palpatine, for some reason we try to append garbage to treecode!")
!         endif

!         ! calculate treecode for neighbor on same level (virtual level)
!         if (dim==3) then
!             call adjacent_block_3D( virt_treecode, neighbor_treecode, direction_names_3d(direction), &
!             level+1, max_treelevel)
!         else
!             call adjacent_block_2D( virt_treecode, neighbor_treecode, direction_names_2d(direction), &
!             level+1, max_treelevel)
!         endif

!     endif
!     ! required: list of neighbors for all leveldiffs -1, 0, +1
!     !
!     ! same level. just a call to adjacent2d/3d very simple
!     ! coarser level. call adjacent and remove last treecode digit -> easy
!     ! finer level. call adjacent2d/3d for a virtual treecode, which is XXXXXY where the new appended
!     ! Y depends on direction
! end subroutine get_neighbor_treecode
