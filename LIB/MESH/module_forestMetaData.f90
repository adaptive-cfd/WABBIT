module module_forestMetaData

    use module_globals

    ! This module contains the memory for the forest meta data.
    !
    ! Until 08/2022, those arrays were passed as arguments to subroutines,
    ! which is generally a good idea.
    ! However, this strategy makes the subroutine calls very long, and if the descriptors
    ! change, they need to be changed everywhere in the code. For example, on 21 Aug 2022,
    ! TE commited a 4000 lines commit because we changed some array structures.
    !
    ! Therefore, we created this "global" module. It is not the nicest style to use
    ! such globals, and their usage should be kept to an absolute minimum. At this time,
    ! it seems they're the lesser of many evils (for example derived datatypes, which include
    ! a perfomance penalty).
    !
    ! This is especially true as we need more descriptors in the future, for example
    ! leaves and virtual nodes of a tree.
    implicit none

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    ! Note how neither HVY_BLOCK (usually the state vector)
    ! nor HVY_TMP (a work array)
    ! are included here !! 

    integer(kind=ik), allocatable, public    :: lgt_block(:, :)

    integer(kind=ik), allocatable, public    :: hvy_neighbor(:,:)
    ! list of active blocks (light data) for each tree
    integer(kind=ik), allocatable, public    :: lgt_active(:,:)
    ! number of active blocks (light data) for each tree
    integer(kind=ik), allocatable, public    :: lgt_n(:)

    ! list of active blocks (heavy data) for each tree
    integer(kind=ik), allocatable, public    :: hvy_active(:,:)
    ! number of active blocks (heavy data) for each tree
    integer(kind=ik), allocatable, public    :: hvy_n(:)

    ! if only a subset of blocks is synchronized, this flag is used to determine which
    logical, allocatable, public :: lgt_BlocksToSync(:)

    ! The following list contains the numerical treecode and the lightID for the active blocks
    ! in a sorted fashion. this is very important for finding blocks. usually, in the rest of the code,
    ! a treecode is an array and this is handy. for finding a block however, this is not true,
    ! here, having a single, unique number is a lot faster. these numbers are called numerical treecodes.
    integer(kind=tsize), allocatable, public :: lgt_sortednumlist(:,:,:)

    ! number of active trees
    integer(kind=ik), public    :: tree_n

end module
