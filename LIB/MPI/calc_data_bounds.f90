! This code does something very simple in a very complicated way.
! In fact, if you had two blocks treecodes for all neighboring relations, everything would be nice
! and smooth.
! Why you ask? because from that you can direcly compute their origin and spacing in a straightforward way,
! and you hence have all coordinates of all point on both blocks.
! Now you can be clever: synchronizing ghosts means copying points with THE SAME COORDINATES. Makes sense, doesnt't it?
! So all you would have to do is compute the sender bounds (from the given, manually set recver bounds)

subroutine compute_sender_buffer_bounds(params, ijkrecv, ijksend, ijkbuffer, dir, leveldiff, gminus, gplus, output_to_file )
    implicit none
    type (type_params), intent(in) :: params
    integer(kind=ik), intent(in) :: ijkrecv(2,3)
    integer(kind=ik), intent(out) :: ijkbuffer(2,3)
    integer(kind=ik), intent(out) :: ijksend(2,3)
    ! leveldiff = 1 ! -1: interpolation, +1: coarsening
    integer(kind=ik), intent(in) :: dir, leveldiff, gminus, gplus
    logical, intent(in) :: output_to_file

    integer(kind=ik), parameter :: Jmax = 6
    integer(kind=ik) :: send_treecode(1:Jmax)
    integer(kind=ik) :: recv_treecode(1:Jmax)
    integer(kind=ik) :: J=4, ineighbor, ishift, ileveldiff
    integer(kind=ik) :: i1, i2, nDirs, min_size, Nsender, i
    integer(kind=ik), dimension(3) :: Bs
    real(kind=rk) :: x0_send(1:3), dx_send(1:3), x0_recv(1:3), dx_recv(1:3), x0_buffer(1:3), dx_buffer(1:3)
    real(kind=rk) :: r1, r2, q1, q2
    ! List of sender and recver pairs, treecodes. They depend on the direction (1:16 in 2D and 1:74 in 3D)
    ! and the level difference, space dimension (2d/3d). They are helper qtys: not related to the treecodes
    ! of the actual computational grid.
    ! These treecodes are easily computed when the sender is on the same level or finer than the recver, by
    ! using the adjacent block compuation. If the neighbor (recver) is coarser, then this is not trivial, as
    ! our grid definition does not allow all blocks to have coarser neighbors. In some cases, the neighbor has to be
    ! on the same level. Consider five blocks:
    !
    ! a b E E
    ! c d E E
    !
    ! then block "c" cannot have a coarser neighbor in the direction to the right! only b,d can.
    !
    ! We computed this list simply by choosind another recver block until we found a valid combination.
    ! For ease of coding, it is now (>03/2020) a hard coded list, see git history for the code.
    integer(kind=ik), dimension(1:74, -1:1, 2:3, 1:Jmax) :: senders, recvers

! write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+"
! write(*,*) "computing sender bounds for dir=", dir, "leveldiff=", leveldiff
! write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+"

    senders = -1
    recvers = -1

    ! the if clause is for code folding only.
    if (.true.) then
        senders( 1, 0, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 1, 0, 2, :) = (/ 0, 3, 3, 1,-1,-1/)
        senders( 2, 0, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 2, 0, 2, :) = (/ 1, 2, 2, 2,-1,-1/)
        senders( 3, 0, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 3, 0, 2, :) = (/ 2, 1, 1, 1,-1,-1/)
        senders( 4, 0, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 4, 0, 2, :) = (/ 0, 3, 3, 2,-1,-1/)
        senders( 5,-1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 5,-1, 2, :) = (/ 1, 2, 2, 0, 2,-1/)
        senders( 5, 0, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 5, 0, 2, :) = (/ 1, 2, 2, 0,-1,-1/)
        senders( 5, 1, 2, :) = (/ 2, 1, 1, 1,-1,-1/)
        recvers( 5, 1, 2, :) = (/ 1, 2, 2,-1,-1,-1/)
        senders( 6,-1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 6,-1, 2, :) = (/ 0, 3, 3, 0, 3,-1/)
        senders( 6, 0, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 6, 0, 2, :) = (/ 0, 3, 3, 0,-1,-1/)
        senders( 6, 1, 2, :) = (/ 3, 0, 0, 0,-1,-1/)
        recvers( 6, 1, 2, :) = (/ 0, 3, 3,-1,-1,-1/)
        senders( 7,-1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 7,-1, 2, :) = (/ 3, 0, 0, 0, 0,-1/)
        senders( 7, 0, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 7, 0, 2, :) = (/ 3, 0, 0, 0,-1,-1/)
        senders( 7, 1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 7, 1, 2, :) = (/ 3, 0, 0,-1,-1,-1/)
        senders( 8,-1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 8,-1, 2, :) = (/ 2, 1, 1, 0, 1,-1/)
        senders( 8, 0, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 8, 0, 2, :) = (/ 2, 1, 1, 0,-1,-1/)
        senders( 8, 1, 2, :) = (/ 1, 2, 2, 2,-1,-1/)
        recvers( 8, 1, 2, :) = (/ 2, 1, 1,-1,-1,-1/)
        senders( 9,-1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers( 9,-1, 2, :) = (/ 0, 3, 3, 1, 3,-1/)
        senders( 9, 1, 2, :) = (/ 2, 1, 1, 1,-1,-1/)
        recvers( 9, 1, 2, :) = (/ 0, 3, 3,-1,-1,-1/)
        senders(10,-1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers(10,-1, 2, :) = (/ 0, 3, 3, 1, 2,-1/)
        senders(10, 1, 2, :) = (/ 3, 0, 0, 0,-1,-1/)
        recvers(10, 1, 2, :) = (/ 1, 2, 2,-1,-1,-1/)
        senders(11,-1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers(11,-1, 2, :) = (/ 2, 1, 1, 1, 1,-1/)
        senders(11, 1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers(11, 1, 2, :) = (/ 2, 1, 1,-1,-1,-1/)
        senders(12,-1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers(12,-1, 2, :) = (/ 2, 1, 1, 1, 0,-1/)
        senders(12, 1, 2, :) = (/ 1, 2, 2, 2,-1,-1/)
        recvers(12, 1, 2, :) = (/ 3, 0, 0,-1,-1,-1/)
        senders(13,-1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers(13,-1, 2, :) = (/ 1, 2, 2, 2, 0,-1/)
        senders(13, 1, 2, :) = (/ 2, 1, 1, 1,-1,-1/)
        recvers(13, 1, 2, :) = (/ 3, 0, 0,-1,-1,-1/)
        senders(14,-1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers(14,-1, 2, :) = (/ 1, 2, 2, 2, 2,-1/)
        senders(14, 1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers(14, 1, 2, :) = (/ 1, 2, 2,-1,-1,-1/)
        senders(15,-1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers(15,-1, 2, :) = (/ 0, 3, 3, 2, 1,-1/)
        senders(15, 1, 2, :) = (/ 3, 0, 0, 0,-1,-1/)
        recvers(15, 1, 2, :) = (/ 2, 1, 1,-1,-1,-1/)
        senders(16,-1, 2, :) = (/ 0, 3, 3, 3,-1,-1/)
        recvers(16,-1, 2, :) = (/ 0, 3, 3, 2, 3,-1/)
        senders(16, 1, 2, :) = (/ 1, 2, 2, 2,-1,-1/)
        recvers(16, 1, 2, :) = (/ 0, 3, 3,-1,-1,-1/)

        senders( 1, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers( 1, 0, 3, :) = (/ 4, 3, 3, 3,-1,-1/)
        senders( 2, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers( 2, 0, 3, :) = (/ 0, 7, 7, 6,-1,-1/)
        senders( 3, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers( 3, 0, 3, :) = (/ 2, 5, 5, 5,-1,-1/)
        senders( 4, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers( 4, 0, 3, :) = (/ 1, 6, 6, 6,-1,-1/)
        senders( 5, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers( 5, 0, 3, :) = (/ 0, 7, 7, 5,-1,-1/)
        senders( 6, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers( 6, 0, 3, :) = (/ 0, 7, 7, 3,-1,-1/)
        senders( 7, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers( 7, 0, 3, :) = (/ 4, 3, 3, 2,-1,-1/)
        senders( 8, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers( 8, 0, 3, :) = (/ 6, 1, 1, 1,-1,-1/)
        senders( 9, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers( 9, 0, 3, :) = (/ 5, 2, 2, 2,-1,-1/)
        senders(10, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(10, 0, 3, :) = (/ 4, 3, 3, 1,-1,-1/)
        senders(11, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(11, 0, 3, :) = (/ 0, 7, 7, 2,-1,-1/)
        senders(12, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(12, 0, 3, :) = (/ 2, 5, 5, 1,-1,-1/)
        senders(13, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(13, 0, 3, :) = (/ 1, 6, 6, 2,-1,-1/)
        senders(14, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(14, 0, 3, :) = (/ 0, 7, 7, 1,-1,-1/)
        senders(15, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(15, 0, 3, :) = (/ 2, 5, 5, 4,-1,-1/)
        senders(16, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(16, 0, 3, :) = (/ 0, 7, 7, 4,-1,-1/)
        senders(17, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(17, 0, 3, :) = (/ 3, 4, 4, 4,-1,-1/)
        senders(18, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(18, 0, 3, :) = (/ 1, 6, 6, 4,-1,-1/)
        senders(19,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(19,-1, 3, :) = (/ 6, 1, 1, 0, 1,-1/)
        senders(19, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(19, 0, 3, :) = (/ 6, 1, 1, 0,-1,-1/)
        senders(19, 1, 3, :) = (/ 1, 6, 6, 6,-1,-1/)
        recvers(19, 1, 3, :) = (/ 6, 1, 1,-1,-1,-1/)
        senders(20,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(20,-1, 3, :) = (/ 7, 0, 0, 0, 0,-1/)
        senders(20, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(20, 0, 3, :) = (/ 7, 0, 0, 0,-1,-1/)
        senders(20, 1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(20, 1, 3, :) = (/ 7, 0, 0,-1,-1,-1/)
        senders(21,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(21,-1, 3, :) = (/ 5, 2, 2, 0, 2,-1/)
        senders(21, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(21, 0, 3, :) = (/ 5, 2, 2, 0,-1,-1/)
        senders(21, 1, 3, :) = (/ 2, 5, 5, 5,-1,-1/)
        recvers(21, 1, 3, :) = (/ 5, 2, 2,-1,-1,-1/)
        senders(22,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(22,-1, 3, :) = (/ 4, 3, 3, 0, 3,-1/)
        senders(22, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(22, 0, 3, :) = (/ 4, 3, 3, 0,-1,-1/)
        senders(22, 1, 3, :) = (/ 3, 4, 4, 4,-1,-1/)
        recvers(22, 1, 3, :) = (/ 4, 3, 3,-1,-1,-1/)
        senders(23,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(23,-1, 3, :) = (/ 2, 5, 5, 0, 5,-1/)
        senders(23, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(23, 0, 3, :) = (/ 2, 5, 5, 0,-1,-1/)
        senders(23, 1, 3, :) = (/ 5, 2, 2, 2,-1,-1/)
        recvers(23, 1, 3, :) = (/ 2, 5, 5,-1,-1,-1/)
        senders(24,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(24,-1, 3, :) = (/ 3, 4, 4, 0, 4,-1/)
        senders(24, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(24, 0, 3, :) = (/ 3, 4, 4, 0,-1,-1/)
        senders(24, 1, 3, :) = (/ 4, 3, 3, 3,-1,-1/)
        recvers(24, 1, 3, :) = (/ 3, 4, 4,-1,-1,-1/)
        senders(25,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(25,-1, 3, :) = (/ 1, 6, 6, 0, 6,-1/)
        senders(25, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(25, 0, 3, :) = (/ 1, 6, 6, 0,-1,-1/)
        senders(25, 1, 3, :) = (/ 6, 1, 1, 1,-1,-1/)
        recvers(25, 1, 3, :) = (/ 1, 6, 6,-1,-1,-1/)
        senders(26,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(26,-1, 3, :) = (/ 0, 7, 7, 0, 7,-1/)
        senders(26, 0, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(26, 0, 3, :) = (/ 0, 7, 7, 0,-1,-1/)
        senders(26, 1, 3, :) = (/ 7, 0, 0, 0,-1,-1/)
        recvers(26, 1, 3, :) = (/ 0, 7, 7,-1,-1,-1/)
        senders(27,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(27,-1, 3, :) = (/ 4, 3, 3, 3, 2,-1/)
        senders(27, 1, 3, :) = (/ 1, 6, 6, 6,-1,-1/)
        recvers(27, 1, 3, :) = (/ 5, 2, 2,-1,-1,-1/)
        senders(28,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(28,-1, 3, :) = (/ 4, 3, 3, 3, 3,-1/)
        senders(28, 1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(28, 1, 3, :) = (/ 4, 3, 3,-1,-1,-1/)
        senders(29,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(29,-1, 3, :) = (/ 4, 3, 3, 3, 1,-1/)
        senders(29, 1, 3, :) = (/ 2, 5, 5, 5,-1,-1/)
        recvers(29, 1, 3, :) = (/ 6, 1, 1,-1,-1,-1/)
        senders(30,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(30,-1, 3, :) = (/ 4, 3, 3, 3, 0,-1/)
        senders(30, 1, 3, :) = (/ 3, 4, 4, 4,-1,-1/)
        recvers(30, 1, 3, :) = (/ 7, 0, 0,-1,-1,-1/)
        senders(31,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(31,-1, 3, :) = (/ 0, 7, 7, 6, 7,-1/)
        senders(31, 1, 3, :) = (/ 1, 6, 6, 6,-1,-1/)
        recvers(31, 1, 3, :) = (/ 0, 7, 7,-1,-1,-1/)
        senders(32,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(32,-1, 3, :) = (/ 0, 7, 7, 6, 3,-1/)
        senders(32, 1, 3, :) = (/ 5, 2, 2, 2,-1,-1/)
        recvers(32, 1, 3, :) = (/ 4, 3, 3,-1,-1,-1/)
        senders(33,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(33,-1, 3, :) = (/ 0, 7, 7, 6, 5,-1/)
        senders(33, 1, 3, :) = (/ 3, 4, 4, 4,-1,-1/)
        recvers(33, 1, 3, :) = (/ 2, 5, 5,-1,-1,-1/)
        senders(34,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(34,-1, 3, :) = (/ 0, 7, 7, 6, 1,-1/)
        senders(34, 1, 3, :) = (/ 7, 0, 0, 0,-1,-1/)
        recvers(34, 1, 3, :) = (/ 6, 1, 1,-1,-1,-1/)
        senders(35,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(35,-1, 3, :) = (/ 2, 5, 5, 5, 4,-1/)
        senders(35, 1, 3, :) = (/ 1, 6, 6, 6,-1,-1/)
        recvers(35, 1, 3, :) = (/ 3, 4, 4,-1,-1,-1/)
        senders(36,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(36,-1, 3, :) = (/ 2, 5, 5, 5, 0,-1/)
        senders(36, 1, 3, :) = (/ 5, 2, 2, 2,-1,-1/)
        recvers(36, 1, 3, :) = (/ 7, 0, 0,-1,-1,-1/)
        senders(37,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(37,-1, 3, :) = (/ 2, 5, 5, 5, 5,-1/)
        senders(37, 1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(37, 1, 3, :) = (/ 2, 5, 5,-1,-1,-1/)
        senders(38,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(38,-1, 3, :) = (/ 2, 5, 5, 5, 1,-1/)
        senders(38, 1, 3, :) = (/ 4, 3, 3, 3,-1,-1/)
        recvers(38, 1, 3, :) = (/ 6, 1, 1,-1,-1,-1/)
        senders(39,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(39,-1, 3, :) = (/ 1, 6, 6, 6, 6,-1/)
        senders(39, 1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(39, 1, 3, :) = (/ 1, 6, 6,-1,-1,-1/)
        senders(40,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(40,-1, 3, :) = (/ 1, 6, 6, 6, 2,-1/)
        senders(40, 1, 3, :) = (/ 4, 3, 3, 3,-1,-1/)
        recvers(40, 1, 3, :) = (/ 5, 2, 2,-1,-1,-1/)
        senders(41,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(41,-1, 3, :) = (/ 1, 6, 6, 6, 4,-1/)
        senders(41, 1, 3, :) = (/ 2, 5, 5, 5,-1,-1/)
        recvers(41, 1, 3, :) = (/ 3, 4, 4,-1,-1,-1/)
        senders(42,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(42,-1, 3, :) = (/ 1, 6, 6, 6, 0,-1/)
        senders(42, 1, 3, :) = (/ 6, 1, 1, 1,-1,-1/)
        recvers(42, 1, 3, :) = (/ 7, 0, 0,-1,-1,-1/)
        senders(43,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(43,-1, 3, :) = (/ 0, 7, 7, 5, 7,-1/)
        senders(43, 1, 3, :) = (/ 2, 5, 5, 5,-1,-1/)
        recvers(43, 1, 3, :) = (/ 0, 7, 7,-1,-1,-1/)
        senders(44,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(44,-1, 3, :) = (/ 0, 7, 7, 5, 3,-1/)
        senders(44, 1, 3, :) = (/ 6, 1, 1, 1,-1,-1/)
        recvers(44, 1, 3, :) = (/ 4, 3, 3,-1,-1,-1/)
        senders(45,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(45,-1, 3, :) = (/ 0, 7, 7, 5, 6,-1/)
        senders(45, 1, 3, :) = (/ 3, 4, 4, 4,-1,-1/)
        recvers(45, 1, 3, :) = (/ 1, 6, 6,-1,-1,-1/)
        senders(46,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(46,-1, 3, :) = (/ 0, 7, 7, 5, 2,-1/)
        senders(46, 1, 3, :) = (/ 7, 0, 0, 0,-1,-1/)
        recvers(46, 1, 3, :) = (/ 5, 2, 2,-1,-1,-1/)
        senders(47,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(47,-1, 3, :) = (/ 0, 7, 7, 3, 6,-1/)
        senders(47, 1, 3, :) = (/ 5, 2, 2, 2,-1,-1/)
        recvers(47, 1, 3, :) = (/ 1, 6, 6,-1,-1,-1/)
        senders(48,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(48,-1, 3, :) = (/ 0, 7, 7, 3, 7,-1/)
        senders(48, 1, 3, :) = (/ 4, 3, 3, 3,-1,-1/)
        recvers(48, 1, 3, :) = (/ 0, 7, 7,-1,-1,-1/)
        senders(49,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(49,-1, 3, :) = (/ 0, 7, 7, 3, 5,-1/)
        senders(49, 1, 3, :) = (/ 6, 1, 1, 1,-1,-1/)
        recvers(49, 1, 3, :) = (/ 2, 5, 5,-1,-1,-1/)
        senders(50,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(50,-1, 3, :) = (/ 0, 7, 7, 3, 4,-1/)
        senders(50, 1, 3, :) = (/ 7, 0, 0, 0,-1,-1/)
        recvers(50, 1, 3, :) = (/ 3, 4, 4,-1,-1,-1/)
        senders(51,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(51,-1, 3, :) = (/ 4, 3, 3, 2, 3,-1/)
        senders(51, 1, 3, :) = (/ 1, 6, 6, 6,-1,-1/)
        recvers(51, 1, 3, :) = (/ 4, 3, 3,-1,-1,-1/)
        senders(52,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(52,-1, 3, :) = (/ 4, 3, 3, 2, 1,-1/)
        senders(52, 1, 3, :) = (/ 3, 4, 4, 4,-1,-1/)
        recvers(52, 1, 3, :) = (/ 6, 1, 1,-1,-1,-1/)
        senders(53,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(53,-1, 3, :) = (/ 6, 1, 1, 1, 0,-1/)
        senders(53, 1, 3, :) = (/ 1, 6, 6, 6,-1,-1/)
        recvers(53, 1, 3, :) = (/ 7, 0, 0,-1,-1,-1/)
        senders(54,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(54,-1, 3, :) = (/ 6, 1, 1, 1, 1,-1/)
        senders(54, 1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(54, 1, 3, :) = (/ 6, 1, 1,-1,-1,-1/)
        senders(55,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(55,-1, 3, :) = (/ 5, 2, 2, 2, 2,-1/)
        senders(55, 1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(55, 1, 3, :) = (/ 5, 2, 2,-1,-1,-1/)
        senders(56,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(56,-1, 3, :) = (/ 5, 2, 2, 2, 0,-1/)
        senders(56, 1, 3, :) = (/ 2, 5, 5, 5,-1,-1/)
        recvers(56, 1, 3, :) = (/ 7, 0, 0,-1,-1,-1/)
        senders(57,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(57,-1, 3, :) = (/ 4, 3, 3, 1, 3,-1/)
        senders(57, 1, 3, :) = (/ 2, 5, 5, 5,-1,-1/)
        recvers(57, 1, 3, :) = (/ 4, 3, 3,-1,-1,-1/)
        senders(58,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(58,-1, 3, :) = (/ 4, 3, 3, 1, 2,-1/)
        senders(58, 1, 3, :) = (/ 3, 4, 4, 4,-1,-1/)
        recvers(58, 1, 3, :) = (/ 5, 2, 2,-1,-1,-1/)
        senders(59,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(59,-1, 3, :) = (/ 0, 7, 7, 2, 7,-1/)
        senders(59, 1, 3, :) = (/ 5, 2, 2, 2,-1,-1/)
        recvers(59, 1, 3, :) = (/ 0, 7, 7,-1,-1,-1/)
        senders(60,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(60,-1, 3, :) = (/ 0, 7, 7, 2, 5,-1/)
        senders(60, 1, 3, :) = (/ 7, 0, 0, 0,-1,-1/)
        recvers(60, 1, 3, :) = (/ 2, 5, 5,-1,-1,-1/)
        senders(61,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(61,-1, 3, :) = (/ 2, 5, 5, 1, 4,-1/)
        senders(61, 1, 3, :) = (/ 5, 2, 2, 2,-1,-1/)
        recvers(61, 1, 3, :) = (/ 3, 4, 4,-1,-1,-1/)
        senders(62,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(62,-1, 3, :) = (/ 2, 5, 5, 1, 5,-1/)
        senders(62, 1, 3, :) = (/ 4, 3, 3, 3,-1,-1/)
        recvers(62, 1, 3, :) = (/ 2, 5, 5,-1,-1,-1/)
        senders(63,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(63,-1, 3, :) = (/ 1, 6, 6, 2, 6,-1/)
        senders(63, 1, 3, :) = (/ 4, 3, 3, 3,-1,-1/)
        recvers(63, 1, 3, :) = (/ 1, 6, 6,-1,-1,-1/)
        senders(64,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(64,-1, 3, :) = (/ 1, 6, 6, 2, 4,-1/)
        senders(64, 1, 3, :) = (/ 6, 1, 1, 1,-1,-1/)
        recvers(64, 1, 3, :) = (/ 3, 4, 4,-1,-1,-1/)
        senders(65,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(65,-1, 3, :) = (/ 0, 7, 7, 1, 7,-1/)
        senders(65, 1, 3, :) = (/ 6, 1, 1, 1,-1,-1/)
        recvers(65, 1, 3, :) = (/ 0, 7, 7,-1,-1,-1/)
        senders(66,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(66,-1, 3, :) = (/ 0, 7, 7, 1, 6,-1/)
        senders(66, 1, 3, :) = (/ 7, 0, 0, 0,-1,-1/)
        recvers(66, 1, 3, :) = (/ 1, 6, 6,-1,-1,-1/)
        senders(67,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(67,-1, 3, :) = (/ 2, 5, 5, 4, 5,-1/)
        senders(67, 1, 3, :) = (/ 1, 6, 6, 6,-1,-1/)
        recvers(67, 1, 3, :) = (/ 2, 5, 5,-1,-1,-1/)
        senders(68,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(68,-1, 3, :) = (/ 2, 5, 5, 4, 1,-1/)
        senders(68, 1, 3, :) = (/ 5, 2, 2, 2,-1,-1/)
        recvers(68, 1, 3, :) = (/ 6, 1, 1,-1,-1,-1/)
        senders(69,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(69,-1, 3, :) = (/ 0, 7, 7, 4, 7,-1/)
        senders(69, 1, 3, :) = (/ 3, 4, 4, 4,-1,-1/)
        recvers(69, 1, 3, :) = (/ 0, 7, 7,-1,-1,-1/)
        senders(70,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(70,-1, 3, :) = (/ 0, 7, 7, 4, 3,-1/)
        senders(70, 1, 3, :) = (/ 7, 0, 0, 0,-1,-1/)
        recvers(70, 1, 3, :) = (/ 4, 3, 3,-1,-1,-1/)
        senders(71,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(71,-1, 3, :) = (/ 3, 4, 4, 4, 4,-1/)
        senders(71, 1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(71, 1, 3, :) = (/ 3, 4, 4,-1,-1,-1/)
        senders(72,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(72,-1, 3, :) = (/ 3, 4, 4, 4, 0,-1/)
        senders(72, 1, 3, :) = (/ 4, 3, 3, 3,-1,-1/)
        recvers(72, 1, 3, :) = (/ 7, 0, 0,-1,-1,-1/)
        senders(73,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(73,-1, 3, :) = (/ 1, 6, 6, 4, 6,-1/)
        senders(73, 1, 3, :) = (/ 2, 5, 5, 5,-1,-1/)
        recvers(73, 1, 3, :) = (/ 1, 6, 6,-1,-1,-1/)
        senders(74,-1, 3, :) = (/ 0, 7, 7, 7,-1,-1/)
        recvers(74,-1, 3, :) = (/ 1, 6, 6, 4, 2,-1/)
        senders(74, 1, 3, :) = (/ 6, 1, 1, 1,-1,-1/)
        recvers(74, 1, 3, :) = (/ 5, 2, 2,-1,-1,-1/)
    endif

    send_treecode = -1
    recv_treecode = -1
    ijksend       = 1
    ijkbuffer     = 1
    Bs            = params%Bs


    ! this file contains a list of all neighboring relations
    ! used for development in 2D only.
#ifdef DEV
    if (output_to_file) then
        if ((dir == 1) .and. (params%rank == 0)) then
            open(16, file='neighbor_blocks2D.dat', status='replace')
            do ineighbor = 1, 16
                do ileveldiff = -1, 1
                    send_treecode = senders(ineighbor, ileveldiff, 2, :)
                    recv_treecode = recvers(ineighbor, ileveldiff, 2, :)

                    call get_block_spacing_origin_array( send_treecode(1:J), real(Bs*2**J, kind=rk), Bs, params%dim, x0_send, dx_send )
                    call get_block_spacing_origin_array( recv_treecode(1:J-ileveldiff), real(Bs*2**J, kind=rk), Bs, params%dim, x0_recv, dx_recv )

                    write(16,*) ineighbor, ileveldiff, x0_send, dx_send, x0_recv, dx_recv
                enddo
            enddo
            close(16)
        endif
    endif
#endif

    !***************************************************************************
    ! Compute sender bounds from recver bounds
    !***************************************************************************
    send_treecode = senders(dir, leveldiff, params%dim, :)
    recv_treecode = recvers(dir, leveldiff, params%dim, :)

    call get_block_spacing_origin_array( send_treecode(1:J), real(Bs*2**(J), kind=rk), Bs, params%dim, x0_send, dx_send )
    call get_block_spacing_origin_array( recv_treecode(1:J-leveldiff), real(Bs*2**(J), kind=rk), Bs, params%dim, x0_recv, dx_recv )

    do i = 1, params%dim
        ! shift to zero at the origin (which is g+1, actually)
        r1 = real(ijkrecv(1,i) - (gminus+1), kind=rk)
        r2 = real(ijkrecv(2,i) - (gminus+1), kind=rk)

        ! there's a very simple relation between sender and recver boundarys.
        ! we only need to define the recver bounds

        ! there is only two options here: either int or XX.5 floats.
        ! q1 = round_one_digit( (r1*dx_recv(i) + x0_recv(i) - x0_send(i)) / dx_send(i) )
        ! q2 = round_one_digit( (r2*dx_recv(i) + x0_recv(i) - x0_send(i)) / dx_send(i) )
        q1 = ( (r1*dx_recv(i) + x0_recv(i) - x0_send(i)) / dx_send(i) )
        q2 = ( (r2*dx_recv(i) + x0_recv(i) - x0_send(i)) / dx_send(i) )

        i1 = floor(q1) + (gminus+1)
        i2 = ceiling(q2) + (gminus+1)
        ijksend(1:2, i) = (/i1, i2/)

        ! write(*,*) "x0_send", x0_send, "dx_send", dx_send
        ! write(*,*) "x0_recv", x0_recv, "dx_recv", dx_recv
        ! ! NOTE at the moment, we really just computed the upper/lower bounds of the patch
        ! ! (selected on the receiver) on the sender. In the next step, we will extend the sender
        ! ! patch to *contain* what the receiver wants, but also a little more, if required, for interpolation.
        ! write(*,*) "dim=", i, "bounds=", i1, i2, "Dir=", dir, g+1, Bs(i)+g, "~", q1, q2,&
        !  (r1*dx_recv(i) + x0_recv(i) - x0_send(i)) / dx_send(i), (r2*dx_recv(i) + x0_recv(i) - x0_send(i)) / dx_send(i)
    enddo


    !***********************************************************************
    ! extend sender bounds as required by interpolation
    !***********************************************************************
    ! Sender patch extension. required if one-sided interp is to be avoided or the sender patch
    ! simply is too small for one-sided interpolation

    if (leveldiff == -1) then ! senseless for leveldiff /= -1
        ! the number S is how many extra coarse points on the sender side you use to
        ! avoid one-sided interpolation. The actual formula is S = (order-2 )/2 so for
        ! 6th order you would require 2 extra ones.
        ! NOTE: you can set it to S=0 and then the code uses one-sided interpolation stencils.
        ! NOTE: S is symmetric, we add the layer on all sides.
        if (params%order_predictor == "multiresolution_4th" ) then
            S  = 1
            min_size = 4
        elseif (params%order_predictor == "multiresolution_2nd" ) then
            S  = 0
            min_size = 2
        elseif (params%order_predictor == "multiresolution_6th" ) then
            S  = 2
            min_size = 6
        else
            call abort(2875490, "The predictor method is unknown")
        endif

        ! asymmetric shift. If the interpolation domain becomes too small for the stencil,
        ! then we can also extend it to the interior of the block, while still
        ! including only the redundant lines. The shift A (asymmetric shift) does that.
        ! COARSE grid points (sender side, red area in python figs).
        A = 0

        ! first we apply symmetric extension, i.e. make the sender patch larger
        ! by S on all sides. This allows excluding one-sided interpolation stencils,
        ! but requires on the other side to have the ghost node layer on the interpolating
        ! block already filled (two stages!)
        do i = 1, params%dim
            ijksend(1, i) = max(1,         ijksend(1, i) - S)
            ijksend(2, i) = min(Bs(i)+gminus+gplus, ijksend(2, i) + S)
        enddo

        ! then we possibly use asymmetric extension, i.e. make the patch larger in the
        ! direction of the interior. We use this also if we do not have enough points
        ! for the interpolation stencil. This way, one can still set S=0 but ensure having
        ! enough points.
        do i = 1, params%dim
            ! patch at least A (defined above) but possibly the required number to
            ! reach min_size.
            A = max(A, min_size - (ijksend(2,i)-ijksend(1,i)+1))
            if ( ijksend(1, i) <= gminus+1) then
                ijksend(2, i) = ijksend(2, i) + A
            else
                ijksend(1, i) = ijksend(1, i) - A
            endif
        enddo
    endif

    !***********************************************************************
    ! compute buffer (RESPRE) bounds
    !***********************************************************************
    ! set buffer for INTERPOLATION
    if (leveldiff == -1 ) then
        ! origin of buffer is lower point of sender patch
        x0_buffer(:) = real(ijksend(1, :)-(gminus+1), kind=rk) * dx_send(:) + x0_send(:)
        ! buffer spacing is fine (recv)
        dx_buffer(:) = dx_recv(:)

        do i = 1, params%dim
            ! find recv coordinates in buffer coordinates
            r1 = real(ijkrecv(1,i)-(gminus+1), kind=rk)
            r2 = real(ijkrecv(2,i)-(gminus+1), kind=rk)

            ! there's a very simple relation between sender and buffer boundarys.
            ! Buffer: x = x0 + dx*(ix-1)
            ! Recver: x = x0 + dx*(ix-(g+1))
            ! we only need to define the recver bounds (NOTE: RESPRE buffer starts 1,1,1 not g+1,g+1,g+1)
            i1 = nint( (r1*dx_recv(i) + x0_recv(i) - x0_buffer(i)) / dx_buffer(i) ) +1
            i2 = nint( (r2*dx_recv(i) + x0_recv(i) - x0_buffer(i)) / dx_buffer(i) ) +1
            ijkbuffer(:, i) = (/ i1, i2 /)
        enddo

    endif

    ! Set buffer bounds for the COARSENING cases.
    if (leveldiff == +1)  then
        ijkbuffer(1, 1:3) = 1
        do i = 1, params%dim
            Nsender = ijksend(2, i) - ijksend(1, i) + 1
            ijkbuffer(2, i) = (Nsender+1)/2
        enddo
    endif
end subroutine compute_sender_buffer_bounds


logical function patch_crosses_periodic_BC(x0, dx, ijk, dim)
    implicit none
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)
    integer(kind=ik) , intent(in) :: ijk(2,3), dim
    real(kind=rk) :: x(1:8,1:3)

    ! corners of patch
    x(1,1:3) = x0 + (/ dx(1)*ijk(1,1), dx(2)*ijk(1,2), dx(3)*ijk(1,3) /)
    x(2,1:3) = x0 + (/ dx(1)*ijk(2,1), dx(2)*ijk(2,2), dx(3)*ijk(2,3) /)

    x(3,1:3) = x0 + (/ dx(1)*ijk(2,1), dx(2)*ijk(1,2), dx(3)*ijk(1,3) /)
    x(4,1:3) = x0 + (/ dx(1)*ijk(1,1), dx(2)*ijk(2,2), dx(3)*ijk(1,3) /)
    x(5,1:3) = x0 + (/ dx(1)*ijk(1,1), dx(2)*ijk(1,2), dx(3)*ijk(2,3) /)

    x(6,1:3) = x0 + (/ dx(1)*ijk(2,1), dx(2)*ijk(2,2), dx(3)*ijk(1,3) /)
    x(7,1:3) = x0 + (/ dx(1)*ijk(1,1), dx(2)*ijk(2,2), dx(3)*ijk(2,3) /)
    x(8,1:3) = x0 + (/ dx(1)*ijk(2,1), dx(2)*ijk(1,2), dx(3)*ijk(2,3) /)

    patch_crosses_periodic_BC = .false.
    ! note we normalize xl, yl, zl to 1.0 just in these routines here
    if ( minval(x(:,1:dim)) < 0.0_rk .or. maxval(x(:,1:dim)) > 1.0_rk) then
        patch_crosses_periodic_BC = .true.
    endif

end function


! set recv bounds for different neighborhood or patch relations
subroutine set_recv_bounds( params, data_bounds, relation, level_diff, gminus, gplus)
    implicit none

    type (type_params), intent(in)                  :: params
    !> data_bounds
    integer(kind=ik), intent(inout)                 :: data_bounds(2,3)
    !> neighborhood or family relation, id from dirs
    !! -8:-1 is mother/daughter relation
    !! 0 is full block relation, level
    !! 1:74 is neighborhood relation
    integer(kind=ik), intent(in)                    :: relation
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: level_diff, gminus, gplus

    integer(kind=ik) :: Bs(1:3), g

    Bs = params%Bs
    g = gminus

    ! set 1 and not -1 (or anything else), because 2D bounds ignore 3rd dimension
    ! and thus cycle from 1:1
    data_bounds(:,:) = 1

    if ( params%dim == 3 ) then
        !---3D------Neighborhood
        select case(relation)
            ! '__1/___'
        case(1)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = 1
            data_bounds(2,3) = g

            ! '__2/___'
        case(2)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = Bs(2)+g+1
            data_bounds(2,2) = Bs(2)+g+gplus
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

            ! '__3/___'
        case(3)
            if (level_diff /= 0) return
            data_bounds(1,1) = 1
            data_bounds(2,1) = g
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

            ! '__4/___'
        case(4)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = 1
            data_bounds(2,2) = g
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

            ! '__5/___'
        case(5)
            if (level_diff /= 0) return
            data_bounds(1,1) = Bs(1)+g+1
            data_bounds(2,1) = Bs(1)+g+gplus
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

            ! '__6/___'
        case(6)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = Bs(3)+g+1
            data_bounds(2,3) = Bs(3)+g+gplus

            ! '_12/___'
        case(7)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = Bs(2)+g+1
            data_bounds(2,2) = Bs(2)+g+gplus
            data_bounds(1,3) = 1
            data_bounds(2,3) = g

            ! '_13/___'
        case(8)
            if (level_diff /= 0) return
            data_bounds(1,1) = 1
            data_bounds(2,1) = g
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = 1
            data_bounds(2,3) = g

            ! '_14/___'
        case(9)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = 1
            data_bounds(2,2) = g
            data_bounds(1,3) = 1
            data_bounds(2,3) = g

            ! '_15/___'
        case(10)
            if (level_diff /= 0) return
            data_bounds(1,1) = Bs(1)+g+1
            data_bounds(2,1) = Bs(1)+g+gplus
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = 1
            data_bounds(2,3) = g

            ! '_62/___'
        case(11)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = Bs(2)+g+1
            data_bounds(2,2) = Bs(2)+g+gplus
            data_bounds(1,3) = Bs(3)+g+1
            data_bounds(2,3) = Bs(3)+g+gplus

            ! '_63/___'
        case(12)
            if (level_diff /= 0) return
            data_bounds(1,1) = 1
            data_bounds(2,1) = g
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = Bs(3)+g+1
            data_bounds(2,3) = Bs(3)+g+gplus

            ! '_64/___'
        case(13)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = 1
            data_bounds(2,2) = g
            data_bounds(1,3) = Bs(3)+g+1
            data_bounds(2,3) = Bs(3)+g+gplus

            ! '_65/___'
        case(14)
            if (level_diff /= 0) return
            data_bounds(1,1) = Bs(1)+g+1
            data_bounds(2,1) = Bs(1)+g+gplus
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = Bs(3)+g+1
            data_bounds(2,3) = Bs(3)+g+gplus

            ! '_23/___'
        case(15)
            if (level_diff /= 0) return
            data_bounds(1,1) = 1
            data_bounds(2,1) = g
            data_bounds(1,2) = Bs(2)+g+1
            data_bounds(2,2) = Bs(2)+g+gplus
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

            ! '_25/___'
        case(16)
            if (level_diff /= 0) return
            data_bounds(1,1) = Bs(1)+g+1
            data_bounds(2,1) = Bs(1)+g+gplus
            data_bounds(1,2) = Bs(2)+g+1
            data_bounds(2,2) = Bs(2)+g+gplus
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

            ! '_43/___'
        case(17)
            if (level_diff /= 0) return
            data_bounds(1,1) = 1
            data_bounds(2,1) = g
            data_bounds(1,2) = 1
            data_bounds(2,2) = g
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

            ! '_45/___'
        case(18)
            if (level_diff /= 0) return
            data_bounds(1,1) = Bs(1)+g+1
            data_bounds(2,1) = Bs(1)+g+gplus
            data_bounds(1,2) = 1
            data_bounds(2,2) = g
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

        case(19,20,21,22)
            data_bounds(1,3) = 1
            data_bounds(2,3) = g
            select case(relation)
            case(19) ! '123/___'
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus

            case(20) ! '134/___'
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,2) = 1
                data_bounds(2,2) = g

            case(21) ! '145/___'
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = 1
                data_bounds(2,2) = g

            case(22) ! '152/___'
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus

            end select

        case(23,24,25,26)
            data_bounds(1,3) = Bs(3)+g+1
            data_bounds(2,3) = Bs(3)+g+gplus
            select case(relation)
            case(23) ! '623/___'
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus

            case(24) ! '634/___'
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,2) = 1
                data_bounds(2,2) = g

            case(25) ! '645/___'
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = 1
                data_bounds(2,2) = g

            case(26) ! '652/___'
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus

            end select

        case(27,28,29,30)
            if ( level_diff == -1 ) then
                data_bounds(1,3) = 1
                data_bounds(2,3) = g
                select case(relation)
                case(27) ! '__1/123'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+g+gplus

                case(28) ! '__1/134'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(29) ! '__1/145'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+g+gplus
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(30) ! '__1/152'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+g+gplus
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+g+gplus
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,3) = 1
                data_bounds(2,3) = g
                select case(relation)
                case(27) ! '__1/123'
                    data_bounds(1,1) = g+(Bs(1))/2 + 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2))/2

                case(28) ! '__1/134'
                    data_bounds(1,1) = g+(Bs(1))/2 + 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = g+(Bs(2))/2 + 1
                    data_bounds(2,2) = Bs(2)+g

                case(29) ! '__1/145'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1))/2
                    data_bounds(1,2) = g+(Bs(2))/2 + 1
                    data_bounds(2,2) = Bs(2)+g

                case(30) ! '__1/152'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1))/2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2))/2
                end select

            end if

        case(31,32,33,34)
            if ( level_diff == -1 ) then
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus
                ! first, third dimension
                select case(relation)
                case(31) ! '__2/123'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(32) ! '__2/623'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+g+gplus

                case(33) ! '__2/152'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+g+gplus
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(34) ! '__2/652'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+g+gplus
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+g+gplus

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus
                ! first, third dimension
                select case(relation)
                case(31) ! '__2/123'
                    data_bounds(1,1) = g+(Bs(1))/2 + 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = g+(Bs(3))/2 + 1
                    data_bounds(2,3) = Bs(3)+g

                case(32) ! '__2/623'
                    data_bounds(1,1) = g+(Bs(1))/2 + 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3))/2

                case(33) ! '__2/152'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1))/2
                    data_bounds(1,3) = g+(Bs(3))/2 + 1
                    data_bounds(2,3) = Bs(3)+g

                case(34) ! '__2/652'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1))/2
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3))/2

                end select

            end if

        case(35,36,37,38)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                ! second, third dimension
                select case(relation)
                case(35) ! '__3/123'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+g+gplus
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(37) ! '__3/134'
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(38) ! '__3/634'
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+g+gplus

                case(36) ! '__3/623'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+g+gplus
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+g+gplus

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                ! second, third dimension
                select case(relation)
                case(35) ! '__3/123'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2))/2
                    data_bounds(1,3) = g+(Bs(3))/2 + 1
                    data_bounds(2,3) = Bs(3)+g

                case(37) ! '__3/134'
                    data_bounds(1,2) = g+(Bs(2))/2 + 1
                    data_bounds(2,2) = Bs(2)+g
                    data_bounds(1,3) = g+(Bs(3))/2 + 1
                    data_bounds(2,3) = Bs(3)+g

                case(38) ! '__3/634'
                    data_bounds(1,2) = g+(Bs(2))/2 + 1
                    data_bounds(2,2) = Bs(2)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3))/2

                case(36) ! '__3/623'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2))/2
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3))/2

                end select

            end if

        case(39,40,41,42)
            if ( level_diff == -1 ) then
                data_bounds(1,2) = 1
                data_bounds(2,2) = g
                ! first, third dimension
                select case(relation)
                case(40) ! '__4/634'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+g+gplus

                case(39) ! '__4/134'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(41) ! '__4/145'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+g+gplus
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(42) ! '__4/645'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+g+gplus
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+g+gplus

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,2) = 1
                data_bounds(2,2) = g
                ! first, third dimension
                select case(relation)
                case(40) ! '__4/634'
                    data_bounds(1,1) = g+(Bs(1))/2 + 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3))/2

                case(39) ! '__4/134'
                    data_bounds(1,1) = g+(Bs(1))/2 + 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = g+(Bs(3))/2 + 1
                    data_bounds(2,3) = Bs(3)+g

                case(41) ! '__4/145'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1))/2
                    data_bounds(1,3) = g+(Bs(3))/2 + 1
                    data_bounds(2,3) = Bs(3)+g

                case(42) ! '__4/645'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1))/2
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3))/2

                end select

            end if

        case(43,44,45,46)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                ! second, third dimension
                select case(relation)
                case(45) ! '__5/152'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+g+gplus
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(43) ! '__5/145'
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(44) ! '__5/645'
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+g+gplus

                case(46) ! '__5/652'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+g+gplus
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+g+gplus

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                ! second, third dimension
                select case(relation)
                case(45) ! '__5/152'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2))/2
                    data_bounds(1,3) = g+(Bs(3))/2 + 1
                    data_bounds(2,3) = Bs(3)+g

                case(43) ! '__5/145'
                    data_bounds(1,2) = g+(Bs(2))/2 + 1
                    data_bounds(2,2) = Bs(2)+g
                    data_bounds(1,3) = g+(Bs(3))/2 + 1
                    data_bounds(2,3) = Bs(3)+g

                case(44) ! '__5/645'
                    data_bounds(1,2) = g+(Bs(2))/2 + 1
                    data_bounds(2,2) = Bs(2)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3))/2

                case(46) ! '__5/652'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2))/2
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3))/2

                end select

            end if

        case(47,48,49,50)
            if ( level_diff == -1 ) then
                data_bounds(1,3) = Bs(3)+g+1
                data_bounds(2,3) = Bs(3)+g+gplus
                select case(relation)
                case(47) ! '__6/623'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+g+gplus

                case(48) ! '__6/634'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(49) ! '__6/645'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+g+gplus
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(50) ! '__6/652'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+g+gplus
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+g+gplus

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,3) = Bs(3)+g+1
                data_bounds(2,3) = Bs(3)+g+gplus
                select case(relation)
                case(47) ! '__6/623'
                    data_bounds(1,1) = g+(Bs(1))/2 + 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2))/2

                case(48) ! '__6/634'
                    data_bounds(1,1) = g+(Bs(1))/2 + 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = g+(Bs(2))/2 + 1
                    data_bounds(2,2) = Bs(2)+g

                case(49) ! '__6/645'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1))/2
                    data_bounds(1,2) = g+(Bs(2))/2 + 1
                    data_bounds(2,2) = Bs(2)+g

                case(50) ! '__6/652'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1))/2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2))/2

                end select

            end if

        case(51,52)
            if ( level_diff == -1 ) then
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus
                data_bounds(1,3) = 1
                data_bounds(2,3) = g
                select case(relation)
                case(51) ! '_12/123'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g

                case(52) ! '_12/152'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+g+gplus
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus
                data_bounds(1,3) = 1
                data_bounds(2,3) = g
                select case(relation)
                case(51) ! '_12/123'
                    data_bounds(1,1) = g+(Bs(1))/2 + 1
                    data_bounds(2,1) = Bs(1)+g

                case(52) ! '_12/152'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1))/2
                end select

            end if

        case(53,54)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,3) = 1
                data_bounds(2,3) = g
                select case(relation)
                case(54) ! '_13/134'
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(53) ! '_13/123'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+g+gplus
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,3) = 1
                data_bounds(2,3) = g
                select case(relation)
                case(54) ! '_13/134'
                    data_bounds(1,2) = g+(Bs(2))/2 + 1
                    data_bounds(2,2) = Bs(2)+g

                case(53) ! '_13/123'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2))/2
                end select

            end if

        case(55,56)
            if ( level_diff == -1 ) then
                data_bounds(1,2) = 1
                data_bounds(2,2) = g
                data_bounds(1,3) = 1
                data_bounds(2,3) = g
                select case(relation)
                case(55) ! '_14/134'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g

                case(56) ! '_14/145'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+g+gplus

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,2) = 1
                data_bounds(2,2) = g
                data_bounds(1,3) = 1
                data_bounds(2,3) = g
                select case(relation)
                case(55) ! '_14/134'
                    data_bounds(1,1) = g+(Bs(1))/2 + 1
                    data_bounds(2,1) = Bs(1)+g

                case(56) ! '_14/145'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1))/2

                end select

            end if

        case(57,58)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,3) = 1
                data_bounds(2,3) = g
                select case(relation)
                case(57) ! '_15/145'
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(58) ! '_15/152''
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+g+gplus

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,3) = 1
                data_bounds(2,3) = g
                select case(relation)
                case(57) ! '_15/145'
                    data_bounds(1,2) = g+(Bs(2))/2 + 1
                    data_bounds(2,2) = Bs(2)+g

                case(58) ! '_15/152''
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2))/2

                end select

            end if

        case(59,60)
            if ( level_diff == -1 ) then
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus
                data_bounds(1,3) = Bs(3)+g+1
                data_bounds(2,3) = Bs(3)+g+gplus
                select case(relation)
                case(59) ! '_62/623'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g

                case(60) ! '_62/652'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+g+gplus
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus
                data_bounds(1,3) = Bs(3)+g+1
                data_bounds(2,3) = Bs(3)+g+gplus
                select case(relation)
                case(59) ! '_62/623'
                    data_bounds(1,1) = g+(Bs(1))/2 + 1
                    data_bounds(2,1) = Bs(1)+g

                case(60) ! '_62/652'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1))/2
                end select

            end if

        case(61,62)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,3) = Bs(3)+g+1
                data_bounds(2,3) = Bs(3)+g+gplus
                select case(relation)
                case(62) ! '_63/634'
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(61) ! '_63/623'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+g+gplus
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,3) = Bs(3)+g+1
                data_bounds(2,3) = Bs(3)+g+gplus
                select case(relation)
                case(62) ! '_63/634'
                    data_bounds(1,2) = g+(Bs(2))/2 + 1
                    data_bounds(2,2) = Bs(2)+g

                case(61) ! '_63/623'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2))/2
                end select

            end if

        case(63,64)
            if ( level_diff == -1 ) then
                data_bounds(1,2) = 1
                data_bounds(2,2) = g
                data_bounds(1,3) = Bs(3)+g+1
                data_bounds(2,3) = Bs(3)+g+gplus
                select case(relation)
                case(63) ! '_64/634'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g

                case(64) ! '_64/645'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+g+gplus

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,2) = 1
                data_bounds(2,2) = g
                data_bounds(1,3) = Bs(3)+g+1
                data_bounds(2,3) = Bs(3)+g+gplus
                select case(relation)
                case(63) ! '_64/634'
                    data_bounds(1,1) = g+(Bs(1))/2 + 1
                    data_bounds(2,1) = Bs(1)+g

                case(64) ! '_64/645'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1))/2

                end select

            end if

        case(65,66)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,3) = Bs(3)+g+1
                data_bounds(2,3) = Bs(3)+g+gplus
                select case(relation)
                case(65) ! '_65/645'
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(66) ! '_65/652'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+g+gplus

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,3) = Bs(3)+g+1
                data_bounds(2,3) = Bs(3)+g+gplus
                select case(relation)
                case(65) ! '_65/645'
                    data_bounds(1,2) = g+(Bs(2))/2 + 1
                    data_bounds(2,2) = Bs(2)+g

                case(66) ! '_65/652'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2))/2

                end select

            end if

        case(67,68)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus
                select case(relation)
                case(67) ! '_23/123'
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(68) ! '_23/236''
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+g+gplus
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus
                select case(relation)
                case(67) ! '_23/123'
                    data_bounds(1,3) = g+(Bs(3))/2 + 1
                    data_bounds(2,3) = Bs(3)+g

                case(68) ! '_23/236''
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3))/2
                end select

            end if

        case(69,70)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus
                select case(relation)
                case(69) ! '_25/152'
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(70) ! '_25/652''
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+g+gplus
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus
                select case(relation)
                case(69) ! '_25/152'
                    data_bounds(1,3) = g+(Bs(3))/2 + 1
                    data_bounds(2,3) = Bs(3)+g

                case(70) ! '_25/652''
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3))/2
                end select

            end if

        case(71,72)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,2) = 1
                data_bounds(2,2) = g
                select case(relation)
                case(71) ! '_43/134'
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(72) ! '_43/634''
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+g+gplus
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,2) = 1
                data_bounds(2,2) = g
                select case(relation)
                case(71) ! '_43/134'
                    data_bounds(1,3) = g+(Bs(3))/2 + 1
                    data_bounds(2,3) = Bs(3)+g

                case(72) ! '_43/634''
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3))/2
                end select

            end if

        case(73,74)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = 1
                data_bounds(2,2) = g
                select case(relation)
                case(73) ! '_45/145'
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(74) ! '_45/645'
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+g+gplus
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = 1
                data_bounds(2,2) = g
                select case(relation)
                case(73) ! '_45/145'
                    data_bounds(1,3) = g+(Bs(3))/2 + 1
                    data_bounds(2,3) = Bs(3)+g

                case(74) ! '_45/645'
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3))/2
                end select

            end if

        !---3D------Full block relation
        case(0)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(2)+g
        !---3D------family relation, assume values on finer side are already WD in mallat-ordering
        case(-1, -2, -3, -4, -5, -6, -7, -8)
            if (level_diff == 0) then
                return
            ! only transfer SC if this is the finer block
            elseif (level_diff == -1) then
                data_bounds(1,1) = g+1
                data_bounds(2,1) = g+Bs(1)/2
                data_bounds(1,2) = g+1
                data_bounds(2,2) = g+Bs(2)/2
                data_bounds(1,3) = g+1
                data_bounds(2,3) = g+Bs(2)/2
            else
                if (modulo(-relation, 2) == 0) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+Bs(1)/2
                else
                    data_bounds(1,1) = g+1+Bs(1)/2
                    data_bounds(2,1) = g+Bs(1)
                endif
                if (modulo(-relation/2, 2) == 0) then
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+Bs(2)/2
                else
                    data_bounds(1,2) = g+1+Bs(2)/2
                    data_bounds(2,2) = g+Bs(2)
                endif
                if (modulo(-relation/4, 2) == 0) then
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+Bs(2)/2
                else
                    data_bounds(1,3) = g+1+Bs(3)/2
                    data_bounds(2,3) = g+Bs(3)
                endif
            endif

        end select
    else
        !---2D------Neighborhood
        select case(relation)
            ! '__N'
        case(1)
            if (level_diff /= 0) return
            data_bounds(1,1) = Bs(1)+g+1
            data_bounds(2,1) = Bs(1)+g+gplus
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g

            ! '__E'
        case(2)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = 1
            data_bounds(2,2) = g

            ! '__S'
        case(3)
            if (level_diff /= 0) return
            data_bounds(1,1) = 1
            data_bounds(2,1) = g
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g

            ! '__W'
        case(4)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = Bs(2)+g+1
            data_bounds(2,2) = Bs(2)+g+gplus

            ! '_NE'
        case(5)
            data_bounds(1,1) = Bs(1)+g+1
            data_bounds(2,1) = Bs(1)+g+gplus
            data_bounds(1,2) = 1
            data_bounds(2,2) = g

            ! '_NW'
        case(6)
            data_bounds(1,1) = Bs(1)+g+1
            data_bounds(2,1) = Bs(1)+g+gplus
            data_bounds(1,2) = Bs(2)+g+1
            data_bounds(2,2) = Bs(2)+g+gplus

            ! '_SE'
        case(7)
            data_bounds(1,1) = 1
            data_bounds(2,1) = g
            data_bounds(1,2) = 1
            data_bounds(2,2) = g

            ! '_SW'
        case(8)
            data_bounds(1,1) = 1
            data_bounds(2,1) = g
            data_bounds(1,2) = Bs(2)+g+1
            data_bounds(2,2) = Bs(2)+g+gplus

            ! 'NNE'
        case(9)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = 1
                data_bounds(2,2) = Bs(2)+g

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = g+(Bs(2))/2 +1
                data_bounds(2,2) = Bs(2)+g

            end if

            ! 'NNW'
        case(10)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs(2)+g+gplus

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = g+1
                data_bounds(2,2) = g+(Bs(2))/2

            end if

            ! 'SSE'
        case(11)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,2) = 1
                data_bounds(2,2) = Bs(2)+g

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,2) = g+(Bs(2))/2 + 1
                data_bounds(2,2) = Bs(2)+g

            end if

            ! 'SSW'
        case(12)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs(2)+g+gplus

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = g
                data_bounds(1,2) = g+1
                data_bounds(2,2) = g+(Bs(2))/2

            end if

            ! 'ENE'
        case(13)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = 1
                data_bounds(2,2) = g

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = g+1
                data_bounds(2,1) = g+(Bs(1))/2
                data_bounds(1,2) = 1
                data_bounds(2,2) = g

            end if

            ! 'ESE'
        case(14)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = Bs(1)+g
                data_bounds(1,2) = 1
                data_bounds(2,2) = g

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = g+(Bs(1))/2  + 1
                data_bounds(2,1) = Bs(1)+g
                data_bounds(1,2) = 1
                data_bounds(2,2) = g

            end if

            ! 'WNW'
        case(15)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs(1)+g+gplus
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = g+1
                data_bounds(2,1) = g+(Bs(1))/2
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus

            end if

            ! 'WSW'
        case(16)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = Bs(1)+g
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = g+(Bs(1))/2 + 1
                data_bounds(2,1) = Bs(1)+g
                data_bounds(1,2) = Bs(2)+g+1
                data_bounds(2,2) = Bs(2)+g+gplus

            end if

        !---2D------Full block relation
        case(0)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g

        !---2D------family relation, assume values on finer side are already WD in mallat-ordering
        case(-1, -2, -3, -4)
            if (level_diff == 0) then
                return
            ! only transfer SC if this is the finer block
            elseif (level_diff == -1) then
                data_bounds(1,1) = g+1
                data_bounds(2,1) = g+Bs(1)/2
                data_bounds(1,2) = g+1
                data_bounds(2,2) = g+Bs(2)/2
            else
                if (modulo(-relation, 2) == 0) then
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+Bs(1)/2
                else
                    data_bounds(1,1) = g+1+Bs(1)/2
                    data_bounds(2,1) = g+Bs(1)
                endif
                if (modulo(-relation/2, 2) == 0) then
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+Bs(2)/2
                else
                    data_bounds(1,2) = g+1+Bs(2)/2
                    data_bounds(2,2) = g+Bs(2)
                endif
            endif

        end select
    end if



end subroutine set_recv_bounds
