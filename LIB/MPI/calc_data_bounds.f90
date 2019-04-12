! This code does something very simple in a very complicated way.
! In fact, if you had two blocks treecodes for all neighboring relations, everything would be nice
! and smooth.
! Why you ask? because from that you can direcly compute their origin and spacing in a straightforward way,
! and you hence have all coordinates of all point on both blocks.
! Now you can be clever: synchronizing ghosts means copying points with THE SAME COORDINATES. Makes sense, doesnt't it?
! So all you would have to do is compute the sender bounds (from the given, manually set recver bounds)

subroutine compute_sender_buffer_bounds(params, ijkrecv, ijksend, ijkbuffer, dir, leveldiff, TYPE)
    implicit none
    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in) :: ijkrecv(2,3)
    integer(kind=ik), intent(out) :: ijkbuffer(2,3)
    integer(kind=ik), intent(out) :: ijksend(2,3)
    ! leveldiff = 1 ! -1: interpolation, +1: coarsening
    integer(kind=ik), intent(in) :: dir, leveldiff, TYPE

    integer(kind=ik), parameter :: Jmax = 6
    integer(kind=ik) :: send_treecode(1:Jmax)
    integer(kind=ik) :: recv_treecode(1:Jmax)
    integer(kind=ik) :: J=4, i, ishift
    integer(kind=ik) :: i1, i2, g, nDirs, min_size,Nsender
    integer(kind=ik), dimension(3) :: Bs
    integer(kind=ik) :: shifts(1:17, 1:3)
    real(kind=rk) :: x0_send(1:3), dx_send(1:3), x0_recv(1:3), dx_recv(1:3), x0_buffer(1:3), dx_buffer(1:3)
    real(kind=rk) :: r1, r2
    logical :: invalid, coarser_neighbor_possible

    shifts(1,:) = (/0,0,0/)
    shifts(2,:) = (/1,0,0/)
    shifts(3,:) = (/0,1,0/)
    shifts(4,:) = (/0,0,1/)
    shifts(5,:) = (/0,1,1/)
    shifts(6,:) = (/1,1,0/)
    shifts(7,:) = (/1,0,1/)
    shifts(8,:) = (/1,1,1/)
    shifts(9,:) = (/-1,0,1/)
    shifts(10,:) = (/1,-1,-1/)
    shifts(11,:) = (/-1,1,-1/)
    shifts(12,:) = (/-1,-1,-1/)
    shifts(13,:) = (/-1,-1,1/)
    shifts(14,:) = (/0,-1,-1/)
    shifts(15,:) = (/-1,0,0/)
    shifts(16,:) = (/0,-1,0/)
    shifts(17,:) = (/0,0,-1/)

    ijksend = 1
    ijkbuffer = 1
    ishift = 1
    invalid = .true.

    g = params%n_ghosts
    Bs = params%Bs
    nDirs = 74
    if (params%dim == 2) nDirs = 16

    ! check if all dimensions are 1 for this patch, if so, skip it
    if ( ijkrecv(2,1)-ijkrecv(1,1)==0 .and. ijkrecv(2,2)-ijkrecv(1,2)==0 .and. ijkrecv(2,3)-ijkrecv(1,3)==0) then
        ! this neighborhood relation is invalid for the leveldiff, we can skip it.
        return
    endif

    !***************************************************************************
    ! Compute sender bounds from recver bounds
    !***************************************************************************
    do while ( invalid )
        ! just choose any point far away from the boundary, to avoid periodicity
        ! when crossing the periodic boundary, x0 jumps. It may well happen that the block cannot have
        ! a coarser neighbor in the direction. In this case, we use the SHIFT variable to choose another
        ! one.
        call encoding(send_treecode, (/shifts(ishift,1)+(2**J)/2, shifts(ishift,2)+(2**J)/2, &
            shifts(ishift,3)+(2**J)/2/),params%dim, (2**J)**params%dim, J)

        ! fetch the neighbors treecode.
        call get_neighbor_treecode( send_treecode, recv_treecode, dir, &
        leveldiff, J, params%dim, Jmax, coarser_neighbor_possible)

        ! in the coarser neighbor case, finding a valid treecode can be tricky. not all treecodes
        ! have a VALID coarser neighbor in the specified direction. If that happens, we have to choose
        ! another sender coordinate
        if (leveldiff == +1 .and. .not. coarser_neighbor_possible) then
            ! write(*,*) "impossible"
            ! the treecodes have the same coarser ID: this was invalid. they share the same mother block,
            ! which means it cannot be the right coarser neighbor
            ishift = ishift + 1
            if (ishift>size(shifts,1)) call abort(551, "no more shifts")
            cycle
        endif

        call get_block_spacing_origin( params, send_treecode(1:J), x0_send, dx_send )
        call get_block_spacing_origin( params, recv_treecode(1:J-leveldiff), x0_recv, dx_recv )

        ! if the ghost nodes patch on the recver is on the domain boundary, then our algorithm
        ! cannot work: the sender will be at the opposite site. So be sure not to choose such a
        ! patch
        if (patch_crosses_periodic_BC(x0_recv, dx_recv, ijkrecv, params%dim)) then
            call abort(5551, "Yoda, the force is not with us. Bounds pass periodic BC: impossible to auto-generate bounds.")
        endif

        do i = 1, params%dim
            ! shift to zero at the origin (which is g+1, actually)
            r1 = real(ijkrecv(1,i) - (g+1), kind=rk)
            r2 = real(ijkrecv(2,i) - (g+1), kind=rk)

            ! there's a very simple relation between sender and recver boundarys.
            ! we only need to define the recver bounds
            i1 = nint( (r1*dx_recv(i) + x0_recv(i) - x0_send(i)) / dx_send(i) ) + (g+1)
            i2 = nint( (r2*dx_recv(i) + x0_recv(i) - x0_send(i)) / dx_send(i) ) + (g+1)
            ijksend(1:2, i) = (/i1, i2/)

            ! NOTE at the moment, we really just computed the upper/lower bounds of the patch
            ! (selected on the receiver) on the sender. In the next step, we will extend the sender
            ! patch to *contain* what the receiver wants, but also a little more, if required, for interpolation.

            ! check if the bounds we computed are valid, i.e. referring strictly to
            ! the interior of blocks. If not, try again.
            if (i1>=g+1 .and. i1<=Bs(i)+g .and. i2>=g+1 .and. i2<=Bs(i)+g) then
                invalid = .false.
            else
                invalid = .true.
            endif
            ! if one direcion yields unreasonable bounds, we can skip the others of course
            ! thus we exit the i loop
            if (invalid) exit
        enddo

        ishift = ishift + 1

        if (ishift>size(shifts,1)) then
            write(*,*) "Dir=", dir, "leveldiff=", leveldiff
            call abort(386739635, "Master yoda. The bounds computation failed, because we were unable &
            & to find a valid pair of blocks for the given neighbor-code. May the force be with us. Tip: &
            & try changing number_ghost_nodes.")
        endif
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
        ijksend(1, 1:params%dim) = ijksend(1, 1:params%dim) - S
        ijksend(2, 1:params%dim) = ijksend(2, 1:params%dim) + S

        ! then we possibly use asymmetric extension, i.e. make the patch larger in the
        ! direction of the interior. We use this also if we do not have enough points
        ! for the interpolation stencil. This way, one can still set S=0 but ensure having
        ! enough points.

        do i = 1, params%dim
            ! patch at least A (defined above) but possibly the required number to
            ! reach min_size.
            A = max(A, min_size - (ijksend(2,i)-ijksend(1,i)+1))
            if ( ijksend(1, i) <= g+1) then
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
        x0_buffer(:) = real(ijksend(1, :)-(g+1), kind=rk) * dx_send(:) + x0_send(:)
        ! buffer spacing is fine (recv)
        dx_buffer(:) = dx_recv(:)

        do i = 1, params%dim
            ! find recv coordinates in buffer coordinates
            r1 = real(ijkrecv(1,i)-(g+1), kind=rk)
            r2 = real(ijkrecv(2,i)-(g+1), kind=rk)

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

! we cannot use the one in module_mesh because it has a circular dependency (makefile)
subroutine get_block_spacing_origin( params, treecode, x0, dx )

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)             :: params
    integer(kind=ik), intent(in)               :: treecode(1:)
    !> output
    real(kind=rk), dimension(1:3), intent(out) :: x0, dx
    ! loop variables and shortcuts
    integer(kind=ik)                           :: ix,iy,iz,level
    integer(kind=ik), dimension(3)             :: Bs

    bs = params%Bs

    ! fetch this blocks level:
    level = size(treecode)

    ! compute its coordinates in ijk space
    call decoding( treecode, ix, iy, iz, level)

    ! the spacing on a block is the basic spacing Lx/Bs of the coarsest block (if there
    ! is only one block, j=0) divided by 2 for each level, thus the 2^-j factor
    dx = (/1.0_rk, 1.0_rk, 1.0_rk/)
    dx(1:params%dim) = dx(1:params%dim) * 2.0_rk**(-level) / real(Bs(1:params%dim)-1, kind=rk)
    ! note zero based indexing:
    x0 = real( ((/ix,iy,iz/) - 1)*(Bs-1) ,kind=rk) * dx

end subroutine get_block_spacing_origin



subroutine set_recv_bounds( params, data_bounds, neighborhood, level_diff, data_bounds_type, sender_or_receiver)

    !---------------------------------------------------------------------------------------------
    ! modules

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> data_bounds
    integer(kind=ik), intent(inout)                 :: data_bounds(2,3)
    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: level_diff

    ! data_bounds_type
    integer(kind=ik), intent(in)                   :: data_bounds_type
    ! sender or reciver
    character(len=*), intent(in)                   :: sender_or_receiver

    integer(kind=ik) :: g
    integer(kind=ik), dimension(3) :: Bs
    integer(kind=ik) :: sh_start, sh_end

    !---------------------------------------------------------------------------------------------

    ! grid parameter
    Bs    = params%Bs
    g     = params%n_ghosts

    sh_start = 0
    sh_end   = 0

    if ( data_bounds_type == exclude_redundant ) then
        sh_start = 1
    end if
    if ( data_bounds_type == only_redundant ) then
        sh_end = -g
    end if

    ! set 1 and not -1 (or anything else), because 2D bounds ignore 3rd dimension
    ! and thus cycle from 1:1
    data_bounds(:,:) = 1

    if ( params%dim == 3 ) then
        !---3D------3D------3D------3D------3D------3D------3D------3D---
        select case(neighborhood)
            ! '__1/___'
        case(1)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = 1-sh_end
            data_bounds(2,3) = g+1-sh_start

            ! '__2/___'
        case(2)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = Bs(2)+g+sh_start
            data_bounds(2,2) = Bs(2)+g+g+sh_end
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

            ! '__3/___'
        case(3)
            if (level_diff /= 0) return
            data_bounds(1,1) = 1-sh_end
            data_bounds(2,1) = g+1-sh_start
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

            ! '__4/___'
        case(4)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = 1-sh_end
            data_bounds(2,2) = g+1-sh_start
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

            ! '__5/___'
        case(5)
            if (level_diff /= 0) return
            data_bounds(1,1) = Bs(1)+g+sh_start
            data_bounds(2,1) = Bs(1)+g+g+sh_end
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
            data_bounds(1,3) = Bs(3)+g+sh_start
            data_bounds(2,3) = Bs(3)+g+g+sh_end

            ! '_12/___'
        case(7)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = Bs(2)+g+sh_start
            data_bounds(2,2) = Bs(2)+g+g+sh_end
            data_bounds(1,3) = 1-sh_end
            data_bounds(2,3) = g+1-sh_start

            ! '_13/___'
        case(8)
            if (level_diff /= 0) return
            data_bounds(1,1) = 1-sh_end
            data_bounds(2,1) = g+1-sh_start
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = 1-sh_end
            data_bounds(2,3) = g+1-sh_start

            ! '_14/___'
        case(9)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = 1-sh_end
            data_bounds(2,2) = g+1-sh_start
            data_bounds(1,3) = 1-sh_end
            data_bounds(2,3) = g+1-sh_start

            ! '_15/___'
        case(10)
            if (level_diff /= 0) return
            data_bounds(1,1) = Bs(1)+g+sh_start
            data_bounds(2,1) = Bs(1)+g+g+sh_end
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = 1-sh_end
            data_bounds(2,3) = g+1-sh_start

            ! '_62/___'
        case(11)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = Bs(2)+g+sh_start
            data_bounds(2,2) = Bs(2)+g+g+sh_end
            data_bounds(1,3) = Bs(3)+g+sh_start
            data_bounds(2,3) = Bs(3)+g+g+sh_end

            ! '_63/___'
        case(12)
            if (level_diff /= 0) return
            data_bounds(1,1) = 1-sh_end
            data_bounds(2,1) = g+1-sh_start
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = Bs(3)+g+sh_start
            data_bounds(2,3) = Bs(3)+g+g+sh_end

            ! '_64/___'
        case(13)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = 1-sh_end
            data_bounds(2,2) = g+1-sh_start
            data_bounds(1,3) = Bs(3)+g+sh_start
            data_bounds(2,3) = Bs(3)+g+g+sh_end

            ! '_65/___'
        case(14)
            if (level_diff /= 0) return
            data_bounds(1,1) = Bs(1)+g+sh_start
            data_bounds(2,1) = Bs(1)+g+g+sh_end
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g
            data_bounds(1,3) = Bs(3)+g+sh_start
            data_bounds(2,3) = Bs(3)+g+g+sh_end

            ! '_23/___'
        case(15)
            if (level_diff /= 0) return
            data_bounds(1,1) = 1-sh_end
            data_bounds(2,1) = g+1-sh_start
            data_bounds(1,2) = Bs(2)+g+sh_start
            data_bounds(2,2) = Bs(2)+g+g+sh_end
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

            ! '_25/___'
        case(16)
            if (level_diff /= 0) return
            data_bounds(1,1) = Bs(1)+g+sh_start
            data_bounds(2,1) = Bs(1)+g+g+sh_end
            data_bounds(1,2) = Bs(2)+g+sh_start
            data_bounds(2,2) = Bs(2)+g+g+sh_end
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

            ! '_43/___'
        case(17)
            if (level_diff /= 0) return
            data_bounds(1,1) = 1-sh_end
            data_bounds(2,1) = g+1-sh_start
            data_bounds(1,2) = 1-sh_end
            data_bounds(2,2) = g+1-sh_start
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

            ! '_45/___'
        case(18)
            if (level_diff /= 0) return
            data_bounds(1,1) = Bs(1)+g+sh_start
            data_bounds(2,1) = Bs(1)+g+g+sh_end
            data_bounds(1,2) = 1-sh_end
            data_bounds(2,2) = g+1-sh_start
            data_bounds(1,3) = g+1
            data_bounds(2,3) = Bs(3)+g

        case(19,20,21,22)
            data_bounds(1,3) = 1-sh_end
            data_bounds(2,3) = g+1-sh_start
            select case(neighborhood)
            case(19) ! '123/___'
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end

            case(20) ! '134/___'
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start

            case(21) ! '145/___'
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start

            case(22) ! '152/___'
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end

            end select

        case(23,24,25,26)
            data_bounds(1,3) = Bs(3)+g+sh_start
            data_bounds(2,3) = Bs(3)+g+g+sh_end
            select case(neighborhood)
            case(23) ! '623/___'
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end

            case(24) ! '634/___'
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start

            case(25) ! '645/___'
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start

            case(26) ! '652/___'
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end

            end select

        case(27,28,29,30)
            if ( level_diff == -1 ) then
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start
                select case(neighborhood)
                case(27) ! '__1/123'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+2*g

                case(28) ! '__1/134'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(29) ! '__1/145'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+2*g
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(30) ! '__1/152'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+2*g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+2*g
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start
                select case(neighborhood)
                case(27) ! '__1/123'
                    data_bounds(1,1) = g+(Bs(1)+1)/2
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2)+1)/2

                case(28) ! '__1/134'
                    data_bounds(1,1) = g+(Bs(1)+1)/2
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = g+(Bs(2)+1)/2
                    data_bounds(2,2) = Bs(2)+g

                case(29) ! '__1/145'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1)+1)/2
                    data_bounds(1,2) = g+(Bs(2)+1)/2
                    data_bounds(2,2) = Bs(2)+g

                case(30) ! '__1/152'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1)+1)/2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2)+1)/2
                end select

            end if

        case(31,32,33,34)
            if ( level_diff == -1 ) then
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end
                ! first, third dimension
                select case(neighborhood)
                case(31) ! '__2/123'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(32) ! '__2/623'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+2*g

                case(33) ! '__2/152'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+2*g
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(34) ! '__2/652'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+2*g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+2*g

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end
                ! first, third dimension
                select case(neighborhood)
                case(31) ! '__2/123'
                    data_bounds(1,1) = g+(Bs(1)+1)/2
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = g+(Bs(3)+1)/2
                    data_bounds(2,3) = Bs(3)+g

                case(32) ! '__2/623'
                    data_bounds(1,1) = g+(Bs(1)+1)/2
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3)+1)/2

                case(33) ! '__2/152'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1)+1)/2
                    data_bounds(1,3) = g+(Bs(3)+1)/2
                    data_bounds(2,3) = Bs(3)+g

                case(34) ! '__2/652'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1)+1)/2
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3)+1)/2

                end select

            end if

        case(35,36,37,38)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                ! second, third dimension
                select case(neighborhood)
                case(35) ! '__3/123'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+2*g
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
                    data_bounds(2,3) = Bs(3)+2*g

                case(36) ! '__3/623'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+2*g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+2*g

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                ! second, third dimension
                select case(neighborhood)
                case(35) ! '__3/123'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2)+1)/2
                    data_bounds(1,3) = g+(Bs(3)+1)/2
                    data_bounds(2,3) = Bs(3)+g

                case(37) ! '__3/134'
                    data_bounds(1,2) = g+(Bs(2)+1)/2
                    data_bounds(2,2) = Bs(2)+g
                    data_bounds(1,3) = g+(Bs(3)+1)/2
                    data_bounds(2,3) = Bs(3)+g

                case(38) ! '__3/634'
                    data_bounds(1,2) = g+(Bs(2)+1)/2
                    data_bounds(2,2) = Bs(2)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3)+1)/2

                case(36) ! '__3/623'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2)+1)/2
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3)+1)/2

                end select

            end if

        case(39,40,41,42)
            if ( level_diff == -1 ) then
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                ! first, third dimension
                select case(neighborhood)
                case(40) ! '__4/634'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+2*g

                case(39) ! '__4/134'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(41) ! '__4/145'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+2*g
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(42) ! '__4/645'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+2*g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+2*g

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                ! first, third dimension
                select case(neighborhood)
                case(40) ! '__4/634'
                    data_bounds(1,1) = g+(Bs(1)+1)/2
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3)+1)/2

                case(39) ! '__4/134'
                    data_bounds(1,1) = g+(Bs(1)+1)/2
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,3) = g+(Bs(3)+1)/2
                    data_bounds(2,3) = Bs(3)+g

                case(41) ! '__4/145'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1)+1)/2
                    data_bounds(1,3) = g+(Bs(3)+1)/2
                    data_bounds(2,3) = Bs(3)+g

                case(42) ! '__4/645'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1)+1)/2
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3)+1)/2

                end select

            end if

        case(43,44,45,46)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                ! second, third dimension
                select case(neighborhood)
                case(45) ! '__5/152'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+2*g
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
                    data_bounds(2,3) = Bs(3)+2*g

                case(46) ! '__5/652'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+2*g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+2*g

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                ! second, third dimension
                select case(neighborhood)
                case(45) ! '__5/152'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2)+1)/2
                    data_bounds(1,3) = g+(Bs(3)+1)/2
                    data_bounds(2,3) = Bs(3)+g

                case(43) ! '__5/145'
                    data_bounds(1,2) = g+(Bs(2)+1)/2
                    data_bounds(2,2) = Bs(2)+g
                    data_bounds(1,3) = g+(Bs(3)+1)/2
                    data_bounds(2,3) = Bs(3)+g

                case(44) ! '__5/645'
                    data_bounds(1,2) = g+(Bs(2)+1)/2
                    data_bounds(2,2) = Bs(2)+g
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3)+1)/2

                case(46) ! '__5/652'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2)+1)/2
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3)+1)/2

                end select

            end if

        case(47,48,49,50)
            if ( level_diff == -1 ) then
                data_bounds(1,3) = Bs(3)+g+sh_start
                data_bounds(2,3) = Bs(3)+g+g+sh_end
                select case(neighborhood)
                case(47) ! '__6/623'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+2*g

                case(48) ! '__6/634'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(49) ! '__6/645'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+2*g
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(50) ! '__6/652'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+2*g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+2*g

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,3) = Bs(3)+g+sh_start
                data_bounds(2,3) = Bs(3)+g+g+sh_end
                select case(neighborhood)
                case(47) ! '__6/623'
                    data_bounds(1,1) = g+(Bs(1)+1)/2
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2)+1)/2

                case(48) ! '__6/634'
                    data_bounds(1,1) = g+(Bs(1)+1)/2
                    data_bounds(2,1) = Bs(1)+g
                    data_bounds(1,2) = g+(Bs(2)+1)/2
                    data_bounds(2,2) = Bs(2)+g

                case(49) ! '__6/645'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1)+1)/2
                    data_bounds(1,2) = g+(Bs(2)+1)/2
                    data_bounds(2,2) = Bs(2)+g

                case(50) ! '__6/652'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1)+1)/2
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2)+1)/2

                end select

            end if

        case(51,52)
            if ( level_diff == -1 ) then
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start
                select case(neighborhood)
                case(51) ! '_12/123'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g

                case(52) ! '_12/152'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+2*g
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start
                select case(neighborhood)
                case(51) ! '_12/123'
                    data_bounds(1,1) = g+(Bs(1)+1)/2
                    data_bounds(2,1) = Bs(1)+g

                case(52) ! '_12/152'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1)+1)/2
                end select

            end if

        case(53,54)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start
                select case(neighborhood)
                case(54) ! '_13/134'
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(53) ! '_13/123'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+2*g
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start
                select case(neighborhood)
                case(54) ! '_13/134'
                    data_bounds(1,2) = g+(Bs(2)+1)/2
                    data_bounds(2,2) = Bs(2)+g

                case(53) ! '_13/123'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2)+1)/2
                end select

            end if

        case(55,56)
            if ( level_diff == -1 ) then
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start
                select case(neighborhood)
                case(55) ! '_14/134'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g

                case(56) ! '_14/145'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+2*g

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start
                select case(neighborhood)
                case(55) ! '_14/134'
                    data_bounds(1,1) = g+(Bs(1)+1)/2
                    data_bounds(2,1) = Bs(1)+g

                case(56) ! '_14/145'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1)+1)/2

                end select

            end if

        case(57,58)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start
                select case(neighborhood)
                case(57) ! '_15/145'
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(58) ! '_15/152''
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+2*g

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,3) = 1-sh_end
                data_bounds(2,3) = g+1-sh_start
                select case(neighborhood)
                case(57) ! '_15/145'
                    data_bounds(1,2) = g+(Bs(2)+1)/2
                    data_bounds(2,2) = Bs(2)+g

                case(58) ! '_15/152''
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2)+1)/2

                end select

            end if

        case(59,60)
            if ( level_diff == -1 ) then
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end
                data_bounds(1,3) = Bs(3)+g+sh_start
                data_bounds(2,3) = Bs(3)+g+g+sh_end
                select case(neighborhood)
                case(59) ! '_62/623'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g

                case(60) ! '_62/652'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+2*g
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end
                data_bounds(1,3) = Bs(3)+g+sh_start
                data_bounds(2,3) = Bs(3)+g+g+sh_end
                select case(neighborhood)
                case(59) ! '_62/623'
                    data_bounds(1,1) = g+(Bs(1)+1)/2
                    data_bounds(2,1) = Bs(1)+g

                case(60) ! '_62/652'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1)+1)/2
                end select

            end if

        case(61,62)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,3) = Bs(3)+g+sh_start
                data_bounds(2,3) = Bs(3)+g+g+sh_end
                select case(neighborhood)
                case(62) ! '_63/634'
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(61) ! '_63/623'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+2*g
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,3) = Bs(3)+g+sh_start
                data_bounds(2,3) = Bs(3)+g+g+sh_end
                select case(neighborhood)
                case(62) ! '_63/634'
                    data_bounds(1,2) = g+(Bs(2)+1)/2
                    data_bounds(2,2) = Bs(2)+g

                case(61) ! '_63/623'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2)+1)/2
                end select

            end if

        case(63,64)
            if ( level_diff == -1 ) then
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                data_bounds(1,3) = Bs(3)+g+sh_start
                data_bounds(2,3) = Bs(3)+g+g+sh_end
                select case(neighborhood)
                case(63) ! '_64/634'
                    data_bounds(1,1) = 1
                    data_bounds(2,1) = Bs(1)+g

                case(64) ! '_64/645'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = Bs(1)+2*g

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                data_bounds(1,3) = Bs(3)+g+sh_start
                data_bounds(2,3) = Bs(3)+g+g+sh_end
                select case(neighborhood)
                case(63) ! '_64/634'
                    data_bounds(1,1) = g+(Bs(1)+1)/2
                    data_bounds(2,1) = Bs(1)+g

                case(64) ! '_64/645'
                    data_bounds(1,1) = g+1
                    data_bounds(2,1) = g+(Bs(1)+1)/2

                end select

            end if

        case(65,66)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,3) = Bs(3)+g+sh_start
                data_bounds(2,3) = Bs(3)+g+g+sh_end
                select case(neighborhood)
                case(65) ! '_65/645'
                    data_bounds(1,2) = 1
                    data_bounds(2,2) = Bs(2)+g

                case(66) ! '_65/652'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = Bs(2)+2*g

                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,3) = Bs(3)+g+sh_start
                data_bounds(2,3) = Bs(3)+g+g+sh_end
                select case(neighborhood)
                case(65) ! '_65/645'
                    data_bounds(1,2) = g+(Bs(2)+1)/2
                    data_bounds(2,2) = Bs(2)+g

                case(66) ! '_65/652'
                    data_bounds(1,2) = g+1
                    data_bounds(2,2) = g+(Bs(2)+1)/2

                end select

            end if

        case(67,68)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end
                select case(neighborhood)
                case(67) ! '_23/123'
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(68) ! '_23/236''
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+2*g
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end
                select case(neighborhood)
                case(67) ! '_23/123'
                    data_bounds(1,3) = g+(Bs(3)+1)/2
                    data_bounds(2,3) = Bs(3)+g

                case(68) ! '_23/236''
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3)+1)/2
                end select

            end if

        case(69,70)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end
                select case(neighborhood)
                case(69) ! '_25/152'
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(70) ! '_25/652''
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+2*g
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end
                select case(neighborhood)
                case(69) ! '_25/152'
                    data_bounds(1,3) = g+(Bs(3)+1)/2
                    data_bounds(2,3) = Bs(3)+g

                case(70) ! '_25/652''
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3)+1)/2
                end select

            end if

        case(71,72)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                select case(neighborhood)
                case(71) ! '_43/134'
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(72) ! '_43/634''
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+2*g
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                select case(neighborhood)
                case(71) ! '_43/134'
                    data_bounds(1,3) = g+(Bs(3)+1)/2
                    data_bounds(2,3) = Bs(3)+g

                case(72) ! '_43/634''
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3)+1)/2
                end select

            end if

        case(73,74)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                select case(neighborhood)
                case(73) ! '_45/145'
                    data_bounds(1,3) = 1
                    data_bounds(2,3) = Bs(3)+g

                case(74) ! '_45/645'
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = Bs(3)+2*g
                end select

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start
                select case(neighborhood)
                case(73) ! '_45/145'
                    data_bounds(1,3) = g+(Bs(3)+1)/2
                    data_bounds(2,3) = Bs(3)+g

                case(74) ! '_45/645'
                    data_bounds(1,3) = g+1
                    data_bounds(2,3) = g+(Bs(3)+1)/2
                end select

            end if

        end select
    else
        !---2D------2D------2D------2D------2D------2D------2D------2D---
        select case(neighborhood)
            ! '__N'
        case(1)
            if (level_diff /= 0) return
            data_bounds(1,1) = Bs(1)+g+sh_start
            data_bounds(2,1) = Bs(1)+g+g+sh_end
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g

            ! '__E'
        case(2)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = 1-sh_end
            data_bounds(2,2) = g+1-sh_start

            ! '__S'
        case(3)
            if (level_diff /= 0) return
            data_bounds(1,1) = 1-sh_end
            data_bounds(2,1) = g+1-sh_start
            data_bounds(1,2) = g+1
            data_bounds(2,2) = Bs(2)+g

            ! '__W'
        case(4)
            if (level_diff /= 0) return
            data_bounds(1,1) = g+1
            data_bounds(2,1) = Bs(1)+g
            data_bounds(1,2) = Bs(2)+g+sh_start
            data_bounds(2,2) = Bs(2)+g+g+sh_end

            ! '_NE'
        case(5)
            data_bounds(1,1) = Bs(1)+g+sh_start
            data_bounds(2,1) = Bs(1)+g+g+sh_end
            data_bounds(1,2) = 1-sh_end
            data_bounds(2,2) = g+1-sh_start

            ! '_NW'
        case(6)
            data_bounds(1,1) = Bs(1)+g+sh_start
            data_bounds(2,1) = Bs(1)+g+g+sh_end
            data_bounds(1,2) = Bs(2)+g+sh_start
            data_bounds(2,2) = Bs(2)+g+g+sh_end

            ! '_SE'
        case(7)
            data_bounds(1,1) = 1-sh_end
            data_bounds(2,1) = g+1-sh_start
            data_bounds(1,2) = 1-sh_end
            data_bounds(2,2) = g+1-sh_start

            ! '_SW'
        case(8)
            data_bounds(1,1) = 1-sh_end
            data_bounds(2,1) = g+1-sh_start
            data_bounds(1,2) = Bs(2)+g+sh_start
            data_bounds(2,2) = Bs(2)+g+g+sh_end

            ! 'NNE'
        case(9)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,2) = 1
                data_bounds(2,2) = Bs(2)+g

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,2) = g+(Bs(2)+1)/2
                data_bounds(2,2) = Bs(2)+g

            end if

            ! 'NNW'
        case(10)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs(2)+g+g

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = Bs(1)+g+sh_start
                data_bounds(2,1) = Bs(1)+g+g+sh_end
                data_bounds(1,2) = g+1
                data_bounds(2,2) = g+(Bs(2)+1)/2

            end if

            ! 'SSE'
        case(11)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = 1
                data_bounds(2,2) = Bs(2)+g

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = g+(Bs(2)+1)/2
                data_bounds(2,2) = Bs(2)+g

            end if

            ! 'SSW'
        case(12)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = g+1
                data_bounds(2,2) = Bs(2)+g+g

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = 1-sh_end
                data_bounds(2,1) = g+1-sh_start
                data_bounds(1,2) = g+1
                data_bounds(2,2) = g+(Bs(2)+1)/2

            end if

            ! 'ENE'
        case(13)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs(1)+g+g
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = g+1
                data_bounds(2,1) = g+(Bs(1)+1)/2
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start

            end if

            ! 'ESE'
        case(14)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = Bs(1)+g
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = g+(Bs(1)+1)/2
                data_bounds(2,1) = Bs(1)+g
                data_bounds(1,2) = 1-sh_end
                data_bounds(2,2) = g+1-sh_start

            end if

            ! 'WNW'
        case(15)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = g+1
                data_bounds(2,1) = Bs(1)+g+g
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = g+1
                data_bounds(2,1) = g+(Bs(1)+1)/2
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end

            end if

            ! 'WSW'
        case(16)
            if ( level_diff == -1 ) then
                data_bounds(1,1) = 1
                data_bounds(2,1) = Bs(1)+g
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end

            elseif ( level_diff == 1 ) then
                data_bounds(1,1) = g+(Bs(1)+1)/2
                data_bounds(2,1) = Bs(1)+g
                data_bounds(1,2) = Bs(2)+g+sh_start
                data_bounds(2,2) = Bs(2)+g+g+sh_end

            end if

        end select
    end if
end subroutine set_recv_bounds
