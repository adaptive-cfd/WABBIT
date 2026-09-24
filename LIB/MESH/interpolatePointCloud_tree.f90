
! Interpolates the value of a given field (hvy_block) at locations given by a list 
! of points (= a point cloud). The field is stored on a wabbit-type grid
! (block-based adaptive), and the points are given as a simple array.
! The routine broadcasts the result of the interpolation to all mpiranks.
! We do not assume any structure in the point cloud - the points do not need to be 
! ordered nor is there an implied neighborhood relation between them.
! The routine always interpolates all components of the data vector, i.e. if you
! pass the state vector, it interpolates all its components (ux,uy,uz,p in the ACM case),
! even if you need only the pressure. The overhead is not too large.
subroutine interpolatePointCloud_tree( params, hvy_block, tree_ID, xq, fq, interpolationMethod, sync )

    type (type_params), intent(inout) :: params                   !> user defined parameter structure
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)    !> heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID
    ! point cloud - query points xq(1:3, 1:npoints)
    real(kind=rk), intent(in) :: xq(:,:)
    ! result of the interpolation (function values at query points)
    real(kind=rk), intent(inout) :: fq(:,:)
    ! interpolation method
    character(len=*), intent(in) :: interpolationMethod
    ! if the ghost nodes are already synced, we can skip the sync'ing here
    logical, intent(in) :: sync

    integer(kind=ik) :: ipoint, npoints, interpolationMethod_code, N_support, ix0, iy0, iz0, dim
    integer(kind=ik) :: ix, iy, iz, mpierr, g, lgt_id, hvy_id, k
    real(kind=rk) :: x(1:3), xx, yy, zz, delx, delz, dely, x0(1:3), dx(1:3)

    npoints = size(xq, 2)
    dim     = params%dim

    ! nothing to do?
    if (npoints == 0) return

    ! determine how many ghost points we need.
    select case(interpolationMethod)
    case ("floor")
        g = 0
        interpolationMethod_code = 0 ! avoid string comparison
        N_support = 0
    case ("linear")
        g = 1
        interpolationMethod_code = 1 ! avoid string comparison
        N_support = 1
    case ('delta')
        g = 3
        interpolationMethod_code = 2 ! avoid string comparison
        N_support = 3
    case default
        call abort(26094257, "Unknown interpolation method. Time for a coffee, you need it.")
    end select

    ! check if we have enough ghost nodes
    if ( params%g < N_support ) then
        call abort(26094254,"Error: not enough ghostpoints for interpolation kernel. Knock knock. Who's there? A Bug. Haha.")
    endif

    ! BEFORE WE CAN INTERPOLATE THE GHOTS NODES NEED TO BE FILLED
    if ((sync) .and. (interpolationMethod_code > 0)) then
        ! floor interpolation does not require sync'ing by design
        call sync_ghosts_tree( params, hvy_block, tree_ID=tree_ID, g_minus=g, g_plus=g )
    endif

    npoints = size(xq, 2)
    dim     = params%dim


    if (size(fq,2) /= npoints) then
        call abort(26094258, "Output array size is incorrect. You seem tired, take a nap.")
    endif
    if (size(fq,1) /= size(hvy_block,4)) then
        call abort(26094259, "Output array size is incorrect. All roads lead to Rome. This one doesnt.")
    endif

    ! set the array to -Inf. Why? If a point is NOT interpolated by a CPU, then the CPU
    ! will not touch the data. So it remains that large negative number. In other words, the array looks like this:
    ! CPU1 = (/12.2, 13.4, -Inf, -Inf/)
    ! CPU2 = (/-Inf, -Inf, 7.2,  -39.9/)
    ! now I can just take the MAXIMUM of all values across all CPU
    ! and the final result is:
    ! (/12.2, 13.4, 7.2, -39.9/)
    ! on all CPU.
    fq = -huge(1.0_rk)

    ! now we can interpolate each point
    do ipoint = 1, npoints
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k,tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            call get_block_spacing_origin( params, lgt_id, x0, dx )

            ! not all blocks are relevant: only one single block contains the interpolation 
            ! point we are looking at. Find the block! 
            if (pointInBlock_block(params, xq(:, ipoint), x0, dx)) then
                ! convert interpolation point to integer, nearest (lower) integer
                x = xq(1:3, ipoint) - x0(1:3)
                ix0 = floor( x(1) / dx(1)) + (params%g + 1)
                iy0 = floor( x(2) / dx(2)) + (params%g + 1)
                iz0 = floor( x(3) / dx(3)) + (params%g + 1)

                if (dim==2) iz0 = 1
                

                select case(interpolationMethod_code) ! no string comparison here
                case (0)
                    ! nearest neighbor interpolation (floor; not actually the nearest neighbor)
                    fq(:, ipoint) = hvy_block( ix0, iy0, iz0, :, hvy_id)

                case (1)
                    ! linear interpolation
                    
                    ! fq is a convolution of the data with a tensor-product style
                    ! delta function
                    fq(:, ipoint) = 0.0_rk
                    do iz = iz0-merge(0, N_support, dim==2), iz0+merge(0, N_support, dim==2) 
                        zz = real(iz - (params%g + 1), rk) * dx(3)
                        delz = linearInterpolationKernel(abs(zz - x(3)),dx(3))
    
                        do iy = iy0-N_support, iy0+N_support
                            yy = real(iy - (params%g + 1), rk) * dx(2)
                            dely = linearInterpolationKernel(abs(yy - x(2)),dx(2))
    
                            do ix = ix0-N_support, ix0+N_support
                                xx = real(ix - (params%g + 1), rk) * dx(1)
                                delx = linearInterpolationKernel(abs(xx - x(1)),dx(1))
    
                                fq(:, ipoint) = fq(:, ipoint) + delx * dely * delz * hvy_block( ix, iy, iz, :, hvy_id )
                            enddo
                        enddo
                    enddo   
                case (2)
                    ! delta kernel interpolation
                    N_support = 3
                    
                    ! fq is a convolution of the data with a tensor-product style
                    ! delta function
                    fq(:, ipoint) = 0.0_rk
                    do iz = iz0-merge(0, N_support, dim==2), iz0+merge(0, N_support, dim==2) 
                        zz = real(iz - (params%g + 1), rk) * dx(3)
                        delz = deltaInterpolationKernel(abs(zz - x(3)),dx(3))
    
                        do iy = iy0-N_support, iy0+N_support
                            yy = real(iy - (params%g + 1), rk) * dx(2)
                            dely = deltaInterpolationKernel(abs(yy - x(2)),dx(2))
    
                            do ix = ix0-N_support, ix0+N_support
                                xx = real(ix - (params%g + 1), rk) * dx(1)
                                delx = deltaInterpolationKernel(abs(xx - x(1)),dx(1))
    
                                fq(:, ipoint) = fq(:, ipoint) + delx * dely * delz * hvy_block( ix, iy, iz, :, hvy_id )
                            enddo
                        enddo
                    enddo   

                case default
                    call abort(2505208, 'interpolationMethod_code must be 0, 1, or 2. Not a good day for coding, maybe for hiking?')

                end select
            endif
        end do
    end do

    ! maximum across mpiranks
    call MPI_allreduce( MPI_IN_PLACE, fq, size(fq), MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)


end subroutine interpolatePointCloud_tree