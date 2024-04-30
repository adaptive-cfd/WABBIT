subroutine filter_wrapper(time, params, hvy_block, hvy_tmp, hvy_mask, tree_ID)
    implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(inout)   :: params                       !> user defined parameter structure, hvy_active
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !> heavy data array - block data
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> hvy_mask are qty that depend on the grid and not explicitly on time
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID
    real(kind=rk), dimension(3)         :: dx, x0                       !> spacing and origin of a block
    integer(kind=ik)                    :: k, dF, neqn, lgt_id, i       ! loop variables
    integer(kind=ik)                    :: g, a, nx, ny, nz, nc
    integer(kind=ik), dimension(3)      :: Bs
    integer(kind=2)                     :: n_domain(1:3)                ! surface normal
    integer(kind=ik)                    :: level, hvy_id
    real(kind=rk)                       :: stencil(-19:19)          ! stencil array, note: size is fixed
    integer(kind=ik)                    :: stencil_size         ! filter position (array postion of value to filter)
    real(kind=rk), allocatable, save    :: u_filtered(:,:,:,:)

    Bs = params%Bs
    g  = params%g
    n_domain = 0
    nx = size(hvy_block, 1)
    ny = size(hvy_block, 2)
    nz = size(hvy_block, 3)
    nc = size(hvy_block, 4)

    if (.not. allocated(u_filtered)) allocate( u_filtered(1:nx,1:ny,1:nz,1:nc) )

    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)

        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        select case(params%filter_type)
        case('explicit_3pt')
            ! This filter comes from
            ! A general class of commutative filters for LES in complex geometries
            ! OV Vasilyev, TS Lund, P Moin - Journal of computational physics, 1998
            ! 
            ! Table I, case 1

            stencil_size = 3
            a = (stencil_size-1)/2
            stencil(-a:+a) = (/ 1.0_rk/4.0_rk, -1.0_rk/2.0_rk, 1.0_rk/4.0_rk /)

        case('explicit_5pt')
            ! This filter comes from
            ! A general class of commutative filters for LES in complex geometries
            ! OV Vasilyev, TS Lund, P Moin - Journal of computational physics, 1998
            ! 
            ! Table I, case 5
            stencil_size = 5
            a = (stencil_size-1)/2
            stencil(-a:+a) = (/ -1.0_rk/ 16.0_rk, &
            1.0_rk/  4.0_rk, &
            -3.0_rk/  8.0_rk, &
            1.0_rk/  4.0_rk, &
            -1.0_rk/ 16.0_rk/)

        case('explicit_7pt')
            ! same as above case 10
            stencil_size = 7
            a = (stencil_size-1)/2
            stencil(-a:+a) = (/  1.0_rk/ 64.0_rk, &
            -3.0_rk/ 32.0_rk, &
            15.0_rk/ 64.0_rk, &
            -5.0_rk/ 16.0_rk, &
            15.0_rk/ 64.0_rk, &
            -3.0_rk/ 32.0_rk, &
            1.0_rk/ 64.0_rk/)

        case('explicit_9pt')
            ! not included in Vasilyev 1998
            stencil_size = 9
            a = (stencil_size-1)/2
            stencil(-a:+a) = (/ -1.0_rk/256.0_rk, &
            1.0_rk/ 32.0_rk, &
            -7.0_rk/ 64.0_rk, &
            7.0_rk/ 32.0_rk, &
            -35.0_rk/128.0_rk, &
            7.0_rk/ 32.0_rk, &
            -7.0_rk/ 64.0_rk, &
            1.0_rk/ 32.0_rk, &
            -1.0_rk/256.0_rk/)

        case('explicit_11pt')
            ! not included in Vasilyev 1998
            stencil_size = 11
            a = (stencil_size-1)/2
            stencil(-a:+a) = (/  1.0_rk/1024.0_rk, &
            -5.0_rk/ 512.0_rk, &
            45.0_rk/1024.0_rk, &
            -15.0_rk/ 128.0_rk, &
            105.0_rk/ 512.0_rk, &
            -63.0_rk/ 256.0_rk, &
            105.0_rk/ 512.0_rk, &
            -15.0_rk/ 128.0_rk, &
            45.0_rk/1024.0_rk, &
            -5.0_rk/ 512.0_rk, &
            1.0_rk/1024.0_rk/)
        end select

        stencil(0) = stencil(0) + 1.0_rk

        call blockFilterXYZ_vct( params, hvy_block(:,:,:,:, hvy_id), u_filtered, stencil(-a:+a), -a, +a)
        hvy_block(:,:,:,:, hvy_id) = u_filtered

    enddo
end subroutine filter_wrapper
