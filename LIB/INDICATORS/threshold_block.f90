subroutine threshold_block( params, u, thresholding_component, refinement_status, norm, level, input_is_WD, eps)
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data - this routine is called on one block only, not on the entire grid. hence th 4D array
    !! They are expected to be already wavelet decomposed in Mallat-ordering
    real(kind=rk), intent(inout)        :: u(:, :, :, :)
    !> it can be useful not to consider all components for thresholding here.
    !! e.g. to work only on the pressure or vorticity.
    logical, intent(in)                 :: thresholding_component(:)
    !> main output of this routine is the new satus
    integer(kind=ik), intent(out)       :: refinement_status
    !> If we use L2 or H1 normalization, the threshold eps is level-dependent, hence
    !! we pass the level to this routine
    integer(kind=ik), intent(in)        :: level
    logical, intent(in)                 :: input_is_WD                       !< flag if hvy_block is already wavelet decomposed
    real(kind=rk), intent(inout)        :: norm( size(u,4) )
    !> if different from the default eps (params%eps), you can pass a different value here. This is optional
    !! and used for example when thresholding the mask function.
    real(kind=rk), intent(in), optional :: eps

    integer(kind=ik)                    :: dF, i, j, l, p
    real(kind=rk)                       :: detail( size(u,4) )
    integer(kind=ik)                    :: g, dim, Jmax, nx, ny, nz, nc
    integer(kind=ik), dimension(3)      :: Bs
    real(kind=rk)                       :: eps_use
    ! The WC array contains SC (scaling function coeffs) as well as all WC (wavelet coeffs)
    ! Note: the precise naming of SC/WC is not really important. we just apply
    ! the correct decomposition/reconstruction filters - thats it.
    !
    ! INDEX            2D     3D     LABEL (NAME)
    ! -----            --    ---     ---------------------------------
    ! wc(:,:,:,:,1)    HH    HHH     sc scaling function coeffs
    ! wc(:,:,:,:,2)    HG    HGH     wcx wavelet coeffs
    ! wc(:,:,:,:,3)    GH    GHH     wcy wavelet coeffs
    ! wc(:,:,:,:,4)    GG    GGH     wcxy wavelet coeffs
    ! wc(:,:,:,:,5)          HHG     wcz wavelet coeffs
    ! wc(:,:,:,:,6)          HGG     wcxz wavelet coeffs
    ! wc(:,:,:,:,7)          GHG     wcyz wavelet coeffs
    ! wc(:,:,:,:,8)          GGG     wcxyz wavelet coeffs
    !
    real(kind=rk), allocatable, dimension(:,:,:,:,:), save :: wc
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: u_wc

    nx     = size(u, 1)
    ny     = size(u, 2)
    nz     = size(u, 3)
    nc     = size(u, 4)
    Bs     = params%Bs
    g      = params%g
    dim    = params%dim
    Jmax   = params%Jmax
    detail = -1.0_rk

    if (allocated(u_wc)) then
        if (size(u_wc, 4) > nc) deallocate(u_wc)
    endif
    if (allocated(wc)) then
        if (size(wc, 4) > nc) deallocate(wc)
    endif
    if (.not. allocated(wc)) allocate(wc(1:nx, 1:ny, 1:nz, 1:nc, 1:8) )
    if (.not. allocated(u_wc)) allocate(u_wc(1:nx, 1:ny, 1:nz, 1:nc ) )


#ifdef DEV
    if (.not. allocated(params%GD)) call abort(1213149, "The cat is angry: Wavelet-setup not yet called?")
    if (modulo(Bs(1),2) /= 0) call abort(1213150, "The dog is angry: Block size must be even.")
    if (modulo(Bs(2),2) /= 0) call abort(1213150, "The dog is angry: Block size must be even.")
#endif

    if (.not. input_is_WD) then
        u_wc = u
        call waveletDecomposition_block(params, u_wc) ! data on u (WC/SC) now in Spaghetti order
        call Spaghetti2inflatedMallat_block(params, u_wc, wc)
    else
        call Spaghetti2inflatedMallat_block(params, u, wc)
    endif


    if (params%dim == 2) then
        do p = 1, nc
            ! if all details are smaller than C_eps, we can coarsen.
            ! check interior WC only
            detail(p) = maxval( abs(wc(g+1:Bs(1)+g, g+1:Bs(2)+g, :, p, 2:4)) )
        enddo
    else
        do p = 1, nc
            ! if all details are smaller than C_eps, we can coarsen.
            ! check interior WC only
            detail(p) = maxval( abs(wc(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, p, 2:8)) )
        enddo
    endif

    detail(1:nc) = detail(1:nc) / norm(1:nc)

    ! We could disable detail checking for qtys we do not want to consider,
    ! but this is more work and selective thresholding is rarely used
    do p = 1, nc
        if (.not. thresholding_component(p)) detail(p) = 0.0_rk
    enddo

    ! ich habe die wavelet normalization ausgebruetet und aufgeschrieben.
    ! ich schicke dir die notizen gleich (photos).
    !
    ! also wir brauchen einen scale(level)- dependent threshold, d.h. \epsilon_j
    ! zudem ist dieser abhaengig von der raum dimension d.
    !
    ! Fuer die L^2 normalisierung (mit wavelets welche in der L^\infty norm normalisiert sind) haben wir
    !
    ! \epsilon_j = 2^{-jd/2} \epsilon
    !
    ! d.h. der threshold wird kleiner auf kleinen skalen.
    !
    ! Fuer die vorticity (anstatt der velocity) kommt nochmal ein faktor 2^{-j} dazu, d.h.
    !
    ! \epsilon_j = 2^{-j(d+2)/2} \epsilon
    !
    ! Zum testen waere es gut in 1d oder 2d zu pruefen, ob die L^2 norm von u - u_\epsilon
    ! linear mit epsilon abnimmt, das gleiche koennte man auch fuer H^1 (philipp koennte dies doch mal ausprobieren?).
    !
    ! fuer CVS brauchen wir dann noch \epsilon was von Z (der enstrophy) und der feinsten
    ! aufloesung abhaengt. fuer L^2 normalisierte wavelets ist
    ! der threshold:
    !
    ! \epsilon = \sqrt{2/3 \sigma^2 \ln N}
    !
    ! wobei \sigma^2 die varianz (= 2 Z) der incoh. vorticity ist.
    ! typischerweise erhaelt man diese mit 1-3 iterationen.
    ! als ersten schritt koennen wir einfach Z der totalen stroemung nehmen.
    ! N ist die maximale aufloesung, typicherweise 2^{d J}.
    !

    ! default thresholding level is the one in the parameter struct
    eps_use = params%eps
    ! but if we pass another one, use that.
    if (present(eps)) eps_use = eps

    ! write(*, '("Detail ", es8.1, " eps ", es8.1)') detail(1), eps_use

    select case(params%eps_norm)
    case ("Linfty")
        ! do nothing, our wavelets are normalized in L_infty norm by default, hence
        ! a simple threshold controls this norm
        eps_use = eps_use

    case ("L2")
        ! If we want to control the L2 norm (with wavelets that are normalized in Linfty norm)
        ! we have to have a level-dependent threshold
        eps_use = eps_use * ( 2.0_rk**(-dble((level-Jmax)*params%dim)/2.0_rk) )

    case ("H1")
        ! H1 norm mimicks filtering of vorticity
        eps_use = eps_use * ( 2**(-level*(params%dim+2.0_rk)*0.5_rk) )

    case default
        call abort(20022811, "ERROR:threshold_block.f90:Unknown wavelet normalization!")

    end select

    ! evaluate criterion: if this blocks detail is smaller than the prescribed precision,
    ! the block is tagged as "wants to coarsen" by setting the tag -1
    ! note gradedness and completeness may prevent it from actually going through with that
    if ( maxval(detail) < eps_use) then
        ! coarsen block, -1
        refinement_status = -1
    else
        refinement_status = 0

        ! if (level == 4 .and. params%rank == 0) then
        !     call Spaghetti2Mallat_block(params, u, u_wc)
        !     write(*, '("Detail ", es8.1, " eps ", es8.1)') detail(1), eps_use
        !     do ny = 1,26
        !         write(*, '(26(es8.1, 1x))') u_wc(1:26, ny, 1, 1)
        !     enddo
        ! endif
    end if
end subroutine threshold_block
