subroutine threshold_block( params, u, thresholding_component, refinement_status, norm, level, eps )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data - this routine is called on one block only, not on the entire grid. hence th 4D array.
    real(kind=rk), intent(inout)        :: u(:, :, :, :)
    !> it can be useful not to consider all components for thresholding here.
    !! e.g. to work only on the pressure or vorticity.
    logical, intent(in)                 :: thresholding_component(:)
    !> main output of this routine is the new satus
    integer(kind=ik), intent(out)       :: refinement_status
    ! If we use L2 or H1 normalization, the threshold eps is level-dependent, hence
    ! we pass the level to this routine
    integer(kind=ik), intent(in)        :: level
    !
    real(kind=rk), intent(inout)        :: norm( size(u,4) )
    ! if different from the default eps (params%eps), you can pass a different value here. This is optional
    ! and used for example when thresholding the mask function.
    real(kind=rk), intent(in), optional :: eps

    integer(kind=ik)                    :: dF, i, j, l, p
    real(kind=rk)                       :: detail( size(u,4) )
    integer(kind=ik)                    :: g, dim, Jmax, nx, ny, nz, nc
    integer(kind=ik), dimension(3)      :: Bs
    real(kind=rk)                       :: t0, eps2
    real(kind=rk), allocatable, dimension(:,:,:,:) :: u_wc, sc, wcx, wcy, wcxy

    t0 = MPI_Wtime()
    nx = size(u, 1)
    ny = size(u, 2)
    nz = size(u, 3)
    nc = size(u, 4)
    Bs = params%Bs
    g  = params%g
    dim = params%dim
    Jmax = params%Jmax
    detail = -1.0_rk
    allocate(  sc(1:nx, 1:ny, 1:nz, 1:nc) )
    allocate( wcx(1:nx, 1:ny, 1:nz, 1:nc) )
    allocate( wcy(1:nx, 1:ny, 1:nz, 1:nc) )
    allocate(wcxy(1:nx, 1:ny, 1:nz, 1:nc) )
    allocate(u_wc(1:nx, 1:ny, 1:nz, 1:nc) )


#ifdef DEV
    if (.not. allocated(params%GD)) call abort(1213149, "The cat is angry: Wavelet-setup not yet called?")
    if (modulo(Bs(1),2) /= 0) call abort(1213150, "The dog is angry: Block size must be even.")
    if (modulo(Bs(2),2) /= 0) call abort(1213150, "The dog is angry: Block size must be even.")
#endif

    ! perform the wavlet decomposition of the block
    ! Note we could not reonstruct here, because the neighboring WC/SC are not
    ! synced. However, here, we only check the details on a block, so there is no
    ! need for reconstruction.
    u_wc = u
    call waveletDecomposition_block(params, u_wc) ! data on u (WC/SC) now in Spaghetti order

    ! NOTE: if the coarse reconstruction is performed before this routine is called, then
    ! the WC affected by the coarseExtension are automatically zero. There is no need to reset
    ! them again. -> checked in postprocessing that this is indeed the case.
    call spaghetti2mallat_block(params, u_wc, sc, wcx, wcy, wcxy)

    do p = 1, nc
        ! if all details are smaller than C_eps, we can coarsen.
        ! check interior WC only
        detail(p) = maxval(abs(wcx(g+1:Bs(1)+g,g+1:Bs(2)+g,:,p)))
        detail(p) = max( detail(p), maxval(abs(wcy(g+1:Bs(1)+g,g+1:Bs(2)+g,:,p))) )
        detail(p) = max( detail(p), maxval(abs(wcxy(g+1:Bs(1)+g,g+1:Bs(2)+g,:,p))) )

        detail(p) = detail(p) / norm(p)
    enddo

    deallocate(u_wc, sc, wcx, wcy, wcxy)

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
    eps2 = params%eps
    ! but if we pass another one, use that.
    if (present(eps)) eps2 = eps


    select case(params%eps_norm)
    case ("Linfty")
        ! do nothing, our wavelets are normalized in L_infty norm by default, hence
        ! a simple threshold controls this norm
        eps2 = eps2

    case ("L2")
        ! If we want to control the L2 norm (with wavelets that are normalized in Linfty norm)
        ! we have to have a level-dependent threshold
        eps2 = eps2 * ( 2.0_rk**(-dble((level-Jmax)*params%dim)/2.0_rk) )
        ! if (params%dim==2) eps2 = eps2*0.1

    case ("H1")
        ! H1 norm mimicks filtering of vorticity
        eps2 = eps2 * ( 2**(-level*(params%dim+2.0_rk)*0.5_rk) )

    case default
        call abort(20022811, "ERROR:threshold_block.f90:Unknown wavelet normalization!")

    end select

    ! evaluate criterion: if this blocks detail is smaller than the prescribed precision,
    ! the block is tagged as "wants to coarsen" by setting the tag -1
    ! note gradedness and completeness may prevent it from actually going through with that
    if ( maxval(detail) < eps2) then
        ! coarsen block, -1
        refinement_status = -1
    else
        refinement_status = 0
    end if


    ! timings
    call toc( "threshold_block (w/o ghost synch.)", MPI_Wtime() - t0 )
end subroutine threshold_block
