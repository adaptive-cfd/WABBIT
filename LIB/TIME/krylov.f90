subroutine krylov_time_stepper(time, dt, params, lgt_block, hvy_block, hvy_work, &
    hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n)
    ! use module_blas
    implicit none

    !---------------------------------------------------------------------------
    !> time varible
    real(kind=rk), intent(inout)        :: time, dt
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !---------------------------------------------------------------------------
    integer :: M_max ! M is M_krylov number of subspace
    real(kind=rk), allocatable, save :: H(:,:), phiMat(:,:), H_tmp(:,:)
    integer :: M_iter
    integer :: i, j, k, l, iter
    real(kind=rk) :: normv, eps, beta, err_tolerance
    real(kind=rk) :: h_klein, err, t0

    M_max = params%M_krylov

    ! allocate matrices with largest admissible
    if (.not. allocated(H)) then
        allocate( H(M_max+2,M_max+2) )
        allocate( H_tmp(M_max+2,M_max+2) )
        allocate( phiMat(M_max+2,M_max+2) )
    endif

    phiMat = 0.0_rk
    H = 0.0_rk

    ! check if we have enough RHS slots available: they must be
    if (size(hvy_work,6) < M_max+3 ) then
        call abort(18101817, "not enough registers for krylov method")
    endif

    ! synchronize ghost nodes
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    ! calculate time step
    call calculate_time_step(params, time, hvy_block, hvy_active, hvy_n, lgt_block, &
        lgt_active, lgt_n, dt)

    ! compute norm "normv" of input state vector ("hvy_block")
    call wabbit_norm( params, hvy_block, hvy_active, hvy_n, normv )
    eps = normv * sqrt(epsilon(1.0_rk))


    ! the very last slot (M+3) is the "reference right hand side"
    call RHS_wrapper( time, params, hvy_block, hvy_work(:,:,:,:,:,M_max+3), lgt_block, hvy_active, hvy_n)

    ! compute norm "beta", which is the norm of the reference RHS evaluation
    ! NSF: this guy is called normuu
    call wabbit_norm( params, hvy_work(:,:,:,:,:,M_max+3), hvy_active, hvy_n, beta )
    if (beta < epsilon(1.0_rk)) beta = 1.0_rk

    ! start iteration, fill first slot
    do k = 1, hvy_n
        hvy_work(:,:,:,:,hvy_active(k),1) = hvy_work(:,:,:,:,hvy_active(k),M_max+3) / beta
    enddo

    !**************************************!
    !*** begin interations              ***!
    !**************************************!
    ! we loop over the full size of subspace, then exit prematurely
    ! if possible.
    do M_iter = 1, M_max

        ! perturbed right hand side is first-to-last (M+2) slot
        do k = 1, hvy_n
            hvy_work(:,:,:,:,hvy_active(k),M_max+2) = hvy_block(:,:,:,:,hvy_active(k)) &
            + eps * hvy_work(:,:,:,:,hvy_active(k),M_iter)
        enddo

        ! call RHS with perturbed state vector, stored in slot (M_max+1)
        call sync_ghosts( params, lgt_block, hvy_work(:,:,:,:,:,M_max+2), hvy_neighbor, hvy_active, hvy_n )
        call RHS_wrapper( time, params, hvy_work(:,:,:,:,:,M_max+2), hvy_work(:,:,:,:,:,M_max+1), &
        lgt_block, hvy_active, hvy_n)

        ! linearization of RHS slot (M_max+1)
        do k = 1, hvy_n
            hvy_work(:,:,:,:,hvy_active(k),M_max+1) = ( hvy_work(:,:,:,:,hvy_active(k),M_max+1) &
            -hvy_work(:,:,:,:,hvy_active(k),M_max+3) ) / eps
        enddo

        ! --- inner loop ----
        ! --- ARNOLDI ---
        do iter = 1, M_iter
            call scalarproduct(params, hvy_work(:,:,:,:,:,iter), hvy_work(:,:,:,:,:,M_max+1),&
             hvy_active, hvy_n, H(iter, M_iter) )

             do k = 1, hvy_n
                 hvy_work(:,:,:,:,hvy_active(k),M_max+1) = hvy_work(:,:,:,:,hvy_active(k),M_max+1) &
                 - H(iter,M_iter) * hvy_work(:,:,:,:,hvy_active(k),iter)
             enddo
        enddo
        ! end of inner i =1:j loop
        call wabbit_norm(params, hvy_work(:,:,:,:,:,M_max+1), hvy_active, hvy_n, H(M_iter+1,M_iter) )

        do k = 1, hvy_n
            hvy_work(:,:,:,:,hvy_active(k),M_iter+1) = hvy_work(:,:,:,:,hvy_active(k),M_max+1) / H(M_iter+1,M_iter)
        enddo


        if (params%krylov_subspace_dimension == "dynamic" .or. M_iter == M_max) then
            ! if this is the last iteration, we compute the H matrix and the matrix exponential
            ! and use it to get an error estimate. If the error seems okay, we are done and can
            ! compute the new time step, otherwise we increase M by one and retry.
            ! if we use the dynamic method, we evaluate the error after every iteration to see if
            ! we're good to go.
            h_klein    = H(M_iter+1,M_iter)

            ! create a copy of the H matrix with the right dimension
            H_tmp      = 0.0_rk
            H_tmp(1:M_iter, 1:M_iter) = H(1:M_iter, 1:M_iter)
            H_tmp(M_iter+1,M_iter)    = 0.0_rk
            H_tmp(1,M_iter+1)         = 1.0_rk
            H_tmp(M_iter+1,M_iter+2)  = 1.0_rk

            ! compute matrix exponential
            t0 = MPI_wtime()
            phiMat(1:M_iter+2, 1:M_iter+2) = expM_pade( dt*H_tmp(1:M_iter+2, 1:M_iter+2) )
            phiMat(M_iter+1, M_iter+1)     = h_klein*phiMat(M_iter, M_iter+2)
            call toc( params, "Krylov: matrix exponential", MPI_wtime()-t0)


            ! *** Error estimate ***!
            err = abs( beta*phiMat(M_iter+1,M_iter+1) )

            if (params%krylov_subspace_dimension == "dynamic" .and. M_iter == M_max .and. err > params%krylov_err_threshold ) then
                ! we are at the last krylov subspace M and cannot increase the number any more.
                ! But the error is still too large, hance we decrease the time step now.
                do while (err > params%krylov_err_threshold)
                    dt = 0.90_rk * dt

                    ! compute matrix exponential
                    phiMat = 0.0_rk
                    phiMat(1:M_iter+2, 1:M_iter+2) = expM_pade( dt*H_tmp(1:M_iter+2, 1:M_iter+2) )
                    phiMat(M_iter+1, M_iter+1)     = h_klein*phiMat(M_iter, M_iter+2)

                    ! *** Error estimate ***!
                    err = abs( beta*phiMat(M_iter+1,M_iter+1) )
                enddo
            endif

            if (err <= params%krylov_err_threshold .or. M_iter == M_max) then
                exit
            endif
        endif
    enddo
    !**************************************!
    !*** end of iterations             ****!
    !**************************************!

    ! compute final value of new state vector at new time level
    ! result will be in hvy_block again (inout)
    do iter = 1, M_iter+1
        do k = 1, hvy_n
            hvy_block(:,:,:,:,hvy_active(k)) = hvy_block(:,:,:,:,hvy_active(k)) &
            + beta * hvy_work(:,:,:,:,hvy_active(k),iter) * phiMat(iter,M_iter+1)
        enddo
    enddo

    if (params%rank==0) then
        open(14,file='krylov_err.t',status='unknown',position='append')
        write (14,'(3(g15.8,1x),2(i3,1x))') time, dt, err, M_iter, M_max
        close(14)
    endif

end subroutine krylov_time_stepper


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function expM_pade(H) result(E)
    implicit none
    real(kind=rk),dimension(:,:),intent(in)                 :: H
    real(kind=rk),dimension(size(H,1),size(H,2))            :: E

    real(kind=rk),dimension(size(H,1),size(H,2))            :: expH
    integer                                                 :: m,ideg,lwsp,iexph,ldh,i,j,ii,iflag,ns
    real(kind=rk),allocatable,dimension(:)                  :: wsp
    integer,allocatable,dimension(:)                        :: ipiv
    integer::k

    !*** Berechnung des Matrixexponenten **!
    m     = size(H,1)
    ideg  = 6
    lwsp  = 4*m*m+ideg+1
    iexph = 1
    ldh   = m
    allocate(ipiv(m))   ! integer
    allocate(wsp(lwsp)) ! rk

    !                  !1.0_rk!, da dt*H -> H uebergeben wird, sonst steht hier dt
    call  DGPADM(ideg,m,1.0_rk,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag )
    if(iflag.lt.0) stop "error in computing exp(t*H)"

    do j=1,m
        do i=1,m
            ii = (i+iexph-1) + (j-1)*m
            E(i,j) = wsp(ii)
        end do
    end do
    !**************************************!

end function expM_pade


!  *----------------------------------------------------------------------|
subroutine DGPADM( ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag )
    implicit none
    integer ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
    double precision t, H(ldh,m), wsp(lwsp)

    !    *-----Purpose----------------------------------------------------------|
    !    *
    !    *     Computes exp(t*H), the matrix exponential of a general matrix in
    !    *     full, using the irreducible rational Pade approximation to the
    !    *     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),
    !      *     combined with scaling-and-squaring.
    !      *
    !      *-----Arguments--------------------------------------------------------|
    !      *
    !      *     ideg      : (input) the degre of the diagonal Pade to be used.
    !      *                 a value of 6 is generally satisfactory.
    !      *
    !      *     m         : (input) order of H.
    !      *
    !*     H(ldh,m)  : (input) argument matrix.
    !*
    !*     t         : (input) time-scale (can be < 0).
    !*
    !*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
    !*
    !*     ipiv(m)   : (workspace)
    !*
    !*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
    !*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
    !*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*                 NOTE: if the routine was called with wsp(iptr),
    !*                       then exp(tH) will start at wsp(iptr+iexph-1).
    !*
    !*     ns        : (output) number of scaling-squaring used.
    !*
    !*     iflag     : (output) exit flag.
    !*                      0 - no problem
    !*                     <0 - problem
    !*
    !*----------------------------------------------------------------------|
    !*     Roger B. Sidje (rbs@maths.uq.edu.au)
    !*     EXPOKIT: Software Package for Computing Matrix Exponentials.
    !*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
    !*----------------------------------------------------------------------|
    !*
    integer mm,i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,iput,iget
    double precision hnorm,scale,scale2,cp,cq

    intrinsic INT,ABS,DBLE,LOG,MAX

    !*---  check restrictions on input parameters ...
    mm = m*m
    iflag = 0
    if ( ldh.lt.m ) iflag = -1
    if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
    if ( iflag.ne.0 ) stop 'bad sizes (in input of DGPADM)'
    !*
    !*---  initialise pointers ...
    !*
    icoef = 1
    ih2 = icoef + (ideg+1)
    ip  = ih2 + mm
    iq  = ip + mm
    ifree = iq + mm
    !*
    !*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2;
    !*     and set scale = t/2^ns ...
    !*
    do i = 1,m
        wsp(i) = 0.0d0
    enddo
    do j = 1,m
        do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
        enddo
    enddo
    hnorm = 0.0d0
    do i = 1,m
        hnorm = MAX( hnorm,wsp(i) )
    enddo
    hnorm = ABS( t*hnorm )

    ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
    scale = t / DBLE(2**ns)
    scale2 = scale*scale
    !*
    !*---  compute Pade coefficients ...
    !*
    i = ideg+1
    j = 2*ideg+1
    wsp(icoef) = 1.0d0
    do k = 1,ideg
        wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
    enddo
    !*
    !*---  H2 = scale2*H*H ...
    !*
    call DGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,0.0d0,wsp(ih2),m )
    !*
    !*---  initialize p (numerator) and q (denominator) ...
    !*
    cp = wsp(icoef+ideg-1)
    cq = wsp(icoef+ideg)
    do j = 1,m
        do i = 1,m
            wsp(ip + (j-1)*m + i-1) = 0.0d0
            wsp(iq + (j-1)*m + i-1) = 0.0d0
        enddo
        wsp(ip + (j-1)*(m+1)) = cp
        wsp(iq + (j-1)*(m+1)) = cq
    enddo
    !*
    !*---  Apply Horner rule ...
    !*
    iodd = 1
    k = ideg - 1
100  continue
    iused = iodd*iq + (1-iodd)*ip
    call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iused),m, &
    wsp(ih2),m, 0.0d0,wsp(ifree),m )
    do j = 1,m
        wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
    enddo
    ip = (1-iodd)*ifree + iodd*ip
    iq = iodd*ifree + (1-iodd)*iq
    ifree = iused
    iodd = 1-iodd
    k = k-1
    if ( k.gt.0 )  goto 100
    !*
    !*---  Obtain (+/-)(I + 2*(p\q)) ...
    !*
    if ( iodd .eq. 1 ) then
        call DGEMM( 'n','n',m,m,m, scale,wsp(iq),m, &
        H,ldh, 0.0d0,wsp(ifree),m )
        iq = ifree
    else
        call DGEMM( 'n','n',m,m,m, scale,wsp(ip),m, &
        H,ldh, 0.0d0,wsp(ifree),m )
        ip = ifree
    endif
    call DAXPY( mm, -1.0d0,wsp(ip),1, wsp(iq),1 )
    call DGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
    if ( iflag.ne.0 ) stop 'Problem in DGESV (within DGPADM)'
    call DSCAL( mm, 2.0d0, wsp(ip), 1 )
    do j = 1,m
        wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + 1.0d0
    enddo
    iput = ip
    if ( ns.eq.0 .and. iodd.eq.1 ) then
        call DSCAL( mm, -1.0d0, wsp(ip), 1 )
        goto 200
    endif
    !*
    !*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
    !*
    iodd = 1
    do k = 1,ns
        iget = iodd*ip + (1-iodd)*iq
        iput = (1-iodd)*ip + iodd*iq
        call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iget),m, wsp(iget),m, &
        0.0d0,wsp(iput),m )
        iodd = 1-iodd
    enddo
200  continue
    iexph = iput
END subroutine DGPADM
!*----------------------------------------------------------------------|


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine arnoldi(rhs, v, h, j, params, lgt_block, hvy_active, lgt_active, lgt_n, hvy_n )
!     use module_blas
!     implicit none
!     real(kind=rk),dimension(:,:,:,:,:)      :: rhs ! w
!     real(kind=rk),dimension(:,:,:,:,:)      :: v
!     real(kind=rk),dimension(:,:)            :: h
!     integer                                 :: j
!     !> user defined parameter structure
!     type (type_params), intent(in)      :: params
!     !> light data array
!     integer(kind=ik), intent(in)        :: lgt_block(:, :)
!     !> heavy data array - block data
!     real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
!     !> list of active blocks (heavy data)
!     integer(kind=ik), intent(in)        :: hvy_active(:)
!     !> list of active blocks (light data)
!     integer(kind=ik), intent(in)        :: lgt_active(:)
!     !> number of active blocks (heavy data)
!     integer(kind=ik), intent(in)        :: hvy_n
!     !> number of active blocks (light data)
!     integer(kind=ik), intent(in)        :: lgt_n
!
!     ! real(kind=rk),dimension(:,:,:,:)        :: w
!     ! real(kind=rk),dimension(:,:,:,:,:)      :: v
!     ! real(kind=rk),dimension(:,:)            :: h
!     ! integer                                 :: j
!
!     !!! Arnoldi-Zerlegung
!     !!! Parallelisiert
!     !!! modified Gram Smith
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     integer                                 :: i
!     real(kind=rk)                           :: beta,h_all
!     !!! j ist die Anzahl der bereits vorhandenen Basisvektoren v(:,j)
!     !!! der neue Verktor (Av(:,j)) kommt als  v(:,j+1) von Aussen
!
!     !!! Falls wir von 0 starten
!     if (j == 0) then
! ! beta = norm(w)**2
! ! call get_sum_all(beta)
! ! beta = sqrt(beta)
!         call wabbit_norm(params, lgt_block, rhs, hvy_active, lgt_active, lgt_n, hvy_n, beta)
!
!         if (beta .lt. epsilon(1.0_rk)) beta = 1.0_rk
!
!         v(:,:,:,:,1) = rhs / beta
!         j = 1
!
!         ! end of routine
!         return
!     end if
!
!     !!! Sonst:
!     do i=1,j
!         h_all = scalprod(w, v(:,:,:,:,i))
!         call get_sum_all(h_all)
!         h(i,j)=h_all
!     end do
!
!     do i=1,j
!         w=w-h(i,j)*v(:,:,:,:,i)
!     end do
!
!     v(:,:,:,:,j+1)=w
!     h_all = nrm2(v(:,:,:,:,j+1))**2
!     call get_sum_all(h_all)
!     h_all = sqrt(h_all)
!     h(j+1,j)= h_all
!     if(h(j+1,j) .eq. 0.0_rk) then
!         write(*,*)'h(',j+1,j,')  =0.0'
!         stop
!     else
!         v(:,:,:,:,j+1)=v(:,:,:,:,j+1)/h(j+1,j)
!     end if
!     j=j+1
!
! end subroutine arnoldi



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_sum_all(inout,COMM)
    implicit none
    real(kind=rk),intent(inout)       :: inout
    integer,optional                  :: COMM
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                     :: tmp
    integer :: mpierr

    if (present(COMM)) then
        call MPI_ALLREDUCE(inout,tmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM,mpierr)
    else
        call MPI_ALLREDUCE(inout,tmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierr)
    end if
    inout = tmp

end subroutine get_sum_all


subroutine wabbit_norm(params, hvy_block, hvy_active, hvy_n, norm)
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n
    !> this is the output of the function:
    real(kind=rk), intent(out)          :: norm

    integer :: k, mpierr, ix, iy, iz, Bs, g

    norm = 0.0_rk

    Bs = params%Bs
    g = params%n_ghosts

    ! loop over active blocks
    if (params%dim == 3) then
        ! 3D
        do k = 1, hvy_n
            do iz = g+1, Bs+g-1 ! Note: loops skip redundant points
            do iy = g+1, Bs+g-1
            do ix = g+1, Bs+g-1
                norm = norm + sum( hvy_block(ix,iy,iz,:,hvy_active(k))**2 )
            enddo
            enddo
            enddo
        enddo
    else
        ! 2D
        do k = 1, hvy_n
            do iy = g+1, Bs+g-1 ! Note: loops skip redundant points
            do ix = g+1, Bs+g-1
                norm = norm + sum( hvy_block(ix,iy,1,:,hvy_active(k))**2 )
            enddo
            enddo
        enddo
    endif

    call get_sum_all(norm, WABBIT_COMM)
    norm = sqrt(norm)

end subroutine wabbit_norm

subroutine scalarproduct(params, hvy_block1, hvy_block2, hvy_active, hvy_n, result)
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block1(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_block2(:, :, :, :, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n
    !> this is the output of the function:
    real(kind=rk), intent(out)          :: result

    integer :: k, mpierr, ix, iy, iz, Bs, g, ieqn

    result = 0.0_rk

    Bs = params%Bs
    g = params%n_ghosts

    ! loop over active blocks
    if (params%dim == 3) then
        ! 3D
        do k = 1, hvy_n
            do iz = g+1, Bs+g-1 ! Note: loops skip redundant points
            do iy = g+1, Bs+g-1
            do ix = g+1, Bs+g-1
            do ieqn = 1, params%n_eqn
                result = result + hvy_block1(ix,iy,iz,ieqn,hvy_active(k)) * hvy_block2(ix,iy,iz,ieqn,hvy_active(k))
            enddo
            enddo
            enddo
            enddo
        enddo
    else
        ! 2D
        do k = 1, hvy_n
            do iy = g+1, Bs+g-1 ! Note: loops skip redundant points
            do ix = g+1, Bs+g-1
            do ieqn = 1, params%n_eqn
                result = result + hvy_block1(ix,iy,1,ieqn,hvy_active(k)) * hvy_block2(ix,iy,1,ieqn,hvy_active(k))
            enddo
            enddo
            enddo
        enddo
    endif

    call get_sum_all(result, WABBIT_COMM)

end subroutine scalarproduct
