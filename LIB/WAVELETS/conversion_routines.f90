! ! Inflated Mallat ordering is a HACK. I was just easier to code.
! ! That does not mean its wrong, it isn't. But is uses extra memory (although
! ! that is negligible) and does unnecessary copy actions.
! subroutine spaghetti2inflatedMallat_block(params, u, wc, sc_only)
!     implicit none
!     type (type_params), intent(in) :: params
!     real(kind=rk), dimension(:,:,:,:), intent(inout) :: u
!     ! The WC array contains SC (scaling function coeffs) as well as all WC (wavelet coeffs)
!     ! Note: the precise naming of SC/WC is not really important. we just apply
!     ! the correct decomposition/reconstruction filters - thats it.
!     !
!     ! INDEX            2D     3D     LABEL
!     ! -----            --    ---     ---------------------------------
!     ! wc(:,:,:,:,1)    HH    HHH     sc scaling function coeffs
!     ! wc(:,:,:,:,2)    GH    GHH     wcy wavelet coeffs
!     ! wc(:,:,:,:,3)    HG    HGH     wcx wavelet coeffs
!     ! wc(:,:,:,:,4)    GG    GGH     wcxy wavelet coeffs
!     ! wc(:,:,:,:,5)          HHG     wcz wavelet coeffs
!     ! wc(:,:,:,:,6)          GHG     wcyz wavelet coeffs
!     ! wc(:,:,:,:,7)          HGG     wcxz wavelet coeffs
!     ! wc(:,:,:,:,8)          GGG     wcxyz wavelet coeffs
!     !
!     real(kind=rk), dimension(:,:,:,:,:), intent(inout) :: wc
!     !> Option to only convert SC
!     logical, optional, intent(in)  :: sc_only
!     integer(kind=ik) :: nx, ny, nz, nc, i_s
!     logical sc_only_set

!     sc_only_set = .false.
!     if (present(sc_only)) sc_only_set = sc_only

!     nx = size(u, 1)
!     ny = size(u, 2)
!     nz = size(u, 3)
!     nc = size(u, 4)

! #ifdef DEV
!     if (.not.areArraysSameSize(u, wc(:,:,:,:,1))) then
!         call abort(27222119, "Allocated arrays are not compatible?! Time for a drink.")
!     endif
! #endif

!     wc = 0.0_rk

!     ! note that what we call "Mallat ordering" here is in fact the "inflated" Mallat
!     ! in the sense that Nx*Ny data gives 4/8 * Nx*Ny decomposition.

!     ! compute index shift for odd g, ignoring outmost elements
!     i_s = 1 + modulo(g, 2)  ! 1 for even, 2 for odd

!     ! copy from Spaghetti to inflated Mallat ordering
!     if (params%dim == 2) then
!         wc(     i_s:nx  :2, i_s:ny  :2, :, :, 1) = u(   i_s:nx:2,   i_s:ny:2, :, :)
!         if (.not. sc_only_set) then
!             wc( i_s:nx-1:2, i_s:ny  :2, :, :, 2) = u( 1+i_s:nx:2,   i_s:ny:2, :, :)
!             wc( i_s:nx  :2, i_s:ny-1:2, :, :, 3) = u(   i_s:nx:2, 1+i_s:ny:2, :, :)
!             wc( i_s:nx-1:2, i_s:ny-1:2, :, :, 4) = u( 1+i_s:nx:2, 1+i_s:ny:2, :, :)
!         endif
!     else
!         wc(     i_s:nx  :2, i_s:ny  :2, i_s:nz  :2, :, 1) = u(   i_s:nx:2,   i_s:ny:2,   i_s:nz:2, :)
!         if (.not. sc_only_set) then
!             wc( i_s:nx-1:2, i_s:ny  :2, i_s:nz  :2, :, 2) = u( 1+i_s:nx:2,   i_s:ny:2,   i_s:nz:2, :)
!             wc( i_s:nx  :2, i_s:ny-1:2, i_s:nz  :2, :, 3) = u(   i_s:nx:2, 1+i_s:ny:2,   i_s:nz:2, :)
!             wc( i_s:nx-1:2, i_s:ny-1:2, i_s:nz  :2, :, 4) = u( 1+i_s:nx:2, 1+i_s:ny:2,   i_s:nz:2, :)

!             wc( i_s:nx  :2, i_s:ny  :2, i_s:nz-1:2, :, 5) = u(   i_s:nx:2,   i_s:ny:2, 1+i_s:nz:2, :)
!             wc( i_s:nx-1:2, i_s:ny  :2, i_s:nz-1:2, :, 6) = u( 1+i_s:nx:2,   i_s:ny:2, 1+i_s:nz:2, :)
!             wc( i_s:nx  :2, i_s:ny-1:2, i_s:nz-1:2, :, 7) = u(   i_s:nx:2, 1+i_s:ny:2, 1+i_s:nz:2, :)
!             wc( i_s:nx-1:2, i_s:ny-1:2, i_s:nz-1:2, :, 8) = u( 1+i_s:nx:2, 1+i_s:ny:2, 1+i_s:nz:2, :)
!         endif
!     endif
! end subroutine


! ! Inflated Mallat ordering is a HACK. I was just easier to code.
! ! That does not mean its wrong, it isn't. But is uses extra memory (although
! ! that is negligible) and does unnecessary copy actions.
! subroutine inflatedMallat2spaghetti_block(params, wc, u, sc_only)
!     implicit none
!     type (type_params), intent(in) :: params
!     real(kind=rk), dimension(:,:,:,:), intent(inout) :: u
!     ! The WC array contains SC (scaling function coeffs) as well as all WC (wavelet coeffs)
!     ! Note: the precise naming of SC/WC is not really important. we just apply
!     ! the correct decomposition/reconstruction filters - thats it.
!     !
!     ! INDEX            2D     3D     LABEL
!     ! -----            --    ---     ---------------------------------
!     ! wc(:,:,:,:,1)    HH    HHH     sc scaling function coeffs
!     ! wc(:,:,:,:,2)    GH    GHH     wcy wavelet coeffs
!     ! wc(:,:,:,:,3)    HG    HGH     wcx wavelet coeffs
!     ! wc(:,:,:,:,4)    GG    GGH     wcxy wavelet coeffs
!     ! wc(:,:,:,:,5)          HHG     wcz wavelet coeffs
!     ! wc(:,:,:,:,6)          GHG     wcyz wavelet coeffs
!     ! wc(:,:,:,:,7)          HGG     wcxz wavelet coeffs
!     ! wc(:,:,:,:,8)          GGG     wcxyz wavelet coeffs
!     !
!     real(kind=rk), dimension(:,:,:,:,:), intent(inout) :: wc
!     !> Option to only convert SC
!     logical, optional, intent(in)  :: sc_only

!     integer(kind=ik) :: nx, ny, nz, nc, i_s
!     logical sc_only_set

!     sc_only_set = .false.
!     if (present(sc_only)) sc_only_set = sc_only

!     nx = size(u, 1)
!     ny = size(u, 2)
!     nz = size(u, 3)
!     nc = size(u, 4)

! #ifdef DEV
!     if (.not.areArraysSameSize(u, wc(:,:,:,:,1))) then
!         call abort(27222119, "Allocated arrays are not compatible?! Time for a drink. You look handsome today.")
!     endif
! #endif

!     ! note that what we call "Mallat ordering" here is in fact the "inflated" Mallat
!     ! in the sense that Nx*Ny data gives 4/8 * Nx*Ny decomposition.

!     ! compute index shift for odd g, ignoring outmost elements
!     i_s = 1 + modulo(g, 2)  ! 1 for even, 2 for odd

!     ! copy from inflated Mallat to Spaghetti ordering
!     if (params%dim == 2) then
!         u(  i_s:nx:2,   i_s:ny:2, :, :) = wc( i_s:nx  :2, i_s:ny  :2, :, :, 1)
!         if (.not. sc_only_set) then
!             u( 1+i_s:nx:2,   i_s:ny:2, :, :) = wc( i_s:nx-1:2, i_s:ny  :2, :, :, 2)
!             u(   i_s:nx:2, 1+i_s:ny:2, :, :) = wc( i_s:nx  :2, i_s:ny-1:2, :, :, 3)
!             u( 1+i_s:nx:2, 1+i_s:ny:2, :, :) = wc( i_s:nx-1:2, i_s:ny-1:2, :, :, 4)
!         endif
!     else
!         u(  i_s:nx:2,   i_s:ny:2,   i_s:nz:2, :) = wc( i_s:nx  :2, i_s:ny  :2, i_s:nz  :2, :, 1)
!         if (.not. sc_only_set) then
!             u( 1+i_s:nx:2,   i_s:ny:2,   i_s:nz:2, :) = wc( i_s:nx-1:2, i_s:ny  :2, i_s:nz  :2, :, 2)
!             u(   i_s:nx:2, 1+i_s:ny:2,   i_s:nz:2, :) = wc( i_s:nx  :2, i_s:ny-1:2, i_s:nz  :2, :, 3)
!             u( 1+i_s:nx:2, 1+i_s:ny:2,   i_s:nz:2, :) = wc( i_s:nx-1:2, i_s:ny-1:2, i_s:nz  :2, :, 4)

!             u(   i_s:nx:2,   i_s:ny:2, 1+i_s:nz:2, :) = wc( i_s:nx  :2, i_s:ny  :2, i_s:nz-1:2, :, 5)
!             u( 1+i_s:nx:2,   i_s:ny:2, 1+i_s:nz:2, :) = wc( i_s:nx-1:2, i_s:ny  :2, i_s:nz-1:2, :, 6)
!             u(   i_s:nx:2, 1+i_s:ny:2, 1+i_s:nz:2, :) = wc( i_s:nx  :2, i_s:ny-1:2, i_s:nz-1:2, :, 7)
!             u( 1+i_s:nx:2, 1+i_s:ny:2, 1+i_s:nz:2, :) = wc( i_s:nx-1:2, i_s:ny-1:2, i_s:nz-1:2, :, 8)
!         endif
!     endif

! end subroutine


subroutine spaghetti2Mallat_block(params, u, wc)
    implicit none
    type (type_params), intent(in) :: params
    real(kind=rk), dimension(:,:,:,:), intent(inout) :: u
    ! The WC array contains SC (scaling function coeffs) as well as all WC (wavelet coeffs)
    ! Note: the precise naming of SC/WC is not really important. we just apply
    ! the correct decomposition/reconstruction filters - thats it.
    !
    ! INDEX                                        2D     3D     LABEL
    ! -----------------------------------------    --    ---     ---------------------------------
    ! wc(     1:nx/2,     1:ny/2,     1:nz/2,:)    HH    HHH     sc scaling function coeffs
    ! wc(nx/2+1:nx  ,     1:ny/2,     1:nz/2,:)    GH    GHH     wcy wavelet coeffs
    ! wc(     1:nx/2,ny/2+1:ny  ,     1:nz/2,:)    HG    HGH     wcx wavelet coeffs
    ! wc(nx/2+1:nx  ,ny/2+1:ny  ,     1:nz/2,:)    GG    GGH     wcxy wavelet coeffs
    ! wc(     1:nx/2,     1:ny/2,nz/2+1:nz  ,:)          HHG     wcz wavelet coeffs
    ! wc(nx/2+1:nx  ,     1:ny/2,nz/2+1:nz  ,:)          GHG     wcyz wavelet coeffs
    ! wc(     1:nx/2,ny/2+1:ny  ,nz/2+1:nz  ,:)          HGG     wcxz wavelet coeffs
    ! wc(nx/2+1:nx  ,ny/2+1:ny  ,nz/2+1:nz  ,:)          GGG     wcxyz wavelet coeffs
    !
    ! For different BS or g the actual points could look differently, let's have a look at 1D:
    !                      |---- spaghetti ordering --------------|    |----- Mallat odering ---------------------------|
    ! BS even, g even:     GS GW IS IW IS IW GS GW    (2G,4I,2G)    ->  GS IS IS GS GW IW IW GW    (1G,2I,1G | 1G,2I,1G)
    ! BS even, g odd :        GW IS IW IS IW GS       (1G,4I,1G)    ->     IS IS GS GW IW IW       (   2I,1G | 1G,2I   )
    ! We therefore have: G_SC_L = g//2, G_WC_L = g//2 + 1, Start_WC = (BS+1)//2 + (g//2)*2
    ! In the above:
    !     G* : ghost point, I* interior point
    !     *S : scaling function coeff *W: wavelet coeff  ==> GS: ghost nodes, scaling function coeff
    !
    

    real(kind=rk), dimension(:,:,:,:), intent(inout) :: wc
    integer(kind=ik) :: nx, ny, nz, nc, i_s, sc_e(1:3)

    nx = size(u, 1)
    ny = size(u, 2)
    nz = size(u, 3)
    nc = size(u, 4)
    ! how many points are in the scaling function part, bear in mind we do integer division here
    do i_s = 1,params%dim
        sc_e(i_s) = params%Bs(i_s)/2 + params%g
    end do

#ifdef DEV
    if (.not.areArraysSameSize(u, wc)) then
        call abort(27222119, "Allocated arrays are not compatible?! Time for a drink.")
    endif
#endif

    ! compute offset to first SC
    i_s = modulo(params%g, 2)  ! 0 for even, 1 for odd

    ! all values will be set and copied, so we do not need to wipe it before
    ! wc = 0.0_rk

    if (params%dim == 2) then
        wc(         1:sc_e(1),         1:sc_e(2), :, :) = u( 1+i_s:nx:2, 1+i_s:ny:2, :, :)
        wc( sc_e(1)+1:nx     ,         1:sc_e(2), :, :) = u( 2-i_s:nx:2, 1+i_s:ny:2, :, :)
        wc(         1:sc_e(1), sc_e(2)+1:ny     , :, :) = u( 1+i_s:nx:2, 2-i_s:ny:2, :, :)
        wc( sc_e(1)+1:nx     , sc_e(2)+1:ny     , :, :) = u( 2-i_s:nx:2, 2-i_s:ny:2, :, :)
    else
        wc(         1:sc_e(1),         1:sc_e(2),         1:sc_e(3), :) = u( 1+i_s:nx:2, 1+i_s:ny:2, 1+i_s:nz:2, :)
        wc( sc_e(1)+1:nx     ,         1:sc_e(2),         1:sc_e(3), :) = u( 2-i_s:nx:2, 1+i_s:ny:2, 1+i_s:nz:2, :)
        wc(         1:sc_e(1), sc_e(2)+1:ny     ,         1:sc_e(3), :) = u( 1+i_s:nx:2, 2-i_s:ny:2, 1+i_s:nz:2, :)
        wc( sc_e(1)+1:nx     , sc_e(2)+1:ny     ,         1:sc_e(3), :) = u( 2-i_s:nx:2, 2-i_s:ny:2, 1+i_s:nz:2, :)

        wc(         1:sc_e(1),         1:sc_e(2), sc_e(3)+1:nz     , :) = u( 1+i_s:nx:2, 1+i_s:ny:2, 2-i_s:nz:2, :)
        wc( sc_e(1)+1:nx     ,         1:sc_e(2), sc_e(3)+1:nz     , :) = u( 2-i_s:nx:2, 1+i_s:ny:2, 2-i_s:nz:2, :)
        wc(         1:sc_e(1), sc_e(2)+1:ny     , sc_e(3)+1:nz     , :) = u( 1+i_s:nx:2, 2-i_s:ny:2, 2-i_s:nz:2, :)
        wc( sc_e(1)+1:nx     , sc_e(2)+1:ny     , sc_e(3)+1:nz     , :) = u( 2-i_s:nx:2, 2-i_s:ny:2, 2-i_s:nz:2, :)
    endif
end subroutine



subroutine Mallat2Spaghetti_block(params, wc, u)
    implicit none
    type (type_params), intent(in) :: params
    real(kind=rk), dimension(:,:,:,:), intent(inout) :: u
    ! The WC array contains SC (scaling function coeffs) as well as all WC (wavelet coeffs)
    ! Note: the precise naming of SC/WC is not really important. we just apply
    ! the correct decomposition/reconstruction filters - thats it.
    !
    ! INDEX                                        2D     3D     LABEL
    ! -----------------------------------------    --    ---     ---------------------------------
    ! wc(     1:nx/2,     1:ny/2,     1:nz/2,:)    HH    HHH     sc scaling function coeffs
    ! wc(nx/2+1:nx  ,     1:ny/2,     1:nz/2,:)    GH    GHH     wcy wavelet coeffs
    ! wc(     1:nx/2,ny/2+1:ny  ,     1:nz/2,:)    HG    HGH     wcx wavelet coeffs
    ! wc(nx/2+1:nx  ,ny/2+1:ny  ,     1:nz/2,:)    GG    GGH     wcxy wavelet coeffs
    ! wc(     1:nx/2,     1:ny/2,nz/2+1:nz  ,:)          HHG     wcz wavelet coeffs
    ! wc(nx/2+1:nx  ,     1:ny/2,nz/2+1:nz  ,:)          GHG     wcyz wavelet coeffs
    ! wc(     1:nx/2,ny/2+1:ny  ,nz/2+1:nz  ,:)          HGG     wcxz wavelet coeffs
    ! wc(nx/2+1:nx  ,ny/2+1:ny  ,nz/2+1:nz  ,:)          GGG     wcxyz wavelet coeffs
    !
    ! For different BS or g the actual points could look differently, let's have a look at 1D:
    ! BS even, g even:     GS IS IS GS GW IW IW GW    (1G,2I,1G | 1G,2I,1G)   ->  GS GW IS IW IS IW GS GW    (2G,4I,2G)
    ! BS even, g odd :        IS IS GS GW IW IW       (   2I,1G | 1G,2I   )   ->     GW IS IW IS IW GS       (1G,4I,1G)
    ! We therefore have: G_SC_L = g//2, G_WC_L = g//2 + 1, Start_WC = (BS+1)//2 + (g//2)*2
    real(kind=rk), dimension(:,:,:,:), intent(inout) :: wc
    integer(kind=ik) :: nx, ny, nz, nc, i_s, sc_e(1:3)

    nx = size(u, 1)
    ny = size(u, 2)
    nz = size(u, 3)
    nc = size(u, 4)
    ! how many points are in the scaling function part
    do i_s = 1,params%dim
        sc_e(i_s) = params%Bs(i_s)/2 + params%g
    end do
    
#ifdef DEV
    if (.not.areArraysSameSize(u, wc)) then
        call abort(27222119, "Allocated arrays are not compatible?! Time for a drink.")
    endif
#endif

    ! compute offset to first SC
    i_s = modulo(params%g, 2)  ! 0 for even, 1 for odd

    ! all values will be set and copied, so we do not need to wipe it before
    ! u = 0.0_rk

    ! copy from Mallat to Spaghetti ordering
    if (params%dim == 2) then
        u( 1+i_s:nx:2, 1+i_s:ny:2, :, :) = wc(         1:sc_e(1),         1:sc_e(2), :, :)
        u( 2-i_s:nx:2, 1+i_s:ny:2, :, :) = wc( sc_e(1)+1:nx     ,         1:sc_e(2), :, :)
        u( 1+i_s:nx:2, 2-i_s:ny:2, :, :) = wc(         1:sc_e(1), sc_e(2)+1:ny     , :, :)
        u( 2-i_s:nx:2, 2-i_s:ny:2, :, :) = wc( sc_e(1)+1:nx     , sc_e(2)+1:ny     , :, :)
    else
        u( 1+i_s:nx:2, 1+i_s:ny:2, 1+i_s:nz:2, :) = wc(         1:sc_e(1),         1:sc_e(2),         1:sc_e(3), :)
        u( 2-i_s:nx:2, 1+i_s:ny:2, 1+i_s:nz:2, :) = wc( sc_e(1)+1:nx     ,         1:sc_e(2),         1:sc_e(3), :)
        u( 1+i_s:nx:2, 2-i_s:ny:2, 1+i_s:nz:2, :) = wc(         1:sc_e(1), sc_e(2)+1:ny     ,         1:sc_e(3), :)
        u( 2-i_s:nx:2, 2-i_s:ny:2, 1+i_s:nz:2, :) = wc( sc_e(1)+1:nx     , sc_e(2)+1:ny     ,         1:sc_e(3), :)

        u( 1+i_s:nx:2, 1+i_s:ny:2, 2-i_s:nz:2, :) = wc(         1:sc_e(1),         1:sc_e(2), sc_e(3)+1:nz     , :)
        u( 2-i_s:nx:2, 1+i_s:ny:2, 2-i_s:nz:2, :) = wc( sc_e(1)+1:nx     ,         1:sc_e(2), sc_e(3)+1:nz     , :)
        u( 1+i_s:nx:2, 2-i_s:ny:2, 2-i_s:nz:2, :) = wc(         1:sc_e(1), sc_e(2)+1:ny     , sc_e(3)+1:nz     , :)
        u( 2-i_s:nx:2, 2-i_s:ny:2, 2-i_s:nz:2, :) = wc( sc_e(1)+1:nx     , sc_e(2)+1:ny     , sc_e(3)+1:nz     , :)
    endif
end subroutine



! ! Inflated Mallat ordering is a HACK. I was just easier to code.
! ! That does not mean its wrong, it isn't. But is uses extra memory (although
! ! that is negligible) and does unnecessary copy actions.
! subroutine Mallat2inflatedMallat_block(params, wc, u)
!     implicit none
!     type (type_params), intent(in) :: params
!     ! The WC and u array contains SC (scaling function coeffs) as well as all WC (wavelet coeffs)
!     ! Note: the precise naming of SC/WC is not really important. we just apply
!     ! the correct decomposition/reconstruction filters - thats it.
!     !
!     ! INDEX                                        2D     3D     LABEL
!     ! -----------------------------------------    --    ---     ---------------------------------
!     ! wc(     1:nx/2,     1:ny/2,     1:nz/2,:)    HH    HHH     sc scaling function coeffs
!     ! wc(nx/2+1:nx  ,     1:ny/2,     1:nz/2,:)    GH    GHH     wcy wavelet coeffs
!     ! wc(     1:nx/2,ny/2+1:ny  ,     1:nz/2,:)    HG    HGH     wcx wavelet coeffs
!     ! wc(nx/2+1:nx  ,ny/2+1:ny  ,     1:nz/2,:)    GG    GGH     wcxy wavelet coeffs
!     ! wc(     1:nx/2,     1:ny/2,nz/2+1:nz  ,:)          HHG     wcz wavelet coeffs
!     ! wc(nx/2+1:nx  ,     1:ny/2,nz/2+1:nz  ,:)          GHG     wcyz wavelet coeffs
!     ! wc(     1:nx/2,ny/2+1:ny  ,nz/2+1:nz  ,:)          HGG     wcxz wavelet coeffs
!     ! wc(nx/2+1:nx  ,ny/2+1:ny  ,nz/2+1:nz  ,:)          GGG     wcxyz wavelet coeffs
!     real(kind=rk), dimension(:,:,:,:), intent(inout) :: wc
!     ! INDEX            2D     3D     LABEL
!     ! -----            --    ---     ---------------------------------
!     ! u(:,:,:,:,1)    HH    HHH     sc scaling function coeffs
!     ! u(:,:,:,:,2)    GH    GHH     wcy wavelet coeffs
!     ! u(:,:,:,:,3)    HG    HGH     wcx wavelet coeffs
!     ! u(:,:,:,:,4)    GG    GGH     wcxy wavelet coeffs
!     ! u(:,:,:,:,5)          HHG     wcz wavelet coeffs
!     ! u(:,:,:,:,6)          GHG     wcyz wavelet coeffs
!     ! u(:,:,:,:,7)          HGG     wcxz wavelet coeffs
!     ! u(:,:,:,:,8)          GGG     wcxyz wavelet coeffs
!     real(kind=rk), dimension(:,:,:,:,:), intent(inout) :: u

!     integer(kind=ik) :: nx, ny, nz, nc, g, Bs(1:3)

!     nx = size(u, 1)
!     ny = size(u, 2)
!     nz = size(u, 3)
!     nc = size(u, 4)
!     g  = params%g
!     Bs = params%bs

! #ifdef DEV
!     if (.not.areArraysSameSize(wc, u(:,:,:,:,1))) then
!         call abort(27222119, "Allocated arrays are not compatible?! Time for a drink.")
!     endif
! #endif

!     u = 0.0_rk

!     ! copy from Mallat to inflated Mallat ordering
!     if (modulo(g, 2) == 0) then
!         ! even g
!         if (params%dim == 2) then
!             u( 1:nx:2, 1:ny:2, :, :, 1) = wc(      1:nx/2,      1:ny/2, :, :)
!             u( 1:nx:2, 1:ny:2, :, :, 2) = wc(      1:nx/2, ny/2+1:ny  , :, :)
!             u( 1:nx:2, 1:ny:2, :, :, 3) = wc( nx/2+1:nx  ,      1:ny/2, :, :)
!             u( 1:nx:2, 1:ny:2, :, :, 4) = wc( nx/2+1:nx  , ny/2+1:ny  , :, :)
!         else
!             u( 1:nx:2, 1:ny:2, 1:nz:2, :, 1) = wc(      1:nx/2,      1:ny/2,      1:nz/2, :)
!             u( 1:nx:2, 1:ny:2, 1:nz:2, :, 2) = wc(      1:nx/2, ny/2+1:ny  ,      1:nz/2, :)
!             u( 1:nx:2, 1:ny:2, 1:nz:2, :, 3) = wc( nx/2+1:nx  ,      1:ny/2,      1:nz/2, :)
!             u( 1:nx:2, 1:ny:2, 1:nz:2, :, 4) = wc( nx/2+1:nx  , ny/2+1:ny  ,      1:nz/2, :)

!             u( 1:nx:2, 1:ny:2, 1:nz:2, :, 5) = wc(      1:nx/2,      1:ny/2, nz/2+1:nz  , :)
!             u( 1:nx:2, 1:ny:2, 1:nz:2, :, 6) = wc(      1:nx/2, ny/2+1:ny  , nz/2+1:nz  , :)
!             u( 1:nx:2, 1:ny:2, 1:nz:2, :, 7) = wc( nx/2+1:nx  ,      1:ny/2, nz/2+1:nz  , :)
!             u( 1:nx:2, 1:ny:2, 1:nz:2, :, 8) = wc( nx/2+1:nx  , ny/2+1:ny  , nz/2+1:nz  , :)
!         endif
!     else
!         ! odd g
!         if (params%dim == 2) then
!             u( 2:nx-1:2, 2:ny-1:2, :, :, 1) = wc(      2:nx/2,      2:ny/2, :, :)
!             u( 2:nx-1:2, 2:ny-1:2, :, :, 2) = wc(      2:nx/2, ny/2+1:ny-1, :, :)
!             u( 2:nx-1:2, 2:ny-1:2, :, :, 3) = wc( nx/2+1:nx-1,      2:ny/2, :, :)
!             u( 2:nx-1:2, 2:ny-1:2, :, :, 4) = wc( nx/2+1:nx-1, ny/2+1:ny-1, :, :)
!         else
!             u( 2:nx-1:2, 2:ny-1:2, 2:nz-1:2, :, 1) = wc(      2:nx/2,      2:ny/2,      2:nz/2, :)
!             u( 2:nx-1:2, 2:ny-1:2, 2:nz-1:2, :, 2) = wc(      2:nx/2, ny/2+1:ny-1,      2:nz/2, :)
!             u( 2:nx-1:2, 2:ny-1:2, 2:nz-1:2, :, 3) = wc( nx/2+1:nx-1,      2:ny/2,      2:nz/2, :)
!             u( 2:nx-1:2, 2:ny-1:2, 2:nz-1:2, :, 4) = wc( nx/2+1:nx-1, ny/2+1:ny-1,      2:nz/2, :)

!             u( 2:nx-1:2, 2:ny-1:2, 2:nz-1:2, :, 5) = wc(      2:nx/2,      2:ny/2, nz/2+1:nz-1, :)
!             u( 2:nx-1:2, 2:ny-1:2, 2:nz-1:2, :, 6) = wc(      2:nx/2, ny/2+1:ny-1, nz/2+1:nz-1, :)
!             u( 2:nx-1:2, 2:ny-1:2, 2:nz-1:2, :, 7) = wc( nx/2+1:nx-1,      2:ny/2, nz/2+1:nz-1, :)
!             u( 2:nx-1:2, 2:ny-1:2, 2:nz-1:2, :, 8) = wc( nx/2+1:nx-1, ny/2+1:ny-1, nz/2+1:nz-1, :)
!         endif
!     endif
! end subroutine



! ! Inflated Mallat ordering is a HACK. I was just easier to code.
! ! That does not mean its wrong, it isn't. But is uses extra memory (although
! ! that is negligible) and does unnecessary copy actions.
! subroutine inflatedMallat2Mallat_block(params, u, wc, sc_only)
!     implicit none
!     type (type_params), intent(in) :: params
!     ! The WC and u array contains SC (scaling function coeffs) as well as all WC (wavelet coeffs)
!     ! Note: the precise naming of SC/WC is not really important. we just apply
!     ! the correct decomposition/reconstruction filters - thats it.
!     !
!     ! INDEX                                        2D     3D     LABEL
!     ! -----------------------------------------    --    ---     ---------------------------------
!     ! wc(     1:nx/2,     1:ny/2,     1:nz/2,:)    HH    HHH     sc scaling function coeffs
!     ! wc(nx/2+1:nx  ,     1:ny/2,     1:nz/2,:)    GH    GHH     wcy wavelet coeffs
!     ! wc(     1:nx/2,ny/2+1:ny  ,     1:nz/2,:)    HG    HGH     wcx wavelet coeffs
!     ! wc(nx/2+1:nx  ,ny/2+1:ny  ,     1:nz/2,:)    GG    GGH     wcxy wavelet coeffs
!     ! wc(     1:nx/2,     1:ny/2,nz/2+1:nz  ,:)          HHG     wcz wavelet coeffs
!     ! wc(nx/2+1:nx  ,     1:ny/2,nz/2+1:nz  ,:)          GHG     wcyz wavelet coeffs
!     ! wc(     1:nx/2,ny/2+1:ny  ,nz/2+1:nz  ,:)          HGG     wcxz wavelet coeffs
!     ! wc(nx/2+1:nx  ,ny/2+1:ny  ,nz/2+1:nz  ,:)          GGG     wcxyz wavelet coeffs
!     real(kind=rk), dimension(:,:,:,:), intent(inout) :: wc
!     ! INDEX            2D     3D     LABEL
!     ! -----            --    ---     ---------------------------------
!     ! u(:,:,:,:,1)    HH    HHH     sc scaling function coeffs
!     ! u(:,:,:,:,2)    GH    GHH     wcy wavelet coeffs
!     ! u(:,:,:,:,3)    HG    HGH     wcx wavelet coeffs
!     ! u(:,:,:,:,4)    GG    GGH     wcxy wavelet coeffs
!     ! u(:,:,:,:,5)          HHG     wcz wavelet coeffs
!     ! u(:,:,:,:,6)          GHG     wcyz wavelet coeffs
!     ! u(:,:,:,:,7)          HGG     wcxz wavelet coeffs
!     ! u(:,:,:,:,8)          GGG     wcxyz wavelet coeffs
!     real(kind=rk), dimension(:,:,:,:,:), intent(inout) :: u
!     !> Option to only convert SC
!     logical, optional, intent(in)  :: sc_only

!     integer(kind=ik) :: nx, ny, nz, nc, g, Bs(1:3), i_s, i_e
!     logical sc_only_set

!     sc_only_set = .false.
!     if (present(sc_only)) sc_only_set = sc_only

!     nx = size(u, 1)
!     ny = size(u, 2)
!     nz = size(u, 3)
!     nc = size(u, 4)
!     g  = params%g
!     Bs = params%bs

! #ifdef DEV
!     if (.not.areArraysSameSize(wc, u(:,:,:,:,1))) then
!         call abort(27222119, "Allocated arrays are not compatible?! Time for a drink.")
!     endif
! #endif

!     wc = 0.0_rk

!     ! compute index shift for odd g, ignoring outmost elements
!     i_s = 1 + modulo(g, 2)  ! 1 for even, 2 for odd
!     i_e = 0 + modulo(g, 2)  ! 0 for even, 1 for odd

!     ! copy from inflated Mallat to Mallat ordering
!     if (params%dim == 2) then
!         wc(        i_s:nx/2  ,    i_s:ny/2  , :, :) = u( i_s:nx-i_e:2, i_s:ny-i_e:2, :, :, 1)
!         if (.not. sc_only_set) then
!             wc(    i_s:nx/2  , ny/2+1:ny-i_e, :, :) = u( i_s:nx-i_e:2, i_s:ny-i_e:2, :, :, 2)
!             wc( nx/2+1:nx-i_e,    i_s:ny/2  , :, :) = u( i_s:nx-i_e:2, i_s:ny-i_e:2, :, :, 3)
!             wc( nx/2+1:nx-i_e, ny/2+1:ny-i_e, :, :) = u( i_s:nx-i_e:2, i_s:ny-i_e:2, :, :, 4)
!         endif
!     else
!         wc(        i_s:nx/2  ,    i_s:ny/2  ,    i_s:nz/2  , :) = u( i_s:nx-i_e:2, i_s:ny-i_e:2, i_s:nz-i_e:2, :, 1)
!         if (.not. sc_only_set) then
!             wc(    i_s:nx/2  , ny/2+1:ny-i_e,    i_s:nz/2  , :) = u( i_s:nx-i_e:2, i_s:ny-i_e:2, i_s:nz-i_e:2, :, 2)
!             wc( nx/2+1:nx-i_e,    i_s:ny/2  ,    i_s:nz/2  , :) = u( i_s:nx-i_e:2, i_s:ny-i_e:2, i_s:nz-i_e:2, :, 3)
!             wc( nx/2+1:nx-i_e, ny/2+1:ny-i_e,    i_s:nz/2  , :) = u( i_s:nx-i_e:2, i_s:ny-i_e:2, i_s:nz-i_e:2, :, 4)
    
!             wc(    i_s:nx/2  ,    i_s:ny/2  , nz/2+1:nz-i_e, :) = u( i_s:nx-i_e:2, i_s:ny-i_e:2, i_s:nz-i_e:2, :, 5)
!             wc(    i_s:nx/2  , ny/2+1:ny-i_e, nz/2+1:nz-i_e, :) = u( i_s:nx-i_e:2, i_s:ny-i_e:2, i_s:nz-i_e:2, :, 6)
!             wc( nx/2+1:nx-i_e,    i_s:ny/2  , nz/2+1:nz-i_e, :) = u( i_s:nx-i_e:2, i_s:ny-i_e:2, i_s:nz-i_e:2, :, 7)
!             wc( nx/2+1:nx-i_e, ny/2+1:ny-i_e, nz/2+1:nz-i_e, :) = u( i_s:nx-i_e:2, i_s:ny-i_e:2, i_s:nz-i_e:2, :, 8)
!         endif
!     endif
! end subroutine
