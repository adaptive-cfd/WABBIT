
  !-----------------------------------------------------------------------------
  ! save data. Since you might want to save derived data, such as the vorticity,
  ! the divergence etc., which are not in your state vector, this routine has to
  ! copy and compute what you want to save to the work array.
  !
  ! In the main code, save_fields than saves the first N_fields_saved components of the
  ! work array to file.
  !
  ! NOTE that as we have way more work arrays than actual state variables (typically
  ! for a RK4 that would be >= 4*dim), you can compute a lot of stuff, if you want to.
  !-----------------------------------------------------------------------------
  subroutine PREPARE_SAVE_DATA_NStokes( time, u, g, x0, dx, work )
    use module_helpers , only: choose, list_contains_name
    implicit none
    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time
    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(in) :: u(1:,1:,1:,1:)
    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g
    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)
    ! output in work array.
    real(kind=rk), intent(inout) :: work(1:,1:,1:,1:)
    ! output in work array.
    real(kind=rk), allocatable,save :: tmp_u(:,:,:,:), vort(:,:,:,:), sigma(:,:,:,:)
    ! local variables
    integer(kind=ik)             ::  Bs, nvar, k
    ! variable name
    character(len=80)            :: name

    Bs                  = size(u,1)-2*g   ! number of block sides
    nvar                = size(u,4)       ! number of variables describing the state of the fluid

    ! allocate temporary field
    if ( .not. allocated(tmp_u) )  allocate(tmp_u(size(u,1),size(u,2),size(u,3),nvar))
    tmp_u  = u(:,:,:,:)
    call convert_statevector(tmp_u,'pure_variables')

    ! compute vorticity
    if (  list_contains_name(params_ns%names,'vortx')>0 .or. &
          list_contains_name(params_ns%names,'vort') >0 ) then
      if ( .not. allocated(vort) ) allocate(vort(size(u,1),size(u,2),size(u,3),3))
      call compute_vorticity(tmp_u(:,:,:,UxF), tmp_u(:,:,:,UyF), tmp_u(:,:,:,UzF), &
                            dx, Bs, g, params_ns%discretization, vort)
    end if

    ! compute filter strength
    if (params_ns%filter%name=="bogey_shock" .and. params_ns%filter%save_filter_strength ) then
      if ( .not. allocated(sigma) ) allocate(sigma(size(u,1),size(u,2),size(u,3),3))
      work(:,:,:,1:nvar)=u
      call filter_block(params_ns%filter, time, work(:,:,:,:), Bs, g, x0, dx)
      sigma(:,:,:,1:params_ns%dim)=work(:,:,:,nvar+1:nvar+params_ns%dim)
    endif

    ! save pure state variables (rho, u, v, w, p)
    work(:,:,:,1:nvar)=tmp_u(:,:,:,:)

    ! save additional variables
    do k = nvar+1, params_ns%N_fields_saved
        name = params_ns%names(k)
        select case(trim(name))
        case('vortx','vort')
            work(:,:,:,k)=vort(:,:,:,1)
        case('vorty')
            work(:,:,:,k)=vort(:,:,:,2)
        case('vortz')
            work(:,:,:,k)=vort(:,:,:,3)
        case('sigmax')
            work(:,:,:,k)=sigma(:,:,:,1)
        case('sigmay')
            work(:,:,:,k)=sigma(:,:,:,2)
        case('sigmaz')
            work(:,:,:,k)=sigma(:,:,:,3)
        case('mask')
            if ( params_ns%dim==3 ) then
              !
            else
              call get_mask(work(:,:,1,k), x0, dx, Bs, g )
            end if
        end select
      end do

  end subroutine PREPARE_SAVE_DATA_NStokes


  !-----------------------------------------------------------------------------
  ! when savig to disk, WABBIT would like to know how you named you variables.
  ! e.g. u(:,:,:,1) is called "ux"
  !
  ! the main routine save_fields has to know how you label the stuff you want to
  ! store from the work array, and this routine returns those strings
  !-----------------------------------------------------------------------------
  subroutine FIELD_NAMES_NStokes( N, name )
    implicit none
    ! component index
    integer(kind=ik), intent(in) :: N
    ! returns the name
    character(len=80), intent(out) :: name

    if (allocated(params_ns%names)) then
      name = params_ns%names(N)
    else
      call abort(5554,'Something ricked')
    endif

  end subroutine FIELD_NAMES_NStokes
