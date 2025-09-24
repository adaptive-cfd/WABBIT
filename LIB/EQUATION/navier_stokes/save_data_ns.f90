
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
  subroutine PREPARE_SAVE_DATA_NStokes(time, u, g, x0, dx, work, boundary_flag )
    use module_helpers , only: choose, list_contains_name
    use module_navier_stokes_cases, only: get_mask
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
    ! when implementing boundary conditions, it is necessary to know if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The boundary_flag tells you if the local field is adjacent to a domain boundary:
    ! boundary_flag(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    ! currently only acessible in the local stage
    integer(kind=2)          , intent(in):: boundary_flag(3)
    ! output in work array.
    real(kind=rk), allocatable,save :: tmp_u(:,:,:,:), vort(:,:,:,:), sigma(:,:,:,:), mask(:,:,:)
    ! local variables
    integer(kind=ik)             ::  nvar, k
    integer(kind=ik), dimension(3) :: Bs
    ! variable name
    character(len=cshort)            :: name


    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g                 ! number of block sides
    nvar = params_ns%n_eqn  ! number of variables describing the state of the fluid
    ! allocate temporary field
    if ( .not. allocated(tmp_u) ) call allocate_statevector_ns(tmp_u,Bs,g)

    tmp_u  = u(:,:,:,:)
    work(:,:,:,1:nvar)=u
    if ( .not. ALL(params_ns%periodic_BC)) then
      call compute_boundary_2D( time, g, Bs, dx, x0, tmp_u(:,:,1,:), boundary_flag)
      call compute_boundary_2D( time, g, Bs, dx, x0, work(:,:,1,:), boundary_flag)
    end if

    !+++++++++++++++++++
    ! save pure state variables (rho, u, v, w, p)
    !++++++++++++++++++
    call convert_statevector(work(:,:,:,1:nvar),'pure_variables')
    ! save additional variables
    call convert_statevector(tmp_u,'pure_variables')

    ! +++++++++++++++++
    ! compute vorticity
    ! +++++++++++++++++
    if (  list_contains_name(params_ns%names,'vortx')>0 .or. &
          list_contains_name(params_ns%names,'vort') >0 ) then
      if ( .not. allocated(vort) ) allocate(vort(size(u,1),size(u,2),size(u,3),3))
      ! JB comment: compute_vorticity expects a 4D array, so we have to pass the full state vector
      ! I assume here that UyF = UxF+1 and UzF = UxF+2
      if (Uyf /= UxF+1 .or. (params_ns%dim == 3 .and. Uzf /= UxF+2)) then
        call abort(250703, 'ERROR: UyF and UzF do not match the expected indices for vorticity computation, these should be UxF+1 and UxF+2.')
      end if
      call compute_vorticity(tmp_u(:,:,:,UxF:UxF+2), &
                            dx, Bs, g, params_ns%discretization, vort)
    end if

    ! +++++++++++++++++++++++++++++++++
    ! compute mask and reference values
    ! +++++++++++++++++++++++++++++++++
    if (  list_contains_name(params_ns%names,'mask')>0 ) then
      if ( .not. allocated(mask) ) allocate(mask(size(u,1),size(u,2),size(u,3)))
          call get_mask(params_ns, x0, dx, Bs, g, mask, .true.)   ! the true boolean stands for: make mask colored if possible
    end if

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
            work(:,:,:,k)=mask
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
    character(len=cshort), intent(out) :: name

    if (allocated(params_ns%names)) then
      name = params_ns%names(N)
    else
      call abort(5554,'Something ricked')
    endif

  end subroutine FIELD_NAMES_NStokes
