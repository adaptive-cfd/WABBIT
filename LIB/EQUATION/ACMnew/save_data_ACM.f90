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
subroutine PREPARE_SAVE_DATA_ACM( time, u, g, x0, dx, work, mask )
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

    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! On input, the mask array is correctly filled. You cannot create the full mask here.
    real(kind=rk), intent(inout) :: mask(1:,1:,1:,1:)

    ! local variables
    integer(kind=ik)  :: neqn, nwork, k
    integer(kind=ik), dimension(3) :: Bs
    character(len=80) :: name

    ! number of state variables
    neqn = size(u,4)
    ! number of available work array slots
    nwork = size(work,4)

    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    if (params_acm%geometry == "Insect" .and. Insect%time /= time) then
        call Update_Insect(time, Insect)
    endif



    do k = 1, size(params_acm%names,1)
        name = params_acm%names(k)
        select case(name)
        case('ux', 'Ux', 'UX')
            ! copy state vector
            work(:,:,:,k) = u(:,:,:,1)

        case('uy', 'Uy', 'UY')
            ! copy state vector
            work(:,:,:,k) = u(:,:,:,2)

        case('uz', 'Uz', 'UZ')
            ! copy state vector
            work(:,:,:,k) = u(:,:,:,3)

        case('p', 'P')
            ! copy state vector (do not use 4 but rather neq for 2D runs, where p=3rd component)
            work(:,:,:,k) = u(:,:,:,neqn)

        case('vor')
            ! vorticity
            call compute_vorticity(u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), &
            dx, Bs, g, params_acm%discretization, work(:,:,:,k:k+3))

            if (k /= 1) then
                call abort(19101810,"ACM: if you want to store vor, put it at 1st position of field_names. works only in 2D")
            endif

            if (params_acm%dim /= 2) then
                call abort(19101811,"ACM: storing vor is not possible in 3D (known bug)")
            endif

            if (size(params_acm%names,1) < 3) then
                call abort(19101811,"ACM: storing vor requires at least 3 fields to be saved (known bug)")
            endif


        case('div')
            ! div(u)
            call divergence(u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), dx, Bs, g, params_acm%discretization, work(:,:,:,k))

        case('mask')
            work(:,:,:,k) = mask(:,:,:,1)

        case('usx')
            work(:,:,:,k) = mask(:,:,:,2)

        case('usy')
            work(:,:,:,k) = mask(:,:,:,3)

        case('usz')
            work(:,:,:,k) = mask(:,:,:,4)

        case('color')
            work(:,:,:,k) = mask(:,:,:,5)

        case('sponge')
            work(:,:,:,k) = mask(:,:,:,6)

        end select
    end do

end subroutine


!-----------------------------------------------------------------------------
! when savig to disk, WABBIT would like to know how you named you variables.
! e.g. u(:,:,:,1) is called "ux"
!
! the main routine save_fields has to know how you label the stuff you want to
! store from the work array, and this routine returns those strings
!-----------------------------------------------------------------------------
subroutine FIELD_NAMES_ACM( N, name )
    implicit none
    ! component index
    integer(kind=ik), intent(in) :: N
    ! returns the name
    character(len=80), intent(out) :: name

    if (allocated(params_acm%names)) then
        name = params_acm%names(N)
    else
        call abort(5554,'Something ricked')
    endif

end subroutine FIELD_NAMES_ACM
