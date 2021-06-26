!> \brief module for data precision
! ********************************************************************************************

module module_precision

    implicit none
    ! define data precision parameters
    integer, parameter, public      :: sngl_prec=selected_real_kind(4)
    integer, parameter, public      :: dble_prec=selected_real_kind(8)
    ! default length of strings (short and long character)
    integer, parameter, public      :: cshort=80
    integer, parameter, public      :: clong =120
    integer, parameter, public      :: int_prec=selected_int_kind(8)
    integer, parameter, public      :: maxdigits = 16
    integer, parameter, public      :: tsize = selected_int_kind(maxdigits)
    integer, parameter, public      :: rk=dble_prec
    integer, parameter, public      :: ik=int_prec
    real(kind=rk),parameter, public :: pi  = 4.0_rk * atan(1.0_rk)

contains

    !-----------------------------------------------------------------------------
    ! convert degree to radiant
    !-----------------------------------------------------------------------------
    real(kind=rk) function deg2rad(deg)
        implicit none
        real(kind=rk), intent(in) :: deg
        real(kind=rk),parameter :: pi2 = 3.1415926535897932384626433832795028841971693993751058209749445923078164d0
        deg2rad = deg*pi2/180.d0
        return
    end function

    !-----------------------------------------------------------------------------
    ! radiant to degree
    !-----------------------------------------------------------------------------
    real(kind=rk) function rad2deg(deg)
        implicit none
        real(kind=rk), intent(in) :: deg
        real(kind=rk),parameter :: pi2 = 3.1415926535897932384626433832795028841971693993751058209749445923078164d0
        rad2deg = deg*180.d0/pi2
        return
    end function

    !-----------------------------------------------------------------------------
    ! cross product of two vectors
    !-----------------------------------------------------------------------------
    function cross(a,b)
        implicit none
        real(kind=rk),dimension(1:3),intent(in) :: a,b
        real(kind=rk),dimension(1:3) :: cross
        cross(1) = a(2)*b(3)-a(3)*b(2)
        cross(2) = a(3)*b(1)-a(1)*b(3)
        cross(3) = a(1)*b(2)-a(2)*b(1)
    end function

    !-----------------------------------------------------------------------------
    ! angle between two vectors
    !-----------------------------------------------------------------------------
    function angle_between_vectors(a,b)
        implicit none
        real(kind=rk),dimension(1:3),intent(in) :: a,b
        real(kind=rk) :: angle_between_vectors

        angle_between_vectors = acos(sum(a*b) / (norm2(a)*norm2(b)))

    end function

    !-----------------------------------------------------------------------------
    ! turn vector into unit vector
    !-----------------------------------------------------------------------------
    function unit_vector(a)
        implicit none
        real(kind=rk),dimension(1:3),intent(in) :: a
        real(kind=rk),dimension(1:3) :: unit_vector
        unit_vector = a / norm2(a)
    end function

    !-----------------------------------------------------------------------------
    ! check if two vectors are identical
    !-----------------------------------------------------------------------------
    function is_same_vector(a,b)
        implicit none
        real(kind=rk),dimension(1:),intent(in) :: a,b
        logical :: is_same_vector
        integer :: i

        is_same_vector = .true.

        do i = 1, size(a)
            if ( abs(a(i)-b(i)) >= 1.0e-11 ) then
                is_same_vector = .false.
                return
            endif
        enddo
    end function

    !-----------------------------------------------------------------------------
    ! 2-norm length of vectors
    !-----------------------------------------------------------------------------
    function norm2(a)
        implicit none
        real(kind=rk),dimension(1:3),intent(in) :: a
        real(kind=rk) :: norm2
        norm2 = sqrt( a(1)*a(1) + a(2)*a(2) + a(3)*a(3) )
    end function

    !-----------------------------------------------------------------------------
    ! given a point x, check if it lies in the computational domain centered at zero
    ! (note: we assume [-xl/2...+xl/2] size this is useful for insects )
    !-----------------------------------------------------------------------------
    function periodize_coordinate(x_glob, box)
        real(kind=rk),intent(in) :: x_glob(1:3), box(1:3)
        real(kind=rk),dimension(1:3) :: periodize_coordinate

        periodize_coordinate = x_glob

        if (x_glob(1)<-box(1)*0.5_rk) periodize_coordinate(1)=x_glob(1)+box(1)
        if (x_glob(2)<-box(2)*0.5_rk) periodize_coordinate(2)=x_glob(2)+box(2)
        if (x_glob(3)<-box(3)*0.5_rk) periodize_coordinate(3)=x_glob(3)+box(3)

        if (x_glob(1)>box(1)*0.5_rk) periodize_coordinate(1)=x_glob(1)-box(1)
        if (x_glob(2)>box(2)*0.5_rk) periodize_coordinate(2)=x_glob(2)-box(2)
        if (x_glob(3)>box(3)*0.5_rk) periodize_coordinate(3)=x_glob(3)-box(3)

    end function

    !---------------------------------------------------------------------------
    ! wrapper for random number generator (this may be compiler dependent)
    !---------------------------------------------------------------------------
    real(kind=rk) function rand_nbr()
        implicit none
        call random_number( rand_nbr )
    end function

end module module_precision
