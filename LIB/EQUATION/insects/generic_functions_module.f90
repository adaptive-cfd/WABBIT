module generic_functions_module
    use module_globals
    implicit none
    PRIVATE

    PUBLIC :: generic_function1D, eval

    type generic_function1D
        ! wrapper for a generic function. provides coefficients and the type of function (which has to be derivable!)
        character(len=clong) :: function_name

        ! coefficients for the function, their use depends on the funciton at hand.
        real(kind=rk)   :: coeffs(5)
        real(kind=rk)   :: ai(10), a0
        real(kind=rk)   :: bi(10)

    end type
    contains

    subroutine eval(function, value, derivative)
    type(generic_function1D), intent(in) :: function
    real(kind=rk), intent(inout) :: value
    integer, intent(in) :: derivative
    integer :: i, n, j, power
    real(kind=rk) :: coeff_factor
    real(kind=rk) :: output, k
    integer :: nfft

    output = 0.0_rk
    n = size(function%coeffs)

    select case (function%function_name)
    case("polynomial")

        ! polynomial of the form coeffs(0) + coeffs(1) * x + coeffs(2) * x ** 2 + .. + coeffs(5) * x ** 5

        if (derivative < n) then
            do i = n, derivative + 1, -1
                power = i - 1
                coeff_factor = 1.0_rk
                do j = 0, derivative - 1
                    coeff_factor = coeff_factor * real(power - j, rk)
                end do
                output = output * value + function%coeffs(i) * coeff_factor
            end do
        end if

    case("exponential")
        output = function%coeffs(1) * function%coeffs(2)**derivative &
                 * exp(function%coeffs(2) * value)

    case("sin")
        select case (mod(derivative, 4))
        case(0); output =  sin(function%coeffs(2)*(value + function%coeffs(3)))
        case(1); output =  cos(function%coeffs(2)*(value + function%coeffs(3)))
        case(2); output = -sin(function%coeffs(2)*(value + function%coeffs(3)))
        case(3); output = -cos(function%coeffs(2)*(value + function%coeffs(3)))
        end select
        output = output * function%coeffs(1) * function%coeffs(2)**derivative
        if (derivative==0) then
            output = output + function%coeffs(4)
        endif

    case("square")

        ! approximate a square wave using harmonic functions. this makes calculating derivatives easier

        do i = 0, 2

            select case (mod(derivative, 4))
            case(0); output = output + 4.0_rk / 3.14159 / (2.0_rk * real(i, kind=rk) + 1.0_rk) * sin( (2.0_rk * real(i, kind=rk) + 1.0_rk)  * function%coeffs(2)*(value + function%coeffs(3)))* function%coeffs(1) * ( (2.0_rk * real(i, kind=rk) + 1.0_rk)  * function%coeffs(2))**derivative
            case(1); output = output + 4.0_rk / 3.14159 /  (2.0_rk * real(i, kind=rk) + 1.0_rk)  * cos( (2.0_rk * real(i, kind=rk) + 1.0_rk)  * function%coeffs(2)*(value + function%coeffs(3)))* function%coeffs(1) * ( (2.0_rk * real(i, kind=rk) + 1.0_rk)  * function%coeffs(2))**derivative
            case(2); output = output + 4.0_rk / 3.14159 /  (2.0_rk * real(i, kind=rk) + 1.0_rk)  * -sin( (2.0_rk * real(i, kind=rk) + 1.0_rk)  * function%coeffs(2)*(value + function%coeffs(3)))* function%coeffs(1) * ( (2.0_rk * real(i, kind=rk) + 1.0_rk) * function%coeffs(2))**derivative
            case(3); output = output + 4.0_rk / 3.14159 /  (2.0_rk * real(i, kind=rk) + 1.0_rk)  * -cos( (2.0_rk * real(i, kind=rk) + 1.0_rk)  * function%coeffs(2)*(value + function%coeffs(3)))* function%coeffs(1) * ( (2.0_rk * real(i, kind=rk) + 1.0_rk)  * function%coeffs(2))**derivative
            end select

        enddo
        if (derivative==0) then
            output = output + function%coeffs(4)
        endif

    case("fseries")
        nfft = size(function%ai)
        if (derivative == 0) then
            output = function%a0
        end if 
        select case (mod(derivative, 4))


            case(0)
                do i = 1, nfft
                    k = 3.141592654_rk * 2.0_rk * real(i, kind=rk)
                    output = output+ (k ** derivative) * (function%ai(i) * cos(k * value) + function%bi(i) * sin(k * value))
                end do


            case(1)
                do i = 1, nfft
                    k = 3.141592654_rk * 2.0_rk * real(i, kind=rk)
                    output = output+ (k ** derivative) * (-function%ai(i) * sin(k * value) + function%bi(i) * cos(k * value))
                end do

            case(2)
                do i = 1, nfft
                    k = 3.141592654_rk * 2.0_rk * real(i, kind=rk)
                    output = output+ (k ** derivative) * ( - function%ai(i) * cos(k * value) - function%bi(i) * sin(k * value))
                end do

            case(3)
                do i = 1, nfft
                    k = 3.141592654_rk * 2.0_rk * real(i, kind=rk)
                    output = output+ (k ** derivative) * (function%ai(i) * sin(k * value) - function%bi(i) * cos(k * value))
                end do
        end select



    case default
        write(*,*) 'UNKNOWN FUNCTION: >', function%function_name, '<'
        write(*,*) 'TRIMMED:          >', trim(function%function_name), '<'
        error stop
    end select

    value = output
end subroutine

    
end module generic_functions_module