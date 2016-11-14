! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: com_allowed.f90
! version: 0.4
! author: msr
!
! check if communication between two procs is possible -> if more than one .true. in array
!
! input:    - logical array
!           - size of array
! output:   - .true. or .false.
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
! ********************************************************************************************

logical function com_allowed(array, N)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! array size
    integer, intent(in)                   :: N
    ! array with proc status
    logical, intent(in)                   :: array(N)

    ! loop variable
    integer                               :: i
    ! found one .true. in array
    logical                               :: one_com_allowed

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    com_allowed   = .false.
    one_com_allowed = .false.

!---------------------------------------------------------------------------------------------
! main body

    do i = 1, N

        if ( one_com_allowed .eqv. .true. ) then
            if ( array(i) .eqv. .true. ) then
                ! second allowed communication detected
                com_allowed = .true.
            end if
        end if

        if ( one_com_allowed .eqv. .false. ) then
            if ( array(i) .eqv. .true. ) then
                ! first allowed communication detected
                one_com_allowed = .true.
            end if
        end if

    enddo

end function com_allowed
