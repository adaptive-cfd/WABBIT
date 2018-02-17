!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name get_inicond_from_file.f90
!> \version 0.5
!> \author sm
!
!> \brief call subroutines that read mesh and fields as initial condition from files
!
!>
!! input:
!!           - parameter array
!!
!! output:
!!           - light data array
!!           - heavy data array
!!           - number of active blocks (light and heavy)
!!           - time and iteration
!!
!!
!! = log ======================================================================================
!! \n
!! 28/09/17 - create
!
! ********************************************************************************************

subroutine get_inicond_from_file(params, lgt_block, hvy_block, hvy_n, lgt_n, time, iteration)
!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)        :: params
    !> light data array
    integer(kind=ik), intent(inout)       :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)          :: hvy_block(:, :, :, :, :)
    !> number of heavy and light active blocks
    integer(kind=ik), intent(inout)       :: hvy_n, lgt_n
    !> time loop variables
    real(kind=rk), intent(inout)          :: time
    integer(kind=ik), intent(inout)       :: iteration
    real(kind=rk), dimension(3)           :: domain

    ! cpu time variables for running time calculation
    real(kind=rk)                         :: t0
    ! number of files to read from
    integer(kind=ik)                      :: N_files
    ! loop variable
    integer(kind=ik)                      :: k, dF
!---------------------------------------------------------------------------------------------
! variables initialization

    ! number of files to read from
    N_files = params%number_data_fields
    ! start time
    t0 = MPI_wtime()

    if (params%rank==0) write(*,*) "Reading initial condition from file"
!---------------------------------------------------------------------------------------------
! main body

    ! read time, iteration, domain size and total number of blocks from first input file
    call get_attributes(params%input_files(1), lgt_n, time, iteration, domain)
    if (.not. params%threeD_case) domain(3) = 0.0_rk
    ! print time, iteration and domain on screen
    if (params%rank==0) then
        write(*,'(80("_"))')
        write(*,'("READING: Reading from file ",A)') trim(adjustl(params%input_files(1)))
        write(*,'("time=",g12.4," iteration=", i5)') time, iteration
        write(*,'("Lx=",g12.4," Ly=",g12.4," Lz=",g12.4)') domain
        ! if the domain size doesn't match, proceed, but yell.
        if ((abs(params%Lx-domain(1))>1e-12_rk).or.(abs(params%Ly-domain(2))>1e-12_rk) &
            .or.(abs(params%Lz-domain(3))>1e-12_rk)) then
            write (*,'(A)') " WARNING! Domain size mismatch."
            write (*,'("in memory:   Lx=",es12.4,"Ly=",es12.4,"Lz=",es12.4)') params%Lx, params%Ly, params%Lz
            write (*,'("but in file: Lx=",es12.4,"Ly=",es12.4,"Lz=",es12.4)') domain
            write (*,'(A)') "proceed, with fingers crossed."
        end if
    end if

    ! read treecode from first input file
    call read_mesh(params%input_files(1), params, lgt_n, hvy_n, lgt_block)

    ! read datafields from files into hvy_block array
    do dF = 1, N_files
        call read_field(params%input_files(dF), dF, params, hvy_block, hvy_n )
    end do

    ! timing
    call toc( params, "read_data", MPI_wtime()-t0 )
end subroutine get_inicond_from_file
