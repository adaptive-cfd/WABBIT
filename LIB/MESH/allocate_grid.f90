!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name init_data.f90
!> \version 0.5
!> \author msr
!
!> \brief Allocate grid data (light, heavy, neighbors, active lists etc), initialize
!
!>
!! input:
!!           - parameter array
!!           - light data array
!!           - heavy data array
!!           - neighbor data array
!!           - light and heavy active block list
!!
!! output:
!!           - filled user defined data structure for global params
!!           - initialized light and heavy data arrays
!!
!! = log ======================================================================================
!! \n
!! 04/11/16 - switch to v0.4, now run complete initialization within these subroutine and return
!!            initialized block data to main program \n
!! 07/12/16 - now uses heavy work data array \n
!! 25/01/17 - switch to 3D, v0.5
!
! ********************************************************************************************
subroutine allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_sortednumlist, simulation, hvy_work)

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)                   :: params
    !> light data array
    integer(kind=ik), allocatable, intent(out)       :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), allocatable, intent(out)             :: hvy_block(:, :, :, :, :)
    !> heavy work array
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_work(:, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: lgt_active(:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_active(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), allocatable, intent(out)       :: lgt_sortednumlist(:,:)
    ! local shortcuts:
    integer(kind=ik)                                    :: Bs, g, N_dF, number_blocks,&
                                                      rank, number_procs,number_trees
    !> do we have to allocate everything?
    logical, intent(in) :: simulation
    integer(kind=ik)    :: rk_steps
    real(kind=rk)       :: effective_memory
    integer             :: status, nwork

    !---------------------------------------------------------------------------------------------
    ! interfaces

    !---------------------------------------------------------------------------------------------
    ! variables initialization
    ! set parameters for readability
    rank            = params%rank
    number_blocks   = params%number_blocks
    Bs              = params%Bs
    g               = params%n_ghosts
    N_dF            = params%n_eqn
    number_procs    = params%number_procs
    ! the number of work arrays is determined by the requirement for the time stepper
    ! and the number of fields to be saved.
    nwork           = max( N_dF*(size(params%butcher_tableau,1)), params%N_fields_saved )


    !---------------------------------------------------------------------------------------------
    ! main body


    if (rank == 0) then
        write(*,'(80("_"))')
        write(*,'(A)') "INIT: Beginning memory allocation and initialization."
        write(*,'("INIT: mpisize=",i6)') params%number_procs
        write(*,'("INIT: Bs=",i7," blocks-per-rank=",i7," total blocks=", i7)') Bs, number_blocks, number_blocks*number_procs
    endif


    ! allocate memory
    if ( params%threeD_case ) then
        ! 3D:
        if (rank==0) write(*,'("INIT: Allocating a 3D case.")')

        !---------------------------------------------------------------------------
        allocate( hvy_block( Bs+2*g, Bs+2*g, Bs+2*g, N_dF, number_blocks ) )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_block", product(real(shape(hvy_block)))*8.0e-9, shape(hvy_block)
        endif

        !---------------------------------------------------------------------------
        ! work data (Runge-Kutta substeps and old time level)
        if (simulation) then
            allocate( hvy_work( Bs+2*g, Bs+2*g, Bs+2*g, nwork, number_blocks ) )
            if (rank==0) then
                write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
                "hvy_work", product(real(shape(hvy_work)))*8.0e-9, shape(hvy_work)
            endif
        end if

        !---------------------------------------------------------------------------
        ! 3D: maximal 74 neighbors per block
        allocate( hvy_neighbor( params%number_blocks, 74 ) )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_neighbor", product(real(shape(hvy_neighbor)))*8.0e-9, shape(hvy_neighbor)
        endif

    else

        ! 2D:
        if (rank==0) write(*,'("INIT: Allocating a 2D case.")')

        !---------------------------------------------------------------------------
        allocate( hvy_block( Bs+2*g, Bs+2*g, 1, N_dF, number_blocks ) )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_block", product(real(shape(hvy_block)))*8.0e-9, shape(hvy_block)
        endif

        !---------------------------------------------------------------------------
        ! work data (Runge-Kutta substeps and old time level)
        if (simulation) then
            allocate( hvy_work( Bs+2*g, Bs+2*g, 1, nwork, number_blocks ) )
            if (rank==0) then
                write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",5(i9,1x))') &
                "hvy_work", product(real(shape(hvy_work)))*8.0e-9, shape(hvy_work)
            endif
        end if

        !---------------------------------------------------------------------------
        ! 2D: maximal 16 neighbors per block
        allocate( hvy_neighbor( params%number_blocks, 16 ) )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_neighbor", product(real(shape(hvy_neighbor)))*8.0e-9, shape(hvy_neighbor)
        endif

    end if

    !---------------------------------------------------------------------------)
    allocate( lgt_block( number_procs*number_blocks, params%max_treelevel+extra_lgt_fields) )
    call reset_lgt_data(lgt_block,lgt_active,params%max_treelevel)

    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_block", product(real(shape(lgt_block)))*4.0e-9, shape(lgt_block)
    endif

    !---------------------------------------------------------------------------
    allocate( lgt_sortednumlist( size(lgt_block,1), 2) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_sortednumlist", product(real(shape(lgt_sortednumlist)))*4.0e-9, shape(lgt_sortednumlist)
    endif

    ! is not really needed here because it will be reseted in reset_grid!!!
    !call reset_lgt_data(lgt_block, lgt_active, lgt_n, lgt_sortednumlist)

    !---------------------------------------------------------------------------
    allocate( lgt_active( size(lgt_block, 1) ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_active", product(real(shape(lgt_active)))*4.0e-9, shape(lgt_active)
    endif

    !---------------------------------------------------------------------------
    ! note: 5th dimension in heavy data is block id
    allocate( hvy_active( size(hvy_block, 5) ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "hvy_active", product(real(shape(hvy_active)))*4.0e-9, shape(hvy_active)
    endif


    if (rank == 0 .and. simulation) then
        write(*,'("INIT: System is allocating heavy data for ",i7," blocks and ", i3, " fields" )') number_blocks, N_dF
        write(*,'("INIT: System is allocating light data for ",i7," blocks" )') number_procs*number_blocks
        write(*,'("INIT: System is allocating heavy work data for ",i7," blocks " )') number_blocks

        effective_memory = (dble(size(hvy_block)+size(hvy_work)) + & ! real data
        dble(size(lgt_block)+size(lgt_sortednumlist)+size(hvy_neighbor)+size(lgt_active)+size(hvy_active))/2.0 & ! integer (hence /2)
        )*8.0e-9 ! in GB

        ! note we currently use 8byte per real and integer by default, so all the same bytes per point
        write(*,'("INIT: Measured (true) local (on 1 cpu) memory (hvy_block+hvy_work+lgt_block no ghosts!) is ",g15.3,"GB per mpirank")') &
        effective_memory

        write(*,'("INIT-GLOBAL: Measured (true) TOTAL (on all CPU) memory (hvy_block+hvy_work+lgt_block no ghosts!) is ",g15.3,"GB")') &
        effective_memory*dble(params%number_procs)
    end if

end subroutine allocate_grid





subroutine deallocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_sortednumlist, hvy_work )

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)                   :: params
    !> light data array
    integer(kind=ik), allocatable, intent(out)          :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), allocatable, intent(out)             :: hvy_block(:, :, :, :, :)
    !> heavy work array  )
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_work(:, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: lgt_active(:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_active(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), allocatable, intent(out)       :: lgt_sortednumlist(:,:)

    if (params%rank == 0) then
        write(*,'(80("-"))')
        write(*,'(A)') "FREE: Beginning freeying of memory."
    endif

    if (allocated(hvy_block)) deallocate( hvy_block )
    if (allocated(hvy_work)) deallocate( hvy_work )
    if (allocated(hvy_neighbor)) deallocate( hvy_neighbor )
    if (allocated(lgt_block)) deallocate( lgt_block )
    if (allocated(lgt_sortednumlist)) deallocate( lgt_sortednumlist )
    if (allocated(lgt_active)) deallocate( lgt_active )
    if (allocated(hvy_active)) deallocate( hvy_active )

    if (params%rank == 0) then
        write(*,'(A)') "All memory is cleared!"
        write(*,'(80("-"))')
    endif

end subroutine deallocate_grid
