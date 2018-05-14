!===========================================================================
!> this module implements the MPI_bridge between fluid and particle on the side
!> of the fluid
!> @todo
!> 1. it initializes the bridge see init_bridge()
!> 2. it provides a routine which sends the treecode and further information
!>    to the particle side of the bridge. This is needed to map the
!>    particle coordinates onto the right processors subdomain
!> 3. interpolates the velocity field \f$(u_\alpha)_{ijk}\f$ to a requested position \f$x \in R^3\f$
!>    and sends it back to the requesting processor
!-------------------------------------------------------------------------!!!!
module module_bridge_interface

use mpi
use module_bridge
use module_treelib
use module_mesh
use module_ini_files_parser_mpi
! use constants
! use methods
! use geometry
! use parametersParticle
! use initializationParticle
! use commonParticle
! use chargedParticle
! use bridge

use module_params

implicit none
private

! ! Declarations

! !! Variable relatives to the other world (fluid domain)
! !! - array corresponding to the parameters of the grid (number of points in each direction)
! !!   1. total number of grid points
! !!   2. number of grid points in one process
! integer         , dimension(3,2), private  :: nbGridPointsFluid
!! - array corresponding to the number of fluid processes in each direction (1,1,1 if only one process)
integer         , dimension(3)  , private  :: processesParticle
! !! - streching parameter to modify the coordinates - WARNING: only one to begin
! double precision,                 private  :: stretchingParameter

! !! Useful subroutines
! public :: getParticleProcesses
! public :: getProcessWithPosition
! public :: exchangeParametersParticle
public :: init_bridge
public :: initialize_communicator
public :: send_lgt_data
public :: serve_data_request
public :: position_to_lgt_id
! public :: signalEndCalculation
! public :: particleDistributionOverFluidProcesses
! public :: exchangeDataParticle




 contains

!===========================================================================
!> \brief read out from the command line the amount of fluid processes to start
subroutine getParticleProcesses(params)

!> structure containing the parameters
type(type_params), intent(in)         :: params

! (total) number of fluid processes
character(10)                            :: nbParticleProcesses

  ! Read out the total number of particle processes
  !to start from the given command line
  call get_command_argument(3, nbParticleProcesses)

  ! Assign this number in the first value of processesParticle (1 otherwise)
  processesParticle = 1
  read(nbParticleProcesses,'(i10)') processesParticle(1)
  if ( params%rank==0 ) then
    write(*,*) 'Number of Particle Processes: ', processesParticle(1)
  endif
end subroutine getParticleProcesses
!===========================================================================


!===========================================================================
subroutine send_fixed_params (myBridge, params)
!> \brief this routine sends all fixed parameters (e.g. \f$L_x,L_y)\f$ to the particle world
!! Subroutine-declarations
type(bridgeMPI)  , intent(in)                :: myBridge           !< type bridge on the particle side
type(type_params), intent(in)                :: params             !> structure containing the parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer                                      :: ierr               ! MPI communication error
integer         , dimension(MPI_STATUS_SIZE) :: status             ! MPI communication status
double precision, dimension(3,2)             :: Domain     ! received domain parameters to transform
integer         , dimension(4)             :: discretizationParams ! Number of grid points and processes
character(len=80)                            :: geometry


  ! Only the root particle process receive the parameters from the root fluid process
  if (myBridge%myWorldRank == 0) then                              ! check if the process is the root process

    ! - receive the discretization parameters
    discretizationParams=(/ params%number_block_nodes,  &
                            params%number_blocks,       &
                            params%max_treelevel,       &
                            params%dim,                 &
                            params%number_procs/)

    call MPI_send(discretizationParams, 5, MPI_integer, myBridge%minOtherWorldRank, &
                  parameters_delivery, myBridge%commonWorld,  ierr)

    !! - receive the domain parameters
    Domain(:,1)=(/params%Lx,params%Ly,params%Lz/) !> Domain size
    Domain(:,2)=(/0, 0, 0/)                         !> Domain origin
    call MPI_send(Domain, 6, MPI_double_precision, myBridge%minOtherWorldRank, &
                  parameters_delivery, myBridge%commonWorld, ierr)

    !! - receive the coordinate system type parameters (overwrite the one given in the particle parameters)
    geometry='cartesian' ! 'cartesian' / 'cylindrical'
    call MPI_send(geometry, 80, MPI_character, myBridge%minOtherWorldRank, &
                  parameters_delivery, myBridge%commonWorld, ierr)
  end if
end subroutine send_fixed_params
!===========================================================================



!===========================================================================
subroutine send_lgt_data (lgt_block,lgt_active,lgt_n,params)
!> \brief this routine sends all fixed parameters (e.g. \f$L_x,L_y)\f$ to the particle world
! Subroutine-declarations
type(type_params), intent(in)         :: params             !< structure containing the parameters
integer(kind=ik), intent(in)          :: lgt_block(:, :)    !< light data array
integer(kind=ik), intent(in)          :: lgt_active(:)      !< list of active blocks (light data)
integer(kind=ik), intent(in)          :: lgt_n              !< number of active blocks (light data)

!--------------------------------------------------------------------------------
integer                              :: ierr               ! MPI communication error
integer , dimension(MPI_STATUS_SIZE) :: status             ! MPI communication status
integer                              :: n,m

  ! Only the root particle process receive the parameters from the root fluid process
  if ( params%bridge%myWorldRank == 0 ) then                              ! check if the process is the root process
    ! number of blocks
    n      = size(lgt_block,1)
    ! number of columns in matrix
    m        = params%max_treelevel+2

    ! send number of active and maximal number of blocks
    call MPI_send((/lgt_n,n/), 2, MPI_integer, &
                  params%bridge%minOtherWorldRank, parameters_delivery, &
                  params%bridge%commonWorld,  ierr )
    ! send list of active blocks
    call MPI_send(lgt_active, n, MPI_integer, &
                  params%bridge%minOtherWorldRank, parameters_delivery, &
                  params%bridge%commonWorld,  ierr )
    ! send light data
    call MPI_send(lgt_block, n*m, MPI_integer, &
                  params%bridge%minOtherWorldRank, parameters_delivery, &
                  params%bridge%commonWorld, ierr )
  end if
end subroutine send_lgt_data
!===========================================================================




!===========================================================================
subroutine serve_data_request(lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n,params)
!> \brief this routine sends all fixed parameters (e.g. \f$L_x,L_y)\f$ to the particle world
! Subroutine-declarations
type(type_params), intent(in)         :: params             !< structure containing the parameters
integer(kind=ik), intent(in)          :: lgt_block(:, :)    !< light data array
integer(kind=ik), intent(in)          :: lgt_active(:)      !< list of active blocks (light data)
integer(kind=ik), intent(in)          :: lgt_n              !< number of active blocks (light data)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: hvy_n
!--------------------------------------------------------------------------------
integer                              :: ierr ,k,maxpoints            ! MPI communication error
integer , dimension(MPI_STATUS_SIZE) :: status             ! MPI communication status
integer                              :: Ndat,request,sta
integer         , allocatable                   :: statuses(:,:)   ! MPI communication status array
integer         , allocatable                   :: requests(:)   !  array
integer                              :: send_id           !! rank of sender
double precision, allocatable        :: distributedParticles(:,:) ! stores the distributed value
double precision, allocatable        :: u_inter(:,:) ! stores the distributed value

character(1)                                    :: buf ! Message sent to the fluid world

!> \detail
!! \c serve_data_request is
!!    - if \c MPI_TAG is \c data_positionDelivery
!!          -> interpolate fluid data to requestet position
!!          -> send fluid data back to the sender (i.e. particle rank)
!!    - if \c MPI_TAG is \c end_communication stop waiting for requests

maxpoints=600
k=0

allocate(distributedParticles(4,maxpoints))
allocate(requests(params%bridge%otherWorldSize))
do
    call MPI_probe(MPI_ANY_SOURCE,MPI_ANY_TAG,params%bridge%otherWorld,status,ierr)

    !=========================================
    select case(status(MPI_TAG))
    !======= case 1: DATA delivery
    case(data_positionDelivery)
    !=======
        send_id=status(MPI_SOURCE)
      !  write(*,*) 'send_id:',send_id
        k=k+1
        call MPI_irecv(distributedParticles,4*maxpoints , MPI_double_precision, send_id, &
                    data_positionDelivery,params%bridge%otherWorld, request,ierr)
        call MPI_Get_count( status,MPI_double_precision , Ndat , ierr)
        ! divide Ndat to get the particle Number
        Ndat=Ndat/4
       ! write(*,*) "number of positions to interpolate=",Ndat
         call interpolate_data(lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n,params,distributedParticles(:,1:Ndat),u_inter)
         call MPI_isend(u_inter, 6*Ndat, MPI_double_precision, send_id, data_fieldDelivery, &
                     params%bridge%otherWorld,requests(k),  ierr)
        ! write(*,*) 'number of particles:',Ndat/4,'positiondata=', distributedParticles(:,1)
    !======= case 2: end communication
    case (end_communication)
    !=======
        call MPI_recv(buf, 1, MPI_character, params%bridge%minOtherWorldRank,end_communication,params%bridge%otherWorld, status,ierr)
        allocate(statuses(MPI_STATUS_SIZE,k))
       ! write(*,*)"k=",k
        call MPI_waitall(k, requests, statuses, ierr)
        !write(*,*) "fluidrank",params%bridge%myWorldRank,"processed data request"
        exit
    end select
    !=========================================
enddo
end subroutine serve_data_request
!===========================================================================



!===========================================================================
! subroutine send_params (myBridge, givenParams)
! !! Exchange of a message between the particle process and a fluid process, according to the content of the message
! !! Subroutine-declarations
! type(bridgeMPI)  , intent(in)                :: myBridge           ! type bridge on the particle side
! type(type_params), intent(in), optional      :: givenParams        ! given structure containing the parameters
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! type(type_params)                            :: paramsInit         ! parameter structure effectively used
! integer                                      :: ierr               ! MPI communication error
! integer         , dimension(MPI_STATUS_SIZE) :: status             ! MPI communication status
! integer         , dimension(3,2)             :: discretizationParams ! Number of grid points and processes
! logical         , dimension(3,2)             :: periodicityParams  ! Periodicity parameters of the domain
! double precision, dimension(3,2)             :: receivedDomain     ! received domain parameters to transform

!   ! Check if a parameters structure is given
!   if (present(givenParams)) then                                   ! if a parameters structure is given
!     paramsInit = givenParams                                       ! associate the effective parameters structure to the given one
!   else                                                             ! if no parameter structure is given
!     paramsInit = parameters                                        ! the effective parameters structure if the default one of the module
!   end if                                                           ! end condition regarding the present of the parameters structure

!   ! Only the root particle process receive the parameters from the root fluid process
!   if (myBridge%myWorldRank == 0) then                              ! check if the process is the root process
!     ! - receive the discretization parameters
!     call MPI_send(discretizationParams, 6, MPI_integer, myBridge%minOtherWorldRank, &
!                   parameters_delivery, myBridge%commonWorld, status, ierr)
!     !! - receive the domain parameters
!     call MPI_recv(receivedDomain, 6, MPI_double_precision, myBridge%minOtherWorldRank, &
!                   parameters_delivery, myBridge%commonWorld, status, ierr)
!     !! - receive the coordinate system type parameters (overwrite the one given in the particle parameters)
!     call MPI_recv(paramsInit%domainCoordinates, 100, MPI_character, myBridge%minOtherWorldRank, &
!                   parameters_delivery, myBridge%commonWorld, status, ierr)
!     if (trim(paramsInit%domainCoordinates) == 'cylindrical') &
!         paramsInit%domainCoordinates = 'cylindrical2Pi'            ! adapt the value according to the implementation
!     !! - receive the periodicity parameters
!     call MPI_recv(periodicityParams, 6, MPI_logical, myBridge%minOtherWorldRank, &
!                   parameters_delivery, myBridge%commonWorld, status, ierr)
!     !! - receive the stretching parameter
!     call MPI_recv(stretchingParameter, 1, MPI_double_precision, myBridge%minOtherWorldRank, &
!                   parameters_delivery, myBridge%commonWorld, status, ierr)
!   end if                                                           ! end condition regarding the rank

!   ! Broadcast the received parameters to all particle processes
!   !! - discretization parameters
!   call MPI_bcast(discretizationParams, 6, MPI_integer, 0, myBridge%myWorld, ierr)
!   nbGridPointsFluid(:,1) = discretizationParams(:,1)               ! assign the number of grid points in the fluid world
!   processesFluid(:)      = discretizationParams(:,2)               ! assign the number of processes in the fluid world
!   nbGridPointsFluid(:,2) = nbGridPointsFluid(:,1) / processesFluid(:) ! compute the number of grid points in each process
! !   write(*,*) 'Particle side - amount of grid points received: ', nbGridPointsFluid
!   !! - domain lengths
!   call MPI_bcast(receivedDomain, 6, MPI_double_precision, 0, myBridge%myWorld, ierr)
!   !! - coordinate system type
!   call MPI_bcast(paramsInit%domainCoordinates, 100, MPI_character, 0, myBridge%myWorld, ierr)
!   !! - Effective Periodicity parameters
!   call MPI_bcast(periodicityParams, 6, MPI_logical, 0, myBridge%myWorld, ierr)
!   !! - Streching parameters
!   call MPI_bcast(stretchingParameter, 1, MPI_double_precision, 0, myBridge%myWorld, ierr)

!   ! Treat the received parameters appropriately
!   !! Store the domain parameters
!   paramsInit%domainLength = receivedDomain(:,1)                  ! fluid domain length
!   paramsInit%domainOrigin = receivedDomain(:,2)                  ! fluid domain origin
!   !! Coordinate transformation (in case of NSF calculation) (r <-> theta)
!   paramsInit%domainLength(2) = receivedDomain(1,1)               ! fluid domain length
!   paramsInit%domainLength(1) = receivedDomain(2,1)               ! fluid domain length
!   paramsInit%domainOrigin(2) = receivedDomain(1,2)               ! fluid domain origin
!   paramsInit%domainOrigin(2) = receivedDomain(2,2)               ! fluid domain origin
!   !! Store the pariodicity parameters
!   paramsInit%domainEffectivePeriodicity  = periodicityParams(:,1)
!   paramsInit%domainArtificialPeriodicity = periodicityParams(:,2)

! end subroutine send_params
!===========================================================================




!> \brief Initialize global communicator
!> \details There are two possible scenarios:
!!      + Bridge Modus:
!!            - 2 communicators for the two MPI-WORLDS are created
!!      + Default Modus:
!!            - the MPI_COMM_WORLD is duplicated to WABBIT_COMM:
!!            - thus it creates a new communicator that has a new communication
!!              context but contains the same group of processes as the
!!              input communicator
subroutine initialize_communicator(params)

!> given structure containing the parameters
type(type_params), intent(inout), optional  :: params
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer                                     :: ierr

  ! Initialize the particle MPI world (according to the given kind of hierarchy)
  if (params%bridge_exists) then                             ! if there is a single WABBIT_COMM

    if (params%rank==0) then
       write(*,*) ''
       write(*,*) '-----------------'
       write(*,*) 'Modus: BRIDGE'
       write(*,*) '-----------------'
    end if
    params%bridge%myWorld=params%WABBIT_COMM
    call init_bridge(params%bridge,params)                              ! initialization using split
  else                                                              ! if there are different WABBIT_COMM
     call MPI_comm_dup(WABBIT_COMM, params%bridge%myWorld, ierr)
     call MPI_comm_size(params%bridge%myWorld, params%bridge%myWorldSize, ierr)
     call MPI_comm_rank(params%bridge%myWorld, params%bridge%myWorldRank, ierr)
  end if
  ! after communicator is created set it to the global communicator of WABBIT
  call set_mpi_comm_global(params%bridge%myWorld)
  params%rank         = params%bridge%myWorldRank
  params%number_procs = params%bridge%myWorldSize
  params%WABBIT_COMM  = params%bridge%myWorld
end subroutine initialize_communicator

!===========================================================================


!===========================================================================

subroutine init_bridge(myBridge,params)
!! call every necessary subroutine / function to initialize the fluid world
!! Subroutine-declarations
type(bridgeMPI)  , intent(out)              :: myBridge            ! type bridge on the fluid side
type(type_params), intent(inout)            :: params         ! given structure containing the parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integer                                     :: unitNum = 1         ! file handle
integer                                     :: ierr                ! MPI communication error

 !    paramsInit = params                                       ! associate the effective parameters structure to the given one

  ! Get the total number of fluid processes
  call getParticleProcesses(params)
  ! Initialize the particle MPI world (according to the given kind of hierarchy)
  if (params%bridgeCommonMPI) then                             ! if there is a single WABBIT_COMM
    call createMPIWorlds(myBridge, 1)                              ! initialization using split
  else                                                             ! if there are different WABBIT_COMM
    if (params%bridgeFluidMaster) then
!      write(*,*) "write particlecommand: ", params%particleCommand                   ! if the fluid side is the master
      call createMPIWorldsMaster(myBridge, processesParticle(1), params%particleCommand) ! initialization using spawn (slave)
    else                                                           ! if the fluid side is the slave
      call createMPIWorldsSlave(myBridge)                          ! initialization using spawn (master)
    end if                                                         ! end condition regarding the role of the particle side
  end if                                                           ! end condition regarding the WABBIT_COMM structure
!   write(*,*) 'Fluid side - MPI worlds are created'
   call send_fixed_params(myBridge,params)
!   write(*,*) 'fixed params send to particle world'
  ! ! Make sure the given number of process for the particle world match the number of processes in the parameters file
  ! if (myBridge%myWorldSize /= (paramsInit%NumberProcesses)) then   ! condition on the total number of processes
  !   if (myBridge%myWorldRank == 0) then                            ! only display the error message if this is the root particle process
  !     write(*,*) "The number of processes given in mpiexec for the particles does't match the &
  !                 &number of processes given in the parameters file. Stop!"
  !     write(*,*) "Given number of particle processes: ", myBridge%myWorldSize
  !     write(*,*) "Expected number of particle processes: ", paramsInit%NumberProcesses
  !     call MPI_abort(myBridge%commonWorld, 99, ierr)               ! Abort the execution
  !   end if                                                         ! end condition on the rank in the fluid world
  ! end if                                                           ! end condition on the total number of processes

!   ! Get the total fluid domain values from the fluid world
! !   write(*,*) 'Particle side - Try to communicate parameters with fluid'
!   call exchangeParametersParticle(myBridge, paramsInit)
! !   if (myBridge%myWorldRank == 0) then
! !   write(*,*) 'Particle side - I got from NSF'
! !   write(*,*) 'Number of discretization points: ', nbGridPointsFluid
! !   write(*,*) 'Number of NSF processes: ', processesFluid
! !   end if

end subroutine init_bridge

!===========================================================================

!> \brief function returns the lgt_id to the corresponding treecode
integer(kind=tsize) function treecode2lgt_id(lgt_block,lgt_active,lgt_n,maxtreelevel,treecode)

! input:
implicit none
integer(kind=ik), intent(in)    :: treecode(:)      !< treecode
integer(kind=ik), intent(in)    :: lgt_block(:, :)  !< light data array
integer(kind=ik), intent(in)    :: lgt_active(:)      !< list of active blocks (light data)
integer(kind=ik), intent(in)    :: lgt_n              !< number of active blocks (light data)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer         , allocatable           :: meshlevel(:)        ! mesh level of the blocks
integer                                 :: nblock,i,j,k    ! Number of blocks
integer                                 :: maxtreelevel


    nblock= size(lgt_block,1)
    allocate(meshlevel(lgt_n))
    meshlevel(:)=lgt_block(lgt_active(1:lgt_n),maxtreelevel+1)

    do k=1,lgt_n !loop over blocks k
      ! only loop over active blocks
      j = lgt_active(k)
      i = 1
      do while(lgt_block(j,i)==treecode(maxtreelevel-i) .and. i<=meshlevel(k))
            i=i+1
      enddo

      if ((i-1)==meshlevel(k)) then
        treecode2lgt_id=j
       ! write(*,*) 'light_id=',treecode2lgt_id

        !return
      endif

    enddo

end function







!===========================================================================

!> \brief return the ligth_id for the given position
integer(kind=tsize) function position_to_lgt_id (lgt_block,lgt_active,lgt_n,position, params)
! Function-declarations
double precision, dimension(3)  , intent(in)    :: position         !< position vector \f$ 0<position(i)\le 1\f$
type(type_params),                intent(in)    :: params           !< params of infile
integer(kind=ik), intent(in)    :: lgt_block(:, :)  !< light data array
integer(kind=ik), intent(in)    :: lgt_active(:)      !< list of active blocks (light data)
integer(kind=ik), intent(in)    :: lgt_n              !< number of active blocks (light data)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer, dimension(3)           :: blockLength     ! number of grid points in each block
integer, dimension(4,8)         :: procIndex       ! index of the processes in each direction for the 8 possibilities and activity flag
integer                         :: d               ! count integer for dimension
double precision, dimension(3)  :: dx,x0,x1              ! seperation
integer                         :: k
!> \detail

    do k=1,lgt_n ! loop over all active blocks

       call get_block_spacing_origin( params, lgt_active(k), lgt_block, x0, dx )
       x1=x0(:)+dx(:)*(params%number_block_nodes-1)
       ! check if position is inside the block with origin x0 and size dx
       ! if yes:   found the right block with its lgt_id
       ! if not:   go to the next block
       !                    x0(1) -------------------------- x1(1)
       !                         |                           |
       !                         |         lgt_id            |
       !                         |               *            |
       !                         |               (position)  |
       !                         |                           |
       !                    x0(2) -------------------------- x1(2)

       if (min(position(1)-x0(1),position(2)-x0(2))>0 .and. max(position(1)-x1(1),position(2)-x1(2))<0) then
          if (params%threeD_case) then ! in 3d we have to check the 3. komponent as well
          ! 3D case
              if (position(3)>x0(3) .and. position(3)<x1(3)) then
                position_to_lgt_id=lgt_active(k)
                return
              endif
          else
          ! 2D case
            position_to_lgt_id=lgt_active(k)
            return
          endif ! end 3dcase
       endif


    enddo
     write(*,'("[bridgefluid.f90:] No block found for position=", f3.6," STOP!")')position(1)
     write(*,*),position
     call abort(272372)
end function position_to_lgt_id

!===========================================================================







!===========================================================================

subroutine interpolate_data(lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, lgt_active, lgt_n, hvy_n,params,positions,u_inter)
!! Send a message to the other world to signal the end of the communication session
!! Subroutine-declarations
 !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !>
    real(kind=rk), allocatable, intent(inout)   :: u_inter(:,:)

    real(kind=rk), intent(inout)        :: positions(:,:)

    type(type_params),                intent(in)    :: params           !< params of infile

!================================================================

    integer(kind=ik)                    ::k,Nr_particle,lgt_id
    ! block index
    integer(kind=ik)                    ::ibx, iby, ibz

    integer(kind=ik), allocatable       ::particle_id(:)

    real(kind=rk),  allocatable         ::u(:),v(:),w(:),rho(:),p(:)

    real(kind=rk), dimension(1:3)       :: x0, dx

!================================================================

    Nr_particle=size(positions,2)

    ! allocate state vector
    ! state vector consists of 6 components:
    !   - space components of velocity
    !   - preasure
    !   - density
    allocate(u_inter(6,Nr_particle))
    allocate(u(Nr_particle))
    allocate(v(Nr_particle))
    allocate(w(Nr_particle))
    allocate(p(Nr_particle))
    allocate(rho(Nr_particle))

    !allocate particle ids
    allocate(particle_id(Nr_particle))


    ! convert position to hvy_id

    !write(*,*) "dimensions=",size(hvy_block,1),size(hvy_block,2),size(hvy_block,3),size(hvy_block,4  ),size(hvy_block,5)
    do k=1,Nr_particle
        particle_id(k)   = position_to_lgt_id(lgt_block,lgt_active,lgt_n,positions(:,k),params)
        !! lgt_id
        lgt_id=particle_id(k)
        !! hvy_id
        particle_id(k)   = particle_id(k)-params%bridge%myWorldRank*params%number_blocks
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
        !!! calculate grid point ibx, iby, ibz
        ibx   = int( (positions(1,k)-x0(1))/dx(1) ) + params%number_ghost_nodes +1
        iby   = int( (positions(2,k)-x0(2))/dx(2) ) + params%number_ghost_nodes +1
        if (params%threeD_case) then
          ibz = int( (positions(3,k)-x0(3))/dx(3) ) + params%number_ghost_nodes +1
        else
          ibz = 1
        endif
        ! write(*,'("hvy_id =", i6, " x=", f6.3," [xmin,xmax]=[",f6.3,",",f6.3,"]")')particle_id(k),positions(1,k),x0(1), x0(1)+dx(1)*(params%number_block_nodes-1)
        ! remark: the first 3 dimensions of hvy_block project the grid structure into the array
        ! this means:

   !        * *  * * * * * * *  * *
   !        * *  * * * * * * *  * *
   !            --------------
   !        * *| + + + + + + +| * *
   !        * *| + + + + + + +| * *  -->
   !        * *| + + + + + + +| * *  -->  hvy(ibx,iby,ibz,...)
   !        * *| + + + + + + +| * *  -->
   !        * *| + + + + + + +| * *
   !        * *| + + + + + + +| * *      * ... ghost_node
   !        * *| + + + + + + +| * *      + ... inner grid point
   !            --------------           | ... indication of boundary
   !     ^  * *  * * * * * * *  * *
   ! iby |  * *  * * * * * * *  * *
   !        1 2  3 4 5 6 7 8 9  10 11
   !     -->
   !    ibx

  !interpolate:
  !> \todo write interpolation
      if (params%threeD_case) then
        if (params%number_data_fields /= 5) then
          call abort(333990,"[bridgefluid] number of data fields is less then 5, Stop")
        else
        ! write(*,*)"ibx =",ibx,"iby =",iby,"ibz =",ibz
          u_inter(:,k)   = hvy_block(ibx, iby, ibz,1:params%number_data_fields, particle_id(k) )
        endif
      else
          u_inter(1:3,k) = hvy_block(ibx, iby, ibz,1:3, particle_id(k) )
          u_inter(4,k)   = 0
          u_inter(5,k)   = hvy_block(ibx, iby, ibz,4, particle_id(k) )
      endif
          u_inter(6,k)   = positions(4,k)
          !write(*,'("Nr:",i6,"   rho=",f9.2,"  ux=",f9.2)') int(u_inter(6,k)),u_inter(1,k),u_inter(2,k)
    enddo
    deallocate(particle_id,u,v,w,p,rho)

 end subroutine interpolate_data

!===========================================================================

subroutine signalEndCalculation (myBridge)
!! Send a message to the other world to signal the end of the calculation session
!! Subroutine-declarations
type(bridgeMPI), intent(in)                     :: myBridge        ! type bridge on the particle side
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(1)                                    :: endCalculationMessage ! Message sent to the fluid world
integer                                         :: ierr            ! MPI communication error
integer                                         :: nProc           ! loop index (fluid process rank)

  ! Wait for all particle processes
  call MPI_barrier(myBridge%myWorld, ierr)

  ! Only the root particle process sent the message
  ! NOTE: one should try to implement this via an intercommunicator broadcast
  if (myBridge%myWorldRank == 0) then                              ! only the root process sends the message
    !! - assign the message to send (but this value is not important)
    endCalculationMessage = 'b'
    !! - send the message to the fluid processes
    do nProc = 0, myBridge%otherWorldSize - 1                      ! loop over all fluid processes
      call MPI_send(endCalculationMessage, 1, MPI_character, nProc, end_calculation, &
                    myBridge%otherWorld, ierr)
    end do                                                         ! end loop over all fluid processes
  end if                                                           ! end condition regarding the rank

end subroutine signalEndCalculation



end module module_bridge_interface
