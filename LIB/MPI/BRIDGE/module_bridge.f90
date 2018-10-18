!-------------------------------------------------------------!
!> This module contains the basic tools to create 2 MPI worlds
!-------------------------------------------------------------!
module module_bridge

use MPI



implicit none
! private

! Declarations

!! Variables to define a MPI process in a sub-world
!! - new type bridge which contains useful MPI parameters
type, public           :: bridgeMPI
  integer                :: commonWorld                            ! INTRA communicator in the common world
  integer                :: commonWorldSize                        ! size of the common world
  integer                :: commonWorldRank                        ! rank in the common world
  integer                :: myWorld                                ! INTRA communicator in the local world
  integer                :: myWorldSize                            ! size of the local world
  integer                :: myWorldRank                            ! rank in the local world
  integer                :: otherWorld                             ! INTER communicator of the other world
  integer                :: otherWorldSize                         ! size of the other world
  integer                :: minMyWorldRank                         ! minimum common rank in the local world
  integer                :: maxMyWorldRank                         ! maximum common rank in the local world
  integer                :: minOtherWorldRank                      ! minimum common rank in the other world
  integer                :: maxOtherWorldRank                      ! maximum common rank in the other world
end type                   bridgeMPI

!! MPI tag nomenclature
!! - exchange of contants values (communicators, initial parameters, etc.)
integer, parameter, public :: communicator_exchange = 10           ! tag for exchange of communicators
integer, parameter, public :: parameters_request    = 11           ! tag for the request of parameters
integer, parameter, public :: parameters_delivery   = 12           ! tag for the delivery of parameters
!! - exchange of calculation data
integer, parameter, public :: data_positionDelivery = 20           ! tag for the delivery of position data
integer, parameter, public :: data_fieldDelivery    = 21           ! tag for the delivery of field values
!! - exchange of state message
integer, parameter, public :: end_communication     = 30           ! tag to end the communication tunnel
integer, parameter, public :: end_calculation       = 31           ! tag to end the calculation

!! Useful subroutines
public :: createMPIWorlds
public :: createMPIWorldsMaster
public :: createMPIWorldsSlave
public :: getMinMaxCommonWorldRanks

contains



!> \brief      Creates two mpi worlds by splitting the processes of the common
!>             communicator in two subsets
!
!> \param      bridgeToSet  bridge which contains useful MPI parameters like
!>                          commonWorld(,rank,size) or myWorld(,rank,size)
!> \param      worldIndex world index which is a unique tag for the two
!>                          groups (integer)
!
!> \return      bridgeToSet
subroutine createMPIWorlds (bridgeToSet, worldIndex)
!! distinguish between the sub-worlds by splitting the processes
!! Subroutine-declarations
type(bridgeMPI), intent(out)             :: bridgeToSet            ! type bridge to allocate in a world
integer        , intent(in)              :: worldIndex             ! index to chose the world in which the bridge is allocated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer                                  :: ierr                   ! MPI communication error
integer                                  :: commonComm             ! MPI common communicator

  ! split the MPI world according to the given worldIndex
  !! - set the communicator in the common world
  call MPI_comm_dup(MPI_Comm_world, commonComm, ierr)
  bridgeToSet%commonWorld = commonComm

  !! - determine the size and rank in the common world
  call MPI_comm_size(bridgeToSet%commonWorld, bridgeToSet%commonWorldSize, ierr)
  call MPI_comm_rank(bridgeToSet%commonWorld, bridgeToSet%commonWorldRank, ierr)


  !! - according to the value of worldIndex, split the MPI_Comm_world
  call MPI_comm_split(bridgeToSet%commonWorld, worldIndex, bridgeToSet%commonWorldRank, &
                      bridgeToSet%myWorld, ierr)

  !! - determine the size and rank in the local world
  call MPI_comm_size(bridgeToSet%myWorld, bridgeToSet%myWorldSize, ierr)
  call MPI_comm_rank(bridgeToSet%myWorld, bridgeToSet%myWorldRank, ierr)

  !! determine the size of the other world
  bridgeToSet%otherWorldSize = bridgeToSet%commonWorldSize - bridgeToSet%myWorldSize

  !! set the min and max common rank of the two worlds
  call getMinMaxCommonWorldRanks(bridgeToSet)

  !! get the communicator of the other world
  write(*,*) ' '
  write(*,*) 'Fluid side: total world size: ', bridgeToSet%commonWorldSize
  write(*,*) 'Fluid side: total world rank: ', bridgeToSet%commonWorldRank
  write(*,*) 'Fluid side: my    world size: ', bridgeToSet%myWorldSize
  write(*,*) 'Fluid side: my    world rank: ', bridgeToSet%myWorldRank
  write(*,*) 'Fluid side: other world min : ', bridgeToSet%minOtherWorldRank
  write(*,*) 'Fluid side: other world max : ', bridgeToSet%maxOtherWorldRank
  call MPI_intercomm_create(bridgeToSet%myWorld, 0, bridgeToSet%commonWorld, bridgeToSet%minOtherWorldRank, &
                            communicator_exchange, bridgeToSet%otherWorld, ierr)
end subroutine createMPIWorlds



!-------------------------------------------------------------------------!!!!
!> \brief      create a new set of processes (slave program)
!>             called by a master program - master side
!>             Subroutine-declarations
!
!> \param      bridgeToSet   type bridge to allocate in a world
!> \param      slaveSize     index to chose the world in which the bridge is allocated
!> \param      slaveCommand  command line of the slave program
!
!> \return     bridgeToSet
subroutine createMPIWorldsMaster (bridgeToSet, slaveSize, slaveCommand)


type(bridgeMPI)  , intent(out)           :: bridgeToSet            ! type bridge to allocate in a world
integer          , intent(in)            :: slaveSize              ! index to chose the world in which the bridge is allocated
character(len=*) , intent(in)            :: slaveCommand           ! command line of the slave program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer                                  :: ierr                   ! MPI communication error
integer, dimension(slaveSize)            :: ierrs                  ! MPI communication error array

  ! As master, create a new set of processes for the slave side
  !! - assign my world to the master intracommunicator
  call MPI_comm_dup(MPI_comm_world, bridgeToSet%myWorld, ierr)

  !! - according to the value of slaveSize, create a new set of processes
  call MPI_comm_spawn(slaveCommand, MPI_argv_null, slaveSize, MPI_info_null, 0, &
                      bridgeToSet%myWorld, bridgeToSet%otherWorld, ierrs, ierr)

  !! - determine the size and rank in the local world
  call MPI_comm_size(bridgeToSet%myWorld, bridgeToSet%myWorldSize, ierr)
  call MPI_comm_rank(bridgeToSet%myWorld, bridgeToSet%myWorldRank, ierr)

  !! - determine the size of the other world
  call MPI_comm_remote_size(bridgeToSet%otherWorld, bridgeToSet%otherWorldSize, ierr)

  !! - define the new common world (master get the high ranks in the common world)
  call MPI_intercomm_merge(bridgeToSet%otherWorld, .true., bridgeToSet%commonWorld, ierr)

  !! - Determine the ranks/size values in the common world
  call MPI_comm_size(bridgeToSet%commonWorld, bridgeToSet%commonWorldSize, ierr)
  call MPI_comm_rank(bridgeToSet%commonWorld, bridgeToSet%commonWorldRank, ierr)
  call getMinMaxCommonWorldRanks(bridgeToSet)

end subroutine createMPIWorldsMaster

!!!!-------------------------------------------------------------------------!!!!

subroutine createMPIWorldsSlave (bridgeToSet)
!! create a new set of processes (slave program) called by a master program - slave side
!! Subroutine-declarations
type(bridgeMPI), intent(out)             :: bridgeToSet            ! type bridge to allocate in a world
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer                                  :: ierr                   ! MPI communication error

  ! As slave, determine the master who created the processes and assign the resulting values
  !! - assign my world to the master intracommunicator
  call MPI_comm_dup(MPI_comm_world, bridgeToSet%myWorld, ierr)

  !! - determine the size and rank in the local world
  call MPI_comm_size(bridgeToSet%myWorld, bridgeToSet%myWorldSize, ierr)
  call MPI_comm_rank(bridgeToSet%myWorld, bridgeToSet%myWorldRank, ierr)

  !! - get the intercommunicator of the other world
  call MPI_comm_get_parent(bridgeToSet%otherWorld, ierr)

  !! - determine the size of the other world
  call MPI_comm_remote_size(bridgeToSet%otherWorld, bridgeToSet%otherWorldSize, ierr)

  !! - define the new common world (master get the high ranks in the common world)
  call MPI_intercomm_merge(bridgeToSet%otherWorld, .false., bridgeToSet%commonWorld, ierr)

  !! - Determine the ranks/size values in the common world
  call MPI_comm_size(bridgeToSet%commonWorld, bridgeToSet%commonWorldSize, ierr)
  call MPI_comm_rank(bridgeToSet%commonWorld, bridgeToSet%commonWorldRank, ierr)
  call getMinMaxCommonWorldRanks(bridgeToSet)

end subroutine createMPIWorldsSlave

!!!!-------------------------------------------------------------------------!!!!

subroutine getMinMaxCommonWorldRanks (myBridge)
!! compute the minimum and maximum ranks in the common world of each side
!! Subroutine-declarations
type(bridgeMPI), intent(inout)           :: myBridge               ! type bridge on the local side
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer                                  :: ierr                   ! MPI communication error

  ! Local world
  !! get the minimum rank
  call MPI_allreduce(myBridge%commonWorldRank, myBridge%minMyWorldRank, 1, &
                     MPI_integer, MPI_min, myBridge%myWorld, ierr)
  !! get the maximal rank
  call MPI_allreduce(myBridge%commonWorldRank, myBridge%maxMyWorldRank, 1, &
                     MPI_integer, MPI_max, myBridge%myWorld, ierr)

  ! Fluid world and check
  !!  if all local processes have the first ranks in the common world
  if ((myBridge%minMyWorldRank.eq.0).and. &
      (myBridge%maxMyWorldRank.eq.(myBridge%myWorldSize - 1))) then
    myBridge%minOtherWorldRank = myBridge%myWorldSize
    myBridge%maxOtherWorldRank = myBridge%commonWorldSize-1
  !! if all local processes have the last ranks in the common world
  else if ((myBridge%minMyWorldRank.eq.(myBridge%commonWorldSize - myBridge%myWorldSize)).and. &
           (myBridge%maxMyWorldRank.eq.(myBridge%commonWorldSize - 1))) then
    myBridge%minOtherWorldRank = 0
    myBridge%maxOtherWorldRank = myBridge%commonWorldSize - myBridge%myWorldSize - 1
  !! otherwise one doesn't know the process numbers
  else
    if (myBridge%myWorldRank.eq.0) then                            ! only display the error message if this is the root local process
      write(*,*) 'Unable to compute the ranks of both sides in the common world. Stop!'
      call MPI_abort(myBridge%commonWorld, 99, ierr)               ! Abort the execution
    end if                                                         ! end condition on the rank in the local world
  end if

end subroutine getMinMaxCommonWorldRanks

end module
