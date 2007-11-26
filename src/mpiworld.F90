!-*- F90 -*------------------------------------------------------------
!
!  module: mpiworld / max3d
!
!  create mpi world 
!
!  subs:
!
!  InitializeMPIWorld
!  FinalizeMPIWorld
!
!----------------------------------------------------------------------


module mpiworld

  use strings
  implicit none
  save 

#ifdef MPI
  include 'mpif.h'
#endif /* MPI */

  integer :: numproc = 1
  integer :: myrank = 0
  integer :: mpi_started = 0
  character(len=255) :: mpi_sfxin, mpi_sfxout
  character(len=20)   :: ranklbl
  integer :: mpisize

#ifdef MPI
  integer :: mpicomm
  integer :: mpierr
  integer :: mpistatus(MPI_STATUS_SIZE) 
  integer, parameter :: xtag = 1
  integer, parameter :: ytag = 2
  integer :: ireqs(2), ireqr(2)
  integer :: mpitype 
#endif /* MPI */


contains

  subroutine InitializeMPIWorld

    implicit none      

#ifdef MPI

    mpicomm = MPI_COMM_WORLD
    mpitype = MPI_REAL8

    call MPI_INIT(mpierr)
    
    if ( mpierr .eq. 0 ) call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, mpierr)
    if ( mpierr .eq. 0 ) call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mpierr)

#endif /* MPI */

    mpi_sfxin = cat3('.',i2str(myrank),'.in')
    mpi_sfxout = cat3('.',i2str(myrank),'.out')
    ranklbl = cat5('(',i2str(myrank),'/',i2str(numproc),')')

#ifdef MPI

    if ( mpierr .ne. 0 ) then
       write(6,*) "@ mpistart: failed!"
       stop
    end if
    
    mpi_started = 1

#endif /* MPI */


  end subroutine InitializeMPIWorld
 
  
  subroutine FinalizeMPIWorld

    implicit none      

#ifdef MPI

    call MPI_FINALIZE(mpierr)

#endif /* MPI */

  end subroutine FinalizeMPIWorld


end module mpiworld











