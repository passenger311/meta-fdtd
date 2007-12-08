!-*- F90 -*------------------------------------------------------------
!
!  module: mpiworld / meta3
!
!  create mpi world 
!
!  subs:
!
!  InitializeMPIWorld
!  FinalizeMPIWorld
!
!----------------------------------------------------------------------

!======================================================================
!
!


module mpiworld

  use strings
  use constant
  implicit none
  save 

M4_IFELSE_MPI(`include 'mpif.h'',`')

  integer :: numproc = 1
  integer :: myrank = 0
  integer :: mpi_started = 0
  character(len=STRLNG) :: mpi_sfxin, mpi_sfxout
  character(len=20)   :: ranklbl
  integer :: mpisize

M4_IFELSE_MPI(`
  integer :: mpicomm
  integer :: mpierr
  integer :: mpistatus(MPI_STATUS_SIZE) 
  integer, parameter :: xtag = 1
  integer, parameter :: ytag = 2
  integer :: ireqs(2), ireqr(2)
  integer :: mpitype 
',`')


contains

!----------------------------------------------------------------------

  subroutine InitializeMPIWorld

M4_IFELSE_MPI(`

    mpicomm = MPI_COMM_WORLD
    mpitype = MPI_REAL8

    call MPI_INIT(mpierr)
    
    if ( mpierr .eq. 0 ) call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, mpierr)
    if ( mpierr .eq. 0 ) call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mpierr)

',`')

    mpi_sfxin = '.' // TRIM(i2str(myrank)) // '.in'
    mpi_sfxout = '.' // TRIM(i2str(myrank)) // '.out'
    ranklbl = '(' // TRIM(i2str(myrank)) // '/' // TRIM(i2str(numproc)) // ')'

M4_IFELSE_MPI(`

    if ( mpierr .ne. 0 ) then
       write(6,*) "@ mpistart: failed!"
       stop
    end if
    
    mpi_started = 1

',`')

  end subroutine InitializeMPIWorld
 
!----------------------------------------------------------------------

  subroutine FinalizeMPIWorld

M4_IFELSE_MPI(`call MPI_FINALIZE(mpierr)',`')

  end subroutine FinalizeMPIWorld

!----------------------------------------------------------------------

end module mpiworld

!
! Authors:  J.Hamm, A.Klaedtke
! Modified: 4/12/2007
!
!======================================================================









