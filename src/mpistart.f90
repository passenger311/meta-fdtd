!----------------------------------------------------------------------
!
!  module: mpistart / max3d
!
!  initialise and finalise mpi. 
!
!----------------------------------------------------------------------


module mpistart

  use strings
  implicit none
  save 

  include 'mpif.h'

  integer :: numproc
  integer :: myrank
  integer :: mpicomm
  integer :: mpierr
  integer :: mpistatus(MPI_STATUS_SIZE) 
  integer :: mpi_started = 0
  integer, parameter :: xtag = 1
  integer, parameter :: ytag = 2
  integer :: ireqs(2), ireqr(2)
  integer :: mpisize                     ! elements of 8 bytes (double)
  integer :: mpitype 
  character(len=255) :: mpi_sfxin, mpi_sfxout
  character(len=20)   :: ranklbl


contains

  subroutine mpiinit

    implicit none      

    mpicomm = MPI_COMM_WORLD
    mpitype = MPI_REAL8

    call MPI_INIT(mpierr)
    
    if ( mpierr .eq. 0 ) call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, mpierr)
    if ( mpierr .eq. 0 ) call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mpierr)

    mpi_sfxin = cat3('.',i2str(myrank),'.in')
    mpi_sfxout = cat3('.',i2str(myrank),'.out')
    ranklbl = cat5('(',i2str(myrank),'/',i2str(numproc),')')

    if ( mpierr .ne. 0 ) then
       write(6,*) "@ mpistart: failed!"
       stop
    end if
    
    mpi_started = 1

    return

  end subroutine mpiinit
 
  
  subroutine mpifinalise

    implicit none      


    call MPI_FINALIZE(mpierr)


  end subroutine mpifinalise


end module mpistart












