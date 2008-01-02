!-*- F90 -*------------------------------------------------------------
!
!  module: mpicomms / meta
!
!  does the mpi communication part.
!
!----------------------------------------------------------------------

!======================================================================
!
! After H field is calculated, the tangential field components are
! being send to the right neighbour partition. As all fields are 
! allocated from MIN:MAX coordinates but filled from BEG:END, this means
! than the END layer goes to the right neighbours MIN layer.
! After E field is calculated, the tangential field components are 
! being send to the left neighbour, from the BEG layer to the MAX 
! layer.
!
! As dimensionality is reduced the communicated layer (3D) becomes a line
!(2D) and finally a point (1D).


module mpicomms

  use strings
  use constant
  use mpiworld
  use grid
  use fdtd
  implicit none
  private
  save 

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'MPICOMMS'

  ! --- Public Methods

  public :: InitializeMPIComms
  public :: FinalizeMPIComms
  public :: InitiateHMPIComms
  public :: InitiateEMPIComms
  public :: CompleteHMPIComms
  public :: CompleteEMPIComms
  
  ! --- Private Data
  
M4_IFELSE_MPI({

  integer, parameter :: xtag = 1
  integer, parameter :: ytag = 2
  integer :: ireqs(2), ireqr(2)
  integer :: mpitype 
  integer :: mpisize

  integer :: ides, isrc
  real(kind=8) :: tt
  real(kind=8) :: time_waited (2)
  time_waited = 0

})


contains

!----------------------------------------------------------------------

  subroutine InitializeMPIComms
  
    integer :: sz
    
M4_IFELSE_MPI({

    sz = M4_IFELSE_CF({2},{1})
  
    mpitype = MPI_REAL8

	select case ( M4_SDIM ) 
    case(1)
		mpisize=1*sz                            ! set 1D mpisize
    case(2)
		mpisize=(IMAX-IMIN+1)*sz                ! set 2D mpisize	
	case default
		mpisize=(IMAX-IMIN+1)*(JMAX-JMIN+1)*sz  ! set 3D mpisize
	end	select
     
})
  
     
  end subroutine InitializeMPIComms

  
!----------------------------------------------------------------------

  subroutine FinalizeMPIComms
  
M4_IFELSE_MPI({

        write(6,*) ranklbl, "* time for Comms H = ", time_waited(1)
        write(6,*) ranklbl, "* time for Comms E = ", time_waited(2)
})

  end subroutine FinalizeMPIComms

!----------------------------------------------------------------------

  subroutine InitiateHMPIComms

M4_IFELSE_MPI({

! 1. send tangential H fields to right

     if ( myrank .ne. numproc-1 ) then
        ides= myrank+1

   	    call StartMPELog(1)
        call MPI_ISEND( Hx(M4_MPICOORD({END})),mpisize,mpitype,ides, &
             xtag,mpicomm,ireqs(1),mpierr )
        call MPI_ISEND( Hy(M4_MPICOORD({END})),mpisize,mpitype,ides, &
             ytag,mpicomm,ireqs(2),mpierr )
   	    call StopMPELog(1)

     endif

! 2. receive tangential H fields from left

     if ( myrank .ne. 0 ) then
        isrc= myrank-1

   	    call StartMPELog(2)
        call MPI_IRECV( Hx(M4_MPICOORD({MIN})),mpisize,mpitype,isrc, &
             xtag,mpicomm,ireqr(1),mpierr )
        call MPI_IRECV( Hy(M4_MPICOORD({MIN})),mpisize,mpitype,isrc, &
             ytag,mpicomm,ireqr(2),mpierr )
   	    call StopMPELog(2)

     endif

})

  end subroutine InitiateHMPIComms
 
!----------------------------------------------------------------------

  subroutine CompleteHMPIComms

M4_IFELSE_MPI({

     tt = MPI_WTIME()

     if ( myrank .ne. 0 ) then

        call StartMPELog(3)
        call MPI_WAIT(ireqr(2),mpistatus, mpierr )
        call MPI_WAIT(ireqr(1),mpistatus, mpierr )
	    call StopMPELog(3)

     endif

     if ( myrank .ne. numproc-1 ) then

        call StartMPELog(3)
        call MPI_WAIT(ireqs(2),mpistatus, mpierr )
        call MPI_WAIT(ireqs(1),mpistatus, mpierr )
	    call StopMPELog(3)

     endif

     tt = MPI_WTIME() - tt
     time_waited(1) = time_waited(1) + tt

})

  end subroutine CompleteHMPIComms

!----------------------------------------------------------------------

  subroutine InitiateEMPIComms

M4_IFELSE_MPI({ 

! 1. send tangential E fields to the left

     if ( myrank .ne. 0 ) then
        ides = myrank-1

  	    call StartMPELog(1)
        call MPI_ISEND( Ex(M4_MPICOORD({BEG})),mpisize,mpitype,ides, &
             xtag,mpicomm,ireqs(1),mpierr )
        call MPI_ISEND( Ey(M4_MPICOORD({BEG})),mpisize,mpitype,ides, &
             ytag,mpicomm,ireqs(2),mpierr )
  	    call StopMPELog(1)

     endif

! 2. receive tangential E fields from right

     if ( myrank .ne. numproc-1 ) then
        isrc = myrank+1

  	    call StartMPELog(2)
        call MPI_IRECV( Ex(M4_MPICOORD({MAX})),mpisize,mpitype,isrc, &
             xtag,mpicomm,ireqr(1),mpierr )
        call MPI_IRECV( Ey(M4_MPICOORD({MAX})),mpisize,mpitype,isrc, &
             ytag,mpicomm,ireqr(2),mpierr )
  	    call StopMPELog(2)

     endif

})

  end subroutine InitiateEMPIComms
 
!----------------------------------------------------------------------

  subroutine CompleteEMPIComms

M4_IFELSE_MPI({

     tt = MPI_WTIME()

     if ( myrank .ne. numproc-1 ) then

        call StartMPELog(3)
        call MPI_WAIT(ireqr(2),mpistatus, mpierr )
        call MPI_WAIT(ireqr(1),mpistatus, mpierr )
        call StopMPELog(3)

     endif

     if ( myrank .ne. 0 ) then

        call StartMPELog(3)
        call MPI_WAIT(ireqs(2),mpistatus, mpierr )
        call MPI_WAIT(ireqs(1),mpistatus, mpierr )
        call StopMPELog(3)

     endif

     tt = MPI_WTIME() - tt
     time_waited(2) = time_waited(2) + tt

})

  end subroutine CompleteEMPIComms


!----------------------------------------------------------------------

end module mpicomms

!
! Authors:  J.Hamm
! Modified: 20/12/2007
!
!======================================================================









