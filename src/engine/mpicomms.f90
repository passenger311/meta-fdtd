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
! Note: if dimensionality is reduced the communicated layer (3D) 
! becomes a line (2D) and finally a point (1D).
!
! Note: if pbcs are set, the leftmost and rightmost partition need
! to exchange data.
!
! Note: do not put non-mpi code here as all subroutines will be skipped
! if run without mpi. All in here is purely related to mpi communication.


module mpicomms

  use strings
  use constant
  use mpiworld
  use grid
  use fdtd
  use pbc
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
  logical :: lperiodic, rperiodic

  integer :: idesH, isrcH, idesE, isrcE
  real(kind=8) :: tt
  real(kind=8) :: time_waited (2)
  time_waited = 0

})


contains

!----------------------------------------------------------------------

  subroutine InitializeMPIComms
  
    integer :: sz
    
    mpibound = .false.

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

! do the leftmost or the rightmost partition have a periodic boundary?

    lperiodic = .false.
    rperiodic = .false.

    if ( planepbc(MPI_SDIM*2-1) .eq. 1 ) lperiodic = .true.
    if ( planepbc(MPI_SDIM*2) .eq. 1 ) rperiodic = .true.

! source and destination mpi coordinates

    idesH = myrank + 1
    isrcH = myrank - 1

    idesE = myrank - 1
    isrcE = myrank + 1

    if ( myrank .eq. 0 ) then
       isrcH = MPI_PROC_NULL ! no left neighbour
       idesE = MPI_PROC_NULL
       if ( lperiodic ) isrcH = numproc - 1 ! except if pbc
       if ( rperiodic ) idesE = numproc - 1
    endif

    if ( myrank .eq. numproc - 1 ) then
       isrcE = MPI_PROC_NULL ! no right neighbour
       idesH = MPI_PROC_NULL
       if ( lperiodic ) idesH = 0 ! except if pbc
       if ( rperiodic ) isrcE = 0
    endif

! mpibound is true if the boundary is updated by comms -> bound.f90 has to know!

    mpibound(M4_SDIM*2-1) = .true.
    mpibound(M4_SDIM*2) = .true.

    if ( myrank .eq. 0 .and. .not. lperiodic ) mpibound(M4_SDIM*2-1) = .false.
    if ( myrank .eq. numproc-1 .and. .not. rperiodic ) mpibound(M4_SDIM*2) = .false.

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

    if ( numproc .eq. 1 ) return

M4_IFELSE_MPI({

! 1. send tangential H fields to right

      call StartMPELog(1)
      call MPI_ISEND( Hx(M4_MPICOORD({END})),mpisize,mpitype,idesH, &
           xtag,mpicomm,ireqs(1),mpierr )
      call MPI_ISEND( Hy(M4_MPICOORD({END})),mpisize,mpitype,idesH, &
           ytag,mpicomm,ireqs(2),mpierr )
      call StopMPELog(1)

! 2. receive tangential H fields from left

      call StartMPELog(2)
      call MPI_IRECV( Hx(M4_MPICOORD({MIN})),mpisize,mpitype,isrcH, &
           xtag,mpicomm,ireqr(1),mpierr )
      call MPI_IRECV( Hy(M4_MPICOORD({MIN})),mpisize,mpitype,isrcH, &
           ytag,mpicomm,ireqr(2),mpierr )
      call StopMPELog(2)
       

})

  end subroutine InitiateHMPIComms
 
!----------------------------------------------------------------------

  subroutine CompleteHMPIComms

    if ( numproc .eq. 1 ) return

M4_IFELSE_MPI({

    tt = MPI_WTIME()

! 1. wait till receive completed 

    call StartMPELog(3)
    call MPI_WAIT(ireqr(2),mpistatus, mpierr )
    call MPI_WAIT(ireqr(1),mpistatus, mpierr )
    call StopMPELog(3)
    
! 2. wait till send completed 

    call StartMPELog(3)
    call MPI_WAIT(ireqs(2),mpistatus, mpierr )
    call MPI_WAIT(ireqs(1),mpistatus, mpierr )
    call StopMPELog(3)
     
    tt = MPI_WTIME() - tt
    time_waited(1) = time_waited(1) + tt

})

  end subroutine CompleteHMPIComms

!----------------------------------------------------------------------

  subroutine InitiateEMPIComms

    if ( numproc .eq. 1 ) return

M4_IFELSE_MPI({ 

! 1. send tangential E fields to the left

    call StartMPELog(1)
    call MPI_ISEND( Ex(M4_MPICOORD({BEG})),mpisize,mpitype,idesE, &
         xtag,mpicomm,ireqs(1),mpierr )
    call MPI_ISEND( Ey(M4_MPICOORD({BEG})),mpisize,mpitype,idesE, &
         ytag,mpicomm,ireqs(2),mpierr )
    call StopMPELog(1)
       
! 2. receive tangential E fields from right

    call StartMPELog(2)
    call MPI_IRECV( Ex(M4_MPICOORD({MAX})),mpisize,mpitype,isrcE, &
         xtag,mpicomm,ireqr(1),mpierr )
    call MPI_IRECV( Ey(M4_MPICOORD({MAX})),mpisize,mpitype,isrcE, &
         ytag,mpicomm,ireqr(2),mpierr )
    call StopMPELog(2)
    

})

  end subroutine InitiateEMPIComms
 
!----------------------------------------------------------------------

  subroutine CompleteEMPIComms

    if ( numproc .eq. 1 ) return

M4_IFELSE_MPI({

    tt = MPI_WTIME()

! 1. wait till receive completed 

    call StartMPELog(3)
    call MPI_WAIT(ireqr(2),mpistatus, mpierr )
    call MPI_WAIT(ireqr(1),mpistatus, mpierr )
    call StopMPELog(3)
       
! 2. wait till send completed 

    call StartMPELog(3)
    call MPI_WAIT(ireqs(2),mpistatus, mpierr )
    call MPI_WAIT(ireqs(1),mpistatus, mpierr )
    call StopMPELog(3)
    
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









