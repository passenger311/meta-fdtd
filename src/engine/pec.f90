!-*- F90 -*------------------------------------------------------------
!
!  module: pec / meta
!
!  set metallic boundary conditions.
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------
 
!======================================================================
!
!

module pec
 
  use constant
  use mpiworld
  use grid
  use fdtd

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'PEC'

  ! --- Public Methods

  public :: InitializePec
  public :: FinalizePec
  public :: StepHBoundPec
  public :: StepEBoundPec

  ! --- Public Data

  ! --- Constants

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializePec(planebound,num)

    integer :: planebound(6), num
    integer :: i

    M4_WRITE_DBG(". enter InitializePec")

    M4_WRITE_DBG({"got planebound(i): ",  (planebound(i),i=1, M4_SDIM*2)})

    M4_WRITE_DBG(". exit InitializePec")

  end subroutine InitializePec

!----------------------------------------------------------------------

  subroutine FinalizePec

    M4_WRITE_DBG(". enter FinalizePec")
    M4_WRITE_DBG(". exit FinalizePec")

  end subroutine FinalizePec

!----------------------------------------------------------------------

  subroutine StepHBoundPec(i)

    integer :: i

  end subroutine StepHBoundPec

!----------------------------------------------------------------------

  subroutine StepEBoundPec(i)

    integer :: i

    select case (i)

       case ( 1 )
          Ez(M4_COORD(IBEG,JBEG:JMAX,KBEG:KMAX))=0.0
          Ey(M4_COORD(IBEG,JBEG:JMAX,KBEG:KMAX))=0.0
       case ( 2 )
          Ez(M4_COORD(IMAX,JBEG:JMAX,KBEG:KMAX))=0.0
          Ey(M4_COORD(IMAX,JBEG:JMAX,KBEG:KMAX))=0.0
       case ( 3 )
          Ez(M4_COORD(IBEG:IMAX,JBEG,KBEG:KMAX))=0.0
          Ex(M4_COORD(IBEG:IMAX,JBEG,KBEG:KMAX))=0.0
       case ( 4 )
          Ez(M4_COORD(IBEG:IMAX,JMAX,KBEG:KMAX))=0.0
          Ex(M4_COORD(IBEG:IMAX,JMAX,KBEG:KMAX))=0.0
       case ( 5 )
          Ey(M4_COORD(IBEG:IMAX,JBEG:JMAX,KBEG))=0.0
          Ex(M4_COORD(IBEG:IMAX,JBEG:JMAX,KBEG))=0.0
       case ( 6 )
          Ey(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMAX))=0.0
          Ex(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMAX))=0.0

       end select

  end subroutine StepEBoundPec

!----------------------------------------------------------------------

end module pec

!
! Authors:  J.Hamm, S.Scholz, A.Klaedtke
! Modified: 4/12/2007
!
!======================================================================
