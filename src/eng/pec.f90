!-*- F90 -*------------------------------------------------------------
!
!  module: pec / meta3
!
!  set metallic boundary conditions.
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

    M4_WRITE_DBG({"got planebound(i): ",  (planebound(i),i=1, M4_DIM*2)})

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
          Ez(IBEG,JBEG:JMAX,KBEG:KMAX)=0.0
          Ey(IBEG,JBEG:JMAX,KBEG:KMAX)=0.0
       case ( 2 )
          Ez(IMAX,JBEG:JMAX,KBEG:KMAX)=0.0
          Ey(IMAX,JBEG:JMAX,KBEG:KMAX)=0.0
       case ( 3 )
          Ez(IBEG:IMAX,JBEG,KBEG:KMAX)=0.0
          Ex(IBEG:IMAX,JBEG,KBEG:KMAX)=0.0
       case ( 4 )
          Ez(IBEG:IMAX,JMAX,KBEG:KMAX)=0.0
          Ex(IBEG:IMAX,JMAX,KBEG:KMAX)=0.0
       case ( 5 )
          Ey(IBEG:IMAX,JBEG:JMAX,KBEG)=0.0
          Ex(IBEG:IMAX,JBEG:JMAX,KBEG)=0.0
       case ( 6 )
          Ey(IBEG:IMAX,JBEG:JMAX,KMAX)=0.0
          Ex(IBEG:IMAX,JBEG:JMAX,KMAX)=0.0

       end select

  end subroutine StepEBoundPec

!----------------------------------------------------------------------

end module pec

!
! Authors:  J.Hamm, S.Scholz, A.Klaedtke
! Modified: 4/12/2007
!
!======================================================================
