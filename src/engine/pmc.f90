!-*- F90 -*------------------------------------------------------------
!
!  module: pmc / meta
!
!  symmetric boundary conditions.
!
!----------------------------------------------------------------------
 
!======================================================================
!
!

module pmc
 
  use constant
  use mpiworld
  use grid
  use fdtd

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'PMC'

  ! --- Public Methods

  public :: InitializePmc
  public :: FinalizePmc
  public :: StepHBoundPmc
  public :: StepEBoundPmc

  ! --- Public Data

  ! --- Constants

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializePmc(planebound,num)

    integer :: planebound(6), num
    integer :: i

    M4_WRITE_DBG(". enter InitializePmc")

    M4_WRITE_DBG({"got planebound(i): ",  (planebound(i),i=1, M4_SDIM*2)})

    M4_WRITE_DBG(". exit InitializePmc")

  end subroutine InitializePmc

!----------------------------------------------------------------------

  subroutine FinalizePmc

    M4_WRITE_DBG(". enter FinalizePmc")
    M4_WRITE_DBG(". exit FinalizePmc")

  end subroutine FinalizePmc

!----------------------------------------------------------------------

  subroutine StepHBoundPmc(i)

    integer :: i

    select case (i)

       case ( 1 )
          Hz(M4_COORD(IMIN,JBEG:JMAX,KBEG:KMAX))=0.0 
          Hy(M4_COORD(IMIN,JBEG:JMAX,KBEG:KMAX))=0.0 
       case ( 2 )
          Hz(M4_COORD(IEND,JBEG:JMAX,KBEG:KMAX))=0.0
          Hy(M4_COORD(IEND,JBEG:JMAX,KBEG:KMAX))=0.0
       case ( 3 )
          Hz(M4_COORD(IBEG:IMAX,JMIN,KBEG:KMAX))=0.0
          Hx(M4_COORD(IBEG:IMAX,JMIN,KBEG:KMAX))=0.0 
       case ( 4 )
          Hz(M4_COORD(IBEG:IMAX,JEND,KBEG:KMAX))=0.0
          Hx(M4_COORD(IBEG:IMAX,JEND,KBEG:KMAX))=0.0 
       case ( 5 )
          Hy(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMIN))=0.0
          Hx(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMIN))=0.0
       case ( 6 )
          Hy(M4_COORD(IBEG:IMAX,JBEG:JMAX,KEND))=0.0
          Hx(M4_COORD(IBEG:IMAX,JBEG:JMAX,KEND))=0.0

       end select

     end subroutine StepHBoundPmc

!----------------------------------------------------------------------

  subroutine StepEBoundPmc(i)

    integer :: i

  end subroutine StepEBoundPmc


!----------------------------------------------------------------------

end module pmc

!
! Authors:  J.Hamm
! Modified: 24/1/2008
!
!======================================================================
