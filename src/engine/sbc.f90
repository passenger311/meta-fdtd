!-*- F90 -*------------------------------------------------------------
!
!  module: sbc / meta
!
!  symmetric boundary conditions.
!
!----------------------------------------------------------------------
 
!======================================================================
!
!

module sbc
 
  use constant
  use mpiworld
  use grid
  use fdtd

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'SBC'

  ! --- Public Methods

  public :: InitializeSbc
  public :: FinalizeSbc
  public :: StepHBoundSbc
  public :: StepEBoundSbc

  ! --- Public Data

  ! --- Constants

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializeSbc(planebound,num)

    integer :: planebound(6), num
    integer :: i

    M4_WRITE_DBG(". enter InitializeSbc")

    M4_WRITE_DBG({"got planebound(i): ",  (planebound(i),i=1, M4_SDIM*2)})

    M4_WRITE_DBG(". exit InitializeSbc")

  end subroutine InitializeSbc

!----------------------------------------------------------------------

  subroutine FinalizeSbc

    M4_WRITE_DBG(". enter FinalizeSbc")
    M4_WRITE_DBG(". exit FinalizeSbc")

  end subroutine FinalizeSbc

!----------------------------------------------------------------------

  subroutine StepHBoundSbc(i)

    integer :: i

    select case (i)

    case ( 1 )
!       Hz(M4_COORD(IMIN,JBEG:JMAX,KBEG:KMAX))=Hz(M4_COORD(IBEG,JBEG:JMAX,KBEG:KMAX))
!       Hy(M4_COORD(IMIN,JBEG:JMAX,KBEG:KMAX))=Hy(M4_COORD(IBEG,JBEG:JMAX,KBEG:KMAX))
    case ( 2 )
!       Hz(M4_COORD(IMAX,JBEG:JMAX,KBEG:KMAX))=Hz(M4_COORD(IMAX-1,JBEG:JMAX,KBEG:KMAX))
!       Hy(M4_COORD(IMAX,JBEG:JMAX,KBEG:KMAX))=Hy(M4_COORD(IMAX-1,JBEG:JMAX,KBEG:KMAX))
    case ( 3 )
!       Hz(M4_COORD(IBEG:IMAX,JMIN,KBEG:KMAX))=Hz(M4_COORD(IBEG:IMAX,JBEG,KBEG:KMAX))
!       Hx(M4_COORD(IBEG:IMAX,JMIN,KBEG:KMAX))=Hx(M4_COORD(IBEG:IMAX,JBEG,KBEG:KMAX))
    case ( 4 )
!       Hz(M4_COORD(IBEG:IMAX,JMAX,KBEG:KMAX))=Hz(M4_COORD(IBEG:IMAX,JMAX-1,KBEG:KMAX))
!       Hx(M4_COORD(IBEG:IMAX,JMAX,KBEG:KMAX))=Hx(M4_COORD(IBEG:IMAX,JMAX-1,KBEG:KMAX))
    case ( 5 )
!       Hy(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMIN))=Hy(M4_COORD(IBEG:IMAX,JBEG:JMAX,KBEG))
!       Hx(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMIN))=Hx(M4_COORD(IBEG:IMAX,JBEG:JMAX,KBEG))
    case ( 6 )
!       Hy(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMAX))=Hy(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMAX-1))
!       Hx(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMAX))=Hx(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMAX-1))
       
    end select

  end subroutine StepHBoundSbc

!----------------------------------------------------------------------

  subroutine StepEBoundSbc(i)

    integer :: i

    select case (i)

    case ( 1 )
       Ez(M4_COORD(IBEG,JBEG:JMAX,KBEG:KMAX))=Ez(M4_COORD(IBEG+1,JBEG:JMAX,KBEG:KMAX))
       Ey(M4_COORD(IBEG,JBEG:JMAX,KBEG:KMAX))=Ey(M4_COORD(IBEG+1,JBEG:JMAX,KBEG:KMAX))
    case ( 2 )
       Ez(M4_COORD(IMAX,JBEG:JMAX,KBEG:KMAX))=Ez(M4_COORD(IMAX-1,JBEG:JMAX,KBEG:KMAX))
       Ey(M4_COORD(IMAX,JBEG:JMAX,KBEG:KMAX))=Ey(M4_COORD(IMAX-1,JBEG:JMAX,KBEG:KMAX))
    case ( 3 )
       Ez(M4_COORD(IBEG:IMAX,JBEG,KBEG:KMAX))=Ez(M4_COORD(IBEG:IMAX,JBEG+1,KBEG:KMAX))
       Ex(M4_COORD(IBEG:IMAX,JBEG,KBEG:KMAX))=Ex(M4_COORD(IBEG:IMAX,JBEG+1,KBEG:KMAX))
    case ( 4 )
       Ez(M4_COORD(IBEG:IMAX,JMAX,KBEG:KMAX))=Ez(M4_COORD(IBEG:IMAX,JMAX-1,KBEG:KMAX))
       Ex(M4_COORD(IBEG:IMAX,JMAX,KBEG:KMAX))=Ex(M4_COORD(IBEG:IMAX,JMAX-1,KBEG:KMAX))
    case ( 5 )
       Ey(M4_COORD(IBEG:IMAX,JBEG:JMAX,KBEG))=Ey(M4_COORD(IBEG:IMAX,JBEG:JMAX,KBEG+1))
       Ex(M4_COORD(IBEG:IMAX,JBEG:JMAX,KBEG))=Ex(M4_COORD(IBEG:IMAX,JBEG:JMAX,KBEG+1))
    case ( 6 )
       Ey(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMAX))=Ey(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMAX-1))
       Ex(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMAX))=Ex(M4_COORD(IBEG:IMAX,JBEG:JMAX,KMAX-1))

       end select

  end subroutine StepEBoundSbc

!----------------------------------------------------------------------

end module sbc

!
! Authors:  J.Hamm
! Modified: 24/1/2008
!
!======================================================================
