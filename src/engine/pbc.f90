!-*- F90 -*------------------------------------------------------------
!
!  module: pbc / meta
!
!  periodic boundary conditions.
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------

!======================================================================
!
! Note: if numprocs > 1, part of the pbc needs to be handled by 
! mpicomm.
!

module pbc

  use constant
  use strings
  use grid  
  use fdtd

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'PBC'
  logical, private :: modconfigured = .false.

  ! --- Public Methods

  public :: ReadConfigPbc
  public :: InitializePbc
  public :: FinalizePbc
  public :: StepEBoundPbc
  public :: StepHBoundPbc

  ! --- Public Data

  public :: planepbc           ! export for mpicomms

  ! --- Constants

  ! --- Data

  integer :: planepbc(6)       ! 1 if pbc applies 0 if not

  M4_FTYPE :: phasex1 = 1.0, phasex2 = 1.0
  M4_FTYPE :: phasey1 = 1.0, phasey2 = 1.0
  M4_FTYPE :: phasez1 = 1.0, phasez2 = 1.0


contains

!----------------------------------------------------------------------

 subroutine ReadConfigPbc(funit,lcount,string)

    integer :: funit, lcount
    character(len=*) :: string

    character(len=LINELNG) :: line

    M4_WRITE_DBG({". enter ReadConfigPbc"})

    if ( string .ne. "(PBC" ) then
       M4_FATAL_ERROR({"BAD SECTION IDENTIFIER: ReadConfigPbc"})
    endif


    ! TODO: configure phase factors for PBCs here.

    ! call readint(funit, lcount, PMLMAX)
    ! M4_WRITE_DBG({"read PMLMAX: ", PMLMAX})


    ! TODO: add some checks on numerical values

    modconfigured = .true.

    M4_WRITE_DBG({". exit ReadConfigPbc"})

  end subroutine ReadConfigPbc



!----------------------------------------------------------------------

  subroutine InitializePbc(planebound, num)

    integer :: planebound(6), num
    integer :: i

    M4_WRITE_DBG({". enter InitializePbc"})

    planepbc = 0
    do i = 1, 6
       if ( planebound(i) .eq. num ) planepbc(i) = 1
    end do

    M4_WRITE_DBG({"planepml : ",(planepbc(i), i = 1,6)})


    M4_WRITE_DBG({". exit InitializePbc"})

  end subroutine InitializePbc

!----------------------------------------------------------------------

  subroutine FinalizePbc

    M4_WRITE_DBG({". enter FinalizePbc"})

    M4_WRITE_DBG({". exit FinalizePbc"})

  end subroutine FinalizePbc


!----------------------------------------------------------------------

  subroutine StepEBoundPbc(i)
    
    integer :: i
    
    if ( i .gt. 2 .and. M4_IS1D ) return
    if ( i .gt. 4 .and. M4_IS2D ) return

!    if ( numproc .gt. 1 .and. i .eq. M4_SDIM*2-1 .or. i .eq. M4_SDIM*2 ) then

       ! border is about to be exchanged by mpicomms and then calculated by StepHBoundPbc
!       return

!    endif

! TODO: MPI stuff interferes!

    select case ( i )
    case ( 1 ) 
!       Ez(M4_COORD(IMIN,JBEG:JEND,KBEG:KEND))= phasex1 * Ez(M4_COORD(IEND,JBEG:JEND,KBEG:KEND))
!       Ey(M4_COORD(IMIN,JBEG:JEND,KBEG:KEND))= phasex1 * Ey(M4_COORD(IEND,JBEG:JEND,KBEG:KEND))
    case ( 2 ) 
       Ez(M4_COORD(IMAX,JBEG:JEND,KBEG:KEND))= phasex2 * Ez(M4_COORD(IBEG,JBEG:JEND,KBEG:KEND))
       Ey(M4_COORD(IMAX,JBEG:JEND,KBEG:KEND))= phasex2 * Ey(M4_COORD(IBEG,JBEG:JEND,KBEG:KEND))
    case ( 3 ) 
!       Ez(M4_COORD(IBEG:IEND,JMIN,KBEG:KEND))= phasey1 * Ez(M4_COORD(IBEG:IEND,JEND,KBEG:KEND))
!       Ex(M4_COORD(IBEG:IEND,JMIN,KBEG:KEND))= phasey1 * Ex(M4_COORD(IBEG:IEND,JEND,KBEG:KEND))
    case ( 4 )
       Ez(M4_COORD(IBEG:IEND,JMAX,KBEG:KEND))= phasey2 * Ez(M4_COORD(IBEG:IEND,JBEG,KBEG:KEND))
       Ex(M4_COORD(IBEG:IEND,JMAX,KBEG:KEND))= phasey2 * Ex(M4_COORD(IBEG:IEND,JBEG,KBEG:KEND))
    case ( 5 ) 
!       Ey(M4_COORD(IBEG:IEND,JBEG:JEND,KMIN))= phasez1 * Ey(M4_COORD(IBEG:IEND,JBEG:JEND,KEND))
!       Ex(M4_COORD(IBEG:IEND,JBEG:JEND,KMIN))= phasez1 * Ex(M4_COORD(IBEG:IEND,JBEG:JEND,KEND))
    case ( 6 ) 
       Ey(M4_COORD(IBEG:IEND,JBEG:JEND,KMAX))= phasez2 * Ey(M4_COORD(IBEG:IEND,JBEG:JEND,KBEG))
       Ex(M4_COORD(IBEG:IEND,JBEG:JEND,KMAX))= phasez2 * Ex(M4_COORD(IBEG:IEND,JBEG:JEND,KBEG))
    end select

  end subroutine StepEBoundPbc

!----------------------------------------------------------------------

  subroutine StepHBoundPbc(i)
    
    integer :: i
    
    if ( i .gt. 2 .and. M4_IS1D ) return
    if ( i .gt. 4 .and. M4_IS2D ) return

!    if ( numproc .gt. 1 .and. i .eq. M4_SDIM*2-1 .or. i .eq. M4_SDIM*2 ) then

       ! border has been exchanged by mpicomms, now phase multiply

! TODO: MPI stuff interferes!

       select case ( i )
       case ( 1 ) 
          Hz(M4_COORD(IMIN,JBEG:JEND,KBEG:KEND))= phasex1 * Hz(M4_COORD(IEND,JBEG:JEND,KBEG:KEND))
          Hy(M4_COORD(IMIN,JBEG:JEND,KBEG:KEND))= phasex1 * Hy(M4_COORD(IEND,JBEG:JEND,KBEG:KEND))
       case ( 2 ) 
!          Hz(M4_COORD(IMAX,JBEG:JEND,KBEG:KEND))= phasex2 * Hz(M4_COORD(IMAX,JBEG:JEND,KBEG:KEND))
!          Hy(M4_COORD(IMAX,JBEG:JEND,KBEG:KEND))= phasex2 * Hy(M4_COORD(IMAX,JBEG:JEND,KBEG:KEND))
       case ( 3 ) 
          Hz(M4_COORD(IBEG:IEND,JMIN,KBEG:KEND))= phasey1 * Hz(M4_COORD(IBEG:IEND,JEND,KBEG:KEND))
          Hx(M4_COORD(IBEG:IEND,JMIN,KBEG:KEND))= phasey1 * Hx(M4_COORD(IBEG:IEND,JEND,KBEG:KEND))
       case ( 4 )
!          Hz(M4_COORD(IBEG:IEND,JMAX,KBEG:KEND))= phasey2 * Hz(M4_COORD(IBEG:IEND,JMAX,KBEG:KEND))
!          Hx(M4_COORD(IBEG:IEND,JMAX,KBEG:KEND))= phasey2 * Hx(M4_COORD(IBEG:IEND,JMAX,KBEG:KEND))
       case ( 5 ) 
          Hy(M4_COORD(IBEG:IEND,JBEG:JEND,KMIN))= phasez1 * Hy(M4_COORD(IBEG:IEND,JBEG:JEND,KEND))
          Hx(M4_COORD(IBEG:IEND,JBEG:JEND,KMIN))= phasez1 * Hx(M4_COORD(IBEG:IEND,JBEG:JEND,KEND))
       case ( 6 ) 
!          Hy(M4_COORD(IBEG:IEND,JBEG:JEND,KMAX))= phasez2 * Hy(M4_COORD(IBEG:IEND,JBEG:JEND,KMAX))
!          Hx(M4_COORD(IBEG:IEND,JBEG:JEND,KMAX))= phasez2 * Hx(M4_COORD(IBEG:IEND,JBEG:JEND,KMAX))
       end select

!    endif

  end subroutine StepHBoundPbc
  
!----------------------------------------------------------------------

end module pbc

!
! Authors:  J.Hamm
! Modified: 4/1/2008
!
!======================================================================
