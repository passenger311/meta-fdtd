!-*- F90 -*------------------------------------------------------------
!
!  module: matdummy / meta3
!
!  Dummy material module.
!
!  subs:
!
!    InitializeMatDummy
!    FinalizeMatDummy
!    ReadMatDummyObj
!    StepEMatDummy
!    StepHMatDummy
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatDummy module allows to add sources to the electromagnetic 
! field equations


module matdummy

  use constant
  use mpiworld
  use regobj
  use grid
  use fdtd

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'MATDUMMY'

  ! --- Public Methods

  public :: InitializeMatDummy
  public :: FinalizeMatDummy
  public :: StepEMatDummy
  public :: StepHMatDummy

  ! --- Public Data

  ! --- Constants

  integer, parameter :: MAXOBJMATDUMMY = 10

  ! --- Types

  type T_MATDUMMY

     integer :: regidx             ! regobj index


  end type T_MATDUMMY

  ! --- Fields

  type(T_MATDUMMY) :: objmatdummy(MAXOBJMATDUMMY) 
  integer :: numobjmatdummy

contains

!----------------------------------------------------------------------

  subroutine InitializeMatDummy

  end subroutine InitializeMatDummy

!----------------------------------------------------------------------

  subroutine FinalizeMatDummy

  end subroutine FinalizeMatDummy

!----------------------------------------------------------------------

  subroutine ReadMatDummyObj(funit)

    integer:: funit
    character(len=STRLNG) :: file, string
    integer :: ios
    type(T_REGION) :: reg
    type(T_MATDUMMY) :: mat

    numobjmatdummy = numobjmatdummy + 1
    mat = objmatdummy(numobjmatdummy)

! read mat parameters ...

!    read(funit,*) mat%parameter

! read regobj

    read(funit,*) string
    if ( string .eq. "(REGION" ) then
       call ReadRegObj(reg, funit)
       mat%regidx = reg%idx
    else
       write(6,*) "!ERROR NO REGION DEFINED: ReadObjMatDummy/matdummy"
       stop
    end if

    read(funit,*) string
    if ( string .ne. ")" ) then
       write(6,*) "!ERROR BAD TERMINATOR: ReadObjMatDummy/matdummy"
       stop
    end if

! initialize object data

! intialize mat structure ...


  end subroutine ReadMatDummyObj


!----------------------------------------------------------------------


  subroutine StepHMatDummy(ncyc)

    integer :: ncyc
    integer :: n
    type(T_MATDUMMY) :: mat

    M4_REGLOOP_DECL(reg,p,i,j,k,w)

    do n = 1, numobjmatdummy

       mat = objmatdummy(n)
    
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,
       
       ! do some stuff here

       )
              
    end do
    
  end subroutine StepHMatDummy

  
  
!----------------------------------------------------------------------

end module matdummy

! =====================================================================


