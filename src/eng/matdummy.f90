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
  use reglist
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
  public :: ReadMatDummyObj

  ! --- Public Data

  public :: T_MATDUMMY
  public :: matdummyobj
  public :: nummatdummyobj

  ! --- Constants

  integer, parameter :: MAXMATDUMMYOBJ = 10

  ! --- Types

  type T_MATDUMMY

     integer :: regidx             ! regobj index


  end type T_MATDUMMY

  ! --- Fields

  type(T_MATDUMMY) :: matdummyobj(MAXMATDUMMYOBJ) 
  integer :: nummatdummyobj

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
    type(T_REG) :: reg
    type(T_MATDUMMY) :: mat

    nummatdummyobj = nummatdummyobj + 1
    mat = matdummyobj(nummatdummyobj)

! read mat parameters ...

!    read(funit,*) mat%parameter

! read regobj

    read(funit,*) string
    if ( string .eq. "(REGION" ) then
       call ReadRegObj(reg, funit)
       mat%regidx = reg%idx
    else
       write(6,*) "!ERROR NO REGION DEFINED: ReadMatdummyobj/matdummy"
       stop
    end if

    read(funit,*) string
    if ( string .ne. ")" ) then
       write(6,*) "!ERROR BAD TERMINATOR: ReadMatdummyobj/matdummy"
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

    do n = 1, nummatdummyobj

       mat = matdummyobj(n)
       reg = regobj(mat%regidx)
     
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,
       
       ! do some stuff here

       )
              
    end do
    
  end subroutine StepHMatDummy

!----------------------------------------------------------------------


  subroutine StepEMatDummy(ncyc)

    integer :: ncyc
    integer :: n
    type(T_MATDUMMY) :: mat

    M4_REGLOOP_DECL(reg,p,i,j,k,w)

    do n = 1, nummatdummyobj

       mat = matdummyobj(n)
       reg = regobj(mat%regidx)
     
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,
       
       ! do some stuff here

       )
              
    end do
    
  end subroutine StepEMatDummy

  
  
!----------------------------------------------------------------------

end module matdummy

! =====================================================================


