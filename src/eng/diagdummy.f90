!-*- F90 -*------------------------------------------------------------
!
!  module: diagdummy / meta3
!
!  Dummy diagnostics module.
!
!  subs:
!
!    InitializeDiagDummy
!    FinalizeDiagDummy
!    ReadDiagDummyObj
!    StepEDiagDummy
!    StepHDiagDummy
!
!----------------------------------------------------------------------


! =====================================================================
!
! The DiagDummy module allows to add sources to the electromagnetic 
! field equations


module diagdummy

  use constant
  use mpiworld
  use reglist
  use grid
  use fdtd

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'DIAGDUMMY'

  ! --- Public Methods

  public :: InitializeDiagDummy
  public :: FinalizeDiagDummy
  public :: StepEDiagDummy
  public :: StepHDiagDummy
  public :: ReadDiagDummyObj

  ! --- Public Data

  public :: T_DIAGDUMMY
  public :: diagdummyobj
  public :: numdiagdummyobj

  ! --- Constants

  integer, parameter :: MAXDIAGDUMMYOBJ = 10

  ! --- Types

  type T_DIAGDUMMY

     integer :: regidx             ! regobj index


  end type T_DIAGDUMMY

  ! --- Fields

  type(T_DIAGDUMMY) :: diagdummyobj(MAXDIAGDUMMYOBJ) 
  integer :: numdiagdummyobj

contains

!----------------------------------------------------------------------

  subroutine InitializeDiagDummy

  end subroutine InitializeDiagDummy

!----------------------------------------------------------------------

  subroutine FinalizeDiagDummy

  end subroutine FinalizeDiagDummy

!----------------------------------------------------------------------

  subroutine ReadDiagDummyObj(funit)

    integer:: funit
    character(len=STRLNG) :: file, string
    integer :: ios
    type(T_REG) :: reg
    type(T_DIAGDUMMY) :: diag

    numdiagdummyobj = numdiagdummyobj + 1
    diag = diagdummyobj(numdiagdummyobj)

! read diag parameters ...

!    read(funit,*) diag%parameter

! read regobj

    read(funit,*) string
    if ( string .eq. "(REGION" ) then
       call ReadRegObj(reg, funit)
       diag%regidx = reg%idx
    else
       write(6,*) "!ERROR NO REGION DEFINED: ReadDiagdummyobj/diagdummy"
       stop
    end if

    read(funit,*) string
    if ( string .ne. ")" ) then
       write(6,*) "!ERROR BAD TERMINATOR: ReadDiagdummyobj/diagdummy"
       stop
    end if

! initialize object data

! intialize diag structure ...


  end subroutine ReadDiagDummyObj


!----------------------------------------------------------------------


  subroutine StepHDiagDummy(ncyc)

    integer :: ncyc
    integer :: n
    type(T_DIAGDUMMY) :: diag

    M4_REGLOOP_DECL(reg,p,i,j,k,w)

    do n = 1, numdiagdummyobj

       diag = diagdummyobj(n)
       reg = regobj(diag%regidx)
     
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,
       
       ! do some stuff here

       )
              
    end do
    
  end subroutine StepHDiagDummy


!----------------------------------------------------------------------


  subroutine StepEDiagDummy(ncyc)

    integer :: ncyc
    integer :: n
    type(T_DIAGDUMMY) :: diag

    M4_REGLOOP_DECL(reg,p,i,j,k,w)

    do n = 1, numdiagdummyobj

       diag = diagdummyobj(n)
       reg = regobj(diag%regidx)
     
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,
       
       ! do some stuff here

       )
              
    end do
    
  end subroutine StepEDiagDummy

  
  
!----------------------------------------------------------------------

end module diagdummy

! =====================================================================


