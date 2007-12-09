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

! read regions and terminator

   M4_GET_REG_AND_TERMINATOR(diag,"ReadDiagDummyObj/diagdummy")

  end subroutine ReadDiagDummyObj

!----------------------------------------------------------------------

  subroutine InitializeDiagDummy

    integer :: n
    type(T_DIAGDUMMY) :: diag

    do n = 1, numdiagdummyobj

       diag = diagdummyobj(n)

       ! Initialize diagdummy object here

    end do

  end subroutine InitializeDiagDummy

!----------------------------------------------------------------------

  subroutine FinalizeDiagDummy

    integer :: n
    type(T_DIAGDUMMY) :: diag
    
    do n = 1, numdiagdummyobj
       
       diag = diagdummyobj(n)

       ! Finalize diagdummy object here

    end do

  end subroutine FinalizeDiagDummy

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


