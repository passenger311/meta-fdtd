!-*- F90 -*------------------------------------------------------------
!
!  module: diag / meta3
!
!  diagnostics hook.
!
!  subs:
!
!    InitializeMat
!    FinalizeMat
!    StepMatE
!    StepMatH
!
!----------------------------------------------------------------------


! =====================================================================
!
! The Diag hook is structural identical to the Mat modules in providing
! hooks for specific Diag modules. In contrary to Mat modules, Diag
! modules promise not to modify the E and H field values.
!

module diag

! ** add material modules
! 1.
  use diagpdft
! 2.
! **

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'mat'

  ! --- Public Methods

  public :: InitializeDiag
  public :: FinalizeDiag
  public :: StepDiagE
  public :: StepDiagH

  ! --- Public Data

  public :: pfxmat

  ! --- Constants

  character(len=STRLNG), parameter :: pfxmat = 'mat'

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializeDiag

    integer :: err

    call ReadConfig

! ** call material initialize methods
! 1. 
    call InitializeDiagSource
! 2.
! **

  contains
    
    subroutine ReadConfig
      
      character(len=STRLNG) :: file, string
      integer :: ios
      
      file=cat2(pfxmat,sfxin)
        
      open(UNITTMP,FILE=file,STATUS='unknown')
      do
         read(UNITTMP,*, IOSTAT=ios) string
         if(ios .ne. 0) exit
         select case ( string )
! ** read material config sections
! 1.
         case("(MATSOURCE")
            call ReadDiagSourceObj(UNITTMP)
! 2.
! **
         end select
         
      enddo
      close(UNITTMP)

    end subroutine ReadConfig
        
  end subroutine InitializeDiag
    
!----------------------------------------------------------------------

  subroutine FinalizeDiag

! ** call material finalize methods
! 1.
    call FinalizeDiagSource
! 2.
! **

  end subroutine FinalizeDiag

!----------------------------------------------------------------------
 
  subroutine StepDiagE(ncyc)

    integer :: ncyc

! ** call material e-field step methods
! 1.
    call StepDiagSourceE(ncyc)
! 2.
! **

  end subroutine StepDiagE

!----------------------------------------------------------------------

  subroutine StepDiagH(ncyc)

    integer :: ncyc

! ** call material h-field step methods
! 1.
    call StepDiagSourceH(ncyc)
! 2.
! **

  end subroutine StepDiagH

!----------------------------------------------------------------------


end module diag

! =====================================================================






