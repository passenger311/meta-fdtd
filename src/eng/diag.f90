!-*- F90 -*------------------------------------------------------------
!
!  module: diag / meta3
!
!  diagerial hook.
!
!  subs:
!
!    InitializeDiag
!    FinalizeDiag
!    StepDiagE
!    StepDiagH
!
!----------------------------------------------------------------------


! =====================================================================
!
! The Diag module hook provides the hooks to include various diagerial 
! modules into the fdtd algorithm.


module diag

  use strings
  use constant
  M4_FOREACH_DIAG(`use ', `
  ')

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'diag'

  ! --- Public Methods

  public :: InitializeDiag
  public :: FinalizeDiag
  public :: StepDiagE
  public :: StepDiagH

  ! --- Public Data

  public :: pfxdiag

  ! --- Constants

  character(len=STRLNG), parameter :: pfxdiag = 'diag'

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializeDiag

    integer :: err

    call ReadConfig

    M4_FOREACH_DIAG(`call Initialize',`
    ')

  contains
    
    subroutine ReadConfig
      
      character(len=STRLNG) :: file, string
      integer :: ios
      
      file=cat2(pfxdiag,sfxin)
        
      open(UNITTMP,FILE=file,STATUS='unknown')
      do
         read(UNITTMP,*, IOSTAT=ios) string
         if(ios .ne. 0) exit
         select case ( string )

            M4_FOREACH_DIAG2(`case ("(',`")
             call Read',`Obj(UNITTMP)
             ')

         case default
            write(STDERR,*) "!ERROR UNDEFINED DIAG SECTION: ReadConfig/diag"
            stop
         end select
         
      enddo
      close(UNITTMP)

    end subroutine ReadConfig
        
  end subroutine InitializeDiag
    
!----------------------------------------------------------------------

  subroutine FinalizeDiag

    M4_FOREACH_DIAG(`call Finalize',`
    ')

  end subroutine FinalizeDiag

!----------------------------------------------------------------------
 
  subroutine StepDiagE(ncyc)

    integer :: ncyc

    M4_FOREACH_DIAG(`call StepE',`(ncyc)
    ')

  end subroutine StepDiagE

!----------------------------------------------------------------------

  subroutine StepDiagH(ncyc)

    integer :: ncyc

    M4_FOREACH_DIAG(`call StepH',`(ncyc)
    ')

  end subroutine StepDiagH

!----------------------------------------------------------------------


end module diag

! =====================================================================






