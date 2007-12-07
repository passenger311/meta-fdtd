!-*- F90 -*------------------------------------------------------------
!
!  module: mat / meta3
!
!  material hook. 
!
!  subs:
!
!    InitializeMat
!    FinalizeMat
!    StepEMat
!    StepHMat
!
!----------------------------------------------------------------------


! =====================================================================
!
! The Mat module hook provides the hooks to include various material 
! modules into the fdtd algorithm.


module mat

  use strings
  use constant
  M4_FOREACH_MAT(`use ', `
  ')

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'MAT'

  ! --- Public Methods

  public :: InitializeMat
  public :: FinalizeMat
  public :: StepEMat
  public :: StepHMat

  ! --- Public Data

  public :: pfxmat

  ! --- Constants

  character(len=STRLNG), parameter :: pfxmat = 'mat'

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializeMat

    integer :: err

    call ReadConfig

    M4_FOREACH_MAT(`call Initialize',`
    ')

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

            M4_FOREACH_MAT2(`case ("(',`")
             call Read',`Obj(UNITTMP)
             ')

         case default
            write(STDERR,*) "!ERROR UNDEFINED MAT SECTION: ReadConfig/mat"
            stop
            
         end select
         
      enddo
      close(UNITTMP)

    end subroutine ReadConfig
        
  end subroutine InitializeMat
    
!----------------------------------------------------------------------

  subroutine FinalizeMat

    M4_FOREACH_MAT(`call Finalize',`
    ')

  end subroutine FinalizeMat

!----------------------------------------------------------------------
 
  subroutine StepEMat(ncyc)

    integer :: ncyc

    M4_FOREACH_MAT(`call StepE',`(ncyc)
    ')

  end subroutine StepEMat

!----------------------------------------------------------------------

  subroutine StepHMat(ncyc)

    integer :: ncyc

    M4_FOREACH_MAT(`call StepH',`(ncyc)
    ')

  end subroutine StepHMat

!----------------------------------------------------------------------


end module mat

! =====================================================================






