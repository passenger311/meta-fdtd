!-*- F90 -*------------------------------------------------------------
!
!  module: mat / max3d
!
!  material equations.
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
! The Mat module provides the hooks to include various material 
! equations into the fdtd algorithm.


module mat

! ** add material modules
! 1.
  use matsource
! 2.
! **

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'mat'

  ! --- Public Methods

  public :: InitializeMat
  public :: FinalizeMat
  public :: StepMatE
  public :: StepMatH

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

! ** call material initialize methods
! 1. 
    call InitializeMatSource
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
            call ReadMatSourceObj(UNITTMP)
! 2.
! **
         end select
         
      enddo
      close(UNITTMP)

    end subroutine ReadConfig
        
  end subroutine InitializeMat
    
!----------------------------------------------------------------------

  subroutine FinalizeMat

! ** call material finalize methods
! 1.
    call FinalizeMatSource
! 2.
! **

  end subroutine FinalizeMat

!----------------------------------------------------------------------
 
  subroutine StepMatE(ncyc)

    integer :: ncyc

! ** call material e-field step methods
! 1.
    call StepMatSourceE(ncyc)
! 2.
! **

  end subroutine StepMatE

!----------------------------------------------------------------------

  subroutine StepMatH(ncyc)

    integer :: ncyc

! ** call material h-field step methods
! 1.
    call StepMatSourceH(ncyc)
! 2.
! **

  end subroutine StepMatH

!----------------------------------------------------------------------


end module mat

! =====================================================================






