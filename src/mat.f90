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
!    ReadMat
!    StepEMat
!    StepHMat
!
!----------------------------------------------------------------------


! =====================================================================
!
! The Mat module provides the hooks to include various material 
! equations into the fdtd algorithm.


module mat

  use region

! ** add material modules
! 1.
  use matsource
! 2.
! **

  implicit none
  save

  ! --- Constants

  character(len=255), parameter :: pfxmat = 'mat'

contains

!----------------------------------------------------------------------

  subroutine InitializeMat

    integer :: err

! ** call material initialize methods
! 1. 
    call InitializeMatSource
! 2.
! **

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

  subroutine ReadMat

    character(len=STRLNG) :: file, string
    integer :: ios
    
    file=cat2(pfxmat,sfxin)
    
    open(UNITTMP,FILE=file,STATUS='unknown')
    do
       read(UNITTMP,*, IOSTAT=ios) string
       if(ios .ne. 0) exit
       select case ( string )
! ** call material config indentifiers
! 1.
       case("(MATSOURCE")
          call ReadObjMatSource(UNITTMP)
! 2.
! **
       end select

    enddo
    close(UNITTMP)

  end subroutine ReadMat
   
!----------------------------------------------------------------------
 
  subroutine StepEMat(ncyc)

    integer :: ncyc

! ** call material e-field step methods
! 1.
    call StepEMatSource(ncyc)
! 2.
! **

  end subroutine StepEMat

!----------------------------------------------------------------------

  subroutine StepHMat(ncyc)

    integer :: ncyc

! ** call material h-field step methods
! 1.
    call StepHMatSource(ncyc)
! 2.
! **

  end subroutine StepHMat

!----------------------------------------------------------------------


end module mat

! =====================================================================






