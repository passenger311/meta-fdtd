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

!
!
! Material manages the hooks to include various material equations
! into the fdtd algorithm. These material modules must register with 
! Material
!

module mathook

  use region

! 1.
  use mat_source
! 2.

  implicit none
  save

contains

  subroutine InitializeMat

    integer :: err

! 1. 
    call InitializeMatSource
! 2.

  end subroutine InitializeMat
  

  subroutine FinalizeMat

! 1.
    call FinalizeMatSource
! 2.

  end subroutine FinalizeMat


  subroutine ReadMat

    select case ( string )
! 1.
       case("(SOURCE"):
          call ReadMatSource
! 2.
    end select
    call ReadRegion
    

  end subroutine ReadMat
    
  subroutine StepEMat

! 1.
    call StepEMatSource
! 2.

  end subroutine StepEMat


  subroutine StepHMat

! 1.
    call StepHMatSource
! 2.

  end subroutine StepHMat



end module mat






