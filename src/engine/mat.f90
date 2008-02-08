!-*- F90 -*------------------------------------------------------------
!
!  module: mat / meta
!
!  material hook. 
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------


! =====================================================================
!
! The Mat module hook provides the hooks to include various material 
! modules into the fdtd algorithm.


module mat

  use strings
  use constant
  use reglist
  use grid

  M4_FOREACH_MAT({use },{
  })

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'MAT'

  ! --- Public Methods

  public :: ReadConfigMat
  public :: InitializeMat
  public :: FinalizeMat
  public :: StepEMat
  public :: StepHMat
  public :: SumJEMat
  public :: SumKHMat

  ! --- Public Data

  ! --- Constants

  ! --- Data

contains

!----------------------------------------------------------------------

    subroutine ReadConfigMat(funit,string)

      integer :: funit
      character(len=*) :: string
      
      select case ( string )
         
         M4_FOREACH_MAT2({case ("(},{")
         call Read},{Obj(funit)
         })
         
      case default
         M4_FATAL_ERROR({"RECEIVED BAD TOKEN ", TRIM(string) ,": ReadConfig/mat"})
      end select
      
    end subroutine ReadConfigMat

!----------------------------------------------------------------------

  subroutine InitializeMat

    M4_FOREACH_MAT({call Initialize}, {
    })
        
  end subroutine InitializeMat
    
!----------------------------------------------------------------------

  subroutine FinalizeMat

    M4_FOREACH_MAT({call Finalize},{
    })

  end subroutine FinalizeMat

!----------------------------------------------------------------------
 
  subroutine StepEMat(ncyc)

    integer :: ncyc

    M4_FOREACH_MAT({call StepE},{(ncyc)
    })

  end subroutine StepEMat

!----------------------------------------------------------------------

  subroutine StepHMat(ncyc)

    integer :: ncyc

    M4_FOREACH_MAT({call StepH},{(ncyc)
    })

  end subroutine StepHMat

!----------------------------------------------------------------------

  real(kind=8 ) function SumJEMat(mask,ncyc)
    
    logical, dimension(IBEG:IEND,JBEG:JEND,KBEG:KEND) :: mask
    integer :: ncyc
    real(kind=8) :: sum

    sum = 0.

    M4_FOREACH_MAT({sum = sum + SumJE},{(mask,ncyc)
    })

    SumJEMat = sum

  end function SumJEMat

!----------------------------------------------------------------------

  real(kind=8 ) function SumKHMat(mask,ncyc)
    
    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc
    real(kind=8) :: sum

    sum = 0.

    M4_FOREACH_MAT({sum = sum + SumKH},{(mask,ncyc)
    })

    SumKHMat = sum

  end function SumKHMat

!----------------------------------------------------------------------


end module mat

! =====================================================================






