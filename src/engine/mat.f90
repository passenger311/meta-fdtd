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

  use checkpoint
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

    subroutine ReadConfigMat(funit,lcount,string)

      integer :: funit, lcount
      character(len=*) :: string
      
      select case ( string )
         
         M4_FOREACH_MAT2({case ("(},{")
         call Read},{Obj(funit,lcount)
         })
         
      case default
         M4_BADTOKEN_ERROR({string(1:1) .ne. "!"},lcount,string)
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

  subroutine SumJEMat(mask,ncyc,sum,idx,mode)
    
    logical, dimension(IBEG:IEND,JBEG:JEND,KBEG:KEND) :: mask
    logical :: mode
    integer :: ncyc, idx
    real(kind=8) :: sum(MAXEBALCH)

    idx = 1
    M4_FOREACH_MAT({
    call SumJE},{(mask,ncyc,sum,idx,mode)
    })

  end subroutine SumJEMat

!----------------------------------------------------------------------

  subroutine SumKHMat(mask,ncyc,sum,idx,mode)
    
    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    logical :: mode
    integer :: ncyc, idx
    real(kind=8) :: sum(MAXEBALCH)

    idx = 1
    M4_FOREACH_MAT({
    call SumKH},{(mask,ncyc,sum,idx,mode)
    })

  end subroutine SumKHMat

!----------------------------------------------------------------------


end module mat

! =====================================================================






