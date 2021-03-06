!-*- F90 -*------------------------------------------------------------
!
!  module: diag / meta
!
!  diagnostics hook.
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------


! =====================================================================
!
! The Diag module hook provides the hooks to include various diagnostic
! modules into the fdtd algorithm.


module diag

  use strings
  use constant

  M4_FOREACH_DIAG({use }, {
  })

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'DIAG'

  ! --- Public Methods

  public :: ReadConfigDiag
  public :: InitializeDiag
  public :: FinalizeDiag
  public :: StepEDiag
  public :: StepHDiag

  ! --- Public Data

  ! --- Constants

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine ReadConfigDiag(funit,lcount,string)

    integer :: funit,lcount
    character(len=*) :: string

    select case ( string )
       
       M4_FOREACH_DIAG2({case ("(},{")
       call Read},{Obj(funit,lcount)
       })
       
    case default
       M4_BADTOKEN_ERROR({string(1:1) .ne. "!"},lcount,string)
    end select
    
  end subroutine ReadConfigDiag
  
!----------------------------------------------------------------------
  
  subroutine InitializeDiag

    M4_FOREACH_DIAG({call Initialize},{
    })

  end subroutine InitializeDiag
 
!----------------------------------------------------------------------

  subroutine FinalizeDiag

    M4_FOREACH_DIAG({call Finalize},{
    })

  end subroutine FinalizeDiag

!----------------------------------------------------------------------
 
  subroutine StepEDiag(ncyc)

    integer :: ncyc

    M4_FOREACH_DIAG({call StepE},{(ncyc)
    })

  end subroutine StepEDiag

!----------------------------------------------------------------------

  subroutine StepHDiag(ncyc)

    integer :: ncyc

    M4_FOREACH_DIAG({call StepH},{(ncyc)
    })

  end subroutine StepHDiag

!----------------------------------------------------------------------


end module diag

! =====================================================================






