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

  subroutine ReadConfigDiag(funit,string)

    integer :: funit
    character(len=*) :: string

    select case ( string )
       
       M4_FOREACH_DIAG2({case ("(},{")
       call Read},{Obj(UNITTMP)
       })
       
    case default
       if ( string(1:1) .ne. "!" ) then
          M4_FATAL_ERROR({"RECEIVED BAD TOKEN: ReadConfig/diag"})
       endif
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






