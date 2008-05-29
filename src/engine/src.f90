!-*- F90 -*------------------------------------------------------------
!
!  module: src / meta
!
!  source hook. 
!
!----------------------------------------------------------------------


! =====================================================================
!
! The Src module hook provides the hooks to include various source 
! modules into the fdtd algorithm.


module src

  use strings
  use constant
  use reglist
  use grid

  M4_FOREACH_SRC({use },{
  })

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'SRC'

  ! --- Public Methods

  public :: ReadConfigSrc
  public :: InitializeSrc
  public :: FinalizeSrc
  public :: StepESrc
  public :: StepHSrc
  public :: SumJESrc
  public :: SumKHSrc

  ! --- Public Data

  ! --- Constants

  ! --- Data

contains

!----------------------------------------------------------------------

    subroutine ReadConfigSrc(funit,lcount,string)

      integer :: funit, lcount
      character(len=*) :: string
      
      select case ( string )
         
         M4_FOREACH_SRC2({case ("(},{")
         call Read},{Obj(funit,lcount)
         })
         
      case default
         M4_BADTOKEN_ERROR({string(1:1) .ne. "!"},lcount,string)
      end select
      
    end subroutine ReadConfigSrc

!----------------------------------------------------------------------

  subroutine InitializeSrc

    M4_FOREACH_SRC({call Initialize}, {
    })
        
  end subroutine InitializeSrc
    
!----------------------------------------------------------------------

  subroutine FinalizeSrc

    M4_FOREACH_SRC({call Finalize},{
    })

  end subroutine FinalizeSrc

!----------------------------------------------------------------------
 
  subroutine StepESrc(ncyc)

    integer :: ncyc

    M4_FOREACH_SRC({call StepE},{(ncyc)
    })

  end subroutine StepESrc

!----------------------------------------------------------------------

  subroutine StepHSrc(ncyc)

    integer :: ncyc

    M4_FOREACH_SRC({call StepH},{(ncyc)
    })

  end subroutine StepHSrc

!----------------------------------------------------------------------

  real(kind=8 ) function SumJESrc(mask,ncyc)
    
    logical, dimension(IBEG:IEND,JBEG:JEND,KBEG:KEND) :: mask
    integer :: ncyc
    real(kind=8) :: sum

    sum = 0.

    M4_FOREACH_SRC({sum = sum + SumJE},{(mask,ncyc)
    })

    SumJESrc = sum

  end function SumJESrc

!----------------------------------------------------------------------

  real(kind=8 ) function SumKHSrc(mask,ncyc)
    
    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc
    real(kind=8) :: sum

    sum = 0.

    M4_FOREACH_SRC({sum = sum + SumKH},{(mask,ncyc)
    })

    SumKHSrc = sum

  end function SumKHSrc

!----------------------------------------------------------------------


end module src

! =====================================================================






