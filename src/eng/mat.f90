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

  ! --- Public Data

  ! --- Constants

  ! --- Data

contains

!----------------------------------------------------------------------

    subroutine ReadConfigMat(funit,string)

      integer :: funit
      character(len=*) :: string
      
      integer :: ios

      do
         select case ( string )

            M4_FOREACH_MAT2({case ("(},{")
            call Read},{Obj(UNITTMP)
            })

         case default
            M4_FATAL_ERROR({"RECEIVED BAD TOKEN ", TRIM(string) ,": ReadConfig/mat"})
         end select
        
      enddo

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


end module mat

! =====================================================================






