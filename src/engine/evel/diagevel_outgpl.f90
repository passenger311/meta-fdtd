!-*- F90 -*------------------------------------------------------------
!
!  module: diagevel_outgpl / meta
!
!  this module handles GPL output of data related to the diagevel module.
!
!----------------------------------------------------------------------


module diagevel_outgpl

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use fdtd_calc
  use diagevel
  use out_calc

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'DIAGEVEL_OUTGPL'

 ! --- Public Methods

  public :: InitializeDiagevelOutgplObj
  public :: FinalizeDiagevelOutgplObj
  public :: WriteDataDiagevelOutgplObj

contains

!----------------------------------------------------------------------

  subroutine InitializeDiagevelOutgplObj(out)

    type (T_OUT) :: out
    
    out%numnodes = 1

  end subroutine InitializeDiagevelOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeDiagevelOutgplObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeDiagevelOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataDiagevelOutgplObj(out, mode)

    type (T_OUT) :: out
    logical :: mode

    if ( .not. mode ) return

    M4_WRITE_DBG({"write data ",TRIM(out%filename), " ",TRIM(out%fn)})
 
    call WriteValues(out, 0)
 
  contains

    ! **************************************************************** !

    subroutine WriteValues(out, mode)

      type (T_OUT) :: out
      integer :: mode
      type (T_DIAGEVEL) :: diag

      diag = diagevelobj(out%objidx)
    
      M4_WRITE_DBG({"WriteValues!"})
      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      write(out%funit,"(5E17.8E3)") diag%u, diag%ux, diag%uy,  diag%uz,  diag%ux2,  diag%uy2, diag%uz2 

    end subroutine WriteValues

  end subroutine WriteDataDiagevelOutgplObj


end module diagevel_outgpl

!
! Authors:  J.Hamm
! Modified: 24/11/2011
!
!======================================================================
