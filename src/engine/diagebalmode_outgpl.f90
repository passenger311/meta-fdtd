!-*- F90 -*------------------------------------------------------------
!
!  module: diagebalmode_outgpl / meta
!
!  this module handles GPL output of data related to the diagebalmode
!  module.
!
!----------------------------------------------------------------------


module diagebalmode_outgpl

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use fdtd_calc
  use diagebalmode
  use out_calc

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'DIAGEBALMODE_OUTGPL'

 ! --- Public Methods

  public :: InitializeDiagebalModeOutgplObj
  public :: FinalizeDiagebalModeOutgplObj
  public :: WriteDataDiagebalModeOutgplObj

contains

!----------------------------------------------------------------------

  subroutine InitializeDiagebalModeOutgplObj(out)

    type (T_OUT) :: out
    
    out%numnodes = 1

  end subroutine InitializeDiagebalModeOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeDiagebalModeOutgplObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeDiagebalModeOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataDiagebalModeOutgplObj(out, mode)

    type (T_OUT) :: out
    logical :: mode
    type (T_DIAGEBALMODE) :: diag

    if ( .not. mode ) return

    select case (out%fn)
    case('S') ! sum
       diag = diagebalmodeobj(out%objidx)
       call WriteEnergyDensityDiagebalMode(out, diag, .false.)
    case default
       write(out%funit,*) "OUTPUT FUNCTION", out%fn, "NOT IMPLEMENTED" 
    end select
    
  end subroutine WriteDataDiagebalModeOutgplObj

  subroutine WriteEnergyDensityDiagebalMode(out, diag, mode)

    type (T_OUT) :: out
    type (T_DIAGEBALMODE) :: diag
    logical :: mode
    
    write(out%funit,*) diag%energy_density

  end subroutine WriteEnergyDensityDiagebalMode
        
end module diagebalmode_outgpl

!
! Authors:  J.Hamm
! Modified: 1/02/2008
!
!======================================================================
