!-*- F90 -*------------------------------------------------------------
!
!  module: checkpoint / meta
!
!  checkpoint
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------


! =====================================================================
!
! The Checkpoint module


module checkpoint

  implicit none
  public
  save

  logical :: save_state = .true.
  logical :: load_state = .true.

  character(len=20), parameter :: checkpoint_fn = "checkpoint.in"
  integer, parameter :: UNITCHK = 21

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'CHECKPOINT'
  logical, private :: modconfigured

  ! --- Public Methods

  !public :: ReadConfigFdtd
  !public :: InitializeFdtd
  !public :: FinalizeFdtd
  !public :: StepE
  !public :: StepH

contains

  subroutine InitializeCheckpoint

  
  end subroutine InitializeCheckpoint


end module checkpoint
