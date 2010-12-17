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
  
  use constant
  use parse

  implicit none
  public
  save

  logical :: save_state = .false.
  logical :: load_state = .false.
  integer :: detail_level = 1

  character(len=20), parameter :: checkpoint_fn = "checkpoint.in"
  integer, parameter :: UNITCHK = 21

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'CHECKPOINT'
  logical, private :: modconfigured = .false.

  ! --- Public Methods

  !public :: ReadConfigFdtd
  !public :: InitializeFdtd
  !public :: FinalizeFdtd
  !public :: StepE
  !public :: StepH

contains

 subroutine ReadConfigCheckpoint(funit,lcount,string)

   integer :: funit, lcount
   character(len=*) :: string
    
   character(len=LINELNG) :: line
   integer :: v(2)

   M4_WRITE_DBG({". enter ReadConfigCheckpoint"})

   M4_WRITE_DBG({"received token: ", TRIM(string)})
   if ( string .ne. "(CHECKPOINT" ) then
      M4_FATAL_ERROR({"BAD SECTION IDENTIFIER: ReadConfigCheckpoint"})
   endif

    call readlogical(funit, lcount, load_state)
    M4_WRITE_DBG({"load_state: ", load_state})

    call readlogical(funit, lcount, save_state)
    M4_WRITE_DBG({"save_state: ", save_state})

    call readint(funit, lcount, detail_level)
    M4_WRITE_DBG({"detail_level: ", detail_level})
    
    call readtoken(funit, lcount, ")CHECKPOINT")

    modconfigured = .true.

    M4_WRITE_DBG({". exit ReadConfigCheckpoint"})


 end subroutine ReadConfigCheckpoint
 

 !----------------------------------------------------------------------

 subroutine InitializeCheckpoint





  
 end subroutine InitializeCheckpoint




end module checkpoint
