!-*- F90 -*------------------------------------------------------------
!
!  module:list / meta3
!
!  this generic module initializes / finalizes all list modules
!
!  subs:
!
!  InitializeList
!  FinalizeList
!
!----------------------------------------------------------------------


!======================================================================
!
!

module list

  use constant
  use strings
  use reglist
  use buflist
  use outlist

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), parameter :: modname = 'LIST'

  ! --- Public Methods

  public :: InitializeList
  public :: FinalizeList

  ! --- Public Data

  ! --- Constants

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializeList

    M4_WRITE_DBG(". enter InitializeList")

    call InitializeRegList
    call InitializeBufList
    call InitializeOutList

    M4_WRITE_DBG(". exit InitializeList")

  end subroutine InitializeList

!----------------------------------------------------------------------

  subroutine FinalizeList

    M4_WRITE_DBG(". enter FinalizeList")

    call FinalizeRegList
    call FinalizeBufList
    call FinalizeOutList

    M4_WRITE_DBG(". enter FinalizeList")

  end subroutine FinalizeList

!----------------------------------------------------------------------
 

!----------------------------------------------------------------------

end module list


!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
