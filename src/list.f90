!-*- F90 -*------------------------------------------------------------
!
!  module:list / max3d
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

  character(len=20), parameter :: modname = 'list'

  ! --- Public Methods

  public :: InitializeList
  public :: FinalizeList

  ! --- Public Data

  ! --- Constants

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializeList

    call InitializeRegList
    call InitializeBufList
    call InitializeOutList

  end subroutine InitializeList

!----------------------------------------------------------------------

  subroutine FinalizeList

    call FinalizeRegList
    call FinalizeBufList
    call FinalizeOutList

  end subroutine FinalizeList

!----------------------------------------------------------------------
 

!----------------------------------------------------------------------

end module list


!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
