!-*- F90 -*------------------------------------------------------------
! 
!  module: out / meta3
!
!  this module manages general output functionality independent of the
!  specific output format.
!
!  subs:
!
!  InitializeOut
!  FinalizeOut
!  WriteHeaderOut
!  WriteDataOut
!
!----------------------------------------------------------------------


!======================================================================
!
!

module out

  use constant
  use strings
  use reglist
  use outlist
  use mpiworld
  use grid

! ** add output modules
! 1.
  use outgpl
! 2.
! **

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'OUT'

  ! --- Public Methods

  public :: InitializeOut
  public :: FinalizeOut
  public :: WriteHeaderOut
  public :: WriteDataOut

  ! --- Public Data


contains

!----------------------------------------------------------------------

  subroutine InitializeOut

    integer :: n

    do n=1, numoutobj
       
       select case ( outobj(n)%fmt ) 
! ** call output buffer preparation
! 1.
       case ( "GPL" ) 
          call InitializeOutgplObj(outobj(n))
! 2.
! **
       case default
          write(STDERR,*) "!ERROR UNDEFINED OUTPUT FORMAT: CloseOut/out"
          stop
       end select

    enddo

  end subroutine InitializeOut


!----------------------------------------------------------------------

  subroutine FinalizeOut

    integer :: n

    do n=1, numoutobj
       
       select case ( outobj(n)%fmt ) 
! ** call output buffer preparation
! 1.
       case ( "GPL" ) 
          call FinalizeOutgplObj(outobj(n))
! 2.
! **
       case default
          write(STDERR,*) "!ERROR UNDEFINED OUTPUT FORMAT: CloseOut/out"
          stop
       end select

    enddo

  end subroutine FinalizeOut


!----------------------------------------------------------------------

  subroutine WriteHeaderOut
    
    integer :: n

    do n=1, numoutobj
       
       select case ( outobj(n)%fmt ) 
! ** call output write header methods
! 1.
       case ( "GPL" ) 
          call WriteHeaderOutgplObj(outobj(n))
! 2.
! **
       case default
          write(STDERR,*) "!ERROR UNDEFINED OUTPUT FORMAT: CloseOut/out"
          stop
       end select

    enddo
    
  end subroutine WriteHeaderOut
  

!----------------------------------------------------------------------

  subroutine WriteDataOut(ncyc, mode)
    
    integer :: n, ncyc
    logical :: mode

    do n=1, numoutobj
      
       select case ( outobj(n)%fmt ) 
! ** call output write data methods
! 1.
       case ( "GPL" ) 
          call WriteDataOutgplObj(outobj(n), ncyc, mode)
! 2.
! **
       case default
          write(STDERR,*) "!ERROR UNDEFINED OUTPUT FORMAT: CloseOut/out"
          stop
       end select

    enddo
    
  end subroutine WriteDataOut
  
!----------------------------------------------------------------------

end module out

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
