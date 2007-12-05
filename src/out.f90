!-*- F90 -*------------------------------------------------------------
! 
!  module: out / max3d
!
!  this module manages general output functionality independent of the
!  specific output format.
!
!  subs:
!
!  InitializeOut
!  FinalizeOut
!  WriteHeaderOut
!  PrepDataOut
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
  save

contains

!----------------------------------------------------------------------

  subroutine InitializeOut

    numoutlist = 0

! ** call output initialize methods
! 1.
    call InitializeOutgpl
! 2.
! **

  end subroutine InitializeOut

!----------------------------------------------------------------------

  subroutine FinalizeOut

! ** call output finalize methods
! 1.
    call FinalizeOutgpl
! 2.
! **

  end subroutine FinalizeOut


!----------------------------------------------------------------------

  subroutine WriteHeaderOut
    
    integer :: n

    do n=1, numoutlist
       
  
       select case ( outlist(n)%format ) 
! ** call output write header methods
! 1.
       case ( "GPL" ) 
          call WriteHeaderOutgplObj(outlist(n))
! 2.
! **
       case default
          write(STDERR,*) "!ERROR UNDEFINED OUTPUT FORMAT: CloseOut/out"
          stop
       end select

    enddo
    
  end subroutine WriteHeaderOut
  
!----------------------------------------------------------------------

  subroutine PrepDataOut(ncyc)
    
    integer :: n, ncyc

    do n=1, numoutlist
      
       select case ( outlist(n)%format ) 
! ** call output write data methods
! 1.
       case ( "GPL" ) 
          call PrepDataOutgplObj(outlist(n), ncyc)
! 2.
! **
       case default
          write(STDERR,*) "!ERROR UNDEFINED OUTPUT FORMAT: CloseOut/out"
          stop
       end select

    enddo
    
  end subroutine PrepDataOut

!----------------------------------------------------------------------

  subroutine WriteDataOut(ncyc)
    
    integer :: n, ncyc

    do n=1, numoutlist
      
       select case ( outlist(n)%format ) 
! ** call output write data methods
! 1.
       case ( "GPL" ) 
          call WriteDataOutgplObj(outlist(n), ncyc)
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
