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
!  OpenOut
!  CloseOut
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
  use mpiworld
  use grid
  use regobj
  use outobj
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

    numobjout = 0

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

  subroutine OpenOut

    integer :: n

    do n=1, numobjout
       
       select case ( objout(n)%format ) 
! ** call output open methods
! 1.
       case ( "GPL" ) 
          call OpenObjOutgpl(objout(n))
! 2.
! **
       case default
          write(STDERR,*) "!ERROR UNDEFINED OUTPUT FORMAT: OpenOut/out"
          stop
       end select

    end do

  end subroutine OpenOut

!----------------------------------------------------------------------

  subroutine CloseOut

    integer :: n

   do n=1, numobjout
       
       select case ( objout(n)%format ) 
! ** call output close methods
! 1.
       case ( "GPL" ) 
          call CloseObjOutgpl(objout(n))
! 2.
! **
       case default
          write(STDERR,*) "!ERROR UNDEFINED OUTPUT FORMAT: CloseOut/out"
          stop
       end select

    end do

  end subroutine CloseOut

!----------------------------------------------------------------------

  subroutine WriteHeaderOut
    
    integer :: n

    do n=1, numobjout
       
  
       select case ( objout(n)%format ) 
! ** call output write header methods
! 1.
       case ( "GPL" ) 
          call WriteHeaderObjOutgpl(objout(n))
! 2.
! **
       case default
          write(STDERR,*) "!ERROR UNDEFINED OUTPUT FORMAT: CloseOut/out"
          stop
       end select

    enddo
    
  end subroutine WriteHeaderOut
  
!----------------------------------------------------------------------

  subroutine WriteDataOut(ncyc)
    
    integer :: n, ncyc

    do n=1, numobjout
      
       select case ( objout(n)%format ) 
! ** call output write data methods
! 1.
       case ( "GPL" ) 
          call WriteDataObjOutgpl(objout(n), ncyc)
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
