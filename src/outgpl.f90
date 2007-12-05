!----------------------------------------------------------------------
!
!  module: outgpl
!
!  ascii (gnuplot) output module
!
!  subs:
!
!  InitializeOutgpl
!  FinalizeOutgpl
!  OpenOutObjgpl
!  CloseOutObjgpl
!  WriteHeaderOutObjgpl
!  PrepDataOutObjgpl
!  WriteDataOutObjgpl
!  
!----------------------------------------------------------------------


!======================================================================
!
!

module outgpl

  use constant
  use strings
  use mpiworld
  use grid  
  use reglist
  use outlist
! ** add output modules
! 1.
  use outgpl_fdtd
! 2.
! **

  implicit none
  save

contains

!----------------------------------------------------------------------

  subroutine InitializeOutgpl

! ** call output initialize methods
! 1.
    call InitializeOutgplFdtd
! 2.
! **

  end subroutine InitializeOutgpl

!----------------------------------------------------------------------

  subroutine FinalizeOutgpl

! ** call output finalize methods
! 1.
    call InitializeOutgplFdtd
! 2.
! **

  end subroutine FinalizeOutgpl

!----------------------------------------------------------------------

  subroutine OpenOutObjgpl(out)

    type(T_OUT) :: out

    open(UNITTMP,FILE=out%filename,STATUS='unknown')

  end subroutine OpenOutObjgpl

!----------------------------------------------------------------------

  subroutine CloseOutObjgpl(out)

    type(T_OUT) :: out

    close(UNITTMP)

  end subroutine CloseOutObjgpl

!----------------------------------------------------------------------

  subroutine WriteHeaderOutObjgpl(out)

    type(T_OUT) :: out

    type(T_REGION) :: reg
    reg = reglistobj(out%regidx)
    
    call OpenOutObjgpl(out)
    
    write(UNITTMP,*) '# ',out%fmt               ! format
    write(UNITTMP,*) '# ',out%modl              ! module
    write(UNITTMP,*) '# ',out%fn                ! function
    write(UNITTMP,*) '# ',out%mode              ! mode
    write(UNITTMP,*) '# ',out%ns,out%ne,out%dn  ! time frame
    write(UNITTMP,*) '# ',reg%is,reg%ie,reg%di  ! space box
    write(UNITTMP,*) '# ',reg%js,reg%je,reg%dj
    write(UNITTMP,*) '# ',reg%ks,reg%ke,reg%dk

    call CloseOutObjgpl(out)

  end subroutine WriteHeaderOutObjgpl

!----------------------------------------------------------------------

  subroutine PrepDataOutObjgpl(out, ncyc)

    type(T_OUT) :: out
    integer :: ncyc
    
    if ( ncyc .lt. out%ns .or. ncyc .gt. out%ne .or. 
       mod(ncyc - out%ns, out%dn) .ne. 0 ) then
       return
    end if
    
    select case ( out%modl ) 
! ** call output methods
! 1.
    case ("fdtd")
       call PrepDataOutObjgplFdtd(out,ncyc)
! 2.
! **
    end select

  end subroutine PrepDataOutObjgpl

!----------------------------------------------------------------------

  subroutine WriteDataOutObjgpl(out, ncyc)

    type(T_OUT) :: out
    integer :: ncyc
    
    if ( ncyc .lt. out%ns .or. ncyc .gt. out%ne .or. 
       mod(ncyc - out%ns, out%dn) .ne. 0 ) then
       return
    end if
    
    call OpenOutObjgpl(out)
    
    select case ( out%modl ) 
! ** call output methods
! 1.
    case ("fdtd")
       call WriteDataOutObjgplFdtd(out,ncyc)
! 2.
! **
    end select

    call CloseOutObjgpl(out)


  end subroutine WriteDataOutObjgpl
  
!----------------------------------------------------------------------

end module outgpl

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
