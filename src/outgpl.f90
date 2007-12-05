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
!  OpenOutgplObj
!  CloseOutgplObj
!  WriteHeaderOutgplObj
!  PrepDataOutgplObj
!  WriteDataOutgplObj
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

  subroutine OpenOutgplObj(out)

    type(T_OUT) :: out

    open(UNITTMP,FILE=out%filename,STATUS='unknown')

  end subroutine OpenOutgplObj

!----------------------------------------------------------------------

  subroutine CloseOutgplObj(out)

    type(T_OUT) :: out

    close(UNITTMP)

  end subroutine CloseOutgplObj

!----------------------------------------------------------------------

  subroutine WriteHeaderOutgplObj(out)

    type(T_OUT) :: out

    type(T_REGION) :: reg
    reg = reglistobj(out%regidx)
    
    call OpenOutgplObj(out)
    
    write(UNITTMP,*) '# ',out%fmt               ! format
    write(UNITTMP,*) '# ',out%modl              ! module
    write(UNITTMP,*) '# ',out%fn                ! function
    write(UNITTMP,*) '# ',out%mode              ! mode
    write(UNITTMP,*) '# ',out%ns,out%ne,out%dn  ! time frame
    write(UNITTMP,*) '# ',reg%is,reg%ie,reg%di  ! space box
    write(UNITTMP,*) '# ',reg%js,reg%je,reg%dj
    write(UNITTMP,*) '# ',reg%ks,reg%ke,reg%dk

    call CloseOutgplObj(out)

  end subroutine WriteHeaderOutgplObj

!----------------------------------------------------------------------

  subroutine PrepDataOutgplObj(out, ncyc)

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
       call PrepDataOutgplFdtdObj(out,ncyc)
! 2.
! **
    end select

  end subroutine PrepDataOutgplObj

!----------------------------------------------------------------------

  subroutine WriteDataOutgplObj(out, ncyc)

    type(T_OUT) :: out
    integer :: ncyc
    
    if ( ncyc .lt. out%ns .or. ncyc .gt. out%ne .or. 
       mod(ncyc - out%ns, out%dn) .ne. 0 ) then
       return
    end if
    
    call OpenOutgplObj(out)
    
    select case ( out%modl ) 
! ** call output methods
! 1.
    case ("fdtd")
       call WriteDataOutgplFdtdObj(out,ncyc)
! 2.
! **
    end select

    call CloseOutgplObj(out)


  end subroutine WriteDataOutgplObj
  
!----------------------------------------------------------------------

end module outgpl

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
