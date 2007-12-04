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
!  OpenObjOutgpl
!  CloseObjOutgpl
!  WriteHeaderObjOutgpl
!  WriteDataObjOutgpl
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
  use region
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

  subroutine OpenObjOutgpl(out)

    type(T_OUT) :: out

    open(UNITTMP,FILE=out%filename,STATUS='unknown')

  end subroutine OpenObjOutgpl

!----------------------------------------------------------------------

  subroutine CloseObjOutgpl(out)

    type(T_OUT) :: out

    close(UNITTMP)

  end subroutine CloseObjOutgpl

!----------------------------------------------------------------------

  subroutine WriteHeaderObjOutgpl(out)

    type(T_OUT) :: out

    type(T_REGION) :: reg
    reg = objregion(out%regidx)
    
    call OpenObjOutgpl(out)
    
    write(UNITTMP,*) '# ',out%fmt               ! format
    write(UNITTMP,*) '# ',out%modl              ! module
    write(UNITTMP,*) '# ',out%fn                ! function
    write(UNITTMP,*) '# ',out%mode              ! mode
    write(UNITTMP,*) '# ',out%ns,out%ne,out%dn  ! time frame
    write(UNITTMP,*) '# ',reg%is,reg%ie,reg%di  ! space box
    write(UNITTMP,*) '# ',reg%js,reg%je,reg%dj
    write(UNITTMP,*) '# ',reg%ks,reg%ke,reg%dk

    call CloseObjOutgpl(out)

  end subroutine WriteHeaderObjOutgpl

!----------------------------------------------------------------------

  subroutine WriteDataObjOutgpl(out, ncyc)

    type(T_OUT) :: out
    integer :: ncyc
    
    if ( ncyc .lt. out%ns .or. ncyc .gt. out%ne .or. 
       mod(ncyc - out%ns, out%dn) .ne. 0 ) then
       return
    end if
    
    call OpenObjOutgpl(out)
    
    select case ( out%modl ) 
! ** call output methods
! 1.
    case ("fdtd")
       call WriteDataObjOutgplFdtd(out,ncyc)
! 2.
! **
    end select

    call CloseObjOutgpl(out)


  end subroutine WriteDataObjOutgpl
  
!----------------------------------------------------------------------

end module outgpl

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
