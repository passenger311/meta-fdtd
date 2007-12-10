!----------------------------------------------------------------------
!
!  module: outgpl
!
!  ascii (gnuplot) output module
!
!  subs:
!
!  InitializeOutgplObj
!  FinalizeOutgplObj
!  OpenOutgplObj
!  CloseOutgplObj
!  WriteHeaderOutgplObj
!  WriteDataOutgplObj
!  
!----------------------------------------------------------------------


!======================================================================
!
!

module outgpl

  use constant
  use strings
  use reglist
  use outlist

! outgpl modules 
  use fdtd_outgpl
M4_FOREACH_OUTGPL({use }, {
  })

  implicit none
  private
  save
  
  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'OUTGPL'

  ! --- Public Methods

  public :: InitializeOutgplObj
  public :: FinalizeOutgplObj
  public :: OpenOutgplObj
  public :: CloseOutgplObj
  public :: WriteHeaderOutgplObj
  public :: WriteDataOutgplObj

  ! --- Public Data

  ! --- Constants

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializeOutgplObj(out)

    type(T_OUT) :: out

    


! call various finalization methods
    select case ( out%modl ) 
    case ("FDTD")
       call  InitializeFdtdOutgplObj(out)
       M4_FOREACH_OUTGPL2({case ("},{")
             call Initialize},{OutgplObj(UNITTMP)
             })
    end select

  end subroutine InitializeOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeOutgplObj(out)

    type(T_OUT) :: out

    select case ( out%modl ) 

! call various finalization methods
    case ("FDTD")
       call FinalizeFdtdOutgplObj(out)
       M4_FOREACH_OUTGPL2({case ("},{")
             call Finalize},{OutgplObj
             })
    end select

  end subroutine FinalizeOutgplObj


!----------------------------------------------------------------------

  subroutine OpenOutgplObj(out)

    type(T_OUT) :: out

    open(UNITTMP,FILE=out%filename,STATUS="unknown")

  end subroutine OpenOutgplObj

!----------------------------------------------------------------------

  subroutine CloseOutgplObj(out)

    type(T_OUT) :: out

    close(UNITTMP)

  end subroutine CloseOutgplObj

!----------------------------------------------------------------------

  subroutine WriteHeaderOutgplObj(out)

    type(T_OUT) :: out
    type(T_REG) :: reg
    reg = regobj(out%regidx)
    
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

  subroutine WriteDataOutgplObj(out, ncyc, mode)

    type(T_OUT) :: out
    integer :: ncyc
    logical :: mode
    
    if ( ncyc .lt. out%ns .or. ncyc .gt. out%ne .or. &
       mod(ncyc - out%ns, out%dn) .ne. 0 ) then
       return
    end if
    
    if ( mode ) call OpenOutgplObj(out)
    
    select case ( out%modl ) 

! call various output methods
    case ("FDTD")
       call WriteDataFdtdOutgplObj(out,ncyc,mode)
       M4_FOREACH_OUTGPL2({case ("},{")
             call WriteData},{OutgplObj(out,ncyc,mode)
             })
    end select

     if ( mode ) call CloseOutgplObj(out)


  end subroutine WriteDataOutgplObj
  
!----------------------------------------------------------------------

end module outgpl

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
