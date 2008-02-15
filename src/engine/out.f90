!-*- F90 -*------------------------------------------------------------
! 
!  module: out / meta3
!
!  this module manages general output functionality independent of the
!  specific output format.
!
!----------------------------------------------------------------------


!======================================================================
!
! Note: currently supports GPL, VTK and HD5 formats.
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
  use outset
  use outgpl
  use outvtk
M4_IFELSE_HD5({
  use outhd5
})
! ! 2.
! **

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'OUT'

  ! --- Public Methods

  public :: InitializeOut
  public :: FinalizeOut
  public :: WriteDataOut

  ! --- Public Data


contains

!----------------------------------------------------------------------

  subroutine InitializeOut

    integer :: n,m 


    do n=1, numoutobj
       
       do m = n+1, numoutobj
          if ( outobj(n)%filename .eq. outobj(m)%filename .and.&
             outobj(n)%fmt .eq. outobj(m)%fmt ) then
             M4_FATAL_ERROR({"DUPLICATE FILENAME : ",TRIM(outobj(n)%filename)})
          endif
       end do

       M4_IFELSE_DBG({call EchoOutObj(outobj(n))})

       select case ( outobj(n)%fmt ) 
! ** call output buffer preparation
! 1.
       case ( "SET" ) 
          call InitializeOutsetObj(outobj(n))
       case ( "GPL" ) 
          call InitializeOutgplObj(outobj(n))
       case ( "VTK" ) 
          call InitializeOutvtkObj(outobj(n))
M4_IFELSE_HD5({
       case ( "HD5" ) 
          call InitializeOuthd5Obj(outobj(n))
})
! 2.
! **
       case default
          M4_FATAL_ERROR({"UNDEFINED OUTPUT FORMAT ", TRIM(outobj(n)%fmt)})
       end select

    enddo

    M4_WRITE_DBG({". exit InitializeOut/out"})

  end subroutine InitializeOut


!----------------------------------------------------------------------

  subroutine FinalizeOut

    integer :: n

    M4_WRITE_DBG({". enter FinalizeOut/out"})

    do n=1, numoutobj
       
       select case ( outobj(n)%fmt ) 
! ** call output buffer preparation
! 1.
       case ( "SET" ) 
          call FinalizeOutsetObj(outobj(n))
       case ( "GPL" ) 
          call FinalizeOutgplObj(outobj(n))
       case ( "VTK" ) 
          call FinalizeOutvtkObj(outobj(n))
M4_IFELSE_HD5({
       case ( "HD5" ) 
          call FinalizeOuthd5Obj(outobj(n))
})
! 2.
! **
       case default
          M4_FATAL_ERROR({"UNDEFINED OUTPUT FORMAT ", TRIM(outobj(n)%fmt)})
       end select

    enddo

    M4_WRITE_DBG({". exit FinalizeOut/out"})

  end subroutine FinalizeOut

!----------------------------------------------------------------------

  subroutine WriteDataOut(ncyc, mode)
    
    integer :: n, ncyc
    logical :: mode

    do n=1, numoutobj

     
       if ( ncyc .lt. outobj(n)%ns .or. ncyc .gt. outobj(n)%ne .or. &
            mod(ncyc - outobj(n)%ns, outobj(n)%dn) .ne. 0 ) then
          cycle
       end if

       M4_WRITE_DBG({"ncyc=",TRIM(i2str(ncyc)),", write out # ", TRIM(i2str(n)), " : ", &
            TRIM(outobj(n)%fmt), " ", TRIM(outobj(n)%filename)})
   
       select case ( outobj(n)%fmt ) 
! ** call output write data methods
! 1.
       case ( "SET" ) 
          call WriteDataOutsetObj(outobj(n), ncyc, mode)
       case ( "GPL" ) 
          call WriteDataOutgplObj(outobj(n), ncyc, mode)
       case ( "VTK" ) 
          call WriteDataOutvtkObj(outobj(n), ncyc, mode)
M4_IFELSE_HD5({
       case ( "HD5" ) 
          call WriteDataOuthd5Obj(outobj(n), ncyc, mode)
})
! 2.
! **
       case default
          M4_FATAL_ERROR({"UNDEFINED OUTPUT FORMAT ", TRIM(outobj(n)%fmt)})
       end select

    enddo
    
  end subroutine WriteDataOut
  
!----------------------------------------------------------------------

end module out

!
! Authors:  J.Hamm 
! Modified: 6/1/2008
!
!======================================================================
