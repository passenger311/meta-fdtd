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
!  CreateObjOut
!  SetObjOut
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
  use region
! ** add output modules
! 1.
  use outgpl
! 2.
! **

  implicit none
  save

  integer, parameter :: MAXOBJOUT = 500

  ! --- Types

  type T_OUT 

     ! output format, eg "GPL"
     character(len=STRLNG) :: fmt

     ! filename
     character(len=STRLNG) :: filename

     ! output module
     character(len=STRLNG) :: modl

     ! output function and mode
     character(len=STRLNG) :: fn
     character(len=10) :: mode

     ! time slices
     integer ns, ne, dn

     ! index to spatial region
     integer :: regidx
     
     ! index
     integer :: idx

  end type T_OUTBAS

  ! --- Variables

  type (T_OUT) :: objgpl(MAXOBJOUT)
  integer :: numobjout = 0


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

  subroutine ReadObjOut(out, funit)

    integer :: funit
    type(T_OUT) :: out
    type(T_OUT), external :: CreateObjOut
    
    character (len=STRLNG) :: fmt, modl, fn, mode
    integer :: ns, ne, dn
    type(T_REGION) :: reg
    

    out = CreateObjOut()
    ! read output information
    read(funit,*) fmt
    read(funit,*) modl
    read(funit,*) fn, mode
    read(funit,*) ns, ne, dn
    read(funit,*) string
    ! consume region start string
    if ( string .ne. "(region" ) then
       write(STDERR,*) "!ERROR NO REGION DEFINED: ReadObjOut/out"
       stop
    end if
    ! read the region information
    call ReadObjRegion(reg, funit)
    ! write output parameters into out object
    call SetObjOut(out, fn, mode, reg, ns, ne, dn)
    ! consume the closing ")"
    read(funit,*) string
    if ( string .ne. "(region" ) then
       write(STDERR,*) "!ERROR (OUT LACKS ) TERMINATOR: ReadObjOut/out"
       stop
    end if

  end subroutine ReadObjOut

!----------------------------------------------------------------------

  type(T_OUT) function CreateObjOut

    numobjout = numobjout + 1
    objout(numobjout)%idx = numobjout
    
    objout(numobjout)%fn = "none"
    objout(numobjout)%ns = 0
    objout(numobjout)%ne = 0
    objout(numobjout)%dn = 0

    CreateObjOut = objout(numobjout)
   
  end function CreateObjOut

!----------------------------------------------------------------------

  subroutine SetObjOut(out, fn, mode, reg, ns, ne, dn) 

    type (T_OUT) :: out

    out%fn = fn
    out%mode = mode
    out%regidx = reg%idx
    out%ns = ns
    out%ne = ne
    out%dn = dn

  end subroutine SetObjOut

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
