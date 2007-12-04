!-*- F90 -*------------------------------------------------------------
! 
!  module: outobj / max3d
!
!  subs:
!
!  ReadObjOut
!  CreateObjOut
!  SetObjOut
!
!----------------------------------------------------------------------


!======================================================================
!
!

module outobj

  use constant
  use strings
  use mpiworld
  use grid
  use regobj

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
     character(len=10) :: fn
     character(len=10) :: mode

     ! time slices
     integer ns, ne, dn

     ! index to spatial regobj
     integer :: regidx

     ! number of points / nodes
     integer :: numnodes
     
     ! index
     integer :: idx

  end type T_OUT

  ! --- Variables

  type (T_OUT) :: objout(MAXOBJOUT)
  integer :: numobjout = 0


contains

!----------------------------------------------------------------------

  subroutine ReadObjOut(out, funit)

    integer :: funit
    type(T_OUT) :: out
    type(T_OUT), external :: CreateObjOut
    
    character (len=STRLNG) :: fmt, modl, fn, mode, filename, string
    integer :: ns, ne, dn
    type(T_REG) :: reg
    

    out = CreateObjOut()
    ! read output information
    read(funit,*) fmt, modl        ! format and module
    read(funit,*) filename         ! filename
    read(funit,*) fn, mode         ! function and mode
    read(funit,*) ns, ne, dn       ! time frame
    read(funit,*) string
    ! consume regobj start string
    if ( string .ne. "(regobj" ) then
       write(STDERR,*) "!ERROR NO REGION DEFINED: ReadObjOut/out"
       stop
    end if
    ! read the regobj information
    call ReadObjReg(reg, funit) ! spatial regobj
    ! write output parameters into out object
    call SetObjOut(out, fmt, modl, filename, fn, mode, ns, ne, dn, reg)
    ! consume the closing ")"
    read(funit,*) string
    if ( string .ne. "(regobj" ) then
       write(STDERR,*) "!ERROR (OUT LACKS ) TERMINATOR: ReadObjOut/out"
       stop
    end if

  end subroutine ReadObjOut

!----------------------------------------------------------------------

  type(T_OUT) function CreateObjOut

    numobjout = numobjout + 1
    objout(numobjout)%idx = numobjout
    
    objout(numobjout)%fmt  = "none"
    objout(numobjout)%modl = "none"
    objout(numobjout)%fn   = "none"
    objout(numobjout)%ns = 0
    objout(numobjout)%ne = 0
    objout(numobjout)%dn = 0

    CreateObjOut = objout(numobjout)
   
  end function CreateObjOut

!----------------------------------------------------------------------

  subroutine SetObjOut(out, fmt, modl, filename, fn, mode, ns, ne, dn, reg)

    type(T_OUT) :: out
    character(len=STRLNG) :: fmt, modl, filename, fn, mode
    type(T_REG) :: reg
    integer:: ns, ne, dn

    out%fmt = fmt
    out%modl = modl
    out%filename = filename
    out%fn = fn
    out%mode = mode
    out%ns = ns
    out%ne = ne
    out%dn = dn
    out%regidx = reg%idx
    out%numnodes = reg%numnodes

  end subroutine SetObjOut
  
!----------------------------------------------------------------------

end module outobj

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
