!-*- F90 -*------------------------------------------------------------
! 
!  module: outlist / meta3
!
!  subs:
!
!  InitializeOutList
!  FinalizeOutList
!  ReadOutObj
!  CreateOutObj
!  SetOutObj
!
!----------------------------------------------------------------------


!======================================================================
!
!

module outlist

  use constant
  use strings
  use reglist

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), parameter :: modname = 'OUTLIST'

  ! --- Public Methods

  public :: InitializeOutList
  public :: FinalizeOutList
  public :: CreateOutObj
  public :: DestroyOutObj

  ! --- Public Data

  public :: outobj
  public :: numoutobj
  public :: T_OUT

  ! --- Constants

  integer, parameter :: MAXOUTOBJ = 500

  ! --- Types

  type T_OUT 

     character(len=20) :: fmt          ! output format, eg "GPL"
     character(len=STRLNG) :: filename ! filename
     character(len=20) :: modl         ! output module
     character(len=20) :: fn           ! output function
     character(len=20) :: mode         ! output mode
     integer :: ns, ne, dn             ! time frame
     integer :: numnodes = 0           ! number of points / nodes
     integer :: regidx = 0             ! region index
     integer :: bufidx = 0             ! buffer index
     integer :: idx = 0                ! this objects index

  end type T_OUT

  ! --- Data

  type (T_OUT) :: outobj(MAXOUTOBJ)
  integer :: numoutobj = 0

contains

!----------------------------------------------------------------------

  subroutine InitializeOutList
  
    numoutobj = 0

  end subroutine InitializeOutList

!----------------------------------------------------------------------

  subroutine FinalizeOutList
  
    integer :: i
    
    do i = 1, numoutobj 
       call DestroyOutObj(outobj(i))
    end do

    numoutobj = 0

  end subroutine FinalizeOutList

!----------------------------------------------------------------------

  subroutine ReadOutObj(out, regdef, modl, funit)

    integer :: funit
    type(T_REG) :: regdef ! default region
    type(T_OUT) :: out
    type(T_OUT), external :: CreateOutObj
    
    character (len=STRLNG) :: fmt, modl, fn, mode, filename, string
    integer :: ns, ne, dn
    type(T_REG) :: reg
    

    out = CreateOutObj()
    ! read output information
    read(funit,*) filename         ! filename
    read(funit,*) fmt, fn, mode    ! format, function and mode
    read(funit,*) ns, ne, dn       ! time frame
    read(funit,*) string
    ! consume regobj start string
    if ( string .eq. "(REG" ) then
       call ReadRegObj(reg, funit) ! spatial regobj
       read(funit,*) string
       call SetOutObj(out, fmt, modl, filename, fn, mode, ns, ne, dn, reg)
    else
       call SetOutObj(out, fmt, modl, filename, fn, mode, ns, ne, dn, regdef)  
    end if

    if ( string .ne. ")" ) then
       write(STDERR,*) "!ERROR (OUT LACKS ) TERMINATOR: ReadOutObj/out"
       stop
    end if

  end subroutine ReadOutObj

!----------------------------------------------------------------------

  type(T_OUT) function CreateOutObj

    numoutobj = numoutobj + 1
    outobj(numoutobj)%idx = numoutobj
    
    outobj(numoutobj)%fmt  = "none"
    outobj(numoutobj)%modl = "none"
    outobj(numoutobj)%fn   = "none"
    outobj(numoutobj)%ns = 0
    outobj(numoutobj)%ne = 0
    outobj(numoutobj)%dn = 0

    CreateOutObj = outobj(numoutobj)
   
  end function CreateOutObj

!----------------------------------------------------------------------

  subroutine DestroyOutObj(out)

    type(T_OUT) :: out

    out%numnodes = 0
    out%regidx = 0  
    out%bufidx = 0
    out%idx = 0
   
  end subroutine DestroyOutObj

!----------------------------------------------------------------------

  subroutine SetOutObj(out, fmt, modl, filename, fn, mode, ns, ne, dn, reg)

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

  end subroutine SetOutObj
  
!----------------------------------------------------------------------

end module outlist

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
