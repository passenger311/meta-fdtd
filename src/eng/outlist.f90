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
  public :: ReadOutObj
  public :: WriteDbgOutObj

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
     integer :: funit
     logical :: snap                   ! snapshot or continuous mode
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
    logical :: snap
    integer :: ns, ne, dn
    type(T_REG) :: reg
    
    out = CreateOutObj()
    
    M4_WRITE_DBG({". enter ReadOutObj num = ",out%idx})

    ! read output information
    read(funit,*) fmt, snap, filename    ! format and filename
    M4_WRITE_DBG({"fmt snap filename: ",TRIM(fmt)," ",snap," ", TRIM(filename)})
    read(funit,*) fn, mode   ! function and mode
    M4_WRITE_DBG({"fn mode: ", TRIM(fn)," ",TRIM(mode) })
    read(funit,*) ns, ne, dn       ! time frame
    M4_WRITE_DBG({"ns ne dn: ",ns, ne, dn })
    read(funit,*) string
    M4_WRITE_DBG({"got token ",TRIM(string)})

    ! consume regobj start string
    if ( string .eq. "(REG" ) then
       M4_WRITE_DBG({"-> ReadRegObj"})
       call ReadRegObj(reg, funit) ! spatial regobj
       read(funit,*) string
       M4_WRITE_DBG({"got token ",TRIM(string)})
       call SetOutObj(out, fmt, modl, filename, fn, mode, ns, ne, dn, reg)
    else
       M4_WRITE_DBG({"using default region!"})
       call SetOutObj(out, fmt, modl, filename, fn, mode, ns, ne, dn, regdef)  
    end if

    if ( string(1:1) .ne. ")" ) then
       M4_FATAL_ERROR({"BAD TERMINATOR: ReadOutObj/out"})
    end if

    outobj(numoutobj) = out

    M4_WRITE_DBG(". exit ReadOutObj")

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

    outobj(numoutobj)%funit = UNITTMP
    outobj(numoutobj)%snap = .false.

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

  subroutine WriteDbgOutObj(out)

    type(T_OUT) :: out

    M4_WRITE_DBG({"out # ", TRIM(i2str(out%idx)) })
    M4_WRITE_DBG({"  fmt modl filename : ", TRIM(out%fmt), " ", TRIM(out%modl), " ", TRIM(out%filename) })
    M4_WRITE_DBG({"  fn mode : ", TRIM(out%fn), " ", TRIM(out%mode) })
    M4_WRITE_DBG({"  ns ne dn : ", out%ns, out%ne, out%dn })
    M4_WRITE_DBG({"defined over"})
    call WriteDbgRegObj(regobj(out%regidx))

  end subroutine WriteDbgOutObj
  
!----------------------------------------------------------------------


end module outlist

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
