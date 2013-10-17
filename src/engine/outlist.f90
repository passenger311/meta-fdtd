!-*- F90 -*------------------------------------------------------------
! 
!  module: outlist / meta
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------


!======================================================================
!
!

module outlist

  use constant
  use strings
  use parse
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
  public :: EchoOutObj

  ! --- Public Data

  public :: outobj
  public :: numoutobj
  public :: T_OUT

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
     integer :: numsteps =0            ! number of time steps 
     integer :: regidx = 0             ! region index
     integer :: bufidx = 0             ! buffer index
     integer :: objidx = 0             ! associated object index
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

  subroutine ReadOutObj(out, regdef, funit, lcount, modl, isref)

    integer :: funit, lcount
    type(T_REG) :: regdef ! default region
    type(T_OUT) :: out
    character(len=*) :: modl
    logical :: err, eof, isref
    character(len=LINELNG) :: line

    character (len=STRLNG) :: fmt, fn, mode, filename, string, line2
    logical :: snap
    integer :: ns, ne, dn, val(3)
    type(T_REG) :: reg
    
    out = CreateOutObj()
    
    M4_WRITE_DBG({". enter ReadOutObj num = ",out%idx})

    err = .false.

    call readline(funit,lcount,eof,line)
    M4_EOF_ERROR(eof,lcount)
    call getstring(line,fmt,err)
    call getstring(line,filename,err)
    M4_SYNTAX_ERROR({line .ne. ""},lcount,"FMT FILENAME")
    M4_WRITE_DBG({"fmt filename: ",TRIM(fmt)," ", TRIM(filename)})

    call readline(funit,lcount,eof,line)
    M4_EOF_ERROR(eof,lcount)
    ! first try to read timestep information
    call getints(line, val, 3, err)
    ns = val(1)
    ne = val(2)
    dn = val(3)
    ! if this failed then get fn/mode info
    if ( err .or. line .ne. "" ) then
       err = .false.
       call getstring(line,fn,err)
       M4_SYNTAX_ERROR({err},lcount,"FN [MODE [SNAP]]")
       if ( line .eq. "" ) then
          mode = "N"
       else
          call getstring(line,mode,err)
       end if
       M4_SYNTAX_ERROR({err},lcount,"FN [MODE [SNAP]]")
       if ( line .eq. "" ) then
          snap = .true.
       else
          call getlogical(line,snap,err)
       end if
       M4_SYNTAX_ERROR({err .or. line .ne. ""},lcount,"FN [MODE [SNAP]]")
       call readline(funit,lcount,eof,line)
       M4_EOF_ERROR(eof,lcount)
       call getints(line, val, 3, err)
       ns = val(1)
       ne = val(2)
       dn = val(3)
       M4_SYNTAX_ERROR({err .or. line .ne. ""},lcount,"[INTEGERS]")
       M4_WRITE_DBG({"ns ne dn: ",ns, ne, dn })
    else
       M4_WRITE_DBG({"ns ne dn: ",ns, ne, dn })
       fn = "UNKNOWN"
       mode = "X"
       snap = .true.
    end if

    M4_WRITE_DBG({"fn mode snap: ", TRIM(fn)," ",TRIM(mode), snap })

    call readline(funit,lcount,eof,line)
    call getstring(line,string,err)

    if ( string .eq. "(REG" ) then
       M4_WRITE_DBG({"r-> ReadRegObj"})
       call ReadRegObj(reg, regdef, funit, lcount, 0, isref) ! spatial regobj
       call readtoken(funit, lcount, ")OUT")
       call SetOutObj(out, fmt, snap, modl, filename, fn, mode, ns, ne, dn, reg)
    else
       M4_SYNTAX_ERROR({string .ne. ")OUT"},lcount,{")OUT"})
       M4_WRITE_DBG({"using default region!"})
       call SetOutObj(out, fmt, snap, modl, filename, fn, mode, ns, ne, dn, regdef)  
    end if

    outobj(numoutobj) = out

    M4_WRITE_DBG(". exit ReadOutObj")

  end subroutine ReadOutObj

!----------------------------------------------------------------------

  type(T_OUT) function CreateOutObj()

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

  subroutine SetOutObj(out, fmt, snap, modl, filename, fn, mode, ns, ne, dn, reg)

    type(T_OUT) :: out
    character(len=*) :: fmt, modl, filename, fn, mode
    logical :: snap
    type(T_REG) :: reg
    integer:: ns, ne, dn

    out%fmt = fmt
    out%snap = snap
    out%modl = modl
    out%filename = filename
    out%fn = fn
    out%mode = mode
    out%ns = ns
    out%ne = ne
    out%dn = dn
    out%regidx = reg%idx
    out%numnodes = reg%numnodes
    out%numsteps = (ne - ns) / dn + 1

  end subroutine SetOutObj

!----------------------------------------------------------------------

  subroutine EchoOutObj(out)

    type(T_OUT) :: out

    M4_WRITE_INFO({"--- out # ", TRIM(i2str(out%idx)) })
    M4_WRITE_INFO({"fmt modl filename = ", TRIM(out%fmt),&
         " ", TRIM(out%modl), " ", TRIM(out%filename) })
    M4_WRITE_INFO({"fn mode = ", TRIM(out%fn), " ", TRIM(out%mode) })
    M4_WRITE_INFO({"ns ne dn = ", out%ns, out%ne, out%dn })
    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(out%regidx))

  end subroutine EchoOutObj
  
!----------------------------------------------------------------------


end module outlist

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
