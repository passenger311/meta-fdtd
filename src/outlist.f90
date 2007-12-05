!-*- F90 -*------------------------------------------------------------
! 
!  module: outlist / max3d
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
  use mpiworld
  use grid
  use reglist

  implicit none
  save

  integer, parameter :: MAXOUTOBJ = 500

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

  type (T_OUT) :: outobj(MAXOUTOBJ)
  integer :: numoutobj = 0


contains

!----------------------------------------------------------------------

  subroutine ReadOutObj(out, funit)

    integer :: funit
    type(T_OUT) :: out
    type(T_OUT), external :: CreateOutObj
    
    character (len=STRLNG) :: fmt, modl, fn, mode, filename, string
    integer :: ns, ne, dn
    type(T_REG) :: reg
    

    out = CreateOutObj()
    ! read output information
    read(funit,*) fmt, modl        ! format and module
    read(funit,*) filename         ! filename
    read(funit,*) fn, mode         ! function and mode
    read(funit,*) ns, ne, dn       ! time frame
    read(funit,*) string
    ! consume regobj start string
    if ( string .ne. "(regobj" ) then
       write(STDERR,*) "!ERROR NO REGION DEFINED: ReadOutObj/out"
       stop
    end if
    ! read the regobj information
    call ReadRegObj(reg, funit) ! spatial regobj
    ! write output parameters into out object
    call SetOutObj(out, fmt, modl, filename, fn, mode, ns, ne, dn, reg)
    ! consume the closing ")"
    read(funit,*) string
    if ( string .ne. "(regobj" ) then
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
