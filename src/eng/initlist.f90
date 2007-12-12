!-*- F90 -*------------------------------------------------------------
! 
!  module: initlist / meta3
!
!  subs:
!
!  InitializeInitList
!  FinalizeInitList
!  ReadInitObj
!  CreateInitObj
!  SetInitObj
!
!----------------------------------------------------------------------


!======================================================================
!
!

module initlist

  use constant
  use strings
  use reglist

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), parameter :: modname = 'INITLIST'

  ! --- Public Methods

  public :: InitializeInitList
  public :: FinalizeInitList
  public :: CreateInitObj
  public :: DestroyInitObj
  public :: ReadInitObj

  ! --- Public Data

  public :: initobj
  public :: numinitobj
  public :: T_INIT

  ! --- Types

  type T_INIT 

    integer :: numval
    character(len=STRLNG) :: type
    real(kind=8) :: val(27)
    integer :: regidx = 0             ! region index
    integer :: idx = 0                ! this objects index

  end type T_INIT

  ! --- Data

  type (T_INIT) :: initobj(MAXINITOBJ)
  integer :: numinitobj = 0

contains

!----------------------------------------------------------------------

  subroutine InitializeInitList
  
    numinitobj = 0

  end subroutine InitializeInitList

!----------------------------------------------------------------------

  subroutine FinalizeInitList
  
    integer :: i
    
    do i = 1, numinitobj 
       call DestroyInitObj(initobj(i))
    end do

    numinitobj = 0

  end subroutine FinalizeInitList

!----------------------------------------------------------------------

  subroutine ReadInitObj(init, regdef, num, type, funit)
  
    type(T_INIT) :: init
    type(T_REG) :: regdef ! default region
    integer :: num
    character(len=*) :: type
    integer :: funit
    type(T_INIT), external :: CreateInitObj
    type(T_REG) :: reg
    character(len=STRLNG) :: string
    integer :: i
    
    init = CreateInitObj()
    init%type = type

    M4_WRITE_DBG({". enter ReadInitObj num = ",init%idx})

    ! read init information
    read(funit,*) (init%val(i),i=1, num)
    M4_WRITE_DBG({"values: ",(init%val(i),i=1, num)})
    read(funit,*) string
    M4_WRITE_DBG({"got token ",TRIM(string)})

    ! consume regobj start string
    if ( string .eq. "(REG" ) then
       M4_WRITE_DBG({"-> ReadRegObj"})
       call ReadRegObj(reg, regdef, funit) ! spatial regobj
       read(funit,*) string
       M4_WRITE_DBG({"got token ",TRIM(string)})
       init%regidx = reg%idx
       init%numval = num
    else
       M4_WRITE_DBG({"using default region!"})
       init%regidx = regdef%idx
       init%numval = num
    end if

    if ( string(1:1) .ne. ")" ) then
       M4_FATAL_ERROR({"BAD TERMINATOR: ReadInitObj/initlist"})
    end if

    initobj(numinitobj) = init

    M4_WRITE_DBG(". exit ReadInitObj")

  end subroutine ReadInitObj

!----------------------------------------------------------------------

  type(T_INIT) function CreateInitObj

    numinitobj = numinitobj + 1
    initobj(numinitobj)%idx = numinitobj
    
    initobj(numinitobj)%val = 0.

    CreateInitObj = initobj(numinitobj)
   
  end function CreateInitObj

!----------------------------------------------------------------------

  subroutine DestroyInitObj(init)

    type(T_INIT) :: init

    init%regidx = 0  
    init%idx = 0
   
  end subroutine DestroyInitObj

!----------------------------------------------------------------------

end module initlist

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
