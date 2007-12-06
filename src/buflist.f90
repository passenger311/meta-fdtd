!-*- F90 -*------------------------------------------------------------
!
!  module: buflist / max3d
!
!  this generic module allows to temporarily buffer fields.
!
!  subs:
!
!  InitializeBufList
!  FinalizeBufList
!  CreateBufObj  
!  DestroyBufObj
!  FillRealBufObj
!  FillComplexBufObj
!  SetRealBufObj
!  SetComplexBufObj
!  AddRealBufObj
!  AddComplexBufObj
!
!----------------------------------------------------------------------


!======================================================================
!
! m4 macro-preprocessor runs over this file and replaces
! M4 REGLOOP_DECL -> variable declarations for loop
! M4 REGLOOP_EXPR -> loop over region, do expression
!

module buflist

  use constant
  use strings
  use reglist

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), parameter :: modname = 'buflist'

  ! --- Public Methods

  public :: InitializeBufList
  public :: FinalizeBufList
  public :: CreateBufObj
  public :: DestroyBufObj
  public :: FillRealBufObj
  public :: FillComplexBufObj
  public :: SetRealBufObj
  public :: SetComplexBufObj
  public :: AddRealBufObj
  public :: AddComplexBufObj

  ! --- Public Data

  public :: bufobj
  public :: numbufobj
  public :: T_BUF
  public :: REAL_BUF, COMPLEX_BUF

  ! --- Constants

  integer, parameter :: MAXBUFOBJ = 1000
  integer, parameter :: REAL_BUF = 2
  integer, parameter :: COMPLEX_BUF = 3

  ! --- Types

  type T_BUF

     integer :: type = 0                              ! buffer data type
     integer, pointer, dimension(:) :: idata          ! integer buffer
     real(kind=8), pointer, dimension(:) :: rdata     ! real buffer
     complex(kind=8), pointer, dimension(:) :: cdata  ! complex buffer
     integer :: size = 0                              ! size of buffer
     integer :: regidx = 0                            ! region index
     integer :: idx = 0                               ! this objects index

  end type T_BUF

  ! --- Data

  type(T_BUF) :: bufobj(MAXBUFOBJ) 
  integer :: numbufobj

contains

!----------------------------------------------------------------------

  subroutine InitializeBufList

    numbufobj = 0

  end subroutine InitializeBufList

!----------------------------------------------------------------------

  subroutine FinalizeBufList

    integer :: i

    do i = 1, numbufobj 
       
       call DestroyBufObj(bufobj(i))

    end do

    numbufobj = 0

  end subroutine FinalizeBufList

!----------------------------------------------------------------------
 
  type(T_BUF) function CreateBufObj(type,reg)

    integer :: type
    type(T_REG) :: reg
    type(T_BUF) :: buf
    integer :: err

    buf%regidx = reg%idx
    buf%size = reg%numnodes
    buf%type = type
    numbufobj = numbufobj + 1
    buf = bufobj(numbufobj)
    select case (type) 
!    case(INTEGER_BUF) 
!       allocate(buf%idata(buf%size), stat=err)
    case(REAL_BUF) 
       allocate(buf%rdata(buf%size), stat=err)
    case(COMPLEX_BUF) 
       allocate(buf%cdata(buf%size), stat=err)
    end select
    if ( err .ne. 0 ) then
       write(STDERR,*) "!ERROR OUT OF MEMORY: CreateBufObj/bufobj"
       stop
    endif
    CreateBufObj = buf
    
  end function CreateBufObj

!----------------------------------------------------------------------
 
  subroutine DestroyBufObj(buf)

    type(T_BUF) :: buf

    select case (buf%type) 
 !   case(INTEGER_BUF) 
 !      deallocate(buf%idata)
    case(REAL_BUF) 
       deallocate(buf%rdata)
    case(COMPLEX_BUF) 
       deallocate(buf%cdata)
    end select
    
    buf%size = 0
    buf%regidx = 0
    buf%idx = 0 
    buf%type = 0 

  end subroutine DestroyBufObj

!----------------------------------------------------------------------

  subroutine FillRealBufObj(buf, field)

    type(T_BUF) :: buf
    real(kind=8), dimension(:,:,:) :: field

    M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    reg = regobj(buf%regidx)
    M4_REGLOOP_EXPR(reg,p,i,j,k,w,
      buf%rdata(p) = field(i,j,k)
    )

  end subroutine FillRealBufObj

!----------------------------------------------------------------------

  subroutine FillComplexBufObj(buf, field)

    type(T_BUF) :: buf
    complex(kind=8), dimension(:,:,:) :: field

    M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    reg = regobj(buf%regidx)
    M4_REGLOOP_EXPR(reg,p,i,j,k,w,
      buf%cdata(p) = field(i,j,k)
    )

  end subroutine FillComplexBufObj

!----------------------------------------------------------------------

  subroutine AddRealBufObj(buf, field)

    type(T_BUF) :: buf
    real(kind=8), dimension(:,:,:) :: field

    M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    reg = regobj(buf%regidx)
    M4_REGLOOP_EXPR(reg,p,i,j,k,w,
      buf%rdata(p) = buf%rdata(p) + field(i,j,k)
    )

  end subroutine AddRealBufObj

!----------------------------------------------------------------------

  subroutine AddComplexBufObj(buf, field)

    type(T_BUF) :: buf
    complex(kind=8), dimension(:,:,:) :: field

    M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    reg = regobj(buf%regidx)
    M4_REGLOOP_EXPR(reg,p,i,j,k,w,
      buf%cdata(p) = buf%cdata(p) + field(i,j,k)
    )

  end subroutine AddComplexBufObj

!----------------------------------------------------------------------

  subroutine SetRealBufObj(buf, val)

    type(T_BUF) :: buf
    real(kind=8) :: val

    M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    reg = regobj(buf%regidx)
    M4_REGLOOP_EXPR(reg,p,i,j,k,w,
      buf%rdata(p) = val
    )

  end subroutine SetRealBufObj

!----------------------------------------------------------------------

  subroutine SetComplexBufObj(buf, val)

    type(T_BUF) :: buf
    complex(kind=8) :: val

    M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    reg = regobj(buf%regidx)
    M4_REGLOOP_EXPR(reg,p,i,j,k,w,
      buf%cdata(p) = val
    )

  end subroutine SetComplexBufObj


!----------------------------------------------------------------------

end module buflist


!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
