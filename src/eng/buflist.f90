!-*- F90 -*------------------------------------------------------------
!
!  module: buflist / meta3
!
!  this generic module allows to temporarily buffer fields.
!
!  subs:
!
!  InitializeBufList
!  FinalizeBufList
!  CreateBufObj  
!  DestroyBufObj
!  FillBufObj
!  SetBufObj
!  AddBufObj
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

  character(len=20), parameter :: modname = 'BUFLIST'

  ! --- Public Methods

  public :: InitializeBufList
  public :: FinalizeBufList
  public :: CreateBufObj
  public :: DestroyBufObj
  public :: FillBufObj
  public :: SetBufObj
  public :: AddBufObj

  ! --- Public Data

  public :: bufobj
  public :: numbufobj
  public :: T_BUF

  ! --- Types

  type T_BUF

     M4_FTYPE, pointer, dimension(:) :: data          ! field buffer
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
 

  type(T_BUF) function CreateBufObj(reg)

    type(T_REG) :: reg
    integer :: err

    numbufobj = numbufobj + 1
    bufobj(numbufobj)%idx = numbufobj
    bufobj(numbufobj)%regidx = reg%idx
    bufobj(numbufobj)%size = reg%numnodes

    allocate(bufobj(numbufobj)%data(reg%numnodes), stat=err)
    M4_ALLOC_ERROR(err, {"CreatBufObj"})
    CreateBufObj = bufobj(numbufobj)
    
  end function CreateBufObj

!----------------------------------------------------------------------
 
  subroutine DestroyBufObj(buf)

    type(T_BUF) :: buf

    deallocate(buf%data)
    
    buf%size = 0
    buf%regidx = 0
    buf%idx = 0 

  end subroutine DestroyBufObj

!----------------------------------------------------------------------

  subroutine FillBufObj(buf, field)

    type(T_BUF) :: buf
    M4_FTYPE, dimension(:,:,:) :: field

    M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    reg = regobj(buf%regidx)
    M4_REGLOOP_EXPR(reg,p,i,j,k,w,
      buf%data(p) = field(i,j,k)
    )
    
  end subroutine FillBufObj

!----------------------------------------------------------------------

  subroutine AddBufObj(buf, field)

    type(T_BUF) :: buf
    M4_FTYPE, dimension(:,:,:) :: field

    M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    reg = regobj(buf%regidx)
    M4_REGLOOP_EXPR(reg,p,i,j,k,w,
      buf%data(p) = buf%data(p) + field(i,j,k)
    )
    
  end subroutine AddBufObj


!----------------------------------------------------------------------

  subroutine SetBufObj(buf, val)

    type(T_BUF) :: buf
    M4_FTYPE :: val

    M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    reg = regobj(buf%regidx)
    M4_REGLOOP_EXPR(reg,p,i,j,k,w,
      buf%data(p) = val
    )

  end subroutine SetBufObj

!----------------------------------------------------------------------

end module buflist


!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
