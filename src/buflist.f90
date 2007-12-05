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
!  FillIntegerBufObj
!  FillRealBufObj
!  FillComplexBufObj
!
!----------------------------------------------------------------------


!======================================================================
!
!

module buflist

  use constant
  use strings
  use reglist

  implicit none
  save

  ! --- Constants

  integer, parameter :: MAXBUFOBJ = 1000

  integer, parameter :: INTEGER_BUF = 1
  integer, parameter :: REAL_BUF = 2
  integer, parameter :: COMPLEX_BUF = 3


  ! --- Types


  type T_BUF

     integer :: type

     integer, pointer, dimension(:) :: idata
     real(kind=8), pointer, dimension(:) :: rdata
     complex(kind=8), pointer, dimension(:) :: cdata

     integer :: size
     integer :: regidx

     integer :: idx

  end type T_BUF


  ! --- Fields

  type(T_BUF) :: bufobj(MAXBUFOBJ) 
  integer :: numbufobj

contains

!----------------------------------------------------------------------

  subroutine InitializeBufList

    numbufobj = 0

  end subroutine InitializeBufList

!----------------------------------------------------------------------

  subroutine FinalizeBufList

    integer :: n

    do n = 1, numbufobj 
       
       call DestroyBufObj(bufobj(n))

    end do


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
    case(INTEGER_BUF) 
       allocate(buf%idata(buf%size), stat=err)
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
    case(INTEGER_BUF) 
       deallocate(buf%idata)
    case(REAL_BUF) 
       deallocate(buf%rdata)
    case(COMPLEX_BUF) 
       deallocate(buf%cdata)
    end select
    
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
