!-*- F90 -*------------------------------------------------------------
!
!  module: buflist / meta
!
!  this generic module allows to temporarily buffer fields.
!
!  CF,1D,2D,3D
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

     real(kind=8), pointer, dimension(:,:) :: rdata    ! real field buffer
     complex(kind=8), pointer, dimension(:,:) :: cdata ! complex field buffer
     M4_FTYPE, pointer, dimension(:,:) :: data          
     logical :: iscf = .false.						   ! complex data?			
     integer :: size = 0                               ! size of buffer
     integer :: regidx = 0                             ! region index
     integer :: idx = 0                                ! this objects index
	 integer :: numslot								   ! number of buf slots
     integer :: numval								   ! number of reg values
	 
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
 

  type(T_BUF) function CreateBufObj(reg, iscf, numslot)

    type(T_REG) :: reg
    logical :: iscf
    integer :: numslot, err
    

    numbufobj = numbufobj + 1
    bufobj(numbufobj)%idx = numbufobj
    bufobj(numbufobj)%regidx = reg%idx
    bufobj(numbufobj)%size = reg%numnodes
    bufobj(numbufobj)%numslot = numslot
    bufobj(numbufobj)%numval = reg%numval
    bufobj(numbufobj)%iscf = iscf
	
    if ( iscf ) then
      allocate(bufobj(numbufobj)%cdata(reg%numnodes,numslot), stat=err)
      M4_ALLOC_ERROR(err, {"CreatBufObj"})
    else
      allocate(bufobj(numbufobj)%rdata(reg%numnodes,numslot), stat=err)
      M4_ALLOC_ERROR(err, {"CreatBufObj"})    
    end if
    
	M4_IFELSE_CF({
		bufobj(numbufobj)%data => bufobj(numbufobj)%cdata
    },{
		bufobj(numbufobj)%data => bufobj(numbufobj)%rdata
    })    	
    
    CreateBufObj = bufobj(numbufobj)
    
  end function CreateBufObj

!----------------------------------------------------------------------
 
  subroutine DestroyBufObj(buf)

    type(T_BUF) :: buf

    if ( buf%idx .eq. -1 ) return
    
    if ( buf%iscf ) then
       deallocate(buf%cdata)
    else
       deallocate(buf%rdata)    
    end if
    
    buf%size = 0
    buf%regidx = 0
    buf%idx = -1

  end subroutine DestroyBufObj

!----------------------------------------------------------------------

  subroutine FillBufObj(buf, slot, field)

    type(T_BUF) :: buf
    integer :: slot
    M4_FTYPE, dimension(:,:,:) :: field
    M4_REGLOOP_DECL(reg,p,i,j,k,w(buf%numval))  

    if ( buf%idx .eq. -1 ) return

    reg = regobj(buf%regidx)
    M4_REGLOOP_EXPR(reg,p,i,j,k,w,
       buf%data(p,slot) = field(i,j,k)
    )
    
  end subroutine FillBufObj

!----------------------------------------------------------------------

  subroutine AddBufObj(buf, slot, field)

    type(T_BUF) :: buf
    integer :: slot
    M4_FTYPE, dimension(:,:,:) :: field
    M4_REGLOOP_DECL(reg,p,i,j,k,w(buf%numval))  

    if ( buf%idx .eq. -1 ) return

    reg = regobj(buf%regidx)
    M4_REGLOOP_EXPR(reg,p,i,j,k,w,
      buf%data(p,slot) = buf%data(p,slot) + field(i,j,k)
    )
    
  end subroutine AddBufObj

!----------------------------------------------------------------------

  subroutine SetBufObj(buf, slot, val)

    type(T_BUF) :: buf
    integer :: slot
    M4_FTYPE :: val
    M4_REGLOOP_DECL(reg,p,i,j,k,w(buf%numval))  
    
    if ( buf%idx .eq. -1 ) return

    reg = regobj(buf%regidx)
    M4_REGLOOP_EXPR(reg,p,i,j,k,w,
      buf%data(p,slot) = val
    )

  end subroutine SetBufObj

!----------------------------------------------------------------------

  integer function MemoryBufList()

    integer :: i, memusage
	type(T_BUF) :: buf
	type(T_REG) :: reg
    integer :: sz 
        
    do i = 1, numbufobj 
       
       buf = bufobj(i)
       if ( buf%iscf ) then 
        	sz = 16 
       else 
        	sz = 8
       end if
       reg = regobj(buf%regidx)
       
       if ( buf%idx .eq. -1 .or. reg%numnodes .eq. 0 ) cycle
       
       memusage = memusage + reg%numnodes * buf%numslot * sz

    end do

    MemoryBufList = memusage

  end function MemoryBufList

!----------------------------------------------------------------------

end module buflist


!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
