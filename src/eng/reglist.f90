!----------------------------------------------------------------------
!
!  module: reglist / meta3
!
!  a spatial box, mask or pointlist 
!
!  subs:
!
!    InitializeRegList
!    FinalizeRegList
!    CreateRegObj
!    DestroyRegObj
!    ReadRegObj
!    SetPointRegObj
!    SetBoxRegObj
!
!----------------------------------------------------------------------

!======================================================================
!
!  (REG 
!    (POINT
!     100 20 30 
!     100 21 33
!    POINT)
!    (VPOINT
!     100 10 32 0.8
!    )
!    (BOX
!     0 10 20 1 6 20 100 1
!    (VBOX
!     0 10 20 1 6 20 100 0.25 1
!     LIST | MASK | AUTO
!    )
!  REG)
!
!  A regobj is a set of points with a weight attached. BOX, POINT set the 
!  weigth to 1.0, while VPOINT and VBOX allow to set a userdefined weight 
!  value. If only one BOX is defined in a regobj, then no point lists or
!  masks will be allocated. Instead the "isbox" value of the regobj is set
!  to true. 
!
!  In MASK mode each regobj allocates a field of integers within a bounding
!  box (is,js,ks)(ie,je,ke). The field values are zero if the point is not
!  set, otherwise p>0 where the value is the point index to the weight 
!  stored in val(p). 
!
!  In LIST mode a list of coordinates i(p),j(p),k(p) and values val(p) are
!  supplied for each point p=1..pe.
!
!  In AUTO mode (default) either LIST or MASK mode is chosen depending on
!  which one has the lower memory usage. 
! 
!  A loop over a region should always utilize the M4 macros as defined in
!  reglist.m4.


module reglist

  use constant
  use strings

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), parameter :: modname = 'REGLIST'

  ! --- Public Methods

  public :: ReadRegObj
  public :: InitializeRegList
  public :: FinalizeRegList
  public :: CreateRegObjStart
  public :: CreateRegObjEnd
  public :: CreateBoxRegObj
  public :: DestroyRegObj
  public :: EchoRegObj
  public :: SetBoxRegObj
  public :: SetPointRegObj

  ! --- Public Data

  public :: regobj
  public :: numregobj
  public :: T_REG

  ! --- Types

  type T_REG

     logical :: isbox = .false.                 ! is it a box?
     logical :: islist = .false.                ! is it a list?
     integer :: idx = 0                         ! this objects index
     integer :: is = 0, ie = 0, di = 0          ! bounding box start,end
     integer :: js = 0, je = 0, dj = 0
     integer :: ks = 0, ke = 0, dk = 0
     integer :: numnodes                        ! number of points in region
     integer, pointer, dimension(:,:,:) :: mask ! mark set points
     integer, pointer, dimension(:) :: i, j, k  ! coordinate fields
     real(kind=8), pointer, dimension(:) :: val ! weight field
     real(kind=8) :: boxval
     integer :: ps = 0, pe = 0                  ! linear list start, end

  end type T_REG

  ! --- Data

  type(T_REG) :: regobj(MAXREGOBJ) 
  integer :: numregobj

  type(T_REG) :: domreg

  integer, allocatable, dimension(:,:,:) :: tmpmask
  real(kind=8), allocatable, dimension(:) :: tmpval
  integer :: numnodes


contains

!----------------------------------------------------------------------

  subroutine InitializeRegList

    integer :: err

    numregobj = 0

  end subroutine InitializeRegList
  
!----------------------------------------------------------------------

  subroutine FinalizeRegList

    integer :: i
    
    do i = 1, numregobj 
       call DestroyRegObj(regobj(i))
    end do

    numregobj = 0

  end subroutine FinalizeRegList

!----------------------------------------------------------------------

  subroutine ReadRegObj(reg, dom, funit)

    integer :: funit
    type(T_REG) :: reg, dom
    type(T_REG), external :: CreateRegObj

    integer :: ios, err
    integer :: i,j,k, i0, i1, di, j0, j1, dj, k0, k1, dk
    real(kind=8) :: val
    character(len=STRLNG) :: string
    logical :: auto

    M4_WRITE_DBG(". enter ReadRegObj")

    reg = CreateRegObjStart(dom)
   
    auto = .true.
    ! read until an unexpected line is encountered, eg ")"
    do

       read(funit,*) string
    
       select case ( string ) 
       case( "(POINT" ) 
          M4_WRITE_DBG({"got token ", TRIM(string)})
          val = 1.0
          do 
             read(funit,*,iostat = ios) i,j,k
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of point list"})
                exit
             end if
             call SetPointRegObj(reg, i,j,k, val)
          end do
       case( "(BOX" ) 
          M4_WRITE_DBG({"got token ", TRIM(string)})
          val = 1.0
          do 
             read(funit,*,iostat = ios)  i0, i1, di, j0, j1, dj, k0, k1, dk
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of box list"})
                exit
             end if
             call SetBoxRegObj(reg, i0, i1, di, j0, j1, dj, k0, k1, dk, val)
          end do
       case( "(VBOX" ) 
          M4_WRITE_DBG({"got token ", TRIM(string)})
          do 
             read(funit,*,iostat = ios)  i0, i1, di, j0, j1, dj, k0, k1, dk, val
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of box list"})
                exit
             end if
             call SetBoxRegObj(reg, i0, i1, di, j0, j1, dj, k0, k1, dk, val)
          end do
       case( "(VPOINT" ) 
          M4_WRITE_DBG({"got token ", TRIM(string)})
          do 
             read(funit,*,iostat = ios) i,j,k, val
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of point list"})
                exit
             end if
             call SetPointRegObj(reg, i,j,k, val)
          end do
       case("LIST") ! force list mode
          M4_WRITE_DBG({"got token ", TRIM(string)})
          reg%islist = .true.
          auto = .false.
          reg%isbox = .false.
       case("MASK") ! force mask mode
          M4_WRITE_DBG({"got token ", TRIM(string)})
          reg%islist = .false.
          reg%isbox = .false.
          auto = .false.
       case("AUTO") 
       case default
          M4_WRITE_DBG({"got terminator: ", TRIM(string)})
          if ( string(1:1) .ne. ")" ) then
             M4_FATAL_ERROR({"BAD TERMINATOR: ReadRegObj"})
          end if
          exit
       end select

    end do

    call CreateRegObjEnd(reg,auto)

    M4_WRITE_DBG(". exit ReadRegObj")

  end subroutine ReadRegObj

!----------------------------------------------------------------------

  type(T_REG) function CreateBoxRegObj(is,ie,di,js,je,dj,ks,ke,dk)

    integer :: is,ie,di,js,je,dj,ks,ke,dk
    
    M4_WRITE_DBG(". enter CreateBoxRegObj")

    numregobj = numregobj + 1
    regobj(numregobj)%idx = numregobj
    regobj(numregobj)%is = is
    regobj(numregobj)%ie = ie
    regobj(numregobj)%di = di
    regobj(numregobj)%js = js
    regobj(numregobj)%je = je
    regobj(numregobj)%dj = dj
    regobj(numregobj)%ks = ks
    regobj(numregobj)%ke = ke
    regobj(numregobj)%dk = dk
    regobj(numregobj)%isbox = .true.
    regobj(numregobj)%islist = .false.
    regobj(numregobj)%numnodes = ( (ie - is)/di + 1 ) * ( (je - js)/dj + 1 ) * &
         ( (ke - ks)/dk + 1 )

    M4_IFELSE_DBG({call EchoRegObj(regobj(numregobj))})

    M4_WRITE_DBG(". exit CreateBoxRegObj")

    CreateBoxRegObj = regobj(numregobj)

  end function CreateBoxRegObj

!----------------------------------------------------------------------

  type(T_REG) function CreateRegObjStart(dom)

    integer :: err
    type (T_REG) :: dom

    M4_WRITE_DBG(". enter CreateRegObjStart")

    domreg = dom 

    M4_WRITE_DBG(". over domain --->")
    M4_IFELSE_DBG({call EchoRegObj(dom)})

    M4_WRITE_DBG({ "allocating mask, size = ", domreg%numnodes})
    allocate(tmpmask(domreg%is:domreg%ie,domreg%js:domreg%je,domreg%ks:domreg%ke), &
         tmpval(1:domreg%numnodes),stat = err)
    M4_ALLOC_ERROR(err,"ReadRegObj")
    tmpmask = 0
    tmpval = 0.0
    numnodes = 0

    numregobj = numregobj + 1
    regobj(numregobj)%idx = numregobj

    regobj(numregobj)%is = domreg%ie
    regobj(numregobj)%ie = domreg%is
    regobj(numregobj)%js = domreg%je
    regobj(numregobj)%je = domreg%js
    regobj(numregobj)%ks = domreg%ke
    regobj(numregobj)%ke = domreg%ks
  
    !no points
    regobj(numregobj)%ps = 1
    regobj(numregobj)%pe = 0
    !no steps
    regobj(numregobj)%di = 0
    regobj(numregobj)%dj = 0
    regobj(numregobj)%dk = 0

    regobj(numregobj)%isbox = .false.
    regobj(numregobj)%islist = .false.
    CreateRegObjStart = regobj(numregobj)

    M4_WRITE_DBG(". exit CreateRegObjStart")
   
  end function CreateRegObjStart

!----------------------------------------------------------------------

  subroutine CreateRegObjEnd(reg,auto)

    type(T_REG) :: reg
    logical :: auto

    integer :: masksz, listsz, num, i,j,k,p, err
    real(kind=8) :: w

    M4_WRITE_DBG(". enter CreateRegObjEnd")

    if ( numnodes .eq. 0 ) then
       reg%isbox = .true.
       deallocate(tmpmask)
       deallocate(tmpval)
       return
    endif

    ! if the domain is a box then nothing needs to be done. If not
    ! then the list of points needs to be filtered!

    if ( .not. domreg%isbox ) then
       M4_REGLOOP_EXPR(domreg,p,i,j,k,w,{
         tmpmask(i,j,k) = - tmpmask(i,j,k) ! flip sign to mark
       })

       ! erase all other points and count again
       num = 0
       do k = reg%ks, reg%ke
          do j = reg%js, reg%je
             do i = reg%is, reg%ie
                p = tmpmask(i,j,k)
                if ( p .lt. 0 ) then
                   tmpmask(i,j,k) = - tmpmask(i,j,k) ! flip sign back
                else
                   tmpmask(i,j,k) = 0 ! erase all non-flipped
                end if
             end do
          end do
       end do
       if ( num .gt. numnodes ) then
          M4_FATAL_ERROR({"BAD NODE COUNT!"})
       endif
       numnodes = num ! new node number
    end if

    ! decide whether a list or a mask is favorable

    if ( auto  ) then
       M4_WRITE_DBG({"auto mode, trying to find best allocation scheme"})
       masksz = (reg%ie - reg%is + 1) * (reg%je - reg%js + 1) * (reg%ke - reg%ks + 1)
       listsz = numnodes * 3 
       M4_WRITE_DBG({"list size = ", listsz, " / mask size = ", masksz})
       if ( listsz .lt. masksz ) then 
          reg%islist = .true.
       else 
          reg%islist = .false.
       endif
    end if
       
    reg%numnodes = numnodes

    ! allocate a private mask or list and fill it with data from tmpmask, tmpval 

    if ( reg%isbox ) then ! A SINGLE BOX
        
       M4_WRITE_DBG({"setting up reg in single box mode", numnodes})

    else
   
       if ( .not. reg%islist ) then ! A MASK

          M4_WRITE_DBG({"setting up reg in mask mode", numnodes})
          allocate(reg%mask(reg%is:reg%ie,reg%js:reg%je,reg%ks:reg%ke),reg%val(numnodes), stat = err )
          M4_ALLOC_ERROR(err, {"ReadRegObj"})
          reg%mask = 0
          ! 1 point
          reg%pe = 1

          num = 0
          do k = reg%ks, reg%ke
             do j = reg%js, reg%je
                do i = reg%is, reg%ie
                   p = tmpmask(i,j,k)
                   if ( p .gt. 0 ) then
                      num = num + 1
                      reg%mask(i,j,k) = p
                      reg%val(p) = tmpval(p)
                   end if
                end do
             end do
          end do
          if ( num .ne. numnodes ) then
             M4_FATAL_ERROR({"BAD NODE COUNT!"})
          endif
          
       else ! A LIST

          M4_WRITE_DBG({"setting up reg in list mode", numnodes})
          reg%pe = numnodes
          allocate(reg%i(numnodes),reg%j(numnodes),reg%k(numnodes),reg%val(numnodes) , stat = err)
          M4_ALLOC_ERROR(err,{"ReadRegObj"})
          ! 1 step stepping!
          reg%di = reg%ie - reg%is + 1
          reg%dj = reg%je - reg%js + 1
          reg%dk = reg%ke - reg%ks + 1
          
          num = 0
          do k = reg%ks, reg%ke
             do j = reg%js, reg%je
                do i = reg%is, reg%ie
                   p = tmpmask(i,j,k)
                   if ( p .gt. 0 ) then
                      num = num + 1
                      reg%i(num) = i
                      reg%j(num) = j
                      reg%k(num) = k
                      reg%val(num) = tmpval(p)
                   endif
                end do
             end do
          end do
          if ( num .ne. numnodes ) then
             M4_FATAL_ERROR({"BAD NODE COUNT!"})
          endif

       end if
       
    end if

    deallocate(tmpmask)
    deallocate(tmpval)


    M4_WRITE_DBG({" created regobj num = ", numregobj})
    M4_IFELSE_DBG({call EchoRegObj(reg)})

    regobj(numregobj) = reg

    M4_WRITE_DBG(". exit CreateRegObjEnd")

  end subroutine CreateRegObjEnd


! ----------------------------------------------------------------------

  subroutine FilterRegObj(reg, dom)

    type(T_REG) :: reg, dom


  end subroutine FilterRegObj

!----------------------------------------------------------------------


  subroutine SetBoxRegObj(reg, i0,i1,di, j0,j1,dj, k0,k1,dk, val)

    type(T_REG) :: reg
    integer :: i0, i1, di, j0, j1, dj, k0, k1, dk
    real*8 :: val
    integer :: i,j,k,p

    M4_WRITE_DBG({"set box: ",i0, i1, di, j0, j1, dj, k0, k1, dk})

    reg%boxval = val
    reg%isbox = numnodes .eq. 0
    
    ! check whether box is inside grid
    if ( i1 .lt. domreg%is .or. j1 .lt. domreg%js .or. k1 .lt. domreg%ks .or. &
         i0 .gt. domreg%ie .or. j0 .gt. domreg%je .or. k0 .gt. domreg%ke ) then
       M4_WRITE_DBG({"box not in given domain --->"})
       M4_IFELSE_DBG({call EchoRegObj(regobj(numregobj))})
       M4_WRITE_DBG({"return!"})
       return
    end if

    ! clip box to domain
    i0 = Max(i0,domreg%is)
    i1 = Min(i1,domreg%ie)
    j0 = Max(j0,domreg%js)
    j1 = Min(j1,domreg%je)
    k0 = Max(k0,domreg%ks)
    k1 = Min(k1,domreg%ke)
    
    ! resize bounding box
    reg%is = Min(i0,reg%is)
    reg%ie = Max(i1,reg%ie)
    reg%js = Min(j0,reg%js)
    reg%je = Max(j1,reg%je)
    reg%ks = Min(k0,reg%ks)
    reg%ke = Max(k1,reg%ke)

    if ( reg%isbox ) then ! the first box has a stride
       reg%di = di
       reg%dj = dj
       reg%dk = dk
    else ! reset stride to 1,1,1 if multiple boxes etc
       reg%di = 1
       reg%dj = 1
       reg%dk = 1
    end if

    M4_WRITE_DBG({"box now: ",i0, i1, di, j0, j1, dj, k0, k1, dk})

    ! set points
    do k = k0, k1, dk
       do j = j0, j1, dj
          do i = i0, i1, di
             if ( tmpmask(i,j,k) .eq. 0 ) then
                numnodes = numnodes + 1
                tmpmask(i,j,k) = numnodes 
                tmpval(numnodes) = val
             endif
          end do
       end do
    end do

  end subroutine SetBoxRegObj

!----------------------------------------------------------------------

  subroutine SetPointRegObj(reg, i,j,k, val)

    type(T_REG) :: reg
    integer :: i,j,k
    real*8 :: val

    M4_WRITE_DBG({"set point: ",i,j,k})

    reg%boxval = val
    reg%isbox = numnodes .eq. 0

    ! check whether point is inside grid
    if ( i .lt. domreg%is .or. j .lt. domreg%js .or. k .lt. domreg%ks .or. &
         i .gt. domreg%ie .or. j .gt. domreg%je .or. k .gt. domreg%ke ) then
       M4_WRITE_DBG({"point not in given domain ==>"})
       M4_IFELSE_DBG({call EchoRegObj(domreg)})
       M4_WRITE_DBG({"return!"})
       return
    end if

    ! resize bounding box
    reg%is = Min(reg%is, i)
    reg%ie = Max(reg%ie, i)
    reg%js = Min(reg%js, j)
    reg%je = Max(reg%je, j)
    reg%ks = Min(reg%ks, k)
    reg%ke = Max(reg%ke, k)
    reg%di = 1
    reg%dj = 1
    reg%dk = 1

    M4_WRITE_DBG({"point now: ",i,j,k})
 
    ! set point
    if ( tmpmask(i,j,k) .eq. 0 ) then
       numnodes = numnodes + 1
       tmpmask(i,j,k) = numnodes
       tmpval(numnodes) = val
    endif

  end subroutine SetPointRegObj

!----------------------------------------------------------------------

  subroutine DestroyRegObj(reg)

    type(T_REG) :: reg

    M4_WRITE_DBG(". enter DestroyRegObjEnd")

    M4_IFELSE_DBG({call EchoRegObj(regobj(numregobj))})
    if ( .not. reg%isbox .and. numnodes .gt. 0) then
       if ( .not. reg%islist ) then
          deallocate(reg%mask, reg%val)
       else
          deallocate(reg%i, reg%j, reg%k, reg%val)
       end if
    end if

    M4_WRITE_DBG(". exit DestroyRegObjEnd")
    
  end subroutine DestroyRegObj

!----------------------------------------------------------------------

  subroutine EchoRegObj(reg)

    type(T_REG) :: reg

    M4_WRITE_INFO({"--- reg # ", TRIM(i2str(reg%idx)) })
    M4_WRITE_INFO({"isbox = ", reg%isbox })
    M4_WRITE_INFO({"islist = ", reg%islist })
    M4_WRITE_INFO({"numnodes = ", reg%numnodes })
    M4_WRITE_INFO({"is ie di = ", reg%is, reg%ie, reg%di })
    M4_WRITE_INFO({"js je dj = ", reg%js, reg%je, reg%dj })
    M4_WRITE_INFO({"ks ke dk = ", reg%ks, reg%ke, reg%dk })

  end subroutine EchoRegObj

!----------------------------------------------------------------------

end module reglist

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
! 
!======================================================================




