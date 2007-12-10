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
  use grid

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), parameter :: modname = 'REGLIST'

  ! --- Public Methods

  public :: ReadRegObj
  public :: InitializeRegList
  public :: FinalizeRegList
  public :: CreateRegObj
  public :: DestroyRegObj

  ! --- Public Data

  public :: regobj
  public :: numregobj
  public :: T_REG

  ! --- Constants

  integer, parameter :: MAXREGOBJ = 1000

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
     integer :: ps = 0, pe = 0                  ! linear list start, end

  end type T_REG

  ! --- Data

  type(T_REG) :: regobj(MAXREGOBJ) 
  integer :: numregobj

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

  subroutine ReadRegObj(reg, funit)

    integer :: funit
    type(T_REG) :: reg
    type(T_REG), external :: CreateRegObj

    integer :: ios, err
    integer :: i,j,k, i0, i1, di, j0, j1, dj, k0, k1, dk, num, p
    integer :: masksz, listsz
    real(kind=8) :: val
    character(len=STRLNG) :: string
    logical :: auto

    M4_WRITE_DBG(". enter ReadRegObj")

    M4_WRITE_DBG("creating new regobj")
    reg = CreateRegObj()

    M4_WRITE_DBG({ "allocating mask, size = ", GRIDSIZE})
    allocate(tmpmask(IBEG:IEND,JBEG:JEND,KBEG:KEND),tmpval(1:GRIDSIZE),stat = err)
    M4_ALLOC_ERROR(err,"ReadRegObj/reglist")
    tmpmask = 0
    tmpval = 0.0
    numnodes = 0
   
    auto = .true.
    ! read until an unexpected line is encountered, eg ")"
    do

       read(funit,*) string
       
       select case ( string ) 
       case( "(POINT" ) 
          M4_WRITE_DBG({"got ", TRIM(string), " token"})
          reg%isbox = .false.
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
          M4_WRITE_DBG({"got ", TRIM(string), " token"})
          val = 1.0
          do 
             read(funit,*,iostat = ios)  i0, i1, di, j0, j1, dj, k0, k1, dk
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of box list"})
                exit
             else
                reg%isbox = numnodes .eq. 0
             end if
             call SetBoxRegObj(reg, i0, i1, di, j0, j1, dj, k0, k1, dk, val)
          end do
       case( "(VBOX" ) 
          reg%isbox = .false.
          do 
             read(funit,*,iostat = ios)  i0, i1, di, j0, j1, dj, k0, k1, dk, val
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of box list"})
                exit
             else 
                M4_WRITE_DBG({"got a ", TRIM(string), " numnodes = ", numnodes})
             end if
             call SetBoxRegObj(reg, i0, i1, di, j0, j1, dj, k0, k1, dk, val)
          end do
       case( "(VPOINT" ) 
          reg%isbox = .false.
          M4_WRITE_DBG({"got ", TRIM(string), " token"})
          do 
             read(funit,*,iostat = ios) i,j,k, val
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of point list"})
                exit
             end if
             call SetPointRegObj(reg, i,j,k, val)
          end do
       case("LIST") ! force list mode
          M4_WRITE_DBG({"got ", TRIM(string), " token"})
          reg%islist = .true.
          auto = .false.
          reg%isbox = .false.
       case("MASK") ! force mask mode
          M4_WRITE_DBG({"got ", TRIM(string), " token"})
          reg%islist = .false.
          reg%isbox = .false.
          auto = .false.
       case("AUTO") 
       case default
          M4_WRITE_DBG({"read terminator: ", TRIM(string)})
          if ( string(1:1) .ne. ")" ) then
             M4_FATAL_ERROR({"BAD TERMINATOR: ReadRegObj/reglist"})
          end if
          exit
       end select

    end do

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
       
    if ( reg%isbox ) then ! A SINGLE BOX
        
       M4_WRITE_DBG({"setting up reg in single box mode", numnodes})
       reg%numnodes = &
            ( reg%ie - reg%is + 1 ) / reg%di * &
            ( reg%je - reg%js + 1 ) / reg%dj * &
            ( reg%ke - reg%ks + 1 ) / reg%dk    

    else

       reg%numnodes = numnodes
       
       if ( .not. reg%islist ) then ! A MASK

          M4_WRITE_DBG({"setting up reg in mask mode", numnodes})
          allocate(reg%mask(reg%is:reg%ie,reg%js:reg%je,reg%ks:reg%ke),reg%val(numnodes), stat = err )
          M4_ALLOC_ERROR(err, {"ReadRegObj/reg"})
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
             M4_FATAL_ERROR({"BAD NODE COUNT last idx = ", p, " .ne. numnodes = ",numnodes})
          endif
          
       else ! A LIST

          M4_WRITE_DBG({"setting up reg in list mode", numnodes})
          reg%pe = numnodes
          allocate(reg%i(numnodes),reg%j(numnodes),reg%k(numnodes),reg%val(numnodes) , stat = err)
          if ( err .ne. 0 ) then
             M4_FATAL_ERROR({"OUT OF MEMORY: ReadRegObj/reg"})
          end if
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
             M4_FATAL_ERROR({"BAD NODE COUNT last idx = ", num, " .ne. numnodes = ",numnodes})
          endif

       end if
       
    end if

    deallocate(tmpmask)
    deallocate(tmpval)

    M4_WRITE_DBG({" created regobj num = ", numregobj})
    M4_WRITE_DBG({" reg%isbox = ", reg%isbox})
    M4_WRITE_DBG({" reg%islist = ", reg%islist})
    M4_WRITE_DBG({" reg%numnodes = ", reg%numnodes})
    M4_WRITE_DBG({" reg%is reg%ie reg%di = ", reg%is, reg%ie, reg%di})
    M4_WRITE_DBG({" reg%js reg%je reg%dj = ", reg%js, reg%je, reg%dj})
    M4_WRITE_DBG({" reg%ks reg%ke reg%dk = ", reg%ks, reg%ke, reg%dk})

    regobj(numregobj) = reg

    M4_WRITE_DBG(". exit ReadRegObj")

  end subroutine ReadRegObj

!----------------------------------------------------------------------

  type(T_REG) function CreateRegObj

    numregobj = numregobj + 1
    regobj(numregobj)%idx = numregobj

    regobj(numregobj)%is = IMAX
    regobj(numregobj)%ie = IMIN
    regobj(numregobj)%js = JMAX
    regobj(numregobj)%je = JMIN
    regobj(numregobj)%ks = KMAX
    regobj(numregobj)%ke = KMIN
  
    !no points
    regobj(numregobj)%ps = 1
    regobj(numregobj)%pe = 0
    !no steps
    regobj(numregobj)%di = 0
    regobj(numregobj)%dj = 0
    regobj(numregobj)%dk = 0

    regobj(numregobj)%isbox = .false.
    regobj(numregobj)%islist = .false.
    CreateRegObj = regobj(numregobj)
   
  end function CreateRegObj

!----------------------------------------------------------------------

  subroutine SetBoxRegObj(reg, i0,i1,di, j0,j1,dj, k0,k1,dk, val)

    type(T_REG) :: reg
    integer :: i0, i1, di, j0, j1, dj, k0, k1, dk
    real*8 :: val
    integer :: i,j,k,p

    M4_WRITE_DBG({"set box: ",i0, i1, di, j0, j1, dj, k0, k1, dk})

    ! check whether box is inside grid
    if ( i1 .lt. IBEG .or. j1 .lt. JBEG .or. k1 .lt. KBEG .or. &
         i0 .gt. IEND .or. j0 .gt. JEND .or. k0 .gt. KEND ) then
       return
    end if

    ! clip box to grid
    i0 = Max(i0,IBEG)
    i1 = Min(i1,IEND)
    j0 = Max(j0,JBEG)
    j1 = Min(j1,JEND)
    k0 = Max(k0,KBEG)
    k1 = Min(k1,KEND)
    
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

    ! check whether point is inside grid
    if ( i .lt. IBEG .or. j .lt. JBEG .or. k .lt. KBEG .or. &
         i .gt. IEND .or. j .gt. JEND .or. k .gt. KEND ) then
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

    if ( .not. reg%isbox ) then
       if ( .not. reg%islist ) then
          deallocate(reg%mask, reg%val)
       else
          deallocate(reg%i, reg%j, reg%k, reg%val)
       end if
    end if
    
  end subroutine DestroyRegObj

!----------------------------------------------------------------------

end module reglist

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
! 
!======================================================================




