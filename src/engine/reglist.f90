!----------------------------------------------------------------------
!
!  module: reglist / meta
!
!  a spatial box, mask or pointlist 
!
!----------------------------------------------------------------------

!======================================================================
!
!  Region definitions:
!
!  (REG 
!    (POINT
!     100 20 30 
!     100 21 33
!    )POINT
!    (VPOINT
!     100 10 32 0.8
!    )VPOINT
!    (BOX
!     0 10 20 1 6 20 100 1
!    )BOX
!    (VBOX
!     0 10 20 1 6 20 100 0.25 1
!    )VBOX
!  )REG
!
!  A regobj is a set of points with float values. BOX, POINT set the 
!  values to 0.0, while VPOINT and VBOX allow to set userdefined 
!  values. If only one BOX is defined in a regobj, then no point lists or
!  masks will be used. Instead the "isbox" value of the regobj is set
!  to true. 
!
!  Each regobj allocates a field of integers within a bounding box 
!  (is,js,ks)(ie,je,ke). The field values are zero if the point is not
!  set, otherwise p>0 where the value is the point index to the weight 
!  stored in val(p). 
!
!  A list of coordinates i(p),j(p),k(p) and values val(p) are supplied for 
!  each point p=1..pe.
!
!  A loop over a region should always utilize the M4 macros as defined in
!  reglist.m4.
!
!  Value Initializers and Load:
!
!  The definitions of points or boxes increases the numnodes value. It is 
!  possible to initialize all values of all defined nodes to the given point
!  to a given set of values using:
!
!  (FILL
!    ... value set ...
!  )FILL
!  
!  A node-by-node initialization can be done with
!
!  (SET
!     ... value set 1 ...  
!     ... value set 2 ...  
!     ... value set 3 ...  
!     ... value set 4 ...  
!     ...
!  )SET
!
!  As lists of values may grow or multiple values need to be initialized,
!  part or all of the object initialization can be done in an external file,
!  which can be loaded with
!
!  (LOAD
!    filename.dat
!  )LOAD
!
!  where filename.dat has the same format as any section in reglist. In 
!  particular:
!
!  INIT
!  (BOX
!   1 100 1 1 100 1 1 100 1
!  )BOX
!  (SET 
!     ... 100 x 100 x 100 lines of values ...
!  )SET
!  INIT
!
!  Init makes sure that the pointer to the next value to be initialized
!  is aligned to the new box object.
!


module reglist

  use constant
  use strings
  use parse

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
  public :: CompressValRegObj
  public :: DestroyRegObj
  public :: EchoRegObj
  public :: DisplayRegObj
  public :: SetBoxRegObj
  public :: SetPointRegObj
  public :: SetMaskRegObj

  ! --- Public Data

  public :: regobj
  public :: numregobj
  public :: T_REG

  ! --- Types

  type T_REG

     logical :: isbox = .false.                 ! is it a box?
     logical :: islist = .false.                ! list loop mode
     integer :: idx = 0                         ! this objects index
     integer :: is = 0, ie = 0, di = 0          ! bounding box start,end
     integer :: js = 0, je = 0, dj = 0
     integer :: ks = 0, ke = 0, dk = 0
     integer :: numnodes                        ! number of points in region
     integer, pointer, dimension(:,:,:) :: mask ! mark set points
     integer, pointer, dimension(:) :: i, j, k  ! coordinate fields
     integer :: numval                          ! number of values
     real(kind=8),pointer,dimension(:,:) :: val ! values field
     integer :: ps = 0, pe = 0                  ! linear list start, end

     logical :: compressval
     integer, pointer,dimension(:) :: valptr

  end type T_REG

  ! --- Data

  type(T_REG) :: regobj(MAXREGOBJ) 
  integer :: numregobj

  type(T_REG) :: domreg

  integer, allocatable, dimension(:,:,:) :: tmpmask
  real(kind=8), allocatable, dimension(:,:) :: tmpval
  integer :: numnodes, numvalues


contains

!----------------------------------------------------------------------

  subroutine InitializeRegList

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

  subroutine ReadRegObj(reg, dom, funit, lcount, numval)

    integer :: funit, lcount, numval
    type(T_REG) :: reg, dom

    logical :: eof, err 
    integer :: ios, unit

    integer :: pvec(3), bvec(9), defpvec(3), defbvec(9)
    real(kind=8) :: fvec(6)
    character(len=LINELNG) :: string, line, loadfn
    character :: bc, ec
    integer :: lcstack(0:10), ldepth = 0
    logical :: auto

    M4_WRITE_DBG(". enter ReadRegObj")

    ! set default point and box coordinates

    defpvec(1) = dom%is
    defpvec(2) = dom%js
    defpvec(3) = dom%ks

    defbvec(1) = dom%is
    defbvec(2) = dom%ie
    defbvec(3) = 1

    defbvec(4) = dom%js
    defbvec(5) = dom%je
    defbvec(6) = 1

    defbvec(7) = dom%ks
    defbvec(8) = dom%ke
    defbvec(9) = 1

    unit = funit
    reg = CreateRegObjStart(dom,numval)
   
    auto = .false.
    reg%islist = .false.

    do

       err = .false.
       call readline(unit,lcount,eof,line)

       if ( eof ) then
          if ( unit .gt. UNITTMP ) then
             M4_WRITE_DBG({"end of file, close unit = ",TRIM(i2str(unit))})
             close(unit)
             unit = unit - 1 ! return from nested load
             lcount = lcstack(unit-UNITTMP)+1 ! restore line count
             M4_WRITE_INFO({"return to config (",TRIM(i2str(unit-UNITTMP)),")"})
             call readline(unit,lcount,eof,line)
             M4_EOF_ERROR({err},lcount)
             call gettoken(line,")LOAD",err)
             M4_SYNTAX_ERROR({err},lcount,{")LOAD"})
             M4_WRITE_DBG({"consumed )LOAD"})          
             cycle
          else
             M4_EOF_ERROR({err},lcount)
          end if
       end if

       call getstring(line,string,err)
       M4_WRITE_DBG({"got token ",TRIM(string)})
       M4_SYNTAX_ERROR({err},lcount,"[STRING]")
           
       select case ( string ) 
       case( "(POINT" ) 
          do 
             pvec = defpvec
             fvec = 1.0
             call readline(unit,lcount,eof,line)
             M4_EOF_ERROR({eof},lcount)
             call getintvec(line, pvec, 3, ":", err)   ! read up to 3 ints till the : 
             call getfloatvec(line, fvec, 6, ":", err) ! then read up to 6 floats
             if ( err ) then
                M4_SYNTAX_ERROR({line .ne. ")POINT"},lcount,{")POINT"})
                M4_WRITE_DBG({"end of point list"})
                exit
             end if
             call SetPointRegObj(reg, pvec, fvec)
          end do
       case( "(BOX" ) 
          do 
             bvec = defbvec
             fvec = 1.0
             call readline(unit,lcount,eof,line)
             M4_EOF_ERROR({eof},lcount)
             call getintvec(line, bvec, 9, ":", err)   ! read up to 9 ints till the : 
             call getfloatvec(line, fvec, 6, ":", err) ! then read up to 6 floats
             if ( err ) then
                M4_WRITE_DBG({"end of box list"})
                M4_SYNTAX_ERROR({line .ne. ")BOX" },lcount,{")BOX"})
                exit
             end if
             call SetBoxRegObj(reg, bvec, fvec)
          end do
       case("(LOAD")    
          read(unit,*,iostat = ios) loadfn
          if ( ios .ne. 0 ) then
             M4_FATAL_ERROR({"EXPECTED FILENAME AFTER LOAD!"})
          end if
          if ( unit .gt. UNITTMP+10 ) then
          	M4_FATAL_ERROR({"TOO MANY NESTED LOADS!"})
          end if
          M4_WRITE_DBG({"trying to open ",TRIM(loadfn)})
          open(unit+1,FILE=loadfn,STATUS="old", IOSTAT=ios)
          if ( ios .eq. 0 ) then
             lcstack(unit-UNITTMP) = lcount
             unit = unit + 1
             lcount = 1
             M4_WRITE_INFO({"opening config (",TRIM(i2str(unit-10)),"): ", TRIM(loadfn)})
             M4_WRITE_DBG({"file unit = ",TRIM(i2str(unit))})
          else
             M4_WRITE_WARN({"COULD NOT OPEN ",TRIM(loadfn)})
             read(unit,*,iostat = ios) string
             M4_WRITE_DBG({"consumed ",TRIM(string)}) ! consume )LOAD
          end if
       case("CLEAR")    
          numvalues = 1 	
       case("SAVE")    
          numvalues = numnodes + 1 	
       case("(FILL") ! a list of values
          fvec = 1.0
          call readline(unit,lcount,eof,line)
          M4_EOF_ERROR({eof},{lcount})
          call getfloatvec(line, fvec, 6, ":", err) ! then read up to 6 floats
          M4_SYNTAX_ERROR({err},lcount,"[FLOATS]")
          call FillValueRegObj(reg, fvec)
          call readtoken(unit,lcount,")FILL")
       case("(SET") ! a list of values
          do 
             fvec = 1.0
             call readline(unit,lcount,eof,line)
             M4_EOF_ERROR({eof},lcount)
             call getfloatvec(line, fvec, 6, ":", err) ! then read up to 6 floats
             if ( err ) then
                M4_WRITE_DBG({"end of set list"})
                M4_SYNTAX_ERROR({line .ne. ")SET"},lcount,{")SET"})
                exit
             end if
             call SetValueRegObj(reg, fvec)
          end do
       case ("LIST")
          auto = .false.
          reg%islist = .true.
       case ("MASK")
          auto = .false.
          reg%islist = .false.
       case ("AUTO")
          auto = .true.
       case(")REG")
          exit
       case default
          M4_BADTOKEN_ERROR({.true.},lcount,string)
       end select

    end do

    call CreateRegObjEnd(reg,auto)
    
    M4_WRITE_DBG(". exit ReadRegObj")

  end subroutine ReadRegObj

!----------------------------------------------------------------------

  type(T_REG) function CreateBoxRegObj(is,ie,di,js,je,dj,ks,ke,dk)

    integer :: is,ie,di,js,je,dj,ks,ke,dk,err
    
    M4_WRITE_DBG(". enter CreateBoxRegObj")
    
    numregobj = numregobj + 1
    regobj(numregobj)%compressval = .false.
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
    regobj(numregobj)%numnodes = ( (ie - is)/di + 1 ) * ( (je - js)/dj + 1 ) * &
         ( (ke - ks)/dk + 1 )
    regobj(numregobj)%numval = 0

    allocate(regobj(numregobj)%val(regobj(numregobj)%numval,regobj(numregobj)%numnodes) , stat = err)
    M4_ALLOC_ERROR(err, {"CreateBoxRegObj"})
          
    M4_IFELSE_DBG({call EchoRegObj(regobj(numregobj))})

    M4_WRITE_DBG(". exit CreateBoxRegObj")

    CreateBoxRegObj = regobj(numregobj)

  end function CreateBoxRegObj

!----------------------------------------------------------------------

  type(T_REG) function CreateRegObjStart(dom,numval)

    integer :: err, numval
    type (T_REG) :: dom

    M4_WRITE_DBG(". enter CreateRegObjStart")

    domreg = dom 

    M4_WRITE_DBG(". over domain --->")
    M4_IFELSE_DBG({call EchoRegObj(dom)})

    M4_WRITE_DBG({ "allocating mask/value field, size = ", domreg%numnodes})
    allocate(tmpmask(domreg%is:domreg%ie,domreg%js:domreg%je,domreg%ks:domreg%ke), &
         tmpval(1:numval,1:domreg%numnodes),stat = err)
    M4_ALLOC_ERROR(err,"CreateRegObjStart")

    ! initialize mask/values to zero
    tmpmask = 0
    tmpval = 0.0
    numnodes = 0
    numvalues = 1

    numregobj = numregobj + 1
    regobj(numregobj)%idx = numregobj

    ! initialize regobj
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

    regobj(numregobj)%numval = numval
    regobj(numregobj)%isbox = .false.
    regobj(numregobj)%compressval = .false.
    CreateRegObjStart = regobj(numregobj)

    M4_WRITE_DBG(". exit CreateRegObjStart")
   
  end function CreateRegObjStart

!----------------------------------------------------------------------

  subroutine CreateRegObjEnd(reg,auto)

    type(T_REG) :: reg
    logical :: auto

    integer :: masksz, listsz, num, i,j,k,p,err
    real(kind=8) :: w(reg%numval)

    M4_WRITE_DBG(". enter CreateRegObjEnd")

    if ( numnodes .eq. 0 ) then

       reg%isbox = .true.
       deallocate(tmpmask)
       deallocate(tmpval)

       M4_WRITE_DBG({" created regobj num = ", numregobj})
       M4_IFELSE_DBG({call EchoRegObj(reg)})
       
       regobj(numregobj) = reg
       
       M4_WRITE_DBG(". exit CreateRegObjEnd")

       return

    endif

    if ( auto ) then ! try to figure out the best loop mode
       masksz =  (reg%ie - reg%is + 1) * (reg%je - reg%js + 1) * (reg%ke - reg%ks + 1)
       listsz = reg%numnodes * 3
       if ( masksz .gt. listsz ) then 
          reg%islist = .false.
       else
          reg%islist = .true.
       end if
    end if


    ! if the domain is a box then nothing needs to be done. If not
    ! then the list of points needs to be filtered!

    if ( .not. domreg%isbox ) then

       do p  = 1, domreg%pe
          i = domreg%i(p)
          j = domreg%j(p)
          k = domreg%k(p)
          tmpmask(i,j,k) = - tmpmask(i,j,k) ! flip sign to mark
       end do

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

    reg%numnodes = numnodes

    allocate(reg%val(reg%numval,numnodes) , stat = err)
    M4_ALLOC_ERROR(err, {"CreateRegObjEnd"})

    ! allocate a private mask or list and fill it with data from tmpmask, tmpval 

    if ( reg%isbox ) then ! A SINGLE BOX
        
       M4_WRITE_DBG({"setting up reg in single box mode", numnodes})
       num = 0
       do k = reg%ks, reg%ke, reg%dk
         do j = reg%js, reg%je, reg%dj
            do i = reg%is, reg%ie, reg%di
               num = num + 1
               reg%val(:,num) = tmpval(:,num)
            end do
         end do
      end do

    else ! A MASK/LIST
   
       M4_WRITE_DBG({"setting up reg in list/mask mode", numnodes})
       allocate(reg%mask(reg%is:reg%ie,reg%js:reg%je,reg%ks:reg%ke),stat = err ) ! allocate mask
       M4_ALLOC_ERROR(err, {"CreateRegObjEnd"})
       reg%pe = numnodes
       allocate(reg%i(numnodes),reg%j(numnodes),reg%k(numnodes),reg%val(reg%numval,numnodes),stat = err) ! allocate coord list
       M4_ALLOC_ERROR(err,{"CreateRegObjEnd"})
       
       reg%mask = 0 ! empty mask

       num = 0
       do k = reg%ks, reg%ke
          do j = reg%js, reg%je
             do i = reg%is, reg%ie
                p = tmpmask(i,j,k)
                if ( p .gt. 0 ) then
                   num = num + 1
                   reg%mask(i,j,k) = p ! set mask value
                   reg%i(num) = i ! set coords
                   reg%j(num) = j
                   reg%k(num) = k
                   reg%val(:,p) = tmpval(:,p) ! set values
                end if
             end do
          end do
       end do

       if ( num .ne. numnodes ) then
          M4_FATAL_ERROR({"BAD NODE COUNT!"})
       endif
       
    end if

    deallocate(tmpmask)
    deallocate(tmpval)

    M4_WRITE_DBG({" created regobj num = ", numregobj})
    M4_IFELSE_DBG({call EchoRegObj(reg)})

    regobj(numregobj) = reg

    M4_WRITE_DBG(". exit CreateRegObjEnd")

  end subroutine CreateRegObjEnd


!----------------------------------------------------------------------


  subroutine SetBoxRegObj(reg, bvec, fvec)

    type(T_REG) :: reg
    integer :: bvec(9)
    real(kind=8) :: fvec(6)
    integer :: i,j,k,p,v

    if ( reg%idx .eq. -1 ) return

    if ( bvec(3) .eq. 0 .or. bvec(6) .eq. 0 .or. bvec(9) .eq. 0 ) then 
       M4_WRITE_DBG({"box has zero stride -> return"})
       return
    end if
    
    M4_WRITE_DBG({"set box: ",bvec(1),bvec(2),bvec(3),bvec(4),bvec(5),bvec(6),bvec(7),bvec(8),bvec(9)})

    reg%isbox = numnodes .eq. 0
    
    ! check whether box is inside grid
    if ( bvec(2) .lt. domreg%is .or. bvec(5) .lt. domreg%js .or. bvec(8) .lt. domreg%ks .or. &
         bvec(1) .gt. domreg%ie .or. bvec(4) .gt. domreg%je .or. bvec(7) .gt. domreg%ke ) then
       M4_WRITE_DBG({"box not in given domain --->"})
       M4_IFELSE_DBG({call EchoRegObj(reg)})
       M4_WRITE_DBG({"return!"})
       return
    end if

    ! clip box to domain
    bvec(1) = Max(bvec(1),domreg%is)
    bvec(2) = Min(bvec(2),domreg%ie)
    bvec(4) = Max(bvec(4),domreg%js)
    bvec(5) = Min(bvec(5),domreg%je)
    bvec(7) = Max(bvec(7),domreg%ks)
    bvec(8) = Min(bvec(8),domreg%ke)
    
    ! resize bounding box
    reg%is = Min(bvec(1),reg%is)
    reg%ie = Max(bvec(2),reg%ie)
    reg%js = Min(bvec(4),reg%js)
    reg%je = Max(bvec(5),reg%je)
    reg%ks = Min(bvec(7),reg%ks)
    reg%ke = Max(bvec(8),reg%ke)

    if ( reg%isbox ) then ! the first box has a stride
       reg%di = bvec(3)
       reg%dj = bvec(6)
       reg%dk = bvec(9)
    else ! reset stride to 1,1,1 if multiple boxes etc
       reg%di = 1
       reg%dj = 1
       reg%dk = 1
    end if

    M4_WRITE_DBG({"box now: ",bvec(1),bvec(2),bvec(3),bvec(4),bvec(5),bvec(6),bvec(7),bvec(8),bvec(9)})

    ! set points
    do k = bvec(7), bvec(8), bvec(9)
       do j = bvec(4), bvec(5), bvec(6)
          do i = bvec(1), bvec(2), bvec(3)
             p = tmpmask(i,j,k)
             if ( p .eq. 0 ) then
                numnodes = numnodes + 1
                tmpmask(i,j,k) = numnodes 
                do v = 1, reg%numval
                   tmpval(v,numnodes) = fvec(v)
                end do
             else
                do v = 1, reg%numval
                   tmpval(v,p) = fvec(v)
                end do
             end if
          end do
       end do
    end do

    M4_WRITE_DBG({"box set!"})

  end subroutine SetBoxRegObj

!----------------------------------------------------------------------

  subroutine SetPointRegObj(reg, pvec, fvec)

    type(T_REG) :: reg
    integer :: pvec(3),p, v
    real(kind=8) :: fvec(6)

    if ( reg%idx .eq. -1 ) return

    M4_WRITE_DBG({"set point: ",pvec(1),pvec(2),pvec(3)})

    reg%isbox = numnodes .eq. 0

    ! check whether point is inside grid
    if ( pvec(1) .lt. domreg%is .or. pvec(2) .lt. domreg%js .or. pvec(3) .lt. domreg%ks .or. &
         pvec(1) .gt. domreg%ie .or. pvec(2) .gt. domreg%je .or. pvec(3) .gt. domreg%ke ) then
       M4_WRITE_DBG({"point not in given domain ==>"})
       M4_IFELSE_DBG({call EchoRegObj(domreg)})
       M4_WRITE_DBG({"return!"})
       return
    end if

    ! resize bounding box
    reg%is = Min(reg%is, pvec(1))
    reg%ie = Max(reg%ie, pvec(1))
    reg%js = Min(reg%js, pvec(2))
    reg%je = Max(reg%je, pvec(2))
    reg%ks = Min(reg%ks, pvec(3))
    reg%ke = Max(reg%ke, pvec(3))
    reg%di = 1
    reg%dj = 1
    reg%dk = 1

    M4_WRITE_DBG({"point now: ",pvec(1),pvec(2),pvec(3)})
 
    ! set point
    p = tmpmask(pvec(1),pvec(2),pvec(3))
    if ( p .eq. 0 ) then
       numnodes = numnodes + 1
       tmpmask(pvec(1),pvec(2),pvec(3)) = numnodes
       do v = 1, reg%numval
          tmpval(v,numnodes) = fvec(v)
       end do
    else ! override weights
       do v = 1, reg%numval
          tmpval(v,numnodes) = fvec(v)
       end do
    end if
    M4_WRITE_DBG({"point set!"})

  end subroutine SetPointRegObj

!----------------------------------------------------------------------

  subroutine SetValueRegObj(reg, fvec)
  
    type(T_REG) :: reg
    real(kind=8) :: fvec(6)
    integer :: v,i
    
    if ( reg%idx .eq. -1 ) return

    M4_WRITE_DBG({"set value at ",numvalues," of ",numnodes," nodes"})
	
    if ( numvalues .gt. numnodes ) return
    
    do v = 1, reg%numval
       tmpval(v,numvalues) = fvec(v)
    end do

    numvalues = numvalues + 1

    	
  end subroutine SetValueRegObj

!----------------------------------------------------------------------
  
  subroutine FillValueRegObj(reg, fvec)
  
    type(T_REG) :: reg
    real(kind=8) :: fvec(6)
    integer :: v,p
    
    M4_WRITE_DBG({"fill values from ",numvalues," to ",numnodes," nodes"})
    
    do p = numvalues, numnodes	
       do v = 1, reg%numval
          tmpval(v,p) = fvec(v)
       end do
    end do
   
   numvalues = numnodes + 1
	
  end subroutine FillValueRegObj

!----------------------------------------------------------------------

  subroutine DestroyRegObj(reg)

    type(T_REG) :: reg

    if ( reg%idx .eq. -1 ) return
    
    if ( reg%numnodes .eq. 0 ) return

    M4_WRITE_DBG(". enter DestroyRegObjEnd")

    M4_IFELSE_DBG({call EchoRegObj(reg)})

    M4_WRITE_DBG("deallocating reg%val")
    if ( reg%numval .gt. 0 ) then
       deallocate(reg%val)
    end if

    if ( reg%compressval ) then
       M4_WRITE_DBG("deallocating reg%valptr")
       deallocate(reg%valptr)
    end if

    if ( .not. reg%isbox ) then
       M4_WRITE_DBG("deallocating reg%mask")
       deallocate(reg%mask)
       M4_WRITE_DBG("deallocating reg%i reg%j reg%k")
       deallocate(reg%i, reg%j, reg%k)
    end if
    
    reg%idx = -1
    
    M4_WRITE_DBG(". exit DestroyRegObjEnd")
    
  end subroutine DestroyRegObj

!----------------------------------------------------------------------

! materials read values sets which are mostly 1.0,1.0,1.0 (cell 
! filling factor). To reduce memory usage indirect addressing with is 
! valptr is used.

  subroutine CompressValRegObj(reg)

    type(T_REG) :: reg
    integer :: p,v, p_new, err
    logical :: isunit
    real(kind=8), pointer, dimension(:,:) :: newval
    
    if ( reg%numval .eq. 0 .or. reg%numnodes .eq. 0 ) return

    allocate(reg%valptr(1:reg%numnodes), stat = err)
    M4_ALLOC_ERROR(err,{"CompressValRegObj"})

    p_new = 1
    do p = 1, reg%numnodes
       
       isunit = .true.
       do v = 1, reg%numval
          if ( reg%val(v,p) .ne. 1.0 ) then
             isunit = .false.
             exit ! break out
          end if
       end do

       if ( isunit ) then 
          reg%valptr(p) = 1 ! unit tuple will be placed at pos 1
       else
          p_new = p_new + 1
          reg%valptr(p) = p_new ! new pos of value in field
       end if

    end do

    ! new value field
    allocate(newval(1:reg%numval,1:p_new), stat = err)
    M4_ALLOC_ERROR(err,{"CompressValRegObj"})

    do p = 1, reg%numnodes
       newval(:,reg%valptr(p)) = reg%val(:,p)
    end do

    ! free old value field
    deallocate(reg%val)

    ! set to new compressed value field

    reg%val => newval

    reg%compressval = .true.

    regobj(reg%idx) = reg

  end subroutine CompressValRegObj

!----------------------------------------------------------------------

  subroutine SetMaskRegObj(reg,mask,ib,ie,jb,je,kb,ke)

    integer:: ib,ie,jb,je,kb,ke
    logical :: mask(ib:ie,jb:je,kb:ke)
    
    type(T_REG) :: reg
    integer :: p,i,j,k

    if ( reg%numnodes .le. 0 ) return

    if ( reg%isbox ) then

       do k = reg%ks, reg%ke, reg%dk
          do j = reg%js, reg%je, reg%dj
             do i = reg%is, reg%ie, reg%di
                
                mask(i,j,k) = .true.

             enddo !i
          enddo !j
       enddo !k
       
    else !isbox

       do p = 1, reg%pe ! -> p

          i = reg%i(p) ! -> i
          j = reg%j(p) ! -> j
          k = reg%k(p) ! -> k

          mask(i,j,k) = .true.

       enddo ! p

    end if

  end subroutine SetMaskRegObj

!----------------------------------------------------------------------

  subroutine DisplayRegObj(reg)

    type(T_REG) :: reg
	character(len=20) :: tstr
    
    if ( reg%idx .eq. -1 ) return

    if ( reg%isbox ) then
       tstr = "B"
    else
       tstr = "P"
    end if
    
    M4_WRITE_INFO({"@",TRIM(i2str(reg%idx)),&
    	" ",TRIM(tstr)," [",TRIM(i2str(reg%is)),":",TRIM(i2str(reg%ie)),"][", &
    	TRIM(i2str(reg%js)),":",TRIM(i2str(reg%je)),"][", &
    	TRIM(i2str(reg%ks)),":",TRIM(i2str(reg%ke)),"] #",TRIM(i2str(reg%numnodes)) })

  end subroutine DisplayRegObj
    
!----------------------------------------------------------------------

  subroutine EchoRegObj(reg)

    type(T_REG) :: reg

    if ( reg%idx .eq. -1 ) return

    M4_WRITE_INFO({"--- reg # ", TRIM(i2str(reg%idx)) })
    M4_WRITE_INFO({"isbox = ", reg%isbox })
    M4_WRITE_INFO({"numnodes = ", reg%numnodes })
    M4_WRITE_INFO({"numval   = ", reg%numval })
    M4_WRITE_INFO({"is ie di = ", reg%is, reg%ie, reg%di })
    M4_WRITE_INFO({"js je dj = ", reg%js, reg%je, reg%dj })
    M4_WRITE_INFO({"ks ke dk = ", reg%ks, reg%ke, reg%dk })
    M4_WRITE_INFO({"compressval = ", reg%compressval })

  end subroutine EchoRegObj

!----------------------------------------------------------------------

  integer function MemoryRegList()
    
    integer :: i, masksz
    integer :: memusage = 0
    type(T_REG) :: reg
    
    do i = 1, numregobj 
       
       reg = regobj(i)
       if ( reg%idx .eq. -1 .or. reg%numnodes .eq. 0 ) cycle
       
       memusage = memusage + reg%numnodes * reg%numval * 8
       
       if ( .not. reg%isbox ) then
          masksz =  (reg%ie - reg%is + 1) * (reg%je - reg%js + 1) * (reg%ke - reg%ks + 1)
          memusage = memusage +  masksz * 4
          memusage = memusage + reg%numnodes * 3 * 4
       end if
       
    end do
    
    MemoryRegList = memusage
    
  end function MemoryRegList
  
  
  
!----------------------------------------------------------------------

end module reglist

!
! Authors:  J.Hamm 
! Modified: 27/12/2007
! 
!======================================================================




