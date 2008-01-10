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
!    LIST | MASK | AUTO
!  )REG
!
!  A regobj is a set of points with float values. BOX, POINT set the 
!  values to 0.0, while VPOINT and VBOX allow to set userdefined 
!  values. If only one BOX is defined in a regobj, then no point lists or
!  masks will be used. Instead the "isbox" value of the regobj is set
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

  subroutine ReadRegObj(reg, dom, funit, numval)

    integer :: funit, numval
    type(T_REG) :: reg, dom

    integer :: ios,err, unit 
    integer :: v, i,j,k, i0, i1, di, j0, j1, dj, k0, k1, dk
    real(kind=8),dimension(:),allocatable :: val,dval
    character(len=STRLNG) :: string, line, loadfn, skiptill
    logical :: auto

    M4_WRITE_DBG(". enter ReadRegObj")

    unit = funit
    reg = CreateRegObjStart(dom,numval)
   
    allocate(val(1:numval),dval(1:numval), stat=err)
    M4_ALLOC_ERROR(err,{"ReadRegObj"})
    
    ! the default value is (1.,0.,0.,0. ...)
    dval = 0.0
    if ( numval .ne. 0 ) dval(1) = 1.0
    
    auto = .true.
    ! read until an unexpected line is encountered, eg ")"
    skiptill = ""
    do
  	   ! defaults for 1D or 2D
       k0 = 0
       k1 = 0
       k =  0
       dk = 1
       j0 = 0
       j1 = 0
       j  = 0
       dj = 1
       read(unit,*,iostat=ios) line
       if ( ios .ne. 0 ) then
          if ( unit .gt. UNITTMP ) then
             M4_WRITE_DBG({"end of file, close unit = ",TRIM(i2str(unit))})
             close(unit)
             unit = unit - 1 ! return from nested load   
             M4_WRITE_DBG({"back to unit = ",TRIM(i2str(unit))})
             read(unit,*) line ! consume )LOAD
             M4_WRITE_DBG({"consumed ",TRIM(string)}) 
             cycle
          else
             M4_FATAL_ERROR("END OF FILE IN SECTION!")		 
          end if
       end if
       string = TRIM(ADJUSTL(line))
       M4_WRITE_DBG({"got token ", TRIM(string)})
           
       if ( skiptill .ne. "" ) then 
         M4_WRITE_DBG({"skipping line ",TRIM(string)})
         if ( string .eq. skiptill ) skiptill = ""  
         cycle              
       end if
      
       select case ( string ) 
       case( "(POINT" ) 
          do 
             read(unit,*,iostat = ios) M4_READCOORD(i,j,k)
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of point list"})
                exit
             end if
             call SetPointRegObj(reg, i,j,k, dval)
          end do
       case( "(VPOINT" ) 
          do 
             read(unit,*,iostat = ios)  M4_READCOORD(i,j,k), (val(v),v=1,numval)
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of point list"})
                exit
             end if
             call SetPointRegObj(reg, i,j,k, val)
          end do
       case( "(BOX" ) 
          do 
             read(unit,*,iostat = ios) M4_READCOORD({i0, i1, di},{j0, j1, dj},{k0, k1, dk})
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of box list"})
                exit
             end if
             call SetBoxRegObj(reg, i0, i1, di, j0, j1, dj, k0, k1, dk, dval)
          end do
       case( "(VBOX" ) 
          do 
             read(unit,*,iostat = ios) M4_READCOORD({i0, i1, di},{j0, j1, dj},{k0, k1, dk}), &
             		(val(v),v=1,numval)
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of box list"})
                exit
             end if
             call SetBoxRegObj(reg, i0, i1, di, j0, j1, dj, k0, k1, dk, val)
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
		 	unit = unit + 1
            M4_WRITE_DBG({"success, unit = ",TRIM(i2str(unit))})
		  else
          	 M4_WRITE_WARN({"COULD NOT OPEN ",TRIM(loadfn)})
			read(unit,*,iostat = ios) string
          	 M4_WRITE_DBG({"consumed ",TRIM(string)}) ! consume )LOAD
          end if
       case("INIT")    
		  numvalues = numnodes + 1 	
       case("(FILL") ! a list of values
          do 
             read(unit,*,iostat = ios) (val(v),v=1,numval)
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of value list"})
                exit
             end if
             call FillValueRegObj(reg, val)
          end do       
       case("(SET") ! a list of values
          do 
             read(unit,*,iostat = ios) (val(v),v=1,numval)
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of value list"})
                exit
             end if
             call SetValueRegObj(reg, val)
          end do
       case("LIST") ! force list mode
          reg%islist = .true.
          auto = .false.
          reg%isbox = .false.
       case("MASK") ! force mask mode
          reg%islist = .false.
          reg%isbox = .false.
          auto = .false.
       case("AUTO") 
       case default
          if ( string(1:2) .eq. "(!" ) then
            skiptill = cat2(")",string(3:))
            M4_WRITE_DBG({"got token (! -> skiptill = ", TRIM(skiptill)})  
            cycle
          end if
          if ( string(1:1) .ne. ")" ) then
             M4_FATAL_ERROR({"BAD TERMINATOR: ", TRIM(string), " in ReadRegObj"})
          end if
          if ( unit .ne. UNITTMP ) then
             M4_FATAL_ERROR({"TERMINATED IN FILE INCLUDED WITH LOAD!"})
          end if
          exit
       end select

    end do

    call CreateRegObjEnd(reg,auto)

    deallocate(val,dval)
    
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
    regobj(numregobj)%islist = .false.
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
    regobj(numregobj)%islist = .false.
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

    else
   
       if ( .not. reg%islist ) then ! A MASK

          M4_WRITE_DBG({"setting up reg in mask mode", numnodes})
          allocate(reg%mask(reg%is:reg%ie,reg%js:reg%je,reg%ks:reg%ke),stat = err )
          M4_ALLOC_ERROR(err, {"CreateRegObjEnd"})
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
                      reg%val(:,p) = tmpval(:,p)
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
          allocate(reg%i(numnodes),reg%j(numnodes),reg%k(numnodes),reg%val(reg%numval,numnodes) , stat = err)
          M4_ALLOC_ERROR(err,{"CreateRegObjEnd"})
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
                      reg%val(:,num) = tmpval(:,p)
                   end if
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


!----------------------------------------------------------------------


  subroutine SetBoxRegObj(reg, i0,i1,di, j0,j1,dj, k0,k1,dk, val)

    type(T_REG) :: reg
    integer :: i0, i1, di, j0, j1, dj, k0, k1, dk
    real(kind=8) :: val(:)
    integer :: i,j,k,p

    if ( reg%idx .eq. -1 ) return

    if ( di .eq. 0 .or. dj .eq. 0 .or. dk .eq. 0 ) then 
       M4_WRITE_DBG({"box has zero stride -> return"})
       return
    end if
    
    M4_WRITE_DBG({"set box: ",i0, i1, di, j0, j1, dj, k0, k1, dk})

    reg%isbox = numnodes .eq. 0
    
    ! check whether box is inside grid
    if ( i1 .lt. domreg%is .or. j1 .lt. domreg%js .or. k1 .lt. domreg%ks .or. &
         i0 .gt. domreg%ie .or. j0 .gt. domreg%je .or. k0 .gt. domreg%ke ) then
       M4_WRITE_DBG({"box not in given domain --->"})
       M4_IFELSE_DBG({call EchoRegObj(reg)})
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
             p = tmpmask(i,j,k)
             if ( p .eq. 0 ) then
                numnodes = numnodes + 1
                tmpmask(i,j,k) = numnodes 
                tmpval(:,numnodes) = val(:)
             else
             	tmpval(:,p) = val(:)
             end if
          end do
       end do
    end do

    M4_WRITE_DBG({"box set!"})

  end subroutine SetBoxRegObj

!----------------------------------------------------------------------

  subroutine SetPointRegObj(reg, i,j,k, val)

    type(T_REG) :: reg
    integer :: i,j,k,p
    real(kind=8) :: val(:)

    if ( reg%idx .eq. -1 ) return

    M4_WRITE_DBG({"set point: ",i,j,k})

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
    p = tmpmask(i,j,k)
    if ( p .eq. 0 ) then
       numnodes = numnodes + 1
       tmpmask(i,j,k) = numnodes
       tmpval(:,numnodes) = val(:)
	else
	   tmpval(:,p) = val(:)
    end if
    M4_WRITE_DBG({"point set!"})

  end subroutine SetPointRegObj

!----------------------------------------------------------------------

  subroutine SetValueRegObj(reg, val)
  
    type(T_REG) :: reg
	real(kind=8) :: val(:)
	integer :: v

    if ( reg%idx .eq. -1 ) return

    M4_WRITE_DBG({"set value at ",numvalues," of ",numnodes," nodes"})
	
	if ( numvalues .gt. numnodes ) return
	
    do v = 1, reg%numval
       tmpval(:,numvalues) = val(:)
    end do

    numvalues = numvalues + 1

    	
  end subroutine SetValueRegObj

!----------------------------------------------------------------------
  
  subroutine FillValueRegObj(reg, val)
  
    type(T_REG) :: reg
	real(kind=8) :: val(:)
	integer :: v,p
	
    M4_WRITE_DBG({"fill values from ",numvalues," to ",numnodes," nodes"})

    do p = numvalues, numnodes	
      do v = 1, reg%numval
         tmpval(v,p) = val(v)
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
       if ( .not. reg%islist ) then
          M4_WRITE_DBG("deallocating reg%mask")
          deallocate(reg%mask)
       else
          M4_WRITE_DBG("deallocating reg%i reg%j reg%k")
          deallocate(reg%i, reg%j, reg%k)
       end if
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
    
    if ( reg%numval .eq. 0 ) return

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

  subroutine DisplayRegObj(reg)

    type(T_REG) :: reg
	character(len=20) :: tstr
    
    if ( reg%idx .eq. -1 ) return

    if ( reg%isbox ) then
    	tstr = "filled"
	else
		tstr = "points"
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
    M4_WRITE_INFO({"islist = ", reg%islist })
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
       if ( .not. reg%islist ) then
          masksz =  (reg%ie - reg%is + 1) * (reg%je - reg%js + 1) * (reg%ke - reg%ks + 1)
          memusage = memusage +  masksz * 4
       else
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




