!----------------------------------------------------------------------
!
!  module: vallist / meta
!
!  a table of values as needed for frequencies etc.
!
!----------------------------------------------------------------------

!======================================================================
!
!  Value definitions:
!
!  (VAL
!    (SET
!       1.0
!       2.0 
!    )SET   
!    (DOSET
!       0.5           ! offset value: voff
!       0.25          ! stepping: dv
!       1 100 10      ! loop: i0, i1, di
!    )DOSET
!  )VAL
!
!  While (SET-) directly sets a number of values, (DOSET-) can be
!  used to generate a whole set of evenspaced values equivalent to
!
!  do i = i0, i1, di
!    n = n + 1
!    val(n) = voff + i * dv
!  end do
!


module vallist

  use constant
  use strings

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), parameter :: modname = 'VALLIST'

  ! --- Public Methods

  public :: ReadValObj
  public :: InitializeValList
  public :: FinalizeValList
  public :: CreateValObjStart
  public :: CreateValObjEnd
  public :: DestroyValObj
  public :: EchoValObj
  public :: DisplayValObj

  ! --- Public Data

  public :: valobj
  public :: numvalobj
  public :: T_VAL

  ! --- Types

  type T_VAL

     integer :: idx = 0                          ! this objects index
     real(kind=8),pointer,dimension(:,:) :: data ! values field
     integer :: numval
     integer :: numdata

  end type T_VAL

  ! --- Data

  type(T_VAL) :: valobj(MAXVALOBJ) 
  integer :: numvalobj
  real(kind=8), allocatable, dimension(:,:) :: tmpdata
  integer :: numdata

contains

!----------------------------------------------------------------------

  subroutine InitializeValList

    numvalobj = 0

  end subroutine InitializeValList
  
!----------------------------------------------------------------------

  subroutine FinalizeValList

    integer :: i
    
    do i = 1, numvalobj 
       call DestroyValObj(valobj(i))
    end do

    numvalobj = 0

  end subroutine FinalizeValList

!----------------------------------------------------------------------

  subroutine ReadValObj(val, funit, numval)

    integer :: funit, numval
    type(T_VAL) :: val

    integer :: ios,err, unit 
    integer :: i, i0, i1, di
    real(kind=8),dimension(:),allocatable :: v,dv
    character(len=STRLNG) :: string, line, loadfn, skiptill

    M4_WRITE_DBG(". enter ReadValObj")

    unit = funit
    val = CreateValObjStart(numval)
   
    allocate(v(1:numval),dv(1:numval),stat=err)
    M4_ALLOC_ERROR(err,{"ReadValObj"})
    
    ! read until an unexpected line is encountered, eg ")"
    skiptill = ""
    do
  	  
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
       case( "(SET" ) 
          do 
             read(unit,*,iostat = ios) (v(i),i=1,numval)
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of value list"})
                exit
             end if
             numdata = numdata + 1
             tmpdata(:,numdata) = v
          end do
       case( "(DOSET" ) 
             read(unit,*,iostat = ios) (v(i),i=1,numval) 
             read(unit,*,iostat = ios) (dv(i),i=1,numval) 
             read(unit,*,iostat = ios) i0,i1,di
             if ( ios .ne. 0 ) then
                M4_WRITE_DBG({"end of value list"})
                exit
             end if
             do i = i0,i1,di
                numdata = numdata + 1
                tmpdata(:,numdata) = v
                v = v + dv * di
             end do
              read(unit,*,iostat = ios) line ! consume )DOSET
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
       end select

    end do

    call CreateValObjEnd(val)

    deallocate(v,dv)
    
    M4_WRITE_DBG(". exit ReadValObj")

  end subroutine ReadValObj

!----------------------------------------------------------------------

  type(T_VAL) function CreateValObjStart(numval)

    integer :: err, numval

    M4_WRITE_DBG(". enter CreateValObjStart")

    M4_WRITE_DBG({ "allocating value field, size = ", numval})
    allocate(tmpdata(1:numval,MAXVALUES),stat = err)
    M4_ALLOC_ERROR(err,"CreateValObjStart")

    numdata = 0
    tmpdata = 0.0

    numvalobj = numvalobj + 1
    valobj(numvalobj)%idx = numvalobj

    ! initialize valobj
    valobj(numvalobj)%numval = numval

    CreateValObjStart = valobj(numvalobj)

    M4_WRITE_DBG(". exit CreateValObjStart")
   
  end function CreateValObjStart

!----------------------------------------------------------------------

  subroutine CreateValObjEnd(val)

    type(T_VAL) :: val
    integer :: i, err

    M4_WRITE_DBG(". enter CreateValObjEnd")

    if ( numdata .eq. 0 ) then
       deallocate(tmpdata)
       return
    endif

    val%numdata = numdata

    allocate(val%data(val%numval,numdata) , stat = err)
    M4_ALLOC_ERROR(err, {"CreateValObjEnd"})

! copy data accross
    do i = 1, numdata
       val%data(:,i) = tmpdata(:,i)
    end do

    deallocate(tmpdata)

    M4_WRITE_DBG({" created valobj num = ", numvalobj})
    M4_IFELSE_DBG({call EchoValObj(val)})

    valobj(numvalobj) = val

    M4_WRITE_DBG(". exit CreateValObjEnd")

  end subroutine CreateValObjEnd

!----------------------------------------------------------------------

  subroutine DestroyValObj(val)

    type(T_VAL) :: val

    if ( val%idx .eq. -1 ) return
    
    M4_WRITE_DBG(". enter DestroyValObjEnd")

    M4_IFELSE_DBG({call EchoValObj(val)})
    deallocate(val%data)
    val%idx = -1
    
    M4_WRITE_DBG(". exit DestroyValObjEnd")
    
  end subroutine DestroyValObj

!----------------------------------------------------------------------

  subroutine DisplayValObj(val)

    type(T_VAL) :: val
    
    if ( val%idx .eq. -1 ) return

    M4_WRITE_INFO({"@",TRIM(i2str(val%idx)),&
    	" [",TRIM(i2str(val%numdata)),"] #",TRIM(i2str(val%numval)) })

  end subroutine DisplayValObj
    
!----------------------------------------------------------------------

  subroutine EchoValObj(val)

    type(T_VAL) :: val

    if ( val%idx .eq. -1 ) return

    M4_WRITE_INFO({"--- val # ", TRIM(i2str(val%idx)) })
    M4_WRITE_INFO({"numdata   = ", val%numdata })
    M4_WRITE_INFO({"numval   = ", val%numval })

  end subroutine EchoValObj

!----------------------------------------------------------------------

  integer function MemoryValList()
    
    integer :: i, masksz
    integer :: memusage = 0
    type(T_VAL) :: val
    
    do i = 1, numvalobj 
       
       val = valobj(i)
       if ( val%idx .eq. -1 .or. val%numdata .eq. 0 ) cycle
       
       memusage = memusage + val%numdata * val%numval * 8
       
    end do
    
    MemoryValList = memusage
    
  end function MemoryValList
  
    
!----------------------------------------------------------------------

end module vallist

!
! Authors:  J.Hamm 
! Modified: 27/12/2007
! 
!======================================================================




