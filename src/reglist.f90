!----------------------------------------------------------------------
!
!  module: reglist / max3d
!
!  a spatial box, mask or pointlist 
!
!  subs:
!
!    InitializeRegList
!    FinalizeRegList
!
!    reg = CreateRegObj
!    DestroyRegObj
!    ReadRegObj
!    SetPointRegObj
!    SetBoxRegObj
!
!----------------------------------------------------------------------

!======================================================================
!
!  (REGION 
!    (POINT
!     100 20 30 
!     100 21 33
!    )
!    (VPOINT
!     100 10 32 0.8
!    )
!    (BOX
!     0 10 20 1 6 20 100 1
!    (VBOX
!     0 10 20 1 6 20 100 0.25 1
!    LIST | MASK | AUTO
!    )
!  )
!
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
  save

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


  ! --- Fields

  type(T_REG) :: regobj(MAXREGOBJ) 
  integer :: numregobj

  integer, allocatable, dimension(:,:) :: tmpregpoints
  real*8, allocatable, dimension(:) :: tmpregvalues
  integer :: numtmpregpoints

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

  end subroutine FinalizeRegList

!----------------------------------------------------------------------

  subroutine ReadRegObj(reg, funit)

    integer :: funit
    type(T_REG) :: reg
    type(T_REG), external :: CreateRegObj

    integer :: ios, err
    integer :: i,j,k, i0, i1, di, j0, j1, dj, k0, k1, dk, num
    integer :: masksz, listsz
    real*8 :: val
    character(len=STRLNG) :: string
    logical :: auto

    reg = CreateRegObj()

    numtmpregpoints = 0
    err = 0
    if ( err .eq. 0 ) then
       allocate(tmpregpoints(GRIDSIZE, 3),stat = err)
    end if
    if ( err .eq. 0 ) then
       allocate(tmpregvalues(GRIDSIZE),stat = err)
    end if
    if ( err .ne. 0 ) then
       write(STDERR,*) "!ERROR OUT OF MEMORY: ReadRegObj/reg"
       stop
    end if
   
    auto = .true.
    ! read until an unexpected line is encountered, eg ")"
    do

       read(funit,*) string

       select case ( string ) 
       case( "(POINTS" ) 
          reg%isbox = .false.
          val = 1.0
          do 
             read(funit,*,iostat = ios) i,j,k
             if ( ios .ne. 0 ) then
                exit
             end if
             call SetPointRegObj(reg, i,j,k, val)
          end do
       case( "(BOX" ) 
          if ( numtmpregpoints .eq. 0 ) then 
             reg%isbox = .true.
          else
             reg%isbox = .false.
          end if
          val = 1.0
          do 
             read(funit,*,iostat = ios)  i0, i1, di, j0, j1, dj, k0, k1, dk
             if ( ios .ne. 0 ) then
                exit
             end if
             call SetBoxRegObj(reg, i0, i1, di, j0, j1, dj, k0, k1, dk, val)
          end do
       case( "(VBOX" ) 
          reg%isbox = .false.
          do 
             read(funit,*,iostat = ios)  i0, i1, di, j0, j1, dj, k0, k1, dk, val
             if ( ios .ne. 0 ) then
                exit
             end if
             call SetBoxRegObj(reg, i0, i1, di, j0, j1, dj, k0, k1, dk, val)
          end do
       case( "(VPOINTS" ) 
          reg%isbox = .false.
          do 
             read(funit,*,iostat = ios) i,j,k, val
             if ( ios .ne. 0 ) then
                exit
             end if
             call SetPointRegObj(reg, i,j,k, val)
          end do
       case("LIST") ! force list mode
          reg%islist = .true.
          auto = .false.
       case("MASK") ! force mask mode
          reg%islist = .false.
          auto = .false.
       case default
          exit
       end select

    end do

    if ( reg%isbox ) then ! A SINGLE BOX
        
       reg%numnodes = &
            ( reg%ie - reg%is + 1 ) / reg%di * &
            ( reg%je - reg%js + 1 ) / reg%dj * &
            ( reg%ke - reg%ks + 1 ) / reg%dk    

    else

       if ( auto ) then
          masksz = (reg%ie - reg%is) * (reg%je - reg%js) * (reg%ke - reg%ks)
          listsz = numtmpregpoints * 3 
          if ( listsz .lt. masksz ) then 
             reg%islist = .true.
          else 
             reg%islist = .false.
          endif
       end if

       reg%numnodes = numtmpregpoints
       
       if ( .not. reg%islist ) then ! A MASK

          allocate(reg%mask(reg%is:reg%ie,reg%js:reg%je,reg%ks:reg%ke),reg%val(numtmpregpoints), stat = err )
          if ( err .ne. 0 ) then
             write(STDERR,*) "!ERROR OUT OF MEMORY: ReadRegObj/reg"
             stop
          end if
          reg%mask = 0
          do num = 1, numtmpregpoints
             reg%mask(tmpregpoints(num,1),tmpregpoints(num,2),tmpregpoints(num,3)) = num
             reg%val(num) = tmpregvalues(num)
          end do
          ! 1 point!
          reg%pe = 1

       else ! A LIST

          reg%pe = numtmpregpoints
          allocate(reg%i(numtmpregpoints),reg%j(numtmpregpoints),reg%k(numtmpregpoints),reg%val(numtmpregpoints) , stat = err)
          if ( err .ne. 0 ) then
             write(STDERR,*) "!ERROR OUT OF MEMORY: ReadRegObj/reg"
             stop
          end if
          do num = 1, numtmpregpoints
             reg%i(num) = tmpregpoints(num,1)
             reg%j(num) = tmpregpoints(num,2)
             reg%k(num) = tmpregpoints(num,3)
             reg%val(num) = tmpregvalues(num)
          end do
          ! 1 step stepping!
          reg%di = reg%ie - reg%is + 1
          reg%dj = reg%je - reg%js + 1
          reg%dk = reg%ke - reg%ks + 1
       end if
       
    end if
 
    deallocate(tmpregpoints)
    deallocate(tmpregvalues)
    
  end subroutine ReadRegObj

!----------------------------------------------------------------------

  type(T_REG) function CreateRegObj

    numregobj = numregobj + 1
    regobj(numregobj)%idx = numregobj
    
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
    integer :: i,j,k

    ! clip box to grid
    i0 = Max(i0,IBEG)
    j0 = Max(j0,JBEG)
    k0 = Max(k0,KBEG)
    i1 = Min(i1,IEND)
    j1 = Min(j1,JEND)
    k1 = Min(k1,KEND)
    
    ! resize bounding box
    reg%is = Min(i0,reg%is)
    reg%js = Min(j0,reg%js)
    reg%ks = Min(k0,reg%ks)
    reg%ie = Max(i1,reg%ie)
    reg%je = Max(j1,reg%je)
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

    ! set points
    do i = i0, i1, di
       do j = j0, j1, dj
          do k = k0, k1, dk
             numtmpregpoints = numtmpregpoints + 1
             tmpregpoints(numtmpregpoints,1) = i
             tmpregpoints(numtmpregpoints,2) = j
             tmpregpoints(numtmpregpoints,3) = k
             tmpregvalues(numtmpregpoints) = val
          end do
       end do
    end do

  end subroutine SetBoxRegObj

!----------------------------------------------------------------------

  subroutine SetPointRegObj(reg, i,j,k, val)

    type(T_REG) :: reg
    integer :: i,j,k
    real*8 :: val

    ! check whether point is inside grid
    if ( i .lt. IBEG .or. j .lt. JBEG .or. k .lt. KBEG .or. &
         i .gt. IEND .or. j .gt. JEND .or. k .gt. KEND ) then
       return
    end if

    ! resize bounding box
    reg%is = Min(reg%is, i)
    reg%js = Min(reg%js, j)
    reg%ks = Min(reg%ks, k)
    reg%ie = Max(reg%ie, i)
    reg%je = Max(reg%je, j)
    reg%ke = Max(reg%ke, k)
    reg%di = 1
    reg%dj = 1
    reg%dk = 1

    ! set point
    numtmpregpoints = numtmpregpoints + 1
    tmpregpoints(numtmpregpoints,1) = i
    tmpregpoints(numtmpregpoints,2) = j
    tmpregpoints(numtmpregpoints,3) = k
    tmpregvalues(numtmpregpoints) = val

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




