!----------------------------------------------------------------------
!
!  module: region / max3d
!
!  a spatial region, box, mask or pointlist 
!
!  subs:
!
!
!----------------------------------------------------------------------

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
!    MASKED    
!    )
!  )
!
!
!  A region is a set of points with a weight attached. BOX, POINT set the 
!  weigth to 1.0, while VPOINT and VBOX allow to set a userdefined weight 
!  value. If only one BOX is defined in a region, then no point lists or
!  masks will be allocated. Instead the "isbox" value of the region is set
!  to true. Normally each region would have a list of indirect addressed 
!  point coordinates, but if MASKED is specified, then a field of values
!  is allocated within a bounding box (is,js,ks)(ie,je,ke).



module region

  use constant
  use strings
  use grid

  implicit none
  save

  ! --- Constants

  integer, parameter :: MAXREGIONS = 1000

  ! --- Types

  type T_REGION

     ! logical
     logical :: isbox, ismask

     ! index
     integer :: idx

     ! bounding box
     integer :: is, ie, di
     integer :: js, je, dj
     integer :: ks, ke, dk

     ! is the box all filled or do we have a point list?

     real*8, pointer, dimension(:,:,:) :: mask
     integer, pointer, dimension(:) :: i, j, k
     real*8, pointer, dimension(:) :: val

  end type T_REGION

  ! --- Variables  

  type(T_REGION) :: regions(MAXREGIONS) 
  integer :: numregions

  ! --- Fields

  integer, allocatable, dimension(:,:) :: tmpregpoints
  real*8, allocatable, dimension(:) :: tmpregvalues
  integer :: numtmpregpoints

contains

  subroutine InitializeRegion 

    integer :: err

    numregions = 0

  end subroutine InitializeRegion
  
  
  subroutine FinalizeRegion

    integer :: i
    
    do i = 1, numregions 
       call RegionDestroy(regions(i))
    end do

  end subroutine FinalizeRegion


  subroutine ReadObjRegion(funit)

    integer :: funit
    type(T_REGION) :: reg
    type(T_REGION), external :: CreateObjRegion

    integer :: ios, err
    integer :: i,j,k, i0, i1, di, j0, j1, dj, k0, k1, dk, num
    real*8 :: val
    character(len=STRLNG) :: string

    reg = CreateObjRegion()

    numtmpregpoints = 0
    err = 0
    if ( err .eq. 0 ) then
       allocate(tmpregpoints(GRIDSIZE, 3),stat = err)
    end if
    if ( err .eq. 0 ) then
       allocate(tmpregvalues(GRIDSIZE),stat = err)
    end if
    if ( err .ne. 0 ) then
       write(6,*) "ALLOCATION ERROR: ReadObjRegion"
       stop
    end if
   
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
             call SetObjRegionPoint(reg, i,j,k, val)
          end do
       case( "(BOX" ) 
          if ( reg%isbox ) then
             reg%isbox = .false.
          else 
             reg%isbox = .true.
          endif
          val = 1.0
          do 
             read(funit,*,iostat = ios)  i0, i1, di, j0, j1, dj, k0, k1, dk
             if ( ios .ne. 0 ) then
                exit
             end if
             call SetObjRegionBox(reg, i0, i1, di, j0, j1, dj, k0, k1, dk, val)
          end do
       case( "(VBOX" ) 
          reg%isbox = .false.
          do 
             read(funit,*,iostat = ios)  i0, i1, di, j0, j1, dj, k0, k1, dk, val
             if ( ios .ne. 0 ) then
                exit
             end if
             call SetObjRegionBox(reg, i0, i1, di, j0, j1, dj, k0, k1, dk, val)
          end do
       case( "(VPOINTS" ) 
          reg%isbox = .false.
          do 
             read(funit,*,iostat = ios) i,j,k, val
             if ( ios .ne. 0 ) then
                exit
             end if
             call SetObjRegionPoint(reg, i,j,k, val)
          end do
       case("MASKED")
          reg%ismask = .true.
       case default
          exit
       end select

    end do

    if ( .not. reg%isbox ) then
        
       if ( reg%ismask ) then

          allocate(reg%mask(reg%is:reg%ie,reg%js:reg%je,reg%ks:reg%ke), stat = err )
          if ( err .ne. 0 ) then
             write(6,*) "ALLOCATION ERROR: ReadObjRegion"
             stop
          end if
          num = 0
          do i = reg%is, reg%ie
             do j = reg%js, reg%je
                do k = reg%ks, reg%ke
                   num = num + 1
                   reg%mask(tmpregpoints(num,1),tmpregpoints(num,2),tmpregpoints(num,3)) = tmpregvalues(num)
                end do
             end do
          end do
       else
          
          allocate(reg%i(numtmpregpoints),reg%j(numtmpregpoints),reg%k(numtmpregpoints),reg%val(numtmpregpoints) , stat = err)
          if ( err .ne. 0 ) then
             write(6,*) "ALLOCATION ERROR: ReadObjRegion"
             stop
          end if
          num = 0
          do num = 1, numtmpregpoints
             reg%i(num) = tmpregpoints(num,1)
             reg%j(num) = tmpregpoints(num,2)
             reg%k(num) = tmpregpoints(num,3)
             reg%val(num) = tmpregvalues(num)
          end do
    
       end if
       
    end if
 
    deallocate(tmpregpoints)
    deallocate(tmpregvalues)
    
  end subroutine ReadObjRegion



  type(T_REGION) function CreateObjRegion

    numregions = numregions + 1
    regions(numregions)%idx = numregions
    
    regions(numregions)%di = 0
    regions(numregions)%dj = 0
    regions(numregions)%dk = 0

    regions(numregions)%isbox = .false.
    regions(numregions)%ismask = .false.

    CreateObjRegion = regions(numregions)
   
  end function CreateObjRegion



  subroutine SetObjRegionBox(reg, i0,i1,di, j0,j1,dj, k0,k1,dk, val)

    type(T_REGION) :: reg
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
    if ( reg%isbox ) then
       reg%di = di
       reg%dj = dj
       reg%dk = dk
    else
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

  end subroutine SetObjRegionBox

  subroutine SetObjRegionPoint(reg, i,j,k, val)

    type(T_REGION) :: reg
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

  end subroutine SetObjRegionPoint


  subroutine DestroyObjRegion(reg)

    type(T_REGION) :: reg

    if ( .not. reg%isbox ) then
       if ( reg%ismask ) then
          deallocate(reg%mask)
       else
          deallocate(reg%i, reg%j, reg%k, reg%val)
       end if
    end if
    
  end subroutine DestroyObjRegion



end module region






