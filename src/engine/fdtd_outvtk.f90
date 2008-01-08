!-*- F90 -*------------------------------------------------------------
!
!  module: fdtd_outvtk / meta
!
!  this module handles VTK output of data related to the fdtd module.
!
!----------------------------------------------------------------------


module fdtd_outvtk

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use fdtd_calc
  use out_calc

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'FDTD_OUTVTK'

 ! --- Public Methods

  public :: InitializeFdtdOutvtkObj
  public :: FinalizeFdtdOutvtkObj
  public :: WriteDataFdtdOutvtkObj

contains

!----------------------------------------------------------------------

  subroutine InitializeFdtdOutvtkObj(out)

    type (T_OUT) :: out
    type (T_BUF) :: buf

    if ( out%fn .eq. 'Px' .or. out%fn .eq. 'Py' .or. out%fn .eq. 'Pz' .or. &
         out%fn .eq. 'En' ) then

       ! allocate buffer for calculated data

       buf = CreateBufObj(regobj(out%regidx),M4_ISCF,1)
       out%bufidx = buf%idx

       M4_WRITE_DBG({"created data buffer # ",&
            TRIM(i2str(out%bufidx))," for out # ",TRIM(i2str(out%idx))})

       return

    endif

    if ( out%fn .eq. 'P' ) then
       
       buf = CreateBufObj(regobj(out%regidx),M4_ISCF,3)

       M4_WRITE_DBG({"created data buffer # ",&
            TRIM(i2str(out%bufidx))," for out # ",TRIM(i2str(out%idx))})

       return

    endif

  end subroutine InitializeFdtdOutvtkObj


!----------------------------------------------------------------------

  subroutine FinalizeFdtdOutvtkObj(out)

    type (T_OUT) :: out


  end subroutine FinalizeFdtdOutvtkObj


!----------------------------------------------------------------------


  subroutine WriteDataFdtdOutvtkObj(out, mode)

    type (T_OUT) :: out
    type (T_BUF) :: buf
    logical :: mode

    M4_WRITE_DBG({"write data ",TRIM(out%filename), " ",TRIM(out%fn)})

    if ( out%bufidx .ge. 1 ) then 
       M4_WRITE_DBG({"from buffer #",TRIM(i2str(out%bufidx))})
       buf = bufobj(out%bufidx)
    endif

    select case (out%fn)
    case('En')
       call FdtdCalcEn(buf, 1, mode)
       call WriteBuffer(out, buf, 0, mode)
    case('P')
       call FdtdCalcPx(buf, 1, mode)
       call FdtdCalcPy(buf, 2, mode)
       call FdtdCalcPz(buf, 3, mode)
       call WriteBuffer(out, buf, 0, mode)
    case('Px')
       call FdtdCalcPx(buf, 1, mode)
       call WriteBuffer(out, buf, 0, mode)
    case('Py')
       call FdtdCalcPy(buf, 1, mode)
       call WriteBuffer(out, buf, 0, mode)
    case('Pz')
       call FdtdCalcPz(buf, 1, mode)
       call WriteBuffer(out, buf, 0, mode)
    case('E')
       call WriteVector(out, Ex,Ey,Ez, 0, mode)
    case('Eamp')
       call WriteVector(out, Ex,Ey,Ez, 1, mode)
    case('Ephase')
       call WriteVector(out, Ex,Ey,Ez, 2, mode)
    case('H')
       call WriteVector(out, Hx,Hy,Hz, 0, mode)
    case('Hamp')
       call WriteVector(out, Hx,Hy,Hz, 1, mode)
    case('Hphase')
       call WriteVector(out, Hx,Hy,Hz, 2, mode)
    case('Ex')
       call WriteScalar(out, Ex, 0, mode)
    case('Ey')
       call WriteScalar(out, Ey, 0, mode)
    case('Ez')
       call WriteScalar(out, Ez, 0, mode)
    case('Hx')
       call WriteScalar(out, Hx, 0, mode)
    case('Hy')
       call WriteScalar(out, Hy, 0, mode)
    case('Hz')
       call WriteScalar(out, Hz, 0, mode)
    case('Di')         
       call WriteEpsilon(out, mode)
    case default
       write(out%funit,*) "OUTPUT FUNCTION NOT IMPLEMENTED" 
    end select
    
  contains

    ! **************************************************************** !

    subroutine WriteScalar(out, fc, pa, mode)

      type (T_OUT) :: out
      logical :: mode
      integer :: pa ! 0 : real part, 1: amplitude, 2: phase

      M4_FTYPE, dimension(M4_RANGE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)) :: fc
      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  
      real(kind=8) :: val, sum

      if ( .not. mode ) return

      M4_WRITE_DBG({"WriteScalar!"})
      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      sum = 0.0

      reg = regobj(out%regidx)

      if ( out%mode .ne. 'S' ) then
         write(out%funit,*) "SCALARS scalar float 1"
         write(out%funit,*) "LOOKUP_TABLE default"
      endif

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      call PaScalar(pa, fc(i,j,k), val)

      if ( out%mode .ne. 'S' ) then
 
         write(out%funit,"(E15.6E3)") val
         
      else
         sum = sum + val
      endif

      },{if ( reg%is .ne. reg%ie ) write(out%funit,*)}, {if ( reg%js .ne. reg%je ) write(out%funit,*)} )
   

      if ( out%mode .eq. 'S' ) then
         write(out%funit,*) "FIELD FieldSum 1"
         write(out%funit,"(E15.6E3)") sum
      endif

    end subroutine WriteScalar

    ! **************************************************************** !

    subroutine WriteVector(out, fx,fy,fz, pa, mode)

      type (T_OUT) :: out
      logical :: mode
      integer :: pa
      M4_FTYPE, dimension(M4_RANGE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)) :: fx,fy,fz
      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  
      real(kind=8) :: vx,vy,vz, sx,sy,sz

      if ( .not. mode ) return

      M4_WRITE_DBG({"WriteVector!"})
      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      sx = 0.0
      sy = 0.0
      sz = 0.0

      reg = regobj(out%regidx)
      
      if ( out%mode .ne. 'S' ) then
         write(out%funit,*) "VECTORS vectors float"
         write(out%funit,*) "LOOKUP_TABLE default"
      end if

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      call PaVector(pa,fx(i,j,k),fy(i,j,k),fz(i,j,k),vx,vy,vz)

      if ( out%mode .ne. 'S' ) then
         write(out%funit,"(3E15.6E3)") vx,vy,vz
      else
         sx = sx + vx
         sy = sy + vy
         sz = sz + vz
      endif

      },{}, {} )

      if ( out%mode .eq. 'S' ) then
         write(out%funit,*) "FIELD FieldSum 3"
         write(out%funit,"(3E15.6E3)") sx, sy, sz
      endif
   
    end subroutine WriteVector


    ! **************************************************************** !

    subroutine WriteBuffer(out, buf, pa, mode)

      type (T_OUT) :: out
      type (T_BUF) :: buf
      logical :: mode
      integer :: pa ! 0 : real part, 1: amplitude, 2: phase
      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  
      real(kind=8), allocatable :: val(:), sum(:)
      integer :: l

      if ( .not. mode ) return

      allocate(val(buf%numslot),sum(buf%numslot))

      M4_WRITE_DBG({"WriteBuffer!"})
      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      sum = 0.0

      reg = regobj(out%regidx)

      write(out%funit,*) "LOOKUP_TABLE default"

      if ( out%mode .ne. 'S' ) then
         if ( buf%numslot .eq. 3 ) then
            write(out%funit,*) "VECTORS vectors float"
         else
            write(out%funit,*) "SCALARS scalars float ",TRIM(i2str(buf%numslot))
         endif
      endif

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      call PaBuffer(pa, buf, p, val)

      if ( out%mode .ne. 'S' ) then
         write(out%funit,*) (val(l), l=1, buf%numslot,1)
      else
         sum = sum + val
      endif

      },{}, {} )
   

      if ( out%mode .eq. 'S' ) then
         write(out%funit,*) "FIELD FieldSum ", buf%numslot
         write(out%funit,*) (sum(l), l=1, buf%numslot,1)
      endif

      deallocate(val,sum)


    end subroutine WriteBuffer

 
    ! **************************************************************** !

    subroutine WriteEpsilon(out,mode)

      type (T_OUT) :: out
      logical :: mode
      real(kind=8) :: val

      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  

      if ( .not. mode ) return

      reg = regobj(out%regidx)
      
      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      val = 1./epsinv(i,j,k)
      if ( reg%isbox ) then
         write(out%funit,"(I5,(E15.6E3))") p,val
      endif
      
      },{}, {} )

    end subroutine WriteEpsilon


    ! *************************************************************** !

  end subroutine WriteDataFdtdOutvtkObj


end module fdtd_outvtk

!
! Authors:  J.Hamm

! Modified: 4/12/2007
!
!======================================================================
