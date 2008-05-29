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

    if ( out%fn .eq. 'Sx' .or. out%fn .eq. 'Sy' .or. out%fn .eq. 'Sz' .or. &
         out%fn .eq. 'En' ) then

       ! allocate buffer for calculated data

       buf = CreateBufObj(regobj(out%regidx),M4_ISCF,1)
       out%bufidx = buf%idx

       M4_WRITE_DBG({"created data buffer # ",&
            TRIM(i2str(out%bufidx))," for out # ",TRIM(i2str(out%idx))})

       return

    endif

    if ( out%fn .eq. 'S' ) then
       
       buf = CreateBufObj(regobj(out%regidx),M4_ISCF,3)
       out%bufidx = buf%idx

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
    case('S')
       call FdtdCalcSx(buf, 1, mode)
       call FdtdCalcSy(buf, 2, mode)
       call FdtdCalcSz(buf, 3, mode)
       call WriteBuffer(out, buf, 0, mode)
    case('Sx')
       call FdtdCalcSx(buf, 1, mode)
       call WriteBuffer(out, buf, 0, mode)
    case('Sy')
       call FdtdCalcSy(buf, 1, mode)
       call WriteBuffer(out, buf, 0, mode)
    case('Sz')
       call FdtdCalcSz(buf, 1, mode)
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
    case('Eps')         
       call WriteEps(out, mode)
    case('Mu')         
       call WriteMu(out, mode)
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
         write(out%funit,"(A)") "SCALARS scalar float 1"
         write(out%funit,"(A)") "LOOKUP_TABLE default"
      endif

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      call PaScalar(pa, fc(i,j,k), val)

      if ( out%mode .ne. 'S' ) then
 
         write(out%funit,"(E14.6E2)") real(val,4)
         
      else
         sum = sum + val
      endif

      },{}, {} )
   

      if ( out%mode .eq. 'S' ) then
         write(out%funit,"(A)") "FIELD FieldSum 1"
         write(out%funit,"(E14.6E2)") real(sum,4)
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
         write(out%funit,"(A)") "VECTORS vectors float"
      end if

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      call PaVector(pa,fx(i,j,k),fy(i,j,k),fz(i,j,k),vx,vy,vz)

      if ( out%mode .ne. 'S' ) then
         write(out%funit,"(3E14.6E2)") real(vx,4),real(vy,4),real(vz,4)
      else
         sx = sx + vx
         sy = sy + vy
         sz = sz + vz
      endif

      },{}, {} )

      if ( out%mode .eq. 'S' ) then
         write(out%funit,*) "FIELD FieldSum 3"
         write(out%funit,"(3E14.6E2)") real(sx,4), real(sy,4), real(sz,4)
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

      if ( out%mode .ne. 'S' ) then
         if ( buf%numslot .eq. 3 ) then
            write(out%funit,"(A)") "VECTORS vectors float"
         else
            write(out%funit,"(A,A)") "SCALARS scalars float ",TRIM(i2str(buf%numslot))
            write(out%funit,*) "LOOKUP_TABLE default"
         endif
      endif

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      call PaBuffer(pa, buf, p, val)

      if ( out%mode .ne. 'S' ) then
         write(out%funit,*) (real(val(l),4), l=1, buf%numslot,1)
      else
         sum = sum + val
      endif

      },{}, {} )
   

      if ( out%mode .eq. 'S' ) then
         write(out%funit,*) "FIELD FieldSum ", buf%numslot
         write(out%funit,*) (real(sum(l),4), l=1, buf%numslot,1)
      endif

      deallocate(val,sum)


    end subroutine WriteBuffer

 
    ! **************************************************************** !

    subroutine WriteEps(out,mode)

      type (T_OUT) :: out
      logical :: mode
      real(kind=8) :: val

      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  

      if ( .not. mode ) return

      write(out%funit,"(A)") "SCALARS scalar float 3"
      write(out%funit,"(A)") "LOOKUP_TABLE default"

      reg = regobj(out%regidx)
      
      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      write(out%funit,"(E14.6E2)") &
           real(1./epsinvx(i,j,k),4),real(1./epsinvy(i,j,k),4),real(1./epsinvz(i,j,k),4)
      
      },{}, {} )

    end subroutine WriteEps

    ! **************************************************************** !

    subroutine WriteMu(out,mode)

      type (T_OUT) :: out
      logical :: mode
      real(kind=8) :: val

      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  

      if ( .not. mode ) return

      write(out%funit,"(A)") "SCALARS scalar float 3"
      write(out%funit,"(A)") "LOOKUP_TABLE default"

      reg = regobj(out%regidx)
      
      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      write(out%funit,"(3E14.6E2)") &
           real(1./M4_MUINVX(i,j,k),4),real(1./M4_MUINVY(i,j,k),4),real(1./M4_MUINVZ(i,j,k),4)

      },{}, {} )

    end subroutine WriteMu


    ! *************************************************************** !

  end subroutine WriteDataFdtdOutvtkObj


end module fdtd_outvtk

!
! Authors:  J.Hamm

! Modified: 4/12/2007
!
!======================================================================
