!-*- F90 -*------------------------------------------------------------
!
!  module: fdtd_outhd5 / meta
!
!  this module handles HD5 output of data related to the fdtd module.
!
!----------------------------------------------------------------------


module fdtd_outhd5

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use fdtd_calc

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'FDTD_OUTHD5'

 ! --- Public Methods

  public :: InitializeFdtdOuthd5Obj
  public :: FinalizeFdtdOuthd5Obj
  public :: WriteDataFdtdOuthd5Obj

contains

!----------------------------------------------------------------------

  subroutine InitializeFdtdOuthd5Obj(out)

    type (T_OUT) :: out
    type (T_BUF) :: buf

    if ( out%fn .eq. 'Px' .or. out%fn .eq. 'Py' .or. out%fn .eq. 'Pz' .or. &
         out%fn .eq. 'En' ) then

       ! allocate buffer for calculated data

       buf = CreateBufObj(regobj(out%regidx),M4_ISCF,1)
       out%bufidx = buf%idx

       M4_WRITE_DBG({"created data buffer # ",&
            TRIM(i2str(out%bufidx))," for out # ",TRIM(i2str(out%idx))})

    endif

  end subroutine InitializeFdtdOuthd5Obj


!----------------------------------------------------------------------

  subroutine FinalizeFdtdOuthd5Obj(out)

    type (T_OUT) :: out


  end subroutine FinalizeFdtdOuthd5Obj


!----------------------------------------------------------------------


  subroutine WriteDataFdtdOuthd5Obj(out, mode)

    type (T_OUT) :: out
    type (T_BUF) :: buf
    logical :: mode

    M4_WRITE_DBG({"write data ",TRIM(out%filename), " ",TRIM(out%fn)})

    buf = bufobj(out%bufidx)

    select case (out%fn)
    case('Px')
       call FdtdCalcPx(buf, 1, mode)
       call WriteBufData(out, mode)
    case('Py')
       call FdtdCalcPy(buf, 1, mode)
       call WriteBufData(out, mode) 
    case('Pz')
       call FdtdCalcPz(buf, 1, mode)
       call WriteBufData(out, mode) 
    case('En')
       call FdtdCalcEn(buf, 1, mode)
       call WriteBufData(out, mode) 
    case('EH')
       call WriteEH(out, mode)
    case('Ex')
       call WriteField(out, Ex, mode)
    case('Ey')
       call WriteField(out, Ey, mode)
    case('Ez')
       call WriteField(out, Ez, mode)
    case('Hx')
       call WriteField(out, Hx, mode)
    case('Hy')
       call WriteField(out, Hy, mode)
    case('Hz')
       call WriteField(out, Hz, mode)
    case('Di')         
       call WriteEpsilon(out, mode)
    case default
       write(out%funit,*) "OUTPUT FUNCTION NOT IMPLEMENTED" 
    end select
    
  contains

    ! **************************************************************** !

    subroutine WriteEH(out,mode)

      type (T_OUT) :: out
      logical :: mode

      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  

      if ( .not. mode ) return

      reg = regobj(out%regidx)
      
      M4_REGLOOP_WRITE(reg,p,i,j,k,w,   
      out%funit,6E15.6E3, { real(Ex(i,j,k)),real(Ey(i,j,k)),real(Ez(i,j,k)), &
           real(Hx(i,j,k)),real(Hy(i,j,k)),real(Hz(i,j,k)) }
      )
      
    end subroutine WriteEH


    ! **************************************************************** !

    subroutine WriteField(out,field,mode)

      type (T_OUT) :: out
      logical :: mode

      M4_FTYPE, dimension(M4_RANGE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)) :: field
      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  

      if ( .not. mode ) return

      M4_WRITE_DBG({"WriteField!"})
      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})
      reg = regobj(out%regidx)

       M4_REGLOOP_WRITE(reg,p,i,j,k,w,   
       out%funit,M4_IFELSE_CF({2})E15.6E3, { M4_IFELSE_CF({real(field(i,j,k)),aimag(field(i,j,k))},{field(i,j,k)}) }
       )

    end subroutine WriteField


    ! **************************************************************** !

    subroutine WriteEpsilon(out,mode)

      type (T_OUT) :: out
      logical :: mode

      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  

      if ( .not. mode ) return

      reg = regobj(out%regidx)
      
      M4_REGLOOP_WRITE(reg,p,i,j,k,w,
      out%funit,E15.6E3, { 1.0/epsinv(i,j,k) }
      )

    end subroutine WriteEpsilon


    ! **************************************************************** !

    subroutine WriteBufData(out, mode)

      type (T_OUT) :: out
      logical :: mode
      M4_FTYPE :: sum = 0.0
      type (T_BUF) :: buf
      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  

      if ( .not. mode ) return

      reg = regobj(out%regidx)
      buf = bufobj(out%bufidx)
  
      M4_WRITE_DBG({"write data buffer # ",&
           TRIM(i2str(out%bufidx)),&
           " for out # ",TRIM(i2str(out%idx))})

      if ( out%mode .eq. 'S' ) then ! spatial integration
         M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
         sum = sum + buf%data(p,1)
         })
         write(out%funit,"(M4_IFELSE_CF({2})E15.6E3)") M4_IFELSE_CF({real(sum),aimag(sum)},sum)
      else
         M4_REGLOOP_WRITE(reg,p,i,j,k,w,
         out%funit, M4_IFELSE_CF({2})E15.6E3, { M4_IFELSE_CF({real(buf%data(p,1)),aimag(buf%data(p,1))},buf%data(p,1)) }
         )
      endif
      ! clear buffer data after write
      !	  buf%data(:,1) = 0.0	
      
    end subroutine WriteBufData

    ! **************************************************************** !

  end subroutine WriteDataFdtdOuthd5Obj


end module fdtd_outhd5

!
! Authors:  J.Hamm

! Modified: 7/1/2008
!
!======================================================================
