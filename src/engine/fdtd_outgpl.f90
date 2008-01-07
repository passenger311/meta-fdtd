!-*- F90 -*------------------------------------------------------------
!
!  module: fdtd_outgpl / meta
!
!  this module handles GPL output of data related to the fdtd module.
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------


module fdtd_outgpl

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'FDTD_OUTGPL'

 ! --- Public Methods

  public :: InitializeFdtdOutgplObj
  public :: FinalizeFdtdOutgplObj
  public :: WriteDataFdtdOutgplObj

contains

!----------------------------------------------------------------------

  subroutine InitializeFdtdOutgplObj(out)

    type (T_OUT) :: out
    type (T_BUF) :: buf

    if ( out%fn .eq. 'Px' .or. out%fn .eq. 'Py' .or. out%fn .eq. 'Pz' .or. &
         out%fn .eq. 'En' ) then

       ! allocate additional buffer for Load functions here

       buf = CreateBufObj(regobj(out%regidx),M4_ISCF,1)
       out%bufidx = buf%idx

       M4_WRITE_DBG({"created data buffer # ",&
            TRIM(i2str(out%bufidx))," for out # ",TRIM(i2str(out%idx))})

    endif

  end subroutine InitializeFdtdOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeFdtdOutgplObj(out)

    type (T_OUT) :: out


  end subroutine FinalizeFdtdOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataFdtdOutgplObj(out, mode)

    type (T_OUT) :: out
    type (T_BUF) :: buf
    logical :: mode

    M4_WRITE_DBG({"write data ",TRIM(out%filename), " ",TRIM(out%fn)})

    select case (out%fn)
    case('Px')
       call LoadBufPx(out)
       if ( mode ) call WriteBufData(out)
    case('Py')
       call LoadBufPy(out)
       if ( mode ) call WriteBufData(out) 
    case('Pz')
       call LoadBufPz(out)
       if ( mode ) call WriteBufData(out) 
    case('En')
       call LoadBufEn(out)
       if ( mode ) call WriteBufData(out) 
    case('EH')
       if ( mode ) call WriteEH(out)
    case('Ex')
       if ( mode ) call WriteComp(out,Ex)
    case('Ey')
       if ( mode ) call WriteComp(out,Ey)
    case('Ez')
       if ( mode ) call WriteComp(out,Ez)
    case('Hx')
       if ( mode ) call WriteComp(out,Hx)
    case('Hy')
       if ( mode ) call WriteComp(out,Hy)
    case('Hz')
       if ( mode ) call WriteComp(out,Hz)
    case('Di')         
       if ( mode ) call WriteEps(out)
    case default
       write(out%funit,*) "OUTPUT FUNCTION NOT IMPLEMENTED" 
    end select
    
  contains

    subroutine WriteEH(out)

      type (T_OUT) :: out

      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  
      reg = regobj(out%regidx)
      
      M4_REGLOOP_WRITE(reg,p,i,j,k,w,   
      out%funit,6E15.6E3, { real(Ex(i,j,k)),real(Ey(i,j,k)),real(Ez(i,j,k)), &
           real(Hx(i,j,k)),real(Hy(i,j,k)),real(Hz(i,j,k)) }
      )
      
    end subroutine WriteEH


    ! **************************************************************** !

    subroutine WriteComp(out,Comp)

      type (T_OUT) :: out
      M4_FTYPE, dimension(M4_RANGE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)) :: Comp
      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  

      M4_WRITE_DBG({"WriteComp!"})
      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})
      reg = regobj(out%regidx)

       M4_REGLOOP_WRITE(reg,p,i,j,k,w,   
       out%funit,M4_IFELSE_CF({2})E15.6E3, { M4_IFELSE_CF({real(Comp(i,j,k)),aimag(Comp(i,j,k))},{Comp(i,j,k)}) }
       )

    end subroutine WriteComp


    ! **************************************************************** !

    subroutine WriteEps(out)

      type (T_OUT) :: out

      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  
      reg = regobj(out%regidx)
      
      M4_REGLOOP_WRITE(reg,p,i,j,k,w,
      out%funit,E15.6E3, { 1.0/epsinv(i,j,k) }
      )

    end subroutine WriteEps


    ! **************************************************************** !

    subroutine WriteBufData(out)

      implicit none
  
      type (T_OUT) :: out
      M4_FTYPE :: sum = 0.0
      type (T_BUF) :: buf
      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  

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
	  buf%data(:,1) = 0.0	
      
    end subroutine WriteBufData


  ! ***************************************************************** !

  subroutine LoadBufPx(out)

    ! calculates and loads the x-component of the poynting vector
    ! @ Px[i,j,k] = Px(i+1/2,j,k)

    type (T_OUT) :: out
    M4_FTYPE :: val
    type (T_BUF) :: buf
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  
    
    M4_WRITE_DBG({"load data buffer # ",&
         TRIM(i2str(out%bufidx))," for out # ",TRIM(i2str(out%idx))})

    reg = regobj(out%regidx)
    buf = bufobj(out%bufidx)

    M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
    val = 0.125*((Ey(M4_COORD(i,j,k))+Ey(M4_COORD(i+1,j,k)))*Hz(M4_COORD(i,j,k)) &
         + (Ey(M4_COORD(i,j-1,k))+Ey(M4_COORD(i+1,j-1,k)))* Hz(M4_COORD(i,j-1,k)) &
         - (Ez(M4_COORD(i,j,k))+Ez(M4_COORD(i+1,j,k)))*Hy(M4_COORD(i,j,k)) &
         - (Ez(M4_COORD(i,j,k-1))+Ez(M4_COORD(i+1,j,k-1)))* Hy(M4_COORD(i,j,k-1)) &
         ) 
    buf%data(p,1) = buf%data(p,1) + val/(4.0*PI)
    })
      
  end subroutine LoadBufPx

  ! ***************************************************************** !

  subroutine LoadBufPy(out)

    ! calculates and loads the y-component of the poynting vector
    ! @ Py[i,j,k] = Py(i,j+1/2,k)

    type (T_OUT) :: out
    M4_FTYPE :: val
    type (T_BUF) :: buf
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  
    
    M4_WRITE_DBG({"load data buffer # ",&
         TRIM(i2str(out%bufidx))," for out # ",TRIM(i2str(out%idx))})

    reg = regobj(out%regidx)
    buf = bufobj(out%bufidx)
    
    M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
    val = 0.125*( (Ez(M4_COORD(i,j,k))+Ez(M4_COORD(i,j+1,k)))*Hx(M4_COORD(i,j,k)) &
         + (Ez(M4_COORD(i,j,k-1))+Ez(M4_COORD(i,j+1,k-1)))*Hx(M4_COORD(i,j,k-1)) &
         - (Ex(M4_COORD(i,j,k))+Ex(M4_COORD(i,j+1,k)))*Hz(M4_COORD(i,j,k)) &
         - (Ex(M4_COORD(i-1,j,k))+Ex(M4_COORD(i-1,j+1,k)))*Hz(M4_COORD(i-1,j,k)) &
         )
    buf%data(p,1) = buf%data(p,1) + val/(4.0*PI)
    } )

  end subroutine LoadBufPy

  ! ***************************************************************** !

  subroutine LoadBufPz(out)
   
    ! calculates and loads the z-component of the poynting vector
    ! @ Pz[i,j,k] = Pz(i,j,k+1/2)

    type (T_OUT) :: out
    M4_FTYPE :: val
    type (T_BUF) :: buf
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  
    
    M4_WRITE_DBG({"load data buffer # ",TRIM(i2str(out%bufidx)),&
         " for out # ",TRIM(i2str(out%idx))})

    reg = regobj(out%regidx)
    buf = bufobj(out%bufidx)
    
    M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
    
    val = 0.125*( (Ex(M4_COORD(i,j,k))+Ex(M4_COORD(i,j,k+1)))*Hy(M4_COORD(i,j,k)) &
         +(Ex(M4_COORD(i-1,j,k))+Ex(M4_COORD(i-1,j,k+1)))*Hy(M4_COORD(i-1,j,k)) &
         - (Ey(M4_COORD(i,j,k))+Ey(M4_COORD(i,j,k+1)))*Hx(M4_COORD(i,j,k)) &
         - (Ey(M4_COORD(i,j-1,k))+Ey(M4_COORD(i,j-1,k+1)))* Hx(M4_COORD(i,j-1,k)) &
         )
    buf%data(p,1) = buf%data(p,1) + val/(4.0*PI)
    
    } )

  end subroutine LoadBufPz

  ! ***************************************************************** !

  subroutine LoadBufEn(out)

    ! calculates and loads energy density En
    ! @ En[i,j,k] = En(i,j,k)

    type (T_OUT) :: out
    M4_FTYPE :: EEx,EEy,EEz,EHx,EHy,EHz
    real(kind=8) :: eps
    type (T_BUF) :: buf
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  
    
    M4_WRITE_DBG({"load data buffer # ",&
         TRIM(i2str(out%bufidx))," for out # ",TRIM(i2str(out%idx))})

    reg = regobj(out%regidx)
    buf = bufobj(out%bufidx)

    M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
    
    eps = 1.0/epsinv(i,j,k)
    EEx = 0.5*eps*(Ex(M4_COORD(i,j,k))**2+Ex(M4_COORD(i-1,j,k))**2)
    EEy = 0.5*eps*(Ey(M4_COORD(i,j,k))**2+Ey(M4_COORD(i,j-1,k))**2)
    EEz = 0.5*eps*(Ez(M4_COORD(i,j,k))**2+Ez(M4_COORD(i,j,k-1))**2)
    EHx = 0.25*(Hx(M4_COORD(i,j,k))**2+Hx(M4_COORD(i,j,k-1))**2 + &
         Hx(M4_COORD(i,j-1,k))**2+Hx(M4_COORD(i,j-1,k-1))**2)
    EHy = 0.25*(Hy(M4_COORD(i,j,k))**2+Hy(M4_COORD(i-1,j,k))**2 + &
         Hy(M4_COORD(i,j,k-1))**2+Hy(M4_COORD(i-1,j,k-1))**2)
    EHz = 0.25*(Hz(M4_COORD(i,j,k))**2+Hz(M4_COORD(i-1,j,k))**2 + &
         Hz(M4_COORD(i,j-1,k))**2+Hz(M4_COORD(i-1,j-1,k))**2)
    buf%data(p,1) = (EEx+EEy+EEz+EHx+EHy+EHz)/(4.0*PI)
    
    } )
    
  end subroutine LoadBufEn

  ! ***************************************************************** !


  end subroutine WriteDataFdtdOutgplObj


end module fdtd_outgpl

!
! Authors:  S.Scholz,A.Klaedtke,C.Hermann,J.Hamm
! Modified: 4/12/2007
!
!======================================================================
