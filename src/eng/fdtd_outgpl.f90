!-*- F90 -*------------------------------------------------------------
!
!  module: fdtd_outgpl / meta3
!
!  this module handles GPL output of data related to the fdtd module.
!
!  subs:
!
!  InitializeFdtdOutgplObj
!  FinalizeFdtdOutgplObj
!  WriteDataFdtdOutgplObj
!
!----------------------------------------------------------------------


module fdtd_outgpl

  use constant
  use strings
  use outlist
  use reglist
  use buflist
  use mpiworld
  use grid 
  use fdtd

  implicit none
  save

contains

!----------------------------------------------------------------------

  subroutine InitializeFdtdOutgplObj(out)

    type (T_OUT) :: out
    type (T_BUF) :: buf

    if ( out%fn .eq. 'Px' .or. out%fn .eq. 'Py' .or. out%fn .eq. 'Pz' .or. &
         out%fn .eq. 'En' ) then

       ! allocate additional buffer for Load functions here

       buf = CreateBufObj(regobj(out%regidx))
       out%bufidx = buf%idx

    endif

  end subroutine InitializeFdtdOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeFdtdOutgplObj(out)

    type (T_OUT) :: out

    if ( out%fn .eq. 'Px' .or. out%fn .eq. 'Py' .or. out%fn .eq. 'Pz' .or. &
         out%fn .eq. 'En' ) then

       ! destroy additional buffer for Load functions here

       call DestroyBufObj(bufobj(out%bufidx))
       
    endif

  end subroutine FinalizeFdtdOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataFdtdOutgplObj(out, ncyc, mode)

    type (T_OUT) :: out
    type (T_BUF) :: buf
    integer :: ncyc
    logical :: mode
    M4_FTYPE :: clearval = 0.0

    if ( .not. mode .and. out%bufidx .gt. 0 ) then

       buf = bufobj(out%bufidx)
       call SetBufObj(buf,clearval)

    end if

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
    end select
    
  contains

    subroutine WriteEH(out)

      type (T_OUT) :: out

      M4_REGLOOP_DECL(reg,p,i,j,k,w)  
      reg = regobj(out%regidx)
      
      M4_REGLOOP_WRITE(reg,p,i,j,k,w,   
      UNITTMP,6E15.6E3, { Ex(i,j,k),Ey(i,j,k),Ez(i,j,k),Hx(i,j,k),Hy(i,j,k),Hz(i,j,k) }
      )
      
    end subroutine WriteEH


    ! **************************************************************** !

    subroutine WriteComp(out,Comp)

      type (T_OUT) :: out
      M4_FTYPE, dimension(IMIN:KMAX,JMIN:KMAX,KMIN:KMAX) :: Comp
      
      M4_REGLOOP_DECL(reg,p,i,j,k,w)  
      reg = regobj(out%regidx)

       M4_REGLOOP_WRITE(reg,p,i,j,k,w,   
       UNITTMP,E15.6E3, { Comp(i,j,k) }
       )

    end subroutine WriteComp


    ! **************************************************************** !

    subroutine WriteEps(out)

      type (T_OUT) :: out, sum

      M4_REGLOOP_DECL(reg,p,i,j,k,w)  
      reg = regobj(out%regidx)
      
      M4_REGLOOP_WRITE(reg,p,i,j,k,w,
      UNITTMP,E15.6E3, { 1.0/epsinv(i,j,k) }
      )

    end subroutine WriteEps


    ! **************************************************************** !

    subroutine WriteBufData(out)

      implicit none
  
      type (T_OUT) :: out
      M4_FTYPE :: sum = 0.0
      type(T_BUF) :: buf
      M4_REGLOOP_DECL(reg,p,i,j,k,w)  
      reg = regobj(out%regidx)

      buf = bufobj(out%bufidx)
  
      if ( out%mode .ne. 'S' ) then ! spatial integration
         M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
         sum = sum + buf%data(p)
         })
         write(UNITTMP,*) sum
      else
         M4_REGLOOP_WRITE(reg,p,i,j,k,w,
         UNITTMP,E15.6E3, {buf%data(p)}
         )
      endif

    end subroutine WriteBufData


  ! ***************************************************************** !

  subroutine LoadBufPx(out)

    ! calculates and loads the x-component of the poynting vector
    ! @ Px[i,j,k] = Px(i+1/2,j,k)

    type (T_OUT) :: out
    M4_FTYPE :: val
    type(T_BUF) :: buf
    M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    
    reg = regobj(out%regidx)
    buf = bufobj(out%bufidx)

    M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
    val = 0.125*((Ey(i,j,k)+Ey(i+1,j,k))*Hz(i,j,k) &
         + (Ey(i,j-1,k)+Ey(i+1,j-1,k))* Hz(i,j-1,k) &
         - (Ez(i,j,k)+Ez(i+1,j,k))*Hy(i,j,k) &
         - (Ez(i,j,k-1)+Ez(i+1,j,k-1))* Hy(i,j,k-1)) 
    buf%data(p) = buf%data(p) + val/(4.0*PI)
    })
      
  end subroutine LoadBufPx

  ! ***************************************************************** !

  subroutine LoadBufPy(out)

    ! calculates and loads the y-component of the poynting vector
    ! @ Py[i,j,k] = Py(i,j+1/2,k)

    type (T_OUT) :: out
    M4_FTYPE :: val
    type(T_BUF) :: buf
    M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    
    reg = regobj(out%regidx)
    buf = bufobj(out%bufidx)
    
    M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
    val = 0.125*( (Ez(i,j,k)+Ez(i,j+1,k))*Hx(i,j,k) &
         + (Ez(i,j,k-1)+Ez(i,j+1,k-1))*Hx(i,j,k-1) &
         - (Ex(i,j,k)+Ex(i,j+1,k))*Hz(i,j,k) &
         - (Ex(i-1,j,k)+Ex(i-1,j+1,k))*Hz(i-1,j,k))
    buf%data(p) = buf%data(p) + val/(4.0*PI)
    } )

  end subroutine LoadBufPy

  ! ***************************************************************** !

  subroutine LoadBufPz(out)
   
    ! calculates and loads the z-component of the poynting vector
    ! @ Pz[i,j,k] = Pz(i,j,k+1/2)

    type (T_OUT) :: out
    M4_FTYPE :: val
    type(T_BUF) :: buf
    M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    
    reg = regobj(out%regidx)
    buf = bufobj(out%bufidx)
    
    M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
    val = 0.125*( (Ex(i,j,k)+Ex(i,j,k+1))*Hy(i,j,k) &
         +(Ex(i-1,j,k)+Ex(i-1,j,k+1))*Hy(i-1,j,k) &
         - (Ey(i,j,k)+Ey(i,j,k+1))*Hx(i,j,k) &
         - (Ey(i,j-1,k)+Ey(i,j-1,k+1))* Hx(i,j-1,k))
    buf%data(p) = buf%data(p) + val/(4.0*PI)
    } )

  end subroutine LoadBufPz

  ! ***************************************************************** !

  subroutine LoadBufEn(out)

    ! calculates and loads energy density En
    ! @ En[i,j,k] = En(i,j,k)

    type (T_OUT) :: out
    M4_FTYPE :: EEx,EEy,EEz,EHx,EHy,EHz,eps
    type(T_BUF) :: buf
    M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    
    reg = regobj(out%regidx)
    buf = bufobj(out%bufidx)
    
    M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
    eps = 1.0/epsinv(i,j,k)
    EEx = 0.5*eps*(Ex(i,j,k)**2+Ex(i-1,j,k)**2)
    EEy = 0.5*eps*(Ey(i,j,k)**2+Ey(i,j-1,k)**2)
    EEz = 0.5*eps*(Ez(i,j,k)**2+Ez(i,j,k-1)**2)
    EHx = 0.25*(Hx(i,j,k)**2+Hx(i,j,k-1)**2 + &
         Hx(i,j-1,k)**2+Hx(i,j-1,k-1)**2)
    EHy = 0.25*(Hy(i,j,k)**2+Hy(i-1,j,k)**2 + &
         Hy(i,j,k-1)**2+Hy(i-1,j,k-1)**2)
    EHz = 0.25*(Hz(i,j,k)**2+Hz(i-1,j,k)**2 + &
         Hz(i,j-1,k)**2+Hz(i-1,j-1,k)**2)
    buf%data(p) = (EEx+EEy+EEz+EHx+EHy+EHz)/(4.0*PI)
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
