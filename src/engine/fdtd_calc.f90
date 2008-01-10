!-*- F90 -*------------------------------------------------------------
!
!  module: fdtd_calc / meta
!
!  calculate and buffer fields from electromagnetic fields  
!
!----------------------------------------------------------------------


module fdtd_calc

  use constant
  use strings
  use reglist
  use buflist
  use mpiworld
  use grid 
  use fdtd

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'FDTD_CALC'

 ! --- Public Methods

  public :: InitializeFdtdCalc
  public :: FinalizeFdtdCalc
  public :: FdtdCalcEn, FdtdCalcPx, FdtdCalcPy, FdtdCalcPz

contains

!----------------------------------------------------------------------

  subroutine InitializeFdtdCalc


  end subroutine InitializeFdtdCalc

!----------------------------------------------------------------------

  subroutine FinalizeFdtdCalc
   

  end subroutine FinalizeFdtdCalc

!----------------------------------------------------------------------

  subroutine FdtdCalcEn(buf, slot, mode)

    ! Localization: En[i,j,k] = En(i,j,k)

    type (T_BUF) :: buf
    integer :: slot
    logical :: mode
    M4_FTYPE :: EEx,EEy,EEz,EHx,EHy,EHz
    real(kind=8) :: eps
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  
    
    if ( .not. mode ) return

    if ( slot .gt. buf%numslot ) then
       M4_FATAL_ERROR({"NOT ENOUGH BUFFER SLOTS ", &
            TRIM(i2str(slot))," > ", TRIM(i2str(buf%numslot)), " (FdtdCalcEnergy)"})
    endif

    reg = regobj(buf%regidx)

    M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
    
    EEx = 0.5/epsinvx(i,j,k) * (Ex(M4_COORD(i,j,k))**2+Ex(M4_COORD(i-1,j,k))**2)
    EEy = 0.5/epsinvy(i,j,k) * (Ey(M4_COORD(i,j,k))**2+Ey(M4_COORD(i,j-1,k))**2)
    EEz = 0.5/epsinvz(i,j,k) * (Ez(M4_COORD(i,j,k))**2+Ez(M4_COORD(i,j,k-1))**2)
    EHx = 0.25/M4_MUINVX(i,j,k) * (Hx(M4_COORD(i,j,k))**2+Hx(M4_COORD(i,j,k-1))**2 + &
         Hx(M4_COORD(i,j-1,k))**2+Hx(M4_COORD(i,j-1,k-1))**2)
    EHy = 0.25/M4_MUINVY(i,j,k) * (Hy(M4_COORD(i,j,k))**2+Hy(M4_COORD(i-1,j,k))**2 + &
         Hy(M4_COORD(i,j,k-1))**2+Hy(M4_COORD(i-1,j,k-1))**2)
    EHz = 0.25/M4_MUINVZ(i,j,k) * (Hz(M4_COORD(i,j,k))**2+Hz(M4_COORD(i-1,j,k))**2 + &
         Hz(M4_COORD(i,j-1,k))**2+Hz(M4_COORD(i-1,j-1,k))**2)
    buf%data(p,slot) = (EEx+EEy+EEz+EHx+EHy+EHz)/(4.0*PI)
    
    } )
  
    
  end subroutine FdtdCalcEn

!----------------------------------------------------------------------

! Localization: Px[i,j,k] = Px(i+1/2,j,k)

  subroutine FdtdCalcPx(buf, slot, mode)

    type (T_BUF) :: buf
    integer :: slot
    logical :: mode
    M4_FTYPE :: val
    real(kind=8) :: eps
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  
   
    if ( slot .gt. buf%numslot ) then
       M4_FATAL_ERROR({"NOT ENOUGH BUFFER SLOTS ", TRIM(i2str(slot)), &
            " > ", TRIM(i2str(buf%numslot)), " (FdtdCalcPx)"})
    endif
 
    if ( .not. mode ) then
       buf%data(:,slot) = 0.0
    endif

    reg = regobj(buf%regidx)

    M4_REGLOOP_EXPR(reg,p,i,j,k,w, {

     val = 0.125*((Ey(M4_COORD(i,j,k))+Ey(M4_COORD(i+1,j,k)))*M4_CONJ(Hz(M4_COORD(i,j,k))) &
         + (Ey(M4_COORD(i,j-1,k))+Ey(M4_COORD(i+1,j-1,k)))* M4_CONJ(Hz(M4_COORD(i,j-1,k))) &
         - (Ez(M4_COORD(i,j,k))+Ez(M4_COORD(i+1,j,k)))*M4_CONJ(Hy(M4_COORD(i,j,k))) &
         - (Ez(M4_COORD(i,j,k-1))+Ez(M4_COORD(i+1,j,k-1)))* M4_CONJ(Hy(M4_COORD(i,j,k-1))) &
         ) 
    buf%data(p,slot) = buf%data(p,slot) + val/(4.0*PI)
    
    } )
  
    
  end subroutine FdtdCalcPx

!----------------------------------------------------------------------

! Localization: Py[i,j,k] = Py(i-1/2,j+1,k+1/2) where Ey sits at n = n + 1/2

  subroutine FdtdCalcPy(buf, slot, mode)

    type (T_BUF) :: buf
    integer :: slot
    logical :: mode
    M4_FTYPE :: val
    real(kind=8) :: eps
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  
   
    if ( slot .gt. buf%numslot ) then
       M4_FATAL_ERROR({"NOT ENOUGH BUFFER SLOTS ", TRIM(i2str(slot)), &
            " > ", TRIM(i2str(buf%numslot)), " (FdtdCalcPy)"})
    endif
 
    if ( .not. mode ) then
       buf%data(:,slot) = 0.0
    endif

    reg = regobj(buf%regidx)

    M4_REGLOOP_EXPR(reg,p,i,j,k,w, {

    val = 0.125*( (Ez(M4_COORD(i,j,k))+Ez(M4_COORD(i,j+1,k)))*M4_CONJ(Hx(M4_COORD(i,j,k))) &
         + (Ez(M4_COORD(i,j,k-1))+Ez(M4_COORD(i,j+1,k-1)))*M4_CONJ(Hx(M4_COORD(i,j,k-1))) &
         - (Ex(M4_COORD(i,j,k))+Ex(M4_COORD(i,j+1,k)))*M4_CONJ(Hz(M4_COORD(i,j,k))) &
         - (Ex(M4_COORD(i-1,j,k))+Ex(M4_COORD(i-1,j+1,k)))*M4_CONJ(Hz(M4_COORD(i-1,j,k))) &
         )
    buf%data(p,slot) = buf%data(p,slot) + val/(4.0*PI)
    
    } )
  
    
  end subroutine FdtdCalcPy

!----------------------------------------------------------------------

! Localization: Pz[i,j,k] = Pz(i-1/2,j+1/2,k+1).where Ez sits at n = n + 1/2
 
  subroutine FdtdCalcPz(buf, slot, mode)

    type (T_BUF) :: buf
    integer :: slot
    logical :: mode
    M4_FTYPE :: val
    real(kind=8) :: eps
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  
   
    if ( slot .gt. buf%numslot ) then
       M4_FATAL_ERROR({"NOT ENOUGH BUFFER SLOTS ", TRIM(i2str(slot)), &
            " > ", TRIM(i2str(buf%numslot)), " (FdtdCalcPz)"})
    endif
 
    if ( .not. mode ) then
       buf%data(:,slot) = 0.0
    endif

    reg = regobj(buf%regidx)

    M4_REGLOOP_EXPR(reg,p,i,j,k,w, {

    val = 0.125*( (Ex(M4_COORD(i,j,k))+Ex(M4_COORD(i,j,k+1)))*M4_CONJ(Hy(M4_COORD(i,j,k))) &
         +(Ex(M4_COORD(i-1,j,k))+Ex(M4_COORD(i-1,j,k+1)))*M4_CONJ(Hy(M4_COORD(i-1,j,k))) &
         - (Ey(M4_COORD(i,j,k))+Ey(M4_COORD(i,j,k+1)))*M4_CONJ(Hx(M4_COORD(i,j,k))) &
         - (Ey(M4_COORD(i,j-1,k))+Ey(M4_COORD(i,j-1,k+1)))*M4_CONJ(Hx(M4_COORD(i,j-1,k))) &
         )

    buf%data(p,slot) = buf%data(p,slot) + val/(4.0*PI)
    
    } )
  
    
  end subroutine FdtdCalcPz

!----------------------------------------------------------------------

end module fdtd_calc

!
! Authors:  J.Hamm
! Modified: 6/1/2008
!
!======================================================================
