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
    buf%data(p,slot) = (EEx+EEy+EEz+EHx+EHy+EHz)/(4.0*PI)
    
    } )
  
    
  end subroutine FdtdCalcEn

!----------------------------------------------------------------------

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

     val = 0.125*((Ey(M4_COORD(i,j,k))+Ey(M4_COORD(i+1,j,k)))*Hz(M4_COORD(i,j,k)) &
         + (Ey(M4_COORD(i,j-1,k))+Ey(M4_COORD(i+1,j-1,k)))* Hz(M4_COORD(i,j-1,k)) &
         - (Ez(M4_COORD(i,j,k))+Ez(M4_COORD(i+1,j,k)))*Hy(M4_COORD(i,j,k)) &
         - (Ez(M4_COORD(i,j,k-1))+Ez(M4_COORD(i+1,j,k-1)))* Hy(M4_COORD(i,j,k-1)) &
         ) 
    buf%data(p,slot) = buf%data(p,slot) + val/(4.0*PI)
    
    } )
  
    
  end subroutine FdtdCalcPx

!----------------------------------------------------------------------

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

    val = 0.125*( (Ez(M4_COORD(i,j,k))+Ez(M4_COORD(i,j+1,k)))*Hx(M4_COORD(i,j,k)) &
         + (Ez(M4_COORD(i,j,k-1))+Ez(M4_COORD(i,j+1,k-1)))*Hx(M4_COORD(i,j,k-1)) &
         - (Ex(M4_COORD(i,j,k))+Ex(M4_COORD(i,j+1,k)))*Hz(M4_COORD(i,j,k)) &
         - (Ex(M4_COORD(i-1,j,k))+Ex(M4_COORD(i-1,j+1,k)))*Hz(M4_COORD(i-1,j,k)) &
         )
    buf%data(p,slot) = buf%data(p,slot) + val/(4.0*PI)
    
    } )
  
    
  end subroutine FdtdCalcPy

!----------------------------------------------------------------------

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

    val = 0.125*( (Ex(M4_COORD(i,j,k))+Ex(M4_COORD(i,j,k+1)))*Hy(M4_COORD(i,j,k)) &
         +(Ex(M4_COORD(i-1,j,k))+Ex(M4_COORD(i-1,j,k+1)))*Hy(M4_COORD(i-1,j,k)) &
         - (Ey(M4_COORD(i,j,k))+Ey(M4_COORD(i,j,k+1)))*Hx(M4_COORD(i,j,k)) &
         - (Ey(M4_COORD(i,j-1,k))+Ey(M4_COORD(i,j-1,k+1)))* Hx(M4_COORD(i,j-1,k)) &
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
