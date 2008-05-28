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
  use numerics
  use grid 
  use fdtd

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'FDTD_CALC'

  ! --- Private Data

  real(kind=8) :: np_omega, np_kinc(3), np_sx, np_sy, np_sz, np_nref

 ! --- Public Methods

  public :: InitializeFdtdCalc
  public :: FinalizeFdtdCalc
  public :: FdtdCalcEn, FdtdCalcSx, FdtdCalcSy, FdtdCalcSz
  public :: NumericalPhaseVelocity

contains

!----------------------------------------------------------------------

  subroutine InitializeFdtdCalc


  end subroutine InitializeFdtdCalc

!----------------------------------------------------------------------

  subroutine FinalizeFdtdCalc
   

  end subroutine FinalizeFdtdCalc

!----------------------------------------------------------------------


  subroutine NumericalPhaseVelocity(kinc, omega, nref, sx, sy, sz, pvel)

    implicit none

    real(kind=8) :: kinc(3), omega, nref, sx,sy,sz, pvel

    real(kind=8) :: val, xr, xl, xh = 1.0, xacc = 0.00001
    integer :: iter = 100

    xr = 1.5 / nref
    xl = 0.5 / nref

    ! calculate and return numerical phase velocity <pvel> from <omega>
    ! and <kinc>.

    val = 0.

    np_omega = omega
    np_kinc = kinc
    np_nref = nref
    np_sx = sx
    np_sy = sy
    np_sz = sz

    call root_find(xr, val, npfct, xl, xh, xacc, iter)

    pvel = xr

  end subroutine NumericalPhaseVelocity


  real(kind=8) function npfct(vel)
      
    real(kind=8) :: vel ! phase velocity
    
    ! see [tavlov p110]

    npfct = ( 1./np_sx * sin( np_omega * np_kinc(1) * np_sx  / ( 2. * vel ) ))**2 &
         + ( 1./np_sy * sin( np_omega * np_kinc(2) * np_sy   / ( 2. * vel ) ))**2  &
         + ( 1./np_sz * sin( np_omega * np_kinc(3) * np_sz   / ( 2. * vel ) ))**2  &
         - ( np_nref/DT * sin( np_omega * DT / 2. ))**2 
    
  end function npfct



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
    
    EEx = 0.5/epsinvx(i,j,k) * (real(Ex(M4_COORD(i,j,k)))**2+real(Ex(M4_COORD(i-1,j,k)))**2)
    EEy = 0.5/epsinvy(i,j,k) * (real(Ey(M4_COORD(i,j,k)))**2+real(Ey(M4_COORD(i,j-1,k)))**2)
    EEz = 0.5/epsinvz(i,j,k) * (real(Ez(M4_COORD(i,j,k)))**2+real(Ez(M4_COORD(i,j,k-1)))**2)
    EHx = 0.25/M4_MUINVX(i,j,k) * ( real(Hx(M4_COORD(i,j,k)))**2+real(Hx(M4_COORD(i,j,k-1)))**2 + &
         real(Hx(M4_COORD(i,j-1,k)))**2+real(Hx(M4_COORD(i,j-1,k-1)))**2)
    EHy = 0.25/M4_MUINVY(i,j,k) * (real(Hy(M4_COORD(i,j,k)))**2+real(Hy(M4_COORD(i-1,j,k)))**2 + &
         real(Hy(M4_COORD(i,j,k-1)))**2+real(Hy(M4_COORD(i-1,j,k-1)))**2)
    EHz = 0.25/M4_MUINVZ(i,j,k) * (real(Hz(M4_COORD(i,j,k)))**2+real(Hz(M4_COORD(i-1,j,k)))**2 + &
         real(Hz(M4_COORD(i,j-1,k)))**2+real(Hz(M4_COORD(i-1,j-1,k)))**2)
    buf%data(p,slot) = (EEx+EEy+EEz+EHx+EHy+EHz)
    
    } )
  
    
  end subroutine FdtdCalcEn

!----------------------------------------------------------------------

! Localization: Sx[i,j,k] = Sx(i+1/2,j,k)

  subroutine FdtdCalcSx(buf, slot, mode)

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

     val = 0.125*(real(Ey(M4_COORD(i,j,k))+Ey(M4_COORD(i+1,j,k)))*real(Hz(M4_COORD(i,j,k))) &
         + real(Ey(M4_COORD(i,j-1,k))+Ey(M4_COORD(i+1,j-1,k)))* real(Hz(M4_COORD(i,j-1,k))) &
         - real(Ez(M4_COORD(i,j,k))+Ez(M4_COORD(i+1,j,k)))*real(Hy(M4_COORD(i,j,k))) &
         - real(Ez(M4_COORD(i,j,k-1))+Ez(M4_COORD(i+1,j,k-1)))* real(Hy(M4_COORD(i,j,k-1))) &
         ) 
    buf%data(p,slot) = buf%data(p,slot) + val
    
    } )
  
    
  end subroutine FdtdCalcSx

!----------------------------------------------------------------------

! Localization: Sy[i,j,k] = Sy(i-1/2,j+1,k+1/2) where Ey sits at n = n + 1/2

  subroutine FdtdCalcSy(buf, slot, mode)

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

    val = 0.125*( real(Ez(M4_COORD(i,j,k))+Ez(M4_COORD(i,j+1,k)))*real(Hx(M4_COORD(i,j,k))) &
         + real(Ez(M4_COORD(i,j,k-1))+Ez(M4_COORD(i,j+1,k-1)))*real(Hx(M4_COORD(i,j,k-1))) &
         - real(Ex(M4_COORD(i,j,k))+Ex(M4_COORD(i,j+1,k)))*real(Hz(M4_COORD(i,j,k))) &
         - real(Ex(M4_COORD(i-1,j,k))+Ex(M4_COORD(i-1,j+1,k)))*real(Hz(M4_COORD(i-1,j,k))) &
         )
    buf%data(p,slot) = buf%data(p,slot) + val
    
    } )
  
    
  end subroutine FdtdCalcSy

!----------------------------------------------------------------------

! Localization: Sz[i,j,k] = Sz(i-1/2,j+1/2,k+1).where Ez sits at n = n + 1/2
 
  subroutine FdtdCalcSz(buf, slot, mode)

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

    val = 0.125*( real(Ex(M4_COORD(i,j,k))+Ex(M4_COORD(i,j,k+1)))*real(Hy(M4_COORD(i,j,k))) &
         + real(Ex(M4_COORD(i-1,j,k))+Ex(M4_COORD(i-1,j,k+1)))*real(Hy(M4_COORD(i-1,j,k))) &
         - real(Ey(M4_COORD(i,j,k))+Ey(M4_COORD(i,j,k+1)))*real(Hx(M4_COORD(i,j,k))) &
         - real(Ey(M4_COORD(i,j-1,k))+Ey(M4_COORD(i,j-1,k+1)))*real(Hx(M4_COORD(i,j-1,k))) &
         )

    buf%data(p,slot) = buf%data(p,slot) + val
    
    } )
  
    
  end subroutine FdtdCalcSz

!----------------------------------------------------------------------

end module fdtd_calc

!
! Authors:  J.Hamm
! Modified: 6/1/2008
!
!======================================================================
