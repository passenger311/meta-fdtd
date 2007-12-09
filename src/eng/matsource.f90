!-*- F90 -*------------------------------------------------------------
!
!  module: matsource / meta3
!
!  field source equations.
!
!  subs:
!
!    InitializeMatSource
!    FinalizeMatSource
!    ReadMatSourceObj
!    StepEMatSource
!    StepHMatSource
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatSource module allows to add sources to the electromagnetic 
! field equations


module matsource

  use constant
  use mpiworld
  use reglist
  use grid
  use fdtd

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'MATSOURCE'

  ! --- Public Methods

  public :: InitializeMatSource
  public :: FinalizeMatSource
  public :: StepEMatSource
  public :: StepHMatSource
  public :: ReadMatSourceObj

  ! --- Public Data

  public :: T_MATSOURCE
  public :: matsourceobj
  public :: nummatsourceobj


  ! --- Constants

  integer, parameter :: MAXMATSOURCEOBJ = 100

  ! --- Types

  type T_MATSOURCE

     integer :: regidx             ! regobj index

     real(kind=8) :: lambda0       ! vacuum wavelength in units of [dx]
     real(kind=8) :: dlambda0      ! spectral width of vac.wave in units of [dx]
     real(kind=8) :: a0            ! gaussian start value as fraction of peak
     real(kind=8) :: ampl          ! amplitude = 1.
     logical      :: cw            ! go over to cw after peak?
     logical      :: esource       ! e (or h) field source
     real(kind=8) :: vec(3)        ! vector

     real(kind=8) :: gamma         ! some internal values ...
     real(kind=8) :: omega0, omega1 
     real(kind=8) :: npeak, nmax
     real(kind=8) :: es

  end type T_MATSOURCE

  ! --- Fields

  type(T_MATSOURCE) :: matsourceobj(MAXMATSOURCEOBJ) 
  integer :: nummatsourceobj

contains

!----------------------------------------------------------------------

  subroutine ReadMatSourceObj(funit)

    integer:: funit
    character(len=STRLNG) :: file, string
    integer :: ios
    type(T_REG) :: reg
    type(T_MATSOURCE) :: mat

   M4_WRITE_DBG(". enter ReadMatSourceObj/matsource")

    nummatsourceobj = nummatsourceobj + 1
    mat = matsourceobj(nummatsourceobj)

! read parameters

    read(funit,*) mat%lambda0     ! vacuum wavelength in units of [dx]
   M4_WRITE_DBG({"read lambda0: ", mat%lambda0})
    read(funit,*) mat%dlambda0    ! spectral width of vac.wave in units of [dx]
   M4_WRITE_DBG({"read dlambda0: ", mat%dlambda0})
    read(funit,*) mat%a0          ! gaussian start value as fraction of peak
   M4_WRITE_DBG({"read a0: ", mat%a0})
    read(funit,*) mat%ampl        ! amplitude = 1.
   M4_WRITE_DBG({"read ampl: ", mat%ampl})
    read(funit,*) mat%cw          ! go over to cw after peak
   M4_WRITE_DBG({"read cw: ", mat%cw})
    read(funit,*) mat%esource     ! electric/magnetic field source
   M4_WRITE_DBG({"read esource: ", mat%esource})
    read(funit,*) mat%vec(1),mat%vec(2), mat%vec(3) ! vector components
   M4_WRITE_DBG({"read vec(1:3): ", mat%vec(1),mat%vec(2), mat%vec(3) })

! read regions and terminator

   M4_GET_REG_AND_TERMINATOR(mat, "ReadMatSourceObj/matsource")

   M4_WRITE_DBG(". exit ReadMatSourceObj/matsource")

  end subroutine ReadMatSourceObj

!----------------------------------------------------------------------

  subroutine InitializeMatSource

    integer :: n
    type(T_MATSOURCE) :: mat

    do n = 1, nummatsourceobj

       mat = matsourceobj(n)

       ! Initialize matsource object

       mat%omega0 = 2. * PI * 1. / ( mat%lambda0 * DT )
       mat%omega1 = 2. * PI * 1. / ( ( mat%lambda0 + mat%dlambda0 ) * DT )
       mat%gamma = (mat%omega0 - mat%omega1 ) / log(2.0)
       mat%npeak =  sqrt ( - log(mat%a0) / mat%gamma**2 )
       mat%nmax = mat%npeak

    end do

  end subroutine InitializeMatSource

!----------------------------------------------------------------------

  subroutine FinalizeMatSource

    integer :: n
    type(T_MATSOURCE) :: mat

    do n = 1, nummatsourceobj

       mat = matsourceobj(n)

       ! Finalize matdsource object here 

    end do

  end subroutine FinalizeMatSource

!----------------------------------------------------------------------

  subroutine StepEMatSource(ncyc)

    integer :: ncyc
    type(T_MATSOURCE) :: mat
    integer :: n
    real(kind=8) :: es
    M4_REGLOOP_DECL(reg,p,i,j,k,w)

    do n = 1, nummatsourceobj

       mat = matsourceobj(n)

       if ( .not. mat%esource ) return
       
       reg = regobj(mat%regidx)
       
       es = 1.0
       
       if ( mat%cw .and. ncyc .lt. mat%nmax ) then
          es =  exp ( - mat%gamma**2 * ( 1.0 * ncyc - mat%npeak )**2 )
       endif
       
       M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
       
       Ex(i,j,k) = Ex(i,j,k) + mat%vec(1) * es  * mat%ampl * cos(mat%omega0*ncyc) * DT
       Ey(i,j,k) = Ey(i,j,k) + mat%vec(2) * es  * mat%ampl * cos(mat%omega0*ncyc) * DT
       Ez(i,j,k) = Ez(i,j,k) + mat%vec(3) * es  * mat%ampl * cos(mat%omega0*ncyc) * DT
       } )            
       
       
    end do

  end subroutine StepEMatSource

!----------------------------------------------------------------------


  subroutine StepHMatSource(ncyc)

    integer :: ncyc
    integer :: n
    type(T_MATSOURCE) :: mat
    real(kind=8) :: es
    M4_REGLOOP_DECL(reg,p,i,j,k,w)

    do n = 1, nummatsourceobj

       mat = matsourceobj(n)

       if ( mat%esource ) return
       
       reg = regobj(mat%regidx)
       
       es = 1.0

       if ( mat%cw .and. ncyc .lt. mat%nmax ) then
          es =  exp ( - mat%gamma**2 * ( 1.0 * ncyc - mat%npeak )**2 )
       endif
    
       M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
       
       Hx(i,j,k) = Hx(i,j,k) + mat%vec(1) * es  * mat%ampl * cos(mat%omega0*ncyc) * DT
       Hy(i,j,k) = Hy(i,j,k) + mat%vec(2) * es  * mat%ampl * cos(mat%omega0*ncyc) * DT
       Hz(i,j,k) = Hz(i,j,k) + mat%vec(3) * es  * mat%ampl * cos(mat%omega0*ncyc) * DT
       
       } )
       
       
    end do
    
  end subroutine StepHMatSource

  
!----------------------------------------------------------------------

  subroutine EchoMatSourceObj(mat)
  
    type(T_MATSOURCE) :: mat

    write (6,*) "----- Point Source"
       
    write(6,*) "omega0 = ", mat%omega0
    write(6,*) "omega1 = ", mat%omega1
    write(6,*) "gamma = ", mat%gamma
    write(6,*) "npeak = ", mat%npeak
    
 
  end subroutine EchoMatSourceObj
  
!----------------------------------------------------------------------

end module matsource

! =====================================================================


