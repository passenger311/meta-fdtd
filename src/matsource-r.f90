!-*- F90 -*------------------------------------------------------------
!
!  module: matsource / max3d
!
!  field source equations.
!
!  subs:
!
!    InitializeMatSource
!    FinalizeMatSource
!    ReadObjMatSource
!    StepEMatSource
!    StepHMatSource
!    StepEObjMatSource
!    StepHObjMatSource
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatSource module allows to add sources to the electromagnetic 
! field equations


module matsource

  use constant
  use mpiworld
  use regobj
  use grid
  use fdtd

  implicit none
  save

  ! --- Constants

  integer, parameter :: MAXOBJMATSOURCE = 100

  ! --- Types

  type T_MATSOURCE

     integer :: regidx ! regobj index

     real(kind=8) :: lambda0    ! vacuum wavelength in units of [dx]
     real(kind=8) :: dlambda0   ! spectral width of vac.wave in units of [dx]
     real(kind=8) :: a0         ! gaussian start value as fraction of peak
     real(kind=8) :: ampl       ! amplitude = 1.
     logical      :: cw         ! go over to cw after peak?
     logical      :: esource    ! e (or h) field source
     real(kind=8) :: vec(3)     ! vector

! --- internal values

     real(kind=8) :: gamma
     real(kind=8) :: omega0, omega1 
     real(kind=8) :: npeak, nmax
     real(kind=8) :: es

  end type T_MATSOURCE


  ! --- Fields

  type(T_MATSOURCE) :: objmatsource(MAXOBJMATSOURCE) 
  integer :: numobjmatsource

contains

!----------------------------------------------------------------------

  subroutine InitializeMatSource

  end subroutine InitializeMatSource

!----------------------------------------------------------------------

  subroutine FinalizeMatSource

  end subroutine FinalizeMatSource

!----------------------------------------------------------------------

  subroutine ReadObjMatSource(funit)

    integer:: funit
    character(len=STRLNG) :: file, string
    integer :: ios
    type(T_REGION) :: reg
    type(T_MATSOURCE) :: mat

    numobjmatsource = numobjmatsource + 1
    mat = objmatsource(numobjmatsource)

! read parameters

    read(funit,*) mat%lambda0     ! vacuum wavelength in units of [dx]
    read(funit,*) mat%dlambda0    ! spectral width of vac.wave in units of [dx]
    read(funit,*) mat%a0          ! gaussian start value as fraction of peak
    read(funit,*) mat%ampl        ! amplitude = 1.
    read(funit,*) mat%cw          ! go over to cw after peak
    read(funit,*) mat%esource     ! electric field source (or magnetic)
    read(funit,*) mat%vec(1),mat%vec(2), mat%vec(3) ! vector components

! read regobj information
    read(funit,*) string
    if ( string .eq. "(REGION" ) then
       call ReadObjReg(reg, funit)
       mat%regidx = reg%idx
    else
       write(6,*) "!ERROR NO REGION DEFINED: ReadObjMatSource"
       stop
    end if

    read(funit,*) string
    if ( string .ne. ")" ) then
       write(6,*) "!ERROR BAD TERMINATOR: ReadObjMatSource"
       stop
    end if

! initialize object data

    mat%omega0 = 2. * PI * 1. / ( mat%lambda0 * DT )
    mat%omega1 = 2. * PI * 1. / ( ( mat%lambda0 + mat%dlambda0 ) * DT )
    mat%gamma = (mat%omega0 - mat%omega1 ) / log(2.0)
    mat%npeak =  sqrt ( - log(mat%a0) / mat%gamma**2 )
    mat%nmax = mat%npeak

  end subroutine ReadObjMatSource

!----------------------------------------------------------------------
 
  subroutine StepEMatSource(ncyc)

    integer :: ncyc
    integer :: n

    do n = 1, numobjmatsource
       
       call StepEObjMatSource(objmatsource(n),ncyc)

    end do

  end subroutine StepEMatSource

!----------------------------------------------------------------------


  subroutine StepEObjMatSource(mat,ncyc)

    type(T_MATSOURCE) :: mat
    integer :: ncyc
    
    type(T_REGION) :: reg
    integer :: p,i,j,k,il,jl,kl
    real(kind=8) :: es

    reg = objregobj(mat%regidx)

    es = 1.0

    if ( mat%cw .and. ncyc .lt. mat%nmax ) then
       es =  exp ( - mat%gamma**2 * ( 1.0 * ncyc - mat%npeak )**2 )
    endif
    
    ! *** There is a bit of magic in this loop structure. Note, that either 
    ! (di,dj,dk) = (ie-is+1, je-js+1, ke-ks+1) or reg%nump = 1. This means,
    ! that this loop either acts an indirect or a masked loop. ***
    do il = reg%is, reg%ie, reg%di
       do jl = reg%js, reg%je, reg%dj
          do kl = reg%ks, reg%ke, reg%dk
             do p = 1, reg%nump
                if ( reg%islist ) then
                   i = reg%i(p)
                   j = reg%j(p)
                   k = reg%k(p)
                else
                   i = il
                   j = jl
                   k = kl
                end if
                Ex(i,j,k) = Ex(i,j,k) + mat%vec(1) * es  * reg%mask(i,j,k) * mat%ampl * cos(mat%omega0*ncyc) * DT
                Ey(i,j,k) = Ey(i,j,k) + mat%vec(2) * es  * reg%mask(i,j,k) * mat%ampl * cos(mat%omega0*ncyc) * DT
                Ez(i,j,k) = Ez(i,j,k) + mat%vec(3) * es  * reg%mask(i,j,k) * mat%ampl * cos(mat%omega0*ncyc) * DT
                
             end do
          end do
       end do
    end do
       
  end subroutine StepEObjMatSource

!----------------------------------------------------------------------

  subroutine StepHMatSource(ncyc)

    integer :: ncyc
    integer :: n

    do n = 1, numobjmatsource
       
       call StepHObjMatSource(objmatsource(n),ncyc)

    end do

  end subroutine StepHMatSource


!----------------------------------------------------------------------


  subroutine StepHObjMatSource(mat,ncyc)

    type(T_MATSOURCE) :: mat
    integer :: ncyc
    
    type(T_REGION) :: reg
    integer :: p,i,j,k,il,jl,kl
    real(kind=8) :: es

    reg = objregobj(mat%regidx)

    es = 1.0

    if ( mat%cw .and. ncyc .lt. mat%nmax ) then
       es =  exp ( - mat%gamma**2 * ( 1.0 * ncyc - mat%npeak )**2 )
    endif
    
    ! *** There is a bit of magic in this loop structure. Note, that either 
    ! (di,dj,dk) = (ie-is+1, je-js+1, ke-ks+1) or reg%nump = 1. This means,
    ! that this loop either acts an indirect or a masked loop ***
    do il = reg%is, reg%ie, reg%di
       do jl = reg%js, reg%je, reg%dj
          do kl = reg%ks, reg%ke, reg%dk
             do p = 1, reg%nump
                if ( reg%islist ) then
                   i = reg%i(p)
                   j = reg%j(p)
                   k = reg%k(p)
                else
                   i = il
                   j = jl
                   k = kl
                end if
                
                Hx(i,j,k) = Hx(i,j,k) + mat%vec(1) * es  * reg%mask(i,j,k) * mat%ampl * cos(mat%omega0*ncyc) * DT
                Hy(i,j,k) = Hy(i,j,k) + mat%vec(2) * es  * reg%mask(i,j,k) * mat%ampl * cos(mat%omega0*ncyc) * DT
                Hz(i,j,k) = Hz(i,j,k) + mat%vec(3) * es  * reg%mask(i,j,k) * mat%ampl * cos(mat%omega0*ncyc) * DT
                
             end do
          end do
       end do
    end do
       
  end subroutine StepHObjMatSource


!----------------------------------------------------------------------


  subroutine EchoObjMatSource(mat)
  
    type(T_MATSOURCE) :: mat

    write (6,*) "----- Point Source"
       
    write(6,*) "omega0 = ", mat%omega0
    write(6,*) "omega1 = ", mat%omega1
    write(6,*) "gamma = ", mat%gamma
    write(6,*) "npeak = ", mat%npeak
    
 
  end subroutine EchoObjMatSource
  
!----------------------------------------------------------------------

end module matsource

! =====================================================================


