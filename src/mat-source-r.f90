!----------------------------------------------------------------------
!
!   module: mat_source
!
!   electromagnetic source
!
!----------------------------------------------------------------------


module mat_source

  use constant
  use mpiworld
  use grid
  use fdtd

  implicit none
  save

  integer :: ps_i = 24         ! position of point source
  integer :: ps_j = 97
  integer :: ps_k = 49

  real(kind=8) :: lambda0 = 86.59  ! vacuum wavelength in units of [dx]
  real(kind=8) :: dlambda0 = 1    ! spectral width of vac.wave in units of [dx]

  real(kind=8) :: a0 = 0.05       ! gaussian start value as fraction of peak
  real(kind=8) :: ampl = 1.0      ! amplitude = 1.

  logical      :: source_on
  real(kind=8) :: gamma
  real(kind=8) :: omega0, omega1 
  real(kind=8) :: npeak, nmax
  real(kind=8) :: es

contains

  subroutine InitSource

    omega0 = 2. * PI * 1. / ( lambda0 * DT )
    
    omega1 = 2. * PI * 1. / ( ( lambda0 + dlambda0 ) * DT )
    
    gamma = (omega0 - omega1 ) / log(2.0)
    
    npeak =  sqrt ( - log(a0) / gamma**2 )

    source_on = .false.

    nmax = npeak ! nmax = 1e10

    if ( ps_i .ge. IBEG .and. ps_i .le. IEND .and. &
         ps_j .ge. JBEG .and. ps_j .le. JEND .and. &
         ps_k .ge. KBEG .and. ps_k .le. KEND ) then

       source_on = .true.

       write (6,*) "----- Point Source"
       
       write(6,*) "omega0 = ", omega0
       write(6,*) "omega1 = ", omega1
       write(6,*) "gamma = ", gamma
       write(6,*) "npeak = ", npeak
       
    endif
       

  end subroutine InitSource

  subroutine SourceEy(ncyc)

    integer :: ncyc

    if ( .not. source_on ) return

    es = 1.0

    if ( ncyc .lt. nmax ) then
       es =  exp ( - gamma**2 * ( 1.0 * ncyc - npeak )**2 )
    endif

    Ey(ps_i,ps_j,ps_k) = Ey(ps_i,ps_j,ps_k) + es  * ampl * cos(omega0*ncyc) * DT

  end subroutine SourceEy


  subroutine SourceHy(ncyc)

    integer :: ncyc

    if ( .not. source_on ) return

    es = 1.0

    if ( ncyc .lt. nmax ) then
       es =  exp ( - gamma**2 * ( 1.0 * ncyc - npeak )**2 )
    endif

    Hy(ps_i,ps_j,ps_k) = Hy(ps_i,ps_j,ps_k) + es  * ampl * cos(omega0*ncyc) * DT

  end subroutine SourceHy

end module mat_source
