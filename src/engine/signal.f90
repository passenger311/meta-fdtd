
!-*- F90 -*------------------------------------------------------------
!
!  module: signal
!
!  time signal functions
!
!----------------------------------------------------------------------


! =====================================================================
!
! The first 5 parameters of each signal (modulation) function are:
!
! ncyc :  current timestep 
! noff :  offset
! natt :  attack period
! nsus :  sustain period
! ndcy :  decay period
! 
!                    +----------------------+
!                   /|                      |\
!                  / |                      | \
!                 /  |                      |  \
!                /   |                      |   \
!               /    |                      |    \
!              /     |                      |     \
!  0----------/------|----------------------|------+-------------> ncyc
!    noffs    | natt |       nsus           | ndcy |
!       
!
!
! The following parameters are specific parameters for the respective 
! signal function.
!

module signal

  use constant
  use grid
 
  implicit none
  private
  save

  character(len=20), private, parameter :: modname = 'SIGNAL'

  public :: GenericWave


contains

!----------------------------------------------------------------------


real(kind=8) function GenericWave(sigshape ,ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega) 

  character(len=20) :: sigshape
  real(kind=8) :: ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega

  select case ( sigshape ) 

  case ( "Linear" )
     GenericWave = LinearWave(ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega)
  case( "Sech" )
     GenericWave = SechWave(ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega)
  case ( "Ramp" )
     GenericWave = RampWave(ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega )
  case default
     GenericWave = GaussianWave(ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega)

  end select

end function GenericWave

real(kind=8) function StepFunction(n,w,m)
  real(kind=8) :: n, w, m
  real(kind=8) :: x

  x = n - (m*w)

  StepFunction = 1./(exp((0.85*w-x)/(0.02*w))+1.)

end function StepFunction

!---------------------------------------------------------------------

  real(kind=8) function LinearWave(ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega)

    real(kind=8) :: ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega
    real(kind=8) :: gamma
    real(kind=8) :: ncyc0, nend, amp


    ncyc0 = ncyc - noffs
    nend = natt + nsus + ndcy

    LinearWave = 0.

!    gamma = sqrt( log(2.) ) / ( nhwhm * DT )  ! calc gamma from HWHM

    if ( ncyc0 .ge. 0. .and. ncyc0 .le. nend ) then

       amp = 1.0

       ! attack phase
       if ( ncyc0 .le. natt ) then
          amp =  ( natt - ncyc0 ) / natt
       end if

       ! decay phase
       if ( ncyc0 .gt. natt + nsus ) then               
          amp =  1.0 - ( ncyc0 - natt - nsus ) / ndcy
       end if

       LinearWave = amp * cos( (omega+domega*(ncyc0-natt)*DT) * (ncyc0 - natt) * DT + DEG*alpha) ! cos will peak at noffs+natt

    end if

  end function LinearWave



!----------------------------------------------------------------------
  real(kind=8) function RampWave(ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega)

    real(kind=8) :: ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega
    real(kind=8) :: amp, n, nend, awidth, dwidth, stepnr, sheight, m

    ! copied from mathematica
    n = ncyc - noffs
    nend = natt + ndcy + nsus
    awidth = natt/nhwhm
    dwidth = ndcy/nhwhm
    sheight = 1./nhwhm

    if ( ( n .ge. 0 ) .and. ( n .le. ( nend - ( 0.5 * dwidth ) ) ) ) then
       if ( n .le. natt ) then
          stepnr = floor( n/awidth )
          amp = (StepFunction(n,awidth,stepnr)+stepnr)*sheight
       else
          if ( n .ge. ( natt+nsus - dwidth ) ) then
             m = n - ( natt+nsus );
             stepnr = floor( m/dwidth )
             amp = 1. - sheight*(StepFunction(m,dwidth,stepnr)+1.+stepnr)
          else
             amp = 1.
          endif
       endif
    else
       amp = 0.
    endif
   
    RampWave = amp * cos( (omega+domega*(n-natt)*DT) * (n - natt) * DT + DEG*alpha)
    
    
  end function RampWave
    

  real(kind=8) function GaussianWave(ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega)

    real(kind=8) :: ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega
    real(kind=8) :: gamma
    real(kind=8) :: ncyc0, nend, amp


    ncyc0 = ncyc - noffs
    nend = natt + nsus + ndcy
    
    GaussianWave = 0.

    gamma = sqrt( log(2.) ) / ( nhwhm * DT )  ! calc gamma from HWHM

    if ( ncyc0 .ge. 0. .and. ncyc0 .le. nend ) then 

       amp = 1.0

       ! attack phase
       if ( ncyc0 .le. natt ) then
          amp =  exp ( - gamma**2 * ( ( ncyc0 - natt ) * DT )**2 )
       end if
       
       ! decay phase
       if ( ncyc0 .gt. natt + nsus ) then         	
          amp =  exp ( - gamma**2 * ( ( ncyc0 - natt - nsus ) * DT )**2 )
       end if
       
       GaussianWave = amp * cos( (omega+domega*(ncyc0-natt)*DT) * (ncyc0 - natt) * DT + DEG*alpha) ! cos will peak at noffs+natt

    end if

  end function GaussianWave
    
!----------------------------------------------------------------------


  real(kind=8) function SechWave(ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega)

    real(kind=8) :: ncyc, noffs, natt, nsus, ndcy, nhwhm, omega, alpha, domega
    real(kind=8) :: gamma
    real(kind=8) :: ncyc0, nend, amp


    ncyc0 = ncyc - noffs
    nend = natt + nsus + ndcy
    
    SechWave = 0.

    gamma = log( 2. + sqrt(3.) )/ ( nhwhm * DT ) ! calc gamma from HWHM

    if ( ncyc0 .ge. 0. .and. ncyc0 .le. nend ) then 

       amp = 1.0

       ! attack phase
       if ( ncyc0 .le. natt ) then
          amp =  1./cosh( gamma * ( ncyc0 - natt ) * DT )
       end if
       
       ! decay phase
       if ( ncyc0 .gt. natt + nsus ) then         	
          amp =  1./cosh( gamma * ( ncyc0 - natt - nsus ) * DT )
       end if
       
       SechWave = amp * cos( (omega+domega*(ncyc0-natt) * DT ) * (ncyc0 - natt) * DT + DEG*alpha) ! cos will peak at noffs+natt


    end if

  end function SechWave
    
!----------------------------------------------------------------------


end module signal


!
! Authors:  J.Hamm
! Modified: 25/02/2007
!
!======================================================================

