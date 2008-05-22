
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

  public :: GaussianWave


contains

!----------------------------------------------------------------------


  real(kind=8) function GaussianWave(ncyc, noffs, natt, nsus, ndcy, gamma, omega)

    real(kind=8) :: ncyc, noffs, natt, nsus, ndcy
    real(kind=8) :: gamma, omega
    real(kind=8) :: ncyc0, nend, amp


    ncyc0 = ncyc - noffs
    nend = natt + nsus + ndcy
    
    GaussianWave = 0.


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
       
       GaussianWave = amp * cos(omega*(ncyc0 - natt)*DT) ! cos will peak at noffs+natt

    end if

  end function GaussianWave
    
!----------------------------------------------------------------------


end module signal


!
! Authors:  J.Hamm
! Modified: 25/02/2007
!
!======================================================================

