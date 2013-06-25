!-*- F90 -*------------------------------------------------------------
!
!  module: constant / meta
!
!  global parameters.
!
!----------------------------------------------------------------------

!======================================================================
!
!
!

module constant

  implicit none
  public
  save 

  integer, parameter :: MAXOUTOBJ  = 500
  integer, parameter :: MAXREGOBJ  = 500
  integer, parameter :: MAXINITOBJ = 100
  integer, parameter :: MAXVALOBJ  = 100
  integer, parameter :: MAXBUFOBJ  = 100
  integer, parameter :: MAXVALUES  = 10
  integer, parameter :: MAXSRCOBJ  = 50
  integer, parameter :: MAXMATOBJ  = 50
  integer, parameter :: MAXLUMPOBJ = 50
  integer, parameter :: MAXDIAGOBJ = 50
  integer, parameter :: MAXTIMER   = 10
  integer, parameter :: MAXEBALCH  = 100
  integer, parameter :: NUMEBALCH  = 10

  integer, parameter :: STRLNG=80                    ! max. string length
  integer, parameter :: LINELNG=160
  integer, parameter :: STDERR=0
  integer, parameter :: STDOUT=6

  character(len=20), parameter :: sfxin = ".in"
  character(len=20), parameter :: sfxout = ".out"

  integer, parameter :: UNITTMP=10                   ! default file unit number

  ! --- mathematical constants

  real(8), parameter :: PI=3.14159265358979323846    ! PI
  real(8), parameter :: DEG = PI/180.                ! convert deg -> rad 
  complex(kind=8), parameter :: IMAG = complex(0,1)  ! imaginary unit

  ! --- physical constants

  real(8), parameter :: SI_C = 299792458             ! speed of light in SI
  real(8), parameter :: SI_HBAR = 1.054571682364E-34 ! hbar in SI
  real(8), parameter :: SI_EPS0 = 8.85418781762E-12  ! eps0 in SI
  real(8), parameter :: SI_KB = 1.3806503E-23        ! boltzmann constant in SI
  real(8), parameter :: SI_ME = 9.10938188E-31       ! electron mass in SI (kg)
  real(8), parameter :: SI_E = 1.6021773E-19         ! elementary charge in SI
  real(8), parameter :: SI_4PIALPHA = 0.09170123649  ! 4 pi * finestructure constant
  


end module constant

!
! Authors:  J.Hamm
! Modified: 02/09/2009
!
!======================================================================


