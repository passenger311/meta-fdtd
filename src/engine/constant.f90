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

  integer, parameter :: MAXOUTOBJ=1000
  integer, parameter :: MAXREGOBJ=1000
  integer, parameter :: MAXINITOBJ=100
  integer, parameter :: MAXVALOBJ=100
  integer, parameter :: MAXBUFOBJ=500
  integer, parameter :: MAXVALUES=10
 
  integer, parameter :: UNITTMP=10                 ! file unit number
  real(8), parameter :: PI=3.14159265358979323846  ! PI
  integer, parameter :: STRLNG=80                  ! max. string length
  integer, parameter :: STDERR=0                   ! 6 for sr8000k, otherwise 0
  integer, parameter :: STDOUT=6                   ! 

  character(len=20), parameter :: sfxin = ".in"
  character(len=20), parameter :: sfxout = ".out"

contains



end module constant

!
! Authors:  S.Scholz, J.Hamm, A.Klaedtke
! Modified: 4/12/2007
!
!======================================================================


