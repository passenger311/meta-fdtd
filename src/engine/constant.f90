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

  integer, parameter :: MAXOUTOBJ=500
  integer, parameter :: MAXREGOBJ=500
  integer, parameter :: MAXINITOBJ=100
  integer, parameter :: MAXVALOBJ=100
  integer, parameter :: MAXBUFOBJ=100
  integer, parameter :: MAXVALUES=10
  integer, parameter :: MAXMATOBJ=50
  integer, parameter :: MAXDIAGOBJ=50
 
  integer, parameter :: UNITTMP=10                 ! file unit number
  real(8), parameter :: PI=3.14159265358979323846  ! PI
  integer, parameter :: STRLNG=80                  ! max. string length
  integer, parameter :: LINELNG=160
  integer, parameter :: STDERR=0
  integer, parameter :: STDOUT=6

  character(len=20), parameter :: sfxin = ".in"
  character(len=20), parameter :: sfxout = ".out"

end module constant

!
! Authors:  J.Hamm, S.Scholz, A.Klaedtke
! Modified: 4/12/2007
!
!======================================================================


