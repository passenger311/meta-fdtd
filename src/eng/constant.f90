!-*- F90 -*------------------------------------------------------------
!
!  module: constant / meta3
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
 
  integer, parameter :: UNITTMP=10                 ! File number
  real(8), parameter :: PI=3.14159265358979323846  ! PI
  integer, parameter :: STRLNG=80                  ! max. Stringlength
  integer, parameter :: STDERR=0                   ! auf Himiko 6, sonst 0

  character(len=20), parameter :: sfxin = ".in"
  character(len=20), parameter :: sfxout = ".out"

end module constant

!
! Authors:  S.Scholz, J.Hamm, A.Klaedtke
! Modified: 4/12/2007
!
!======================================================================


