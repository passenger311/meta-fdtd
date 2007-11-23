!----------------------------------------------------------------------
!
!  module: constant / max3d
!
!  global parameters.
!
!----------------------------------------------------------------------

module constant

  implicit none
  save 
 
  integer, parameter :: UNITTMP=10                 ! File number
  real(8), parameter :: PI=3.14159265358979323846  ! PI
  integer, parameter :: STRLNG=80                  ! max. Stringlength
  integer, parameter :: STDERR=0                   ! auf Himiko 6, sonst 0

  character(len=255), parameter :: sfxin = '.in'
  character(len=255), parameter :: sfxout = '.out'


end module constant
