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

  integer, parameter :: MAXOUTOBJ=1000
  integer, parameter :: MAXREGOBJ=1000
  integer, parameter :: MAXINITOBJ=100
  integer, parameter :: MAXBUFOBJ=500
 
  integer, parameter :: UNITTMP=10                 ! File number
  real(8), parameter :: PI=3.14159265358979323846  ! PI
  integer, parameter :: STRLNG=80                  ! max. Stringlength
  integer, parameter :: STDERR=0                   ! auf Himiko 6, sonst 0

  character(len=40), parameter :: VERSION_NUMBER = "1.1.1.1"
  character(len=40), parameter :: VERSION_DATE = "2007-2008"
  character(len=40), parameter :: VERSION_AUTHORS = "J.Hamm, A.Klaedtke, S.Scholz, C.Hermann"

  character(len=20), parameter :: sfxin = ".in"
  character(len=20), parameter :: sfxout = ".out"

contains

  subroutine DisplayVersionLine

    write(6,*) "* . version: ", TRIM(VERSION_NUMBER), " (", TRIM(VERSION_DATE),") ."

  end subroutine DisplayVersionLine


end module constant

!
! Authors:  S.Scholz, J.Hamm, A.Klaedtke
! Modified: 4/12/2007
!
!======================================================================


