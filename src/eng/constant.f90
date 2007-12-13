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
 
  integer, parameter :: UNITTMP=10                 ! file unit number
  real(8), parameter :: PI=3.14159265358979323846  ! PI
  integer, parameter :: STRLNG=80                  ! max. string length
  integer, parameter :: STDERR=0                   ! auf Himiko 6, sonst 0
  integer, parameter :: STDOUT=6                   ! 

  character(len=STRLNG),parameter :: VERSION_NUMBER = "1.1.1.91"
  character(len=STRLNG),parameter :: VERSION_DATE = "2007-2008"
  character(len=STRLNG),parameter :: VERSION_AUTHORS = "J.Hamm, A.Klaedtke, S.Scholz, C.Hermann"
  character(len=STRLNG),parameter :: BUILD_ARCH= &
       "M4_ARCH M4_IFELSE_DBG({DBG })M4_IFELSE_OMP({OMP })M4_IFELSE_MPI({MPI })M4_IFELSE_MPELOG({MPELOG })" 
   character(len=STRLNG),parameter :: BUILD_FLAGS= &
        "M4_IFELSE_CF({CF })M4_IFELSE_NG({NG })M4_IFELSE_TE({TE })" 

  character(len=20), parameter :: sfxin = ".in"
  character(len=20), parameter :: sfxout = ".out"

contains

  subroutine DisplayVersion

    write(6,*) "* . VERSION     : ", TRIM(VERSION_NUMBER), " (", TRIM(VERSION_DATE),") ."
    write(6,*) "* . BUILD ARCH  : ",TRIM(BUILD_ARCH) 
    write(6,*) "* . BUILD FLAGS : ",TRIM(BUILD_FLAGS) 
    write(6,*) "* . MAT MODS    : M4_MATLIST"
    write(6,*) "* . DIAG MODS   : M4_DIAGLIST"

  end subroutine DisplayVersion


end module constant

!
! Authors:  S.Scholz, J.Hamm, A.Klaedtke
! Modified: 4/12/2007
!
!======================================================================


