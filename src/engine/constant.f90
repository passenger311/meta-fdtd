!-*- F90 -*------------------------------------------------------------
!
!  module: constant / meta
!
!  global parameters.
!
!  CF,1D,2D,3D
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
  integer, parameter :: STDERR=0                   ! 6 for sr8000k, otherwise 0
  integer, parameter :: STDOUT=6                   ! 

  character(len=STRLNG),parameter :: VERSION_NUMBER = "0.07.12.25"
  character(len=STRLNG),parameter :: VERSION_DATE = "2007-2008"
  character(len=STRLNG),parameter :: VERSION_AUTHORS = &
       "J.Hamm, A.Klaedtke, S.Scholz, C.Hermann"
  character(len=STRLNG),parameter :: BUILD_ARCH= &
       "M4_ARCH() M4_IFELSE_DBG({DBG })M4_IFELSE_OMP({OMP })M4_IFELSE_MPI({MPI })M4_IFELSE_MPELOG({MPELOG })" 
   character(len=STRLNG),parameter :: BUILD_FLAGS= &
        "M4_SDIM({D })M4_IFELSE_CF({CF })M4_IFELSE_NG({NG })M4_IFELSE_TE({TE })" 

  character(len=20), parameter :: sfxin = ".in"
  character(len=20), parameter :: sfxout = ".out"

contains

!----------------------------------------------------------------------

  subroutine DisplayVersion

  
    write(6,*) "* BUILD INFO" 
    write(6,*) "*   %VERSION   : ",TRIM(VERSION_NUMBER), " (", TRIM(VERSION_DATE),")"
    write(6,*) "*   %ARCH      : ",TRIM(BUILD_ARCH) 
    write(6,*) "*   %FLAGS     : ",TRIM(BUILD_FLAGS) 
    write(6,*) "* MODULE INFO"
    write(6,*) "*   %MAT       : M4_MATLIST"
    write(6,*) "*   %DIAG      : M4_DIAGLIST"
    write(6,*) "*   %OUTGPL    : M4_OUTGPLLIST"

  end subroutine DisplayVersion

!----------------------------------------------------------------------

end module constant

!
! Authors:  S.Scholz, J.Hamm, A.Klaedtke
! Modified: 4/12/2007
!
!======================================================================


