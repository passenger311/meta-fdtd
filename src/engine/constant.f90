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
  integer, parameter :: MAXVALOBJ=100
  integer, parameter :: MAXBUFOBJ=500
  integer, parameter :: MAXVALUES=10
 
  integer, parameter :: UNITTMP=10                 ! file unit number
  real(8), parameter :: PI=3.14159265358979323846  ! PI
  integer, parameter :: STRLNG=80                  ! max. string length
  integer, parameter :: STDERR=0                   ! 6 for sr8000k, otherwise 0
  integer, parameter :: STDOUT=6                   ! 

  character(len=STRLNG),parameter :: VERSION_NUMBER = "M4_VERSION"
  character(len=STRLNG),parameter :: VERSION_DATE = "M4_DATE"
  character(len=STRLNG),parameter :: VERSION_AUTHORS = &
       "M4_AUTHORS"
  character(len=STRLNG),parameter :: BUILD= &
       "M4_BUILD"
  character(len=STRLNG),parameter :: FLAVOUR= &
        "M4_FLAVOUR"

  character(len=20), parameter :: sfxin = ".in"
  character(len=20), parameter :: sfxout = ".out"

contains

!----------------------------------------------------------------------

  subroutine DisplayVersion

  
    write(6,*) "* BUILD INFO" 
    write(6,*) "*   %VERSION : ",TRIM(VERSION_NUMBER), " (", TRIM(VERSION_DATE),")"
    write(6,*) "*   %BUILD   : ",TRIM(BUILD) 
    write(6,*) "*   %FLAVOUR : ",TRIM(FLAVOUR) 
    write(6,*) "* MODULE INFO"
    write(6,*) "*   %MAT     : M4_MATLIST"
    write(6,*) "*   %DIAG    : M4_DIAGLIST"
    write(6,*) "*   %OUTGPL  : M4_OUTGPLLIST"
    write(6,*) "*   %OUTVTK  : M4_OUTVTKLIST"

  end subroutine DisplayVersion

!----------------------------------------------------------------------

end module constant

!
! Authors:  S.Scholz, J.Hamm, A.Klaedtke
! Modified: 4/12/2007
!
!======================================================================


