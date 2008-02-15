!-*- F90 -*------------------------------------------------------------
!
!  module: manifest / meta
!
!----------------------------------------------------------------------

!======================================================================
!
!
!

module manifest

  use constant
  use mat
  use diag
  use out
  implicit none
  public
  save 

  character(len=STRLNG),parameter :: VERSION_NUMBER = "M4_VERSION"
  character(len=STRLNG),parameter :: VERSION_DATE = "M4_DATE"
  character(len=STRLNG),parameter :: VERSION_AUTHORS = &
       "M4_AUTHORS"
  character(len=STRLNG),parameter :: BUILD= &
       "M4_BUILD"
  character(len=STRLNG),parameter :: FLAVOUR= &
        "M4_FLAVOUR"

contains

!----------------------------------------------------------------------

  subroutine DisplayManifest

  
    write(6,*) "* BUILD INFO" 
    write(6,*) "*   %VERSION : ",TRIM(VERSION_NUMBER), " (", TRIM(VERSION_DATE),")"
    write(6,*) "*   %BUILD   : ",TRIM(BUILD) 
    write(6,*) "*   %FLAVOUR : ",TRIM(FLAVOUR) 
    write(6,*) "* MODULE INFO"
    write(6,*) "*   %MAT     : M4_MATLIST"
    write(6,*) "*   %DIAG    : M4_DIAGLIST"
    write(6,*) "*   %OUTGPL  : M4_OUTGPLLIST"
    write(6,*) "*   %OUTVTK  : M4_OUTVTKLIST"

  end subroutine DisplayManifest

!----------------------------------------------------------------------

end module manifest

!
! Authors:  J. Hamm
! Modified: 28/1/2007
!
!======================================================================


