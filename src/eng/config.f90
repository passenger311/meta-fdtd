!-*- F90 -*------------------------------------------------------------
!
!  module: config / meta3
!
!  config file reader. 
!
!  subs:
!
!    ReadConfig
!
!----------------------------------------------------------------------


! =====================================================================
!


module config

  use strings
  use constant

  use grid
  use fdtd
  use pml
  use mat 
  use diag

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'CONFIG'

  ! --- Public Methods

  public :: ReadConfig

  ! --- Constants

  character(len=STRLNG), parameter :: pfxconfig = 'config'

  ! --- Data

contains

!----------------------------------------------------------------------


  subroutine ReadConfig
      
      character(len=STRLNG) :: file, string
      integer :: ios
      
      M4_WRITE_DBG({". enter ReadConfig/config"})

      file=cat2(pfxconfig,mpi_sfxin)
      
      M4_WRITE_DBG({"trying to open ", TRIM(file)})

      open(UNITTMP,FILE=file,STATUS="old", IOSTAT=ios)
      M4_OPEN_ERROR(ios,file)

      do
         read(UNITTMP,*, IOSTAT=ios) string

         if(ios .ne. 0) exit
 
         select case ( string )

         case( "(GRID" )
            M4_WRITE_DBG({"got ",TRIM(string),"-> invoking ReadConfigGrid"})
            call ReadConfigGrid(UNITTMP,string)
            
         case( "(FDTD" )
            M4_WRITE_DBG({"got ",TRIM(string),"-> invoking ReadConfigFdtd"})
            call ReadConfigFdtd(UNITTMP,string)

         case( "(PML" )
            M4_WRITE_DBG({"got ",TRIM(string),"-> invoking ReadConfigPml"})
            call ReadConfigPml(UNITTMP,string)   

         case("")
         case default

            if ( string(1:4) .eq. "(MAT" ) then
            M4_WRITE_DBG({"got ",TRIM(string),"-> invoking ReadConfigMat"})
               call ReadConfigMat(UNITTMP,string)
               cycle
            endif
            
            if ( string(1:5) .eq. "(DIAG" ) then
               M4_WRITE_DBG({"got ",TRIM(string),"-> invoking ReadConfigDiag"})
               call ReadConfigDiag(UNITTMP,string)
               cycle
            endif

            if ( string(1:1) .ne. "!" ) then
               M4_FATAL_ERROR({"BAD TOKEN ", TRIM(string) ,": ReadConfig/mat"})
            endif

         end select
         
      enddo
      close(UNITTMP)

      M4_WRITE_DBG({". exit ReadConfig/config"})

    end subroutine ReadConfig
      

!----------------------------------------------------------------------

  end module config

! =====================================================================






